from pyiron_atomistics.lammps.lammps import Lammps
from pyiron_atomistics.lammps.base import LammpsBase
from pyiron_atomistics.atomistics.job.interactive import GenericInteractive
from pyiron_atomistics.lammps.structure import UnfoldingPrism
from pyiron_base import state, Executable
from pyiron_atomistics.atomistics.structure.atoms import ase_to_pyiron, pyiron_to_ase

from pyscal_rdf.rdfsystem import System
from pyscal_rdf import StructureGraph
from rdflib import Graph, Literal, Namespace, XSD, RDF, RDFS, BNode, URIRef, FOAF, SKOS, DCTERMS

import numpy as np
import os
import copy
import ast
import uuid

PROV = Namespace("http://www.w3.org/ns/prov#")
CMSO = Namespace("http://purls.helmholtz-metadaten.de/cmso/")
PODO = Namespace("http://purls.helmholtz-metadaten.de/podo/")
ASO = Namespace("http://purls.helmholtz-metadaten.de/aso/")

class RDFLammps(Lammps):
    def __init__(self, project, job_name):
        super().__init__(project, job_name)
        #self._executable_activate(enforce=True)
        self._method = None
        dbfile = os.path.join(self.project.path, "project.db")
        self.graph = StructureGraph(store="SQLAlchemy", store_file=dbfile)
        self.graph.graph.bind("prov", PROV)    
        self._initial_sample = None
        self._final_sample = None
        self._initial_structure = None    
    
    def _executable_activate(self, enforce = True):
        if self._executable is None or enforce:
            self._executable = Executable(
                codename='lammps',
                module='lammps',
                path_binary_codes=state.settings.resource_paths
            )
            if self.server.cores > 1:
                self._executable.mpi = True

    @property
    def structure(self):
        #struct = GenericInteractive.structure.fget(self)
        #struct = pyiron_to_ase(struct)
        return GenericInteractive.structure.fget(self)


    @structure.setter
    def structure(self, structure):
        if isinstance(structure, System):
            self._initial_sample = structure.sample
            self._initial_structure = structure
            ase_structure = structure.write.ase()
            pyiron_structure = ase_to_pyiron(ase_structure)
        else:
            pyiron_structure = structure
        self._prism = UnfoldingPrism(pyiron_structure.cell)
        GenericInteractive.structure.fset(self, pyiron_structure)
    
    def calc_minimize(
        self,
        ionic_energy_tolerance=0.0,
        ionic_force_tolerance=1e-4,
        e_tol=None,
        f_tol=None,
        max_iter=1000000,
        pressure=None,
        n_print=100,
        style="cg",
    ):
        super(Lammps, self).calc_minimize(
            ionic_energy_tolerance=ionic_energy_tolerance,
            ionic_force_tolerance=ionic_force_tolerance,
            e_tol=e_tol,
            f_tol=f_tol,
            max_iter=max_iter,
            pressure=pressure,
            n_print=n_print,
            style=style,
        )
        self._method = "md_min"
        self._pressure = pressure
    
    def calc_md(
        self,
        temperature=None,
        pressure=None,
        n_ionic_steps=1000,
        time_step=1.0,
        n_print=100,
        temperature_damping_timescale=100.0,
        pressure_damping_timescale=1000.0,
        seed=None,
        tloop=None,
        initial_temperature=None,
        langevin=False,
        delta_temp=None,
        delta_press=None,
    ):
        super(Lammps, self).calc_md(
            temperature=temperature,
            pressure=pressure,
            n_ionic_steps=n_ionic_steps,
            time_step=time_step,
            n_print=n_print,
            temperature_damping_timescale=temperature_damping_timescale,
            pressure_damping_timescale=pressure_damping_timescale,
            seed=seed,
            tloop=tloop,
            initial_temperature=initial_temperature,
            langevin=langevin,
            delta_temp=delta_temp,
            delta_press=delta_press)
        if pressure is None:
            self._method = "md_nvt"
        else:
            self._method = "md_npt"
        self._pressure = pressure
        self._temperature = temperature

    def get_structure_as_system(self):
        fstruct = self.get_structure().to_ase()
        if self._method == "md_npt":
            #rescale volume
            fstruct.cell = np.mean(self.output.cells, axis=0)

        #now change the system
        final_structure = System(fstruct, format="ase")
        #we have to read in the info; what it  is missing is the sysdict
        final_structure.atoms._lattice_constant = self._initial_structure.atoms._lattice_constant
        final_structure.atoms._lattice = self._initial_structure.atoms._lattice
        final_structure._structure_dict = copy.copy(self._initial_structure._structure_dict)
        return final_structure

    def _get_potential_doi(self):
        df = self.potential
        df = df.loc[0]
        if "Citations" in df.keys():
            potdict = ast.literal_eval(df["Citations"][1:-1])
            potential = "https://doi.org/"+potdict[list(potdict.keys())[0]]["doi"]
        else:
            potential = df.Name
        return potential

    def _get_potential_string(self):
        df = self.potential
        df = df.loc[0]
        ps = df["Config"][0]
        return ps.strip().split()

    def _add_initial_structure_to_graph(self):
        #-------------------------------------------------
        # Step 1: Add initial structure if does not exist
        #-------------------------------------------------
        if self._initial_sample not in self.graph.samples:
            #add to graph
            self.graph.add_structure_to_graph(self._initial_structure)
            #reset the sample id
            self._initial_sample = self.graph.sample

    def _add_final_structure_to_graph(self):
        #-------------------------------------------------
        # Step 2: Add final structure
        #-------------------------------------------------        
        fstruct = self.get_structure()
        #now change the system
        final_structure = System(fstruct, format="ase")
        #we have to read in the info; what it  is missing is the sysdict
        final_structure.atoms._lattice_constant = self._initial_structure.atoms._lattice_constant
        final_structure._structure_dict = copy.copy(self._initial_structure._structure_dict)
        #final structure is prepared, now we can simply add it to the system
        #this will make sure there is a new id for the system
        self.graph.add_structure_to_graph(final_structure)
        self._final_sample = self.graph.sample

    def _add_inherited_properties(self):
        #Here we need to add inherited info: CalculatedProperties will be lost
        #Defects will be inherited
        #add vac stuff
        material = list([k[2] for k in self.graph.graph.triples((self._initial_sample, CMSO.hasMaterial, None))])[0]
        defects = list([x[2] for x in self.graph.graph.triples((material, CMSO.hasDefect, None))])
        #now for each defect we copy add this to the final sample
        final_material = list([k[2] for k in self.graph.graph.triples((self._final_sample, CMSO.hasMaterial, None))])[0]
        
        for defect in defects:
            new_defect = URIRef(defect.toPython())
            self.graph.graph.add((final_material, CMSO.hasDefect, new_defect))
            #now fetch all defect based info
            for triple in self.graph.graph.triples((defect, None, None)):
                self.graph.graph.add((new_defect, triple[1], triple[2]))
        
        #now add the special props for vacancy
        initial_simcell = self.graph.graph.value(self._initial_sample, CMSO.hasSimulationCell)
        final_simcell = self.graph.graph.value(self._final_sample, CMSO.hasSimulationCell) 
        for triple in self.graph.graph.triples((initial_simcell, PODO.hasVacancyConcentration, None)):
            self.graph.graph.add((final_simcell, triple[1], triple[2]))
        for triple in self.graph.graph.triples((initial_simcell, PODO.hasNumberOfVacancies, None)):
            self.graph.graph.add((final_simcell, triple[1], triple[2]))

    def _add_provo_mappings(self):
        #-------------------------------------------------
        # Step 3: Add PROV-O mappings
        #-------------------------------------------------        
        #add samples as Entity and their relationship
        self.graph.add((self._initial_sample, RDF.type, PROV.Entity))
        self.graph.add((self._final_sample, RDF.type, PROV.Entity))
        self.graph.add((self._final_sample, PROV.wasDerivedFrom, self._initial_sample))
        
        #add the process that generated the samples
        main_id = uuid.uuid4()
        activity = URIRef(f'activity:{main_id}')
        self.graph.add((activity, RDF.type, PROV.Activity))
        self.graph.add((activity, RDF.type, ASO.StructureOptimization))

        method = URIRef(f'method:{main_id}')
        self.graph.add((method, RDF.type, ASO.MolecularDynamics))
        self.graph.add((activity, ASO.hasMethod, method))

        self.graph.add((activity, ASO.hasRelaxationDOF, ASO.AtomicPosition))
        if self._pressure is None:
            pass
        elif np.isscalar(self._pressure):
            self.graph.add((activity, ASO.hasRelaxationDOF, ASO.CellVolume))
        else: 
            #check if pressure is hydrostatic or not
            axial_all_alike = None not in self._pressure[:3] and np.allclose(
                self._pressure[:3], self._pressure[0]
            )
            shear_all_none = all(p is None for p in self._pressure[3:])
            shear_all_zero = None not in self._pressure[3:] and np.allclose(self._pressure[3:], 0)
            hydrostatic = axial_all_alike and (shear_all_none or shear_all_zero)
            if hydrostatic:
                self.graph.add((activity, ASO.hasRelaxationDOF, ASO.CellVolume))
            else:
                self.graph.add((activity, ASO.hasRelaxationDOF, ASO.CellVolume))
                self.graph.add((activity, ASO.hasRelaxationDOF, ASO.CellShape))

        if self._method == "md_min":
            #depending on the nature of pressure, we have to add different DOF
            pass
        elif self._method == "md_npt":
            self.graph.add((method, ASO.hasStatisticalEnsemble, ASO.NPT))
            temperature = URIRef(f'temperature:{main_id}')
            self.graph.add((temperature, RDF.type, ASO.InputParameter))
            self.graph.add((temperature, RDFS.label, Literal('temperature', datatype=XSD.string)))
            self.graph.add((activity, ASO.hasInputParameter, temperature))
            self.graph.add((temperature, ASO.hasValue, Literal(self._temperature, datatype=XSD.float)))
            self.graph.add((temperature, ASO.hasUnit, URIRef('http://qudt.org/vocab/unit/K')))

            pressure = URIRef(f'pressure:{main_id}')
            self.graph.add((pressure, RDF.type, ASO.InputParameter))
            self.graph.add((pressure, RDFS.label, Literal('pressure', datatype=XSD.string)))
            self.graph.add((activity, ASO.hasInputParameter, pressure))
            self.graph.add((pressure, ASO.hasValue, Literal(self._pressure, datatype=XSD.float)))
            self.graph.add((pressure, ASO.hasUnit, URIRef('http://qudt.org/vocab/unit/GigaPA')))
        
        elif self._method == "md_nvt":
            self.graph.add((method, ASO.hasStatisticalEnsemble, ASO.NVT))
            temperature = URIRef(f'temperature:{main_id}')
            self.graph.add((temperature, RDF.type, ASO.InputParameter))
            self.graph.add((temperature, RDFS.label, Literal('temperature', datatype=XSD.string)))
            self.graph.add((activity, ASO.hasInputParameter, temperature))
            self.graph.add((temperature, ASO.hasValue, Literal(self._temperature, datatype=XSD.float)))
            self.graph.add((temperature, ASO.hasUnit, URIRef('http://qudt.org/vocab/unit/K')))
        
        potential_string = self._get_potential_string()
        potential_string = " ".join(potential_string.split()[1:])
        potential = URIRef(f'potential:{main_id}')

        if 'meam' in potential_string:
            self.graph.add((potential, RDF.type, ASO.MEAM))
        elif 'eam' in potential_string:
            self.graph.add((potential, RDF.type, ASO.EAM))
        elif 'lj' in potential_string:
            self.graph.add((potential, RDF.type, ASO.LennardJones))
        elif 'ace' in potential_string:
            self.graph.add((potential, RDF.type, ASO.MLPotential))
        else:
            self.graph.add((potential, RDF.type, ASO.InteratomicPotential))

        self.graph.add((potential, ASO.hasReference, Literal(self._get_potential_doi())))
        self.graph.add((self._final_sample, PROV.wasGeneratedBy, activity))
        
        #add pyiron
        pyiron_agent = URIRef("http://demo.fiz-karlsruhe.de/matwerk/E457491")
        self.graph.add((pyiron_agent, RDF.type, PROV.SoftwareAgent))
        self.graph.add((activity, PROV.wasAssociatedWith, pyiron_agent))
        self.graph.add((pyiron_agent, RDFS.label, Literal("pyiron")))
        
        #add lammps
        lammps_agent = URIRef("http://demo.fiz-karlsruhe.de/matwerk/E447986")
        self.graph.add((lammps_agent, RDF.type, PROV.SoftwareAgent))
        self.graph.add((pyiron_agent, PROV.actedOnBehalfOf, lammps_agent))
        self.graph.add((lammps_agent, RDFS.label, Literal("LAMMPS")))


    def _add_calculated_properties(self):
        #now add calculated quantity
        #if self._method == "MinimizationMD":
        self.graph.add_calculated_quantity("TotalEnergy", 
            self.output.energy_tot[-1],
            unit='EV', sample=self._final_sample) 

    def collect_rdf(self):
        self._add_initial_structure_to_graph()
        self._add_final_structure_to_graph()
        self._add_inherited_properties()
        self._add_provo_mappings()
        self._add_calculated_properties()
            

    def collect_output(self):
        super(Lammps, self).collect_output()
        #create an updated graph for the final structure, we use rdf_graph directly
        self.collect_rdf()
        

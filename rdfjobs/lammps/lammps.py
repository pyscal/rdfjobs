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

PROV = Namespace("http://www.w3.org/ns/prov#")
CMSO = Namespace("http://purls.helmholtz-metadaten.de/cmso/")
PODO = Namespace("http://purls.helmholtz-metadaten.de/podo/")

#Temporary namespace, which needs to be replaced
MSMO = Namespace("http://purls.helmholtz-metadaten.de/msmo/")

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
        self._method = "MinimizationMD"
    
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
            self._method = "NVTMD"
        else:
            self._method = "NPTMD"
        self._pressure = pressure

    def get_structure_as_system(self):
        fstruct = self.get_structure().to_ase()
        if self._method == "NPTMD":
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
            potdict = ast.literal_eval(df.loc[0]["Citations"][1:-1])
            potential = "https://doi.org/"+potdict[list(potdict.keys())[0]]["doi"]
        else:
            potential = df.Name
        return potential

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
            new_defect = BNode()
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
        method = URIRef(self.name)
        if self._method == "MinimizationMD":
            self.graph.add((method, RDF.type, MSMO.MinimizationMD))
        elif self._method == "NPTMD":
            self.graph.add((method, RDF.type, MSMO.NPTMD))
            self.graph.add((method, MSMO.hasPressure, self._pressure))
            self.graph.add((method, MSMO.hasTemperature, self.output.temperature))

        elif self._method == "NVTMD":
            self.graph.add((method, RDF.type, MSMO.NVTMD))
            self.graph.add((method, MSMO.hasTemperature, self.output.temperature))

        self.graph.add((method, MSMO.usesPotential, self._get_potential_doi()))
        self.graph.add((method, RDF.type, PROV.Activity))
        self.graph.add((self._final_sample, PROV.wasGeneratedBy, method))
        
        #add pyiron
        pyiron_agent = URIRef("http://demo.fiz-karlsruhe.de/matwerk/E457491")
        self.graph.add((pyiron_agent, RDF.type, PROV.SoftwareAgent))
        self.graph.add((method, PROV.wasAssociatedWith, pyiron_agent))
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
        

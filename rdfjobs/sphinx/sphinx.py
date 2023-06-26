from pyiron_atomistics.sphinx.sphinx import Sphinx
from pyiron_atomistics.dft.job.generic import GenericDFTJob
from pyiron_atomistics.vasp.potential import (
    VaspPotentialFile,
    strip_xc_from_potential_name,
    find_potential_file as find_potential_file_vasp,
    VaspPotentialSetter,
)
from pyiron_base import state, Executable
from pyiron_atomistics.atomistics.structure.atoms import ase_to_pyiron, pyiron_to_ase

from pyscal_rdf.rdfsystem import System
from pyscal_rdf import StructureGraph
from rdflib import Graph, Literal, Namespace, XSD, RDF, RDFS, BNode, URIRef, FOAF, SKOS, DCTERMS

import numpy as np
import os
import copy

PROV = Namespace("http://www.w3.org/ns/prov#")
CMSO = Namespace("https://purls.helmholtz-metadaten.de/cmso/")
PODO = Namespace("https://purls.helmholtz-metadaten.de/podo/")

class RDFSphinx(Sphinx):
    def __init__(self, project, job_name, dbfile=None):
        super().__init__(project, job_name)
        self._method = None
        if dbfile is None:
            dbfile = os.path.join(self.project.path, "project.db")
        self.graph = StructureGraph(store="SQLAlchemy", store_file=dbfile)
        self.graph.graph.bind("prov", PROV)    
        self._initial_sample = None
        self._final_sample = None
        self._initial_structure = None    

    def _executable_activate(self, enforce = True):
        if self._executable is None or enforce:
            self._executable = Executable(
                codename='sphinx',
                module='sphinx',
                path_binary_codes=state.settings.resource_paths
            )

    @property
    def structure(self):
        #struct = GenericInteractive.structure.fget(self)
        #struct = pyiron_to_ase(struct)
        return GenericDFTJob.structure.fget(self)


    @structure.setter
    def structure(self, structure):
        if isinstance(structure, System):
            self._initial_sample = structure.sample
            self._initial_structure = structure
            ase_structure = structure.to_ase()
            pyiron_structure = ase_to_pyiron(ase_structure)
        else:
            pyiron_structure = structure
        GenericDFTJob.structure.fset(self, pyiron_structure)
        if pyiron_structure is not None:
            self._potential = VaspPotentialSetter(
                element_lst=pyiron_structure.get_species_symbols().tolist()
            )

    def calc_minimize(
        self,
        electronic_steps=60,
        ionic_steps=None,
        max_iter=None,
        pressure=None,
        algorithm=None,
        retain_charge_density=False,
        retain_electrostatic_potential=False,
        ionic_energy=None,
        ionic_forces=None,
        ionic_energy_tolerance=None,
        ionic_force_tolerance=None,
        volume_only=False,
    ):

        super(Sphinx, self).calc_minimize(
            electronic_steps=electronic_steps,
            ionic_steps=ionic_steps,
            max_iter=max_iter,
            pressure=pressure,
            algorithm=algorithm,
            retain_charge_density=retain_charge_density,
            retain_electrostatic_potential=retain_electrostatic_potential,
            ionic_energy=ionic_energy,
            ionic_forces=ionic_forces,
            ionic_energy_tolerance=ionic_energy_tolerance,
            ionic_force_tolerance=ionic_force_tolerance,
            volume_only=volume_only,
        )
        self._method = "MinimizationDFT"

    def get_structure_as_system(self):
        fstruct = self.get_structure()
        #now change the system
        final_structure = System(fstruct, format="ase")
        #we have to read in the info; what it  is missing is the sysdict
        final_structure.atoms._lattice_constant = self._initial_structure.atoms._lattice_constant
        final_structure._structure_dict = copy.copy(self._initial_structure._structure_dict)
        final_structure.atoms._lattice = self._initial_structure.atoms._lattice
        return final_structure

    def collect_rdf(self):
        #-------------------------------------------------
        # Step 1: Add initial structure if does not exist
        #-------------------------------------------------
        if self._initial_sample not in self.graph.samples:
            #add to graph
            self.graph.add_structure_to_graph(self._initial_structure)
            #reset the sample id
            self._initial_sample = self.graph.sample
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


        #-------------------------------------------------
        # Step 3: Add PROV-O mappings
        #-------------------------------------------------        
        #add samples as Entity and their relationship
        self.graph.add((self._initial_sample, RDF.type, PROV.Entity))
        self.graph.add((self._final_sample, RDF.type, PROV.Entity))
        self.graph.add((self._final_sample, PROV.wasDerivedFrom, self._initial_sample))
        #add the process that generated the samples
        method = URIRef(f'http://example.org/{self._method}')
        self.graph.add((method, RDF.type, PROV.Activity))
        self.graph.add((self._final_sample, PROV.wasGeneratedBy, method))
        #add pyiron
        pyiron_agent = URIRef("http://demo.fiz-karlsruhe.de/matwerk/E457491")
        self.graph.add((pyiron_agent, RDF.type, PROV.SoftwareAgent))
        self.graph.add((method, PROV.wasAssociatedWith, pyiron_agent))
        self.graph.add((pyiron_agent, RDFS.label, Literal("pyiron")))
        #add lammps
        sphinx_agent = URIRef("http://demo.fiz-karlsruhe.de/matwerk/E435779")
        self.graph.add((sphinx_agent, RDF.type, PROV.SoftwareAgent))
        self.graph.add((pyiron_agent, PROV.actedOnBehalfOf, sphinx_agent))
        self.graph.add((sphinx_agent, RDFS.label, Literal("Sphinx")))
        
        #now add calculated quantity
        if self._method == "MinimizationDFT":
            self.graph.add_calculated_quantity("TotalEnergy", 
                self.output.energy_tot[-1],
                unit='EV', sample=self._final_sample)         

    def collect_output(self):
        super(Sphinx, self).collect_output()
        #create an updated graph for the final structure, we use rdf_graph directly
        self.collect_rdf()
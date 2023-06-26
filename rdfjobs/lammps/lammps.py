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

PROV = Namespace("http://www.w3.org/ns/prov#")
CMSO = Namespace("https://purls.helmholtz-metadaten.de/cmso/")
PODO = Namespace("https://purls.helmholtz-metadaten.de/podo/")

class RDFLammps(Lammps):
    def __init__(self, project, job_name, dbfile=None):
        super().__init__(project, job_name)
        #self._executable_activate(enforce=True)
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
                codename='lammps',
                module='lammps',
                path_binary_codes=state.settings.resource_paths
            )

    @property
    def structure(self):
        #struct = GenericInteractive.structure.fget(self)
        #struct = pyiron_to_ase(struct)
        return GenericInteractive.structure.fget(self)


    @structure.setter
    def structure(self, structure):
        if isninstance(structure, System):
            self._initial_sample = structure.sample
            self._initial_structure = structure
            ase_structure = structure.to_ase()
            structure = ase_to_pyiron(ase_structure)
        self._prism = UnfoldingPrism(structure.cell)
        GenericInteractive.structure.fset(self, structure)
    
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
        final_structure = System(atoms, format="ase")
        #we have to read in the info; what it  is missing is the sysdict
        final_structure.atoms._lattice_constant = self._initial_structure._lattice_constant
        final_structure._structure_dict = copy.copy(self._initial_structure._structure_dict)
        #final structure is prepared, now we can simply add it to the system
        #this will make sure there is a new id for the system
        self.graph.add_structure_to_graph(final_structure)
        self._final_sample = self.graph.sample
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
        lammps_agent = URIRef("http://demo.fiz-karlsruhe.de/matwerk/E447986")
        self.graph.add((lammps_agent, RDF.type, PROV.SoftwareAgent))
        self.graph.add((pyiron_agent, PROV.actedOnBehalfOf, lammps_agent))
        self.graph.add((lammps_agent, RDFS.label, Literal("LAMMPS")))
        #

 
            

    def collect_output(self):
        super(Lammps, self).collect_output()
        #create an updated graph for the final structure, we use rdf_graph directly
        self.collect_rdf()
        

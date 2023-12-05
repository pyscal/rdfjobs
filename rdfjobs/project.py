import os
from functools import partial, update_wrapper

#pyiron imports
from pyiron_base import (
    ProjectHDFio,
    JobType,
    JobTypeChoice,
    Creator as CreatorCore,
    deprecate,
    PyironFactory,
)

from pyiron_atomistics.atomistics.generic.object_type import (
    ObjectType,
    ObjectTypeChoice,
)


from pyiron_atomistics import Project as ProjectCore


#creators for structure
from pyscal_rdf import StructureGraph

class Project(ProjectCore):
    def __init__(self, path="", user=None, sql_query=None, default_working_directory=False):
        super(Project, self).__init__(
            path=path,
            user=user,
            sql_query=sql_query,
            default_working_directory=default_working_directory
        )
        self.job_type = JobTypeChoice()
        self._graph_file = os.path.join(self.path, "project.db")
        self._structure_store = os.path.join(self.path, 'rdf_structure_store')
        self.graph = StructureGraph(store="SQLAlchemy", 
            store_file=self._graph_file, 
            structure_store=self._structure_store)
        self._creator = Creator(self)

Project.__doc__ = ProjectCore.__doc__

class Creator(CreatorCore):
    def __init__(self, project):
        super().__init__(project)
        self._structure = StructureFactory(project.graph)

    @property
    def structure(self):
        return self._structure

class StructureFactory(PyironFactory):
    def __init__(self, graph):
        self._graph = graph

    def bulk(self, 
        element,
        repetitions=None, 
        crystalstructure=None,
        a=None,
        covera=None,
        cubic=True,
        add_to_graph=True):

        if crystalstructure is None:
            crystalstructure = self._graph._element_dict[element]['structure']
            if a is None:
                a = self._graph._element_dict[element]['lattice_constant']
        
        if covera is None:
            covera = 1.633

        return self._graph._annotated_make_crystal(crystalstructure,
            repetitions=repetitions,
            lattice_constant=a,
            ca_ratio = covera,
            element = element,
            primitive = not cubic,
            add_to_graph=add_to_graph,
            )

    def grain_boundary(self,
        element,
        axis,
        sigma,
        gb_plane,
        repetitions = (1,1,1),
        crystalstructure=None,
        a=1,
        overlap=0.0,
        add_to_graph=True,
        ):

        return self._graph._annotated_make_grain_boundary(self, 
            axis,
            sigma,
            gb_plane,
            structure = crystalstructure,
            element=element,
            lattice_constant=a,
            repetitions=repetitions,
            overlap=overlap,
            add_to_graph=add_to_graph)

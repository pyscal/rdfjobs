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
        self.graph = StructureGraph(store="SQLAlchemy", store_file=self._graph_file)
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
        crystalstructure=None,
        a=None,
        covera=None,
        cubic=True):

        if crystalstructure is None:
            crystalstructure = self.graph._element_dict[element]['structure']
            if a is None:
                a = self.graph._element_dict[element]['lattice_constant']
        
        if covera = None:
            covera = 1.633

        return self.graph._annotated_make_crystal(crystalstructure,
            lattice_constant=a,
            ca_ratio = covera,
            element = element,
            primitive = not cubic,
            )

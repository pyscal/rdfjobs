from pyiron_base import (
    ProjectHDFio,
    JobType,
    JobTypeChoice,
    Creator as CreatorCore,
    deprecate,
)

from pyiron_atomistics.atomistics.generic.object_type import (
    ObjectType,
    ObjectTypeChoice,
)


from pyiron_atomistics import Project as ProjectCore

class Project(ProjectCore):
    def __init__(self, path="", user=None, sql_query=None, default_working_directory=False):
        super(Project, self).__init__(
            path=path,
            user=user,
            sql_query=sql_query,
            default_working_directory=default_working_directory
        )
        self.job_type = JobTypeChoice()
        #self._creator = Creator(self)

Project.__doc__ = ProjectCore.__doc__

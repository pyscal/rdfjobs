#from rdfjobs.lammps.lammps import RDFLammps
#from rdfjobs.sphinx.sphinx import RDFSphinx

from pyiron_base import JOB_CLASS_DICT
from rdfjobs.project import Project

JOB_CLASS_DICT["RDFLammps"] = "rdfjobs.lammps.lammps"
JOB_CLASS_DICT["RDFSphinx"] = "rdfjobs.sphinx.sphinx"
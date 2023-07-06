# rdfjobs

A collection of [pyiron](https://pyiron.org/) Jobs with RDF elements. Structures and jobs are automatically annotated using [CMSO](https://github.com/Materials-Data-Science-and-Informatics/cmso-ontology) and [PROV-O](https://www.w3.org/TR/2013/REC-prov-o-20130430/), allowing for querying of the generated data. 

## Installation

We strongly recommend creating a conda environment for the installation. [Mamba](https://mamba.readthedocs.io/en/latest/installation.html) can allow the creation of conda environments. If you have existing conda installation, you can do:

```
conda install -c conda-forge mamba
```

Once a conda distribution is available, the following steps will help set up an environment to use `pyscal_rdf`. First step is to clone the repository.

```
git clone https://github.com/pyscal/rdfjobs.git
```

After cloning, an environment can be created from the included file-

```
cd rdfjobs
mamba env create -f environment.yml
```

This will create an environment called `rdfjobs`. It can be activated by,

```
conda activate rdfjobs
```

then, install `rdfjobs` using,

```
pip install .
```

## Using `rdfjobs`

A detailed documentation is coming soon. Please check the [examples](examples/) folder for some examples.


## Acknowledgements

This work is supported by the [NFDI-Matwerk](https://nfdi-matwerk.de/) consortia.





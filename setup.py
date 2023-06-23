from setuptools import setup, find_packages


with open('README.md') as readme_file:
    readme = readme_file.read()

setup(
    name='rdfjobs',
    version='0.0.1',
    author='Abril Azocar Guzman, Sarath Menon',
    author_email='sarath.menon@pyscal.org',
    description='pyiron jobs with RDF elements',
    long_description=readme,
    long_description_content_type='text/markdown',
    packages=find_packages(include=['rdfjobs', 'rdfjobs.*']),
    zip_safe=False,
    download_url = 'https://github.com/pyscal/pyscal_rdf',
    url = 'https://pyscal.org',
    install_requires=['numpy'],
    classifiers=[
        'Programming Language :: Python :: 3'
    ],
    include_package_data=True,
)

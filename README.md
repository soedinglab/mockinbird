# mockinbird: A fully automatic and reproducible PAR-CLIP analysis pipeline

## Installing mockinbird with conda

Clone the repository: `git clone git@github.com:soedinglab/stammp2.git mockinbird`

### Setting up the conda channels
If you have not yet registered the `conda-forge` and `bioconda` channels, please do

```bash
conda config --add channels conda-forge
conda config --add channels bioconda
```

### Installing mockinbird and its dependencies

```bash
cd mockinbird
conda install -f conda.pkgs
pip install .
```

## Building the documentation

### Installing the dependencies

```bash
conda install -f conda_doc.pkg
```

### Building the documentation

```bash
cd doc && make html
```

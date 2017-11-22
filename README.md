# mockinbird: A fully automatic and reproducible PAR-CLIP analysis pipeline

For a detailed description of *mockinbird*, please refer to the [documentation](http://wwwuser.gwdg.de/~compbiol/mockinbird).

## Installing mockinbird with bioconda (recommended)
```bash
conda create -n mockinbird -c bioconda -c conda-forge mockinbird
```

## Installing with conda from the git repository (not recommended)

Clone the repository: `git clone https://github.com/soedinglab/mockinbird.git`

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

## Building the sphinx documentation

### Installing the dependencies

```bash
conda install -f conda_doc.pkg
```

### Building the documentation

```bash
cd doc && make html
```

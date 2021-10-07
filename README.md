# QSAR Pipeline

QSAR Pipeline is an attempt to provide a fast and reproducible pipeline for
 preprocessing and computing molecular descriptors for QSAR modelling.
Molecular descriptors are computed via [RDKit](https://github.com/rdkit/rdkit)
 and [Mordred](https://github.com/mordred-descriptor/mordred).

This library was developed for my own needs, documentation is therefore lacking. 
If you have any questions or problems using it feel free to reach out and open an 
[issue](https://github.com/gfrogat/qsar-pipeline/issues).

## Setup

Install the dependencies via [`conda`](https://docs.conda.io/en/latest/):

```bash
# Option 1: create environment by hand
conda create -n "qsar-pipeline" -c conda-forge python=3.8 rdkit=2021.03.4 openbabel=3.1.1

# Activate environment
conda activate qsar-pipeline
```

then install `qsar_pipeline` via [poetry](https://python-poetry.org/).

```bash
poetry install 
```

## Usage

Usage examples can be found in the [`examples`](./examples) folder.

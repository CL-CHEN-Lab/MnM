[![MnM version](https://img.shields.io/badge/release-1.0.0-blue)](https://github.com/CL-CHEN-Lab/MnM/releases/latest)
[![Python version](https://img.shields.io/badge/Python-3-yellow?logo=python)]
[![Chen Twitter](https://img.shields.io/badge/Share%20it-black?logo=X)](https://twitter.com/TeamChenCurie)


# MnM

<img src="https://xfer.curie.fr/get/nil/DYswdDZesK2/MnM.logo.png" width="150" height="150" />

MnM, a machine learning approach to detect replication states and genomic subpopulations for single-cell DNA replication timing disentanglement from whole-genome sequencing data. It includes a single-cell copy-number imputation method, a replication state classifier and a subpopulation detector.
Requires ≥10 cells.

## Usage

```bash
# Move to the downloaded directory
cd path_to_directory
# Make MnM executable
chmod +x MnM
# See all options and explanations
./MnM -h
# Launch
./MnM -i INPUTFILE [-m] [--sep SEP] [-o OUTPUT] [-n NAME] [-g GENOME] [-w WINDOWSIZE] [--seed SEED] [--maxcells MAXCELLS] [-r] [-s] [--cpu CPU] [--CNcol CNCOL] [--Cellcol CELLCOL] [--groups GROUPS] [-p] [-b] [-v] [-h]
```

## Disclaimer

This software is provided as-is and without warranty. By using this software, you agree to abide by the terms specified in the license.

### Tested versions

This program has been tested with the following software versions. It is recommended you use the same or superior versions.
* Python 3.9.11
* pandas 1.5.0
* numpy 1.24.3
* sklearn 1.1.2
* natsort 8.2.0
* seaborn 0.12.1
* matplotlib 3.6.2
* pybedtools 0.9.0
* multiprocess 0.70.14
* tensorflow 2.13.0
* umap 0.5.3

You can run:
```bash
pip install pandas==1.5.0
pip install numpy==1.24.3
pip install scikit-learn==1.1.2
pip install natsort==8.2.0
pip install seaborn==0.12.1
pip install matplotlib==3.6.2
pip install pybedtools==0.9.0
pip install multiprocess==0.70.14
pip install tensorflow==2.13.0
pip install umap==0.5.3
```

## Contact

For inquiries and authorization requests, please contact [joseph.josephides@curie.fr](mailto:joseph.josephides@curie.fr) and [chunlong.chen@curie.fr](mailto:chunlong.chen@curie.fr).

## MnM releases

MnM v1.0.0

Author: Joseph Josephides.
Institut Curie, Paris.
Last update: 09 Aug 2023.

[![MnM version](https://img.shields.io/badge/release-1.0.0-blue)](https://github.com/CL-CHEN-Lab/MnM/releases/latest)
![Python version](https://img.shields.io/badge/Python-3-yellow?logo=python)
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

## Example with MCF-7 cells

Processed CNV data of MCF-7 cells obtained from [published data](https://www.nature.com/articles/s41467-022-30043-x) with [Kronos scRT](https://github.com/CL-CHEN-Lab/Kronos_scRT) can be used as an example. These data are included in the [scRT/scCNV atlas](https://github.com/CL-CHEN-Lab/MnM/tree/main/scRT_scCNV_Atlas). We can discover replication states and genomic heterogeneity with the following command:

```bash
MnM -i scRT_scCNV_Atlas/Gnan2022_MCF-7/Gnan2022_MCF-7_scCNV_Matrix.tsv.gz -o ~/Documents/MnM_test_MCF-7_MnM_Output -m -g hg38 -n MCF-7 -r -s
```
Where -i is the input file (scCNV matrix), -o the output directory, -m flags that the input file is a matrix and not a BAM file, -g precises the reference genome, -n names the output files, -r reports replication states and -s discovers genomic heterogeneity (subpopulations).

As a result the following figures are produced:

* Distinction of replicating cells from non-replicating cells.
<img src="https://xfer.curie.fr/get/nil/1qIlI8WRNx6/MCF-7_phases_scCNV_heatmap.png" />

* Detection of subpopulations from non-replicating cells visualised:
 *  Genome-wide at the single-cell level.
<img src="https://xfer.curie.fr/get/nil/hwtBDWumJWs/MCF-7_subpopulations_scCNV_heatmap.png" />
 * Chromosome-wide per subpopulation.
<img src="https://xfer.curie.fr/get/nil/y01Cv4WrkQJ/MCF-7_UMAP_Subpopulations.png" />
 * On a UMAP plot.
<img src="https://xfer.curie.fr/get/nil/YJgmnWjjkG5/MCF-7_subpopulations_median_CNs.png" />

And a metadata file lists the replicating state and subpopulation of each cell in a metadata file:
| Cell                             | Phase | Subpopulation |
|----------------------------------|-------|---------------|
| AAACCTGCAACCCAAT-1_First_exp.bam | G1    | 1             |
| AAACCTGGTCATTACG-1_Normal.bam    | S     | 1             |
| AAACGGGTCGGGAAAC-1_Normal.bam    | G1    | 2             |
| AAAGATGAGCTATGCT-1_Normal.bam    | S     | 1             |
| AAAGATGAGTAAGTAC-1_Normal.bam    | S     | 2             |
| AAAGATGGTTTGCCTC-1_Normal.bam    | G1    | 1             |
| CTAGTGAAGTGCTGCC-1_Normal.bam    | S     | 1             |
| ...                              | ...   | ...           |


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
## Citation

If you use MnM please cite the following pre-print:

Josephides, J. M., & Chen, C. L. (2023). MnM: a machine learning approach to detect replication states and genomic subpopulations for single-cell DNA replication timing disentanglement. bioRxiv, 2023-12.

>@article{josephides2023mnm,
>  title={MnM: a machine learning approach to detect replication states and genomic subpopulations for single-cell DNA replication timing disentanglement},
>  author={Josephides, Joseph M and Chen, Chun-Long},
>  journal={bioRxiv},
>  pages={2023--12},
>  year={2023},
>  elocation-id = {2023.12.26.573369},
>  publisher={Cold Spring Harbor Laboratory},
>  doi = {10.1101/2023.12.26.573369},
>  URL = {https://www.biorxiv.org/content/early/2023/12/28/2023.12.26.573369}
>}

## Contact

For inquiries and authorization requests, please contact [joseph.josephides@curie.fr](mailto:joseph.josephides@curie.fr) and [chunlong.chen@curie.fr](mailto:chunlong.chen@curie.fr).

## MnM releases

MnM v1.0.0

Author: Joseph Josephides.
Institut Curie, Paris.
Last update: 09 Aug 2023.

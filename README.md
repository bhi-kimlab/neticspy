# NetICSpy
[![version](https://img.shields.io/pypi/v/neticspy.svg)](https://pypi.org/project/neticspy)

Python implementation of [NetICS](https://doi.org/10.1093/bioinformatics/bty148) (Network-based Integration of Multi-omics data).

## Installation
```shell
$ pip install neticspy
```

## Description


## Usage
On the commandline:
```
$ neticspy --abberant MutationFileName --adj AdjecencyMatrix --beta RestartProb --method RankMethod --network NetworkFileName --deg DEGFileName --dep DEPFileName
```

As a ordinary python package:
```
import neticspy

ranked_list_genes, ranked_scores = neticspy.netics_fun(filenameMu=MutationFileName, filenameAdj=AdjencyMatrix, restart_prob=RestartProb, rank_method_str=RankMethod, filenameNet=NetworkFileName, filenameRNA=DEGFileName, filenamePRDEPFileName)
```

The default value of beta (restart probability) is 0.4.

Available rank aggregation methods are 'SUM', 'MEDIAN' and 'RRA'.

Input files for differentially expressed genes and proteins are optional.

## References
Dimitrakopoulos, Christos, et al. "Network-based integration of multi-omics data for prioritizing cancer genes." *Bioinformatics* 34.14 (2018): 2441-2448.

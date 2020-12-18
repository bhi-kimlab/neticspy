# NetICSpy
[![version](https://img.shields.io/pypi/v/neticspy.svg)](https://pypi.org/project/neticspy)

Python implementation of [NetICS](https://doi.org/10.1093/bioinformatics/bty148) (Network-based Integration of Multi-omics data).

## Installation
```shell
$ pip install neticspy
```

## Description
NetICS performs a per sample bidirectional network diffusion-based method for prioritizing genes based on their proximity to upstream aberration events and to downstream differentially expressed genes and proteins in an interaction network.
The sample-specific gene lists are integrated into an overall ranked list of genes using rank aggregation technique.

## Usage
On the commandline:
```
$ neticspy --aberration ABERRATION \
    --adj ADJACENCY_MATRIX \
    --network NETWORK \
    --degs DEGS \
    --output OUTPUT \
    [--beta RestartProb] \
    [--diffusion-matrix DIFFUSION_MATRIX] \
    [--opposite-diffusion-matrix OPPOSITE_DIFFUSION_MATRIX] \
```

### Arguments
`-a ABERRATION, --aberration ABERRATION`: Input two-column table (without headers) containing genetically aberrant genes for each sample. It contain two columns that map every gene (1st column) to the samples that it it genetically aberrant (2nd column).

`-j ADJACENCY_MATRIX, --adj ADJACENCY_MATRIX`: Adjacency matrix of the directed interaction network.

`-n NETWORK, --network NETWORK`: Input file that contains the list of the genes that are present in the network. They should be in the same order as in the rows of the adjacency matrix given by `--adj`. An example file is given that contains the gene names of the network described in (Wu et al., 2010).

`-d DEGS, --degs DEGS`: List of the names of predefined differentially expressed genes.

`-o OUTPUT, --output OUTPUT`: File to save NetICS result.

`-b BETA, --beta BETA`: Restart probability for the insulated diffusion. Default: 0.4 (for the network from Wu et al., 2010.)

`-f DIFFUSION_MATRIX, --diffusion-matrix DIFFUSION_MATRIX`: To accelerate the diffusion process, you can feed precomputed diffusion matrix into NetICS.

`-g OPPOSITE_DIFFUSION_MATRIX, --opposite-diffusion-matrix OPPOSITE_DIFFUSION_MATRIX`: To accelerate the diffusion process, you can feed precomputed diffusion matrix in the *opposite* direction into NetICS.

As a ordinary python package:
```
import neticspy

ranked_list_genes, ranked_scores = neticspy.netics_fun(filenameMu=MutationFileName, filenameAdj=AdjencyMatrix, restart_prob=RestartProb, rank_method_str=RankMethod, filenameNet=NetworkFileName, filenameRNA=DEGFileName, filenamePRDEPFileName)

or

neticspy.netics_fun(filenameMu=MutationFileName, filenameAdj=AdjencyMatrix, restart_prob=RestartProb, rank_method_str=RankMethod, filenameNet=NetworkFileName, filenameRNA=DEGFileName, filenamePR=DEPFileName, output=OutputFileName)
```

The default value of beta (restart probability) is 0.4.

Available rank aggregation methods are 'SUM', 'MEDIAN' and 'RRA'.

Input files for differentially expressed genes and proteins are optional.

## References
Dimitrakopoulos, Christos, et al. "Network-based integration of multi-omics data for prioritizing cancer genes." *Bioinformatics* 34.14 (2018): 2441-2448.

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
There are two subcommands: `neticspy diffuse` and `neticspy rank`.

### diffuse
Prepare diffusion matrix for given network and restart probability. This subcommand will produce `.npz` file containing precomputed forward/backward diffusion matrices.

```
$ neticspy diffuse --adj ADJ \
    --beta BETA \
    --output OUTPUT 
```
**Arguments**

`-j ADJ, --adj ADJ`: Adjacency matrix of the directed interaction network.

`-b BETA, --beta BETA`: Restart probability for the insulated diffusion. Default: 0.4 (For the network from Wu et al., 2010)

`-o OUTPUT, --output OUTPUT`: Output filename for diffusion matrix in .npz format.

### rank
Run core NetICS algorithm and rank genes by their mediator effect. This subcommand will produce two comma-separated tables: raw result containing sample-gene-diffusion score triplets (raw.txt), and aggregated gene rankings (rank_aggregated.txt) that allow cohort-wise gene prioritization.

```
$ neticspy rank --aberration ABERRATION \
    --diffusion-matrix DIFFUSION_MATRIX \
    --network NETWORK \
    --degs DEGS \
    --output OUTPUT \
```
**Arguments**

`-a ABERRATION, --aberration ABERRATION`: Input two-column table (without headers) containing genetically aberrant genes for each sample. It contain two columns that map every gene (1st column) to the samples that it is genetically aberrant (2nd column).

`-f DIFFUSION_MATRIX, --diffusion-matrix DIFFUSION_MATRIX`: Path to .npz file for diffusion matrix.

`-n NETWORK, --network NETWORK`: Input file that contains the list of the genes that are present in the network.

`-d DEGS, --degs DEGS`: List of the names of predefined differentially expressed genes.

`-o OUTPUT-PREFIX, --output-prefix OUTPUT-PREFIX`: Prefix of the output file to save raw NetICS result and aggregated ranks.

`-v, --verbose`: Increase verbosity.

**As an ordinary python package**
```
import neticspy

netics_result = neticspy.netics_fun(
    filename_aberration,
    filename_network_genes,
    output,
    filename_deg_list,
    restart_prob=0.4,
)
```

## Note
The default value of beta (restart probability) is 0.4.

## References
Dimitrakopoulos, Christos, et al. "Network-based integration of multi-omics data for prioritizing cancer genes." *Bioinformatics* 34.14 (2018): 2441-2448.

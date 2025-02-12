# Transcriptionally-active Riboswitch Tracer Leveraging Edge deTection (TaRTLET)

Finds transcriptional riboswitch activity in RNA-seq data.

## Installation

TaRTLET is currently only compatible on Unix devices. We provide a conda environment to install all the necessary dependencies and packages required for running TaRTLET. From the repository root:

```bash
conda env create --file env.yaml
```

Note: The environment is fairly dense, and we have noticed that the default solver used by conda is unable to practically create the environment. We highly suggest using `conda-libmamba-solver` as the solver instead. To do this:

```bash
conda install -n base conda-libmamba-solver
conda config --set solver libmamba
```

More details may be found on the [Anaconda blog](https://www.anaconda.com/blog/a-faster-conda-for-a-growing-community).

## Using the tool

A walkthrough of the tool's operation using data generated for the manuscript is available in the [tool validation](https://github.com/lepton-7/tartlet-pub) repository.

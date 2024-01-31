# bvbrc_docking

A workflow that runs `fred` or `diffdock` to screen small molecules for receptor proteins. The program can be run as a local script or a `colmena` workflow.

## Setup

We will set up the run environment with Anaconda/Miniconda. 

### Installation

The base environment can be setup with the `env.yml` file in the repo. 

```bash
git clone https://github.com/hengma1001/bvbrc_docking.git
cd bvbrc_docking
conda env create -n bvbrc_docking -f env.yml
pip install .
```

The required `torch` modules are installed using `pip`. It's recommand the user download the correct CUDA version for `pytorch`. (Note: test the pyg with conda install to reduce compilation time. )

```bash
pip install torch==2.1.0 torchvision==0.16.0 torchaudio==2.1.0 
pip install torch-scatter==2.1.2 torch-sparse==0.6.18 torch-cluster==1.6.3 torch-geometric==2.4.0
```

Next we need to setup `DiffDock`, along with `esm`. The Diffdock uses the `esm` package to obtain protein embeddings and structures. 

```bash
git clone https://github.com/gcorso/DiffDock.git
cd DiffDock
git clone https://github.com/facebookresearch/esm
cd esm
git checkout ca8a710
```



## Examples

```bash
python -m bvbrc_docking.run_local -c runs/fred_local.yml
```
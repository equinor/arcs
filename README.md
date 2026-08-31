# ARCS - Automated Reactions for CO<sub>2</sub> Storage
<p align="center">
 <img src="./assets/ARCS_Logo.png" width="300" height="300">
</p>

## Installation

To install ARCS for standard usage:

```
pip install git+https://github.com/equinor/arcs@eqmain
```

ARCS uses the [Astral uv](https://docs.astral.sh/uv/) package manager. Follow
their official instructions on setting up, then do `uv sync`

### Fetching Model Files

Model files are stored using [Git LFS](https://git-lfs.com/), which needs to be
installed when cloning this repository. Install it, then do `git lfs pull` to
fetch the data.

## History / Credits

This is a fork of the original [ARCS](https://github.com/badw/arcs) developed by Benjamin A. D. Williamson. A legacy version of ARCS, based on Bens original version, is deployed [here](https://server-arcs-legacy-dash.radix.equinor.com/) as a demo.

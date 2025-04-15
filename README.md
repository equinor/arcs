# ARCS - Automated Reactions for CO<sub>2</sub> Storage
<p align="center">
 <img src="./assets/ARCS_Logo.png" width="300" height="300">
</p>

## Installation

To install ARCS using the standard `pip` Python package manager, follow these steps after cloning this repository:

```
pip install .
```

ARCS uses the [Poetry](https://python-poetry.org) package and dependency manager. For development and testing, install `poetry` (using `brew install poetry`, `pipx install poetry` or using your preferred package manager), then:

```
poetry sync
```

Refer to the [Poetry documentation](https://python-poetry.org/docs/) for more information.

### Fetching Model Files

To fetch the necessary large model files (pre-computed equations and Gibbs free energy), run:

```
git lfs pull
```

This command assumes you have `git lfs` installed. On macOS you can install it using Homebrew:
```
brew install git-lfs
git lfs install
```
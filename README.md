ehreact
==============================
[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/hesther/ehreact/workflows/CI/badge.svg)](https://github.com/hesther/ehreact/actions?query=workflow%3ACI)


A Python package for extracting and scoring reaction templates based on extended Hasse diagrams. A link to the corresponding manuscript will be available, soon.

**Documentation:** Documentation and a tutorial of EHreact is available at https://hesther.github.io/ehreact/.

## Table of Contents

- [Requirements](#requirements)
- [Installation](#installation)
* [Installation of ehreact](#installation-of-ehreact)
* [Installation of optional packages](#installation-of-optional-packages)
- [Creating a diagram](#creating-a-diagram)
- [Scoring a diagram](#scoring-a-diagram)
- [Contact](#contact)

## Requirements

Ehreact uses Python 3.6+ and a number of python package. We assume you have Conda on your system. If this is not the case yet, install [Miniconda](<https://conda.io/miniconda.html>). However, Conda is optional and you may choose to install packages otherwise.

## Installation

### Installation of ehreact

To install ehreact, make a new conda environment via

1. `conda env create -f environment.yml`
2. `conda activate ehreact`
3. `pip install -e .`

(Note: If you have an old version of conda (< conda 4.4), the second command is `source activate ehreact`). The environment
simply contains python, numpy, rdkit, and typed-argument-parser. If you do not like to use conda, you can instead install these
packages yourself, and then run `pip install -e .`.

To check your installation, try to import ehreact in a Python session::

`import ehreact`

### Installation of optional packages

If you want to automatically atom-map reactions, download ReactionDecoder from [Github](https://github.com/asad/ReactionDecoder/releases) (tested with rdt-2.4.1-jar-with-dependencies.jar), put the jar file into the folder `ReactionDecoder` (within the ehreact package folder) and rename the file to `ReactionDecoder.jar` This is optional, and you can run EHreact without ReactionDecoder.

Plotting Hasse diagrams relies on dot (install for example graphviz) and rsvg-convert. If you do not have these packages on your system already, you can install both via homebrew on Macos::

1. `brew install graphviz`
2. `brew install librsvg`
  
or via apt-get on Linux (Ubuntu, Debian)::

1. `sudo apt-get install graphviz`
2. `sudo apt-get install librsvg2-bin`

There are other options to install, too. Installing dot and rsvg-convert is optional, and you can run EHreact without both, but will not be able to plot diagrams to an image.

## Creating a diagram

Use `train.py` to create a Hasse diagram and specify `--data_path <path>` (the path to the known reactions or substrates), `--train_mode <str>` (`transition_state` for reaction mode if the input is reactions or `single_reactant` for single reactant mode if the input is single substrates), `--save_path <path>` (file to save the created diagram to, optional), and `--save_plot <path>` (file to save an image of the diagram, optional). For example, to create a diagram for the reactions in `ehreact/data/reaction_training.dat`, run

```
python train.py --data_path ehreact/data/reaction_training.dat --save_path test_transition_state.pkl --train_mode transition_state --save_plot test.png
```

which creates an image of the diagram saved under `test.png` and a pickle file of the diagram under `test_transition_state.pkl`. Please refer to the [tutorial](https://hesther.github.io/ehreact/) for more examples and options.

## Scoring on a diagram

Use `predict.py` to score a set of query substrates or reactions. Specify `--test_path <path>` (the path to the query reactions or substrates), `--load_path <path>` (the path to a saved Hasse diagram pickle file), `--predict_mode <str>` (`transition_state` to score reactions, `multi_reactant` to propose products based on reactants and score all reactions, `single_reactant` to propose products and co-substrates based on a reactant and score all reactions). For example, to score the reactions in `ehreact/data/reaction_test.dat`, run

```
python predict.py --test_path ehreact/data/reaction_test.dat --load_path test_transition_state.pkl --predict_mode transition_state
```

Please refer to the [tutorial](https://hesther.github.io/ehreact/) for more examples and options.

## Contact
Feel free to post questions, feedback, bugs or suggestions on github, or email to eheid@mit.edu.

### Copyright

Copyright (c) 2021, Esther Heid


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.5.

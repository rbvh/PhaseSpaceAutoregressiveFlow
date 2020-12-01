# Phase Space Sampling and Inference from Weighted Events with Autoregressive Flows

This repository contains all code used in the paper 'Phase Space Sampling and Inference from Weighted Events with Autoregressive Flows' by Bob Stienen and Rob Verheyen. This code showcases the use of autoregressive flows for the learning of event distributions based on weighted data, including data with negative weights. The interested reader is referred to the paper `2011.13445`.

## What does this repository contain?
This repository contains three distinct things.

1. Slightly modified version of the `nflows` package, which can be found in the folder `nflows-0.13`. When using the code in this repository, please install this version of `nflows` over the official version.
2. A C++ event generator with a custom implementation of VEGAS and an interface to the Madgraph matrix element for ee > ttbar, used to generate the data on which the autoregressive flows can be trained. This generator can be found in the folder `ee_to_ttbar/ME_VEGAS`. That folder also contains the code to calculate the true likelihood of provided data points, as well as code to invert the unit box transformation described in the paper.
3. A whole series of Jupyter Notebooks, containing all code necessary to reproduce the described experiments.

The remainder of this README details how to reproduce the experiments.

## Reproducing the experiments in the paper
### ee -> ttbar
#### 1. Event generation

1. Go to `ee_to_ttbar/ME_VEGAS`
2. Use `make` to compile the code in this folder.
3. Run `generate_vegas` with command-line arguments `-nev n` to indicate the number of events and one of the flags `-weighted`, `-unweighted` or `-flat`
4. For plotting purposes, the generated data can be converted to 4-vectors by running `convert_to_ps.cpp`.

#### 2. Training an autoregressive flow on unweighted events
1. Create unweighted events with instruction set (1) above.
2. Open the notebook `ee_to_ttbar/train_flow_unweighted.ipynb`.
3. Set the hyperparameters in the cell under the HYPERPARAMETERS heading to the settings you wish to use. The default settings are the ones as used in the paper.
4. Direct the notebook to the data created in step 1 by accordingly changing the path in the cell immediately below the LOADING THE TRAINING DATA heading.
5. Run the notebook.

This notebook will create three model files (.pt): `best_ess.pt`, `best_validation.pt` and `final.pt`. We recommend using `best_validation.pt` for any further steps you might take.

#### 3. Training an autoregressive flow on weighted events
1. Create weighted events with instruction set (3) above.
2. Create unweighted events with instruction set (2) above.
3. Open the notebook `ee_to_ttbar/train_flow_unweighted.ipynb`.
4. Define which reference weight you want to use. You can do this in the cell immediately below the heading REWEIGHTING STRATEGY. The setting must be `'min'`, `'mean'`, or `'max'`.
5. Set the hyperparameters in the cell under the HYPERPARAMETERS heading to the settings you wish to use. The default settings are the ones as used in the paper.
6. Direct the notebook to the data and weights created in step 1 by accordingly changing the path in the cell immediately below the LOAD THE DATA AND REWEIGHT heading.
7. Do the same for the test data, by setting the path in the cell LOAD THE TEST DATA to the path of the unweighted data generated in set 2 above.
8. Run the notebook.

This notebook will create three model files (.pt): `flow_model_weighted_{}_best_ess.pt`, `flow_model_weighted_{}_best_validation.pt` and `flow_model_weighted_{}_final.pt` (where {} is replaced by `min`, `mean`, or `max`, depending on your configuration). We recommend using `flow_model_weighted_{}_best_validation.pt` for any further steps you might take.

#### 4. Sampling an autoregressive flow trained on unweighted or weighted events
1. Follow instruction set (4) or (5) above.
2. Open the notebook `ee_to_ttbar/sample_flow.ipynb`
3. In the seventh cell, change the `model_name` variable to the name of the model file (without `.pt`) as it can be found in the `ee_to_ttbar` folder.
4. Run the notebook.

Following these steps will create two files in the `ee_to_ttbar/data` folder, one for the data and one for the corresponding likelihoods. The names of these files depend on the model name defined in step 3.

##### 5. Evaluating unweighting efficiencies
1. Unweighting efficiencies for VEGAS events can be evaluated by running `compute_metrics_from_weights.cpp` and supplying a weights file.
2. Unweighting efficiencies for Flow event samples can be evaluated by running `compute_metrics_from_likelihoods.cpp` and supplying a samples and weights file.

### pp -> ttbar
#### 1. Event generation
1. A Madgraph script and a MadSpin card are included in the `pp_to_ttbar/generation_scripts` folder. 
These may be used by running ./mg5_aMC mg5_script, exiting event generation, replacing the madspin_card.dat in the newly created `ttbar` folder, and running the `ttbar/bin/generate_events` script. 
2. A pythia program and command file are included to shower decayed events, cluster them with `fastjet` and output the samples and (negative) weights directly.

#### 2. Training/sampling an autoregressive flow on negatively weighted events
1. Training of and sampling from the autoregressive flow proceeds similar to the previous experiments
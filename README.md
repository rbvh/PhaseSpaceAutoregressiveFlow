# Phase Space Sampling and Inference from Weighted Events with Autoregressive Flows

This repository contains all code used in the paper 'Phase Space Sampling and Inference from Weighted Events with Autoregressive Flows' by Bob Stienen and Rob Verheyen. This code showcases the use of autoregressive flows for the learning of event distributions based on weighted data, including data with negative weights. The interested reader is referred to the paper [link soon available].

If you used this code in your own scholarly work, don't forget to cite the paper:

    Bob Stienen and Rob Verheyen
    Phase Space Sampling and Inference from Weighted Events with Autoregressive Flows
    [more information soon]

## What does this repository contain?
This repository contains three distinct things.

1. An altered version of the `nflows` package, which can be found in the folder `nflows-0.13`. This version ... . When using the code in this repository, please install this version of `nflows` over the official version.
2. A C++ event generator created by Madgraph, used to generate the data on which the autoregressive flows can be trained. This generator can be found in the folder `ee_to_ttbar/ME_VEGAS`. That folder also contains the code to calculate the true likelihood of provided data points, as well as code to invert the unit box transformation described in the paper. Information on how to use these codes can be found in this folder as well.
3. A whole series of Jupyter Notebooks, containing all code necessary to reproduce the described experiments.

The remainder of this README details how to reproduce the experiments.

## Reproducing the experiments in the paper
### 1. Creating flat events in the space of the ee -> ttbar distribution

1. Go to `ee_to_ttbar/ME_VEGAS`
2. Open the file `generate_flat.cpp` and change the argument of the function in line 79 to the number of data points that you want to generate.
3. Use `make` to compile the code in this folder.
4. Run `generate_flat`.

### 2. Creating unweighted events from ee -> ttbar distribution

1. Go to `ee_to_ttbar/ME_VEGAS`
2. Open the file `generate_vegas.cpp` and at the bottom, uncomment line 490 if it is commented, and comment line 493 if it isn't.
3. On line 487, set the number of data points that you want to generate.
4. Use `make` to compile the code in this folder.
5. Run `generate_vegas`.

### 3. Creating weighted events from ee -> ttbar distribution

1. Go to `ee_to_ttbar/ME_VEGAS`
2. Open the file `generate_vegas.cpp` and at the bottom, uncomment line 493 if it is commented, and comment line 490 if it isn't.
3. On line 487, set the number of data points that you want to generate.
4. Use `make` to compile the code in this folder.
5. Run `generate_vegas`.

### 4. Training an autoregressive flow on unweighted events
1. Create unweighted events with instruction set (2) above.
2. Open the notebook `ee_to_ttbar/train_flow_unweighted.ipynb`.
3. Set the hyperparameters in the cell under the HYPERPARAMETERS heading to the settings you wish to use. The default settings are the ones as used in the paper.
4. Direct the notebook to the data created in step 1 by accordingly changing the path in the cell immediately below the LOADING THE TRAINING DATA heading.
5. Run the notebook.

This notebook will create three model files (.pt): `best_ess.pt`, `best_validation.pt` and `final.pt`. We recommend using `best_validation.pt` for any further steps you might take.

### 5. Training an autoregressive flow on weighted events
1. Create weighted events with instruction set (3) above.
2. Create unweighted events with instruction set (2) above.
3. Open the notebook `ee_to_ttbar/train_flow_unweighted.ipynb`.
4. Define which reference weight you want to use. You can do this in the cell immediately below the heading REWEIGHTING STRATEGY. The setting must be `'min'`, `'mean'`, or `'max'`.
5. Set the hyperparameters in the cell under the HYPERPARAMETERS heading to the settings you wish to use. The default settings are the ones as used in the paper.
6. Direct the notebook to the data and weights created in step 1 by accordingly changing the path in the cell immediately below the LOAD THE DATA AND REWEIGHT heading.
7. Do the same for the test data, by setting the path in the cell LOAD THE TEST DATA to the path of the unweighted data generated in set 2 above.
8. Run the notebook.

This notebook will create three model files (.pt): `flow_model_weighted_{}_best_ess.pt`, `flow_model_weighted_{}_best_validation.pt` and `flow_model_weighted_{}_final.pt` (where {} is replaced by `min`, `mean`, or `max`, depending on your configuration). We recommend using `flow_model_weighted_{}_best_validation.pt` for any further steps you might take.

### 6. Sampling an autoregressive flow trained on unweighted or weighted events
1. Follow instruction set (4) or (5) above.
2. Open the notebook `ee_to_ttbar/sample_flow.ipynb`
3. In the seventh cell, change the `model_name` variable to the name of the model file (without `.pt`) as it can be found in the `ee_to_ttbar` folder.
4. Run the notebook.

Following these steps will create two files in the `ee_to_ttbar/data` folder, one for the data and one for the corresponding likelihoods. The names of these files depend on the model name defined in step 3.
# DLRS_CVR_BAT: Deep-learning-enabled-brain-hemodynamic-mapping-using-resting-state-fMRI

Simple overview of use/purpose.

## Description

An in-depth paragraph about your project and overview of use.

## Getting Started

### Installing

* Create a virtual environment with Python 3.9
```
conda create -n DLRS_CVR_BAT_env python=3.9
conda activate DLRS_CVR_BAT_env
```
* Install packages from requirements.txt
```
pip install -r ./requirements.txt
```
### Executing program
* Preprocess the raw BOLD data
```
matlab ./preprocessing/rs_running.m
```
* Compute the resting-state CVR and BAT maps based on pretrained model
```
python ./src/DLRS_CVR_BAT_inference.py
```


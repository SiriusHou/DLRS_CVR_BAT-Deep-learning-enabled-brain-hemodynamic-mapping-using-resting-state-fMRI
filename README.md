# DLRS_CVR_BAT: Deep-learning-enabled-brain-hemodynamic-mapping-using-resting-state-fMRI

This repository contains the code and pre-trained models for our paper [Deep-learning-enabled brain hemodynamic mapping using resting-state fMRI](https://www.nature.com/articles/s41746-023-00859-y)

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

* Please download SPM12 into ./preprocessing folder. 
* Please download SUIT (https://github.com/jdiedrichsen/suit/releases/tag/3.5) into ./processing/SPM12/toolbox.
### Creating data folder
This pipeline requires each data folder includes both raw BOLD and MPRAGE images. For the code to run automatically, 
it's also essential to generate a parameter_RS.txt and slice_order_RS.txt file for each data folder, as demonstrated 
in the ./template/subject1 folder.

### Downloading pretrained weights
You can obtain the pre-trained weights by submitting a request to Hanzhang Lu. After receiving the weights, they should 
be placed in the ./model directory.

### Executing program
* Preprocess the raw BOLD data using matlab.
```
Open matlab 
Run ./preprocessing/rs_running.m
```
* Compute the resting-state CVR and BAT maps based on pretrained model
```
python ./src/DLRS_CVR_BAT_inference.py
```

### Citation
Please cite our paper if you use DLRS CVR/BAT pipeline in your work:

Hou, X., Guo, P., Wang, P. et al. Deep-learning-enabled brain hemodynamic mapping using resting-state fMRI. npj Digit. Med. 6, 116 (2023). https://doi.org/10.1038/s41746-023-00859-y
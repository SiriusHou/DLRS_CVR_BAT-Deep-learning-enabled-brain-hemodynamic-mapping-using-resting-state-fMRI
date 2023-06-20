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

* Please copy and paste SPM12 into ./preprocessing folder

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


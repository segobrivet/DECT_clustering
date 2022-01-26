# Tumor analysis on spectral images

The aim of this project is to perform an unsupervised clustering of Dual-Energy Computed Tomography (DECT) scans in order to observe separation of tumor tissues or organ anatomical regions as dedicated cluster regions.  

### Environment
This project is developed in Matlab.  
A remote GPU machine is recommanded for heavy computations.

### Data
DECT are 4D data: a 3D body volume over a range of X-ray energy levels.  
*Data are not sharable, so data folders in this repo are empty and need to be filled on local machines.*  


### 1. Image preparation
In `dataset_builder` folder, a large section around tumors is cut from DECT, as well as its ground truth segmentation, and the two are saved as .mat files in `data_tumor` folder.  

### 2. Clustering
In `clustering` folder, `main_clustering.m` script:
- loads a an image file and its associated ground truth segmentation, 
- computes an image clustering,
- at least one cluster region covers the tumor and a Dice score is computed between the cluster region(s) of interest and the ground truth segment.

The unsupervised clustering algorithm is based on mixtures-of-experts models, we integrate spatial context in mixture weights, and construct mixture component densities upon the spectral energy decay curves as functional observations. A dedicated expectation-maximization (EM) algorithm is written to estimate the joint maximum likelihood of the model parameters.

The clustering algorithm is an adaptation from Faicel Chamroukhi repository:  
Curve clustering with the MixFRHLP model and the EM (or a CEM-like) algorithm  
https://github.com/fchamroukhi/mixRHLP_m  


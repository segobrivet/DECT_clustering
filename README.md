# Tumor detection on spectral images

The aim of this project is to segment tumors (or respectively organs) on Dual-Energy Computed Tomography (DECT) scans.  

### Environment
This project is developed in Matlab and uses a remote GPU machine for heavy computations.

### Data
DECT are 4D data: a 3D body volume over a range of X-ray energy levels.  
*Data are owned by AIPHL lab and are not sharable, so data folders in this repo are empty and need to be filled on local machines.*  

### 1. Image preparation
In `dataset_builder` folder, ROIs around tumors (resp. around organs) are cut from DECT, as well as their ground truth segmentations, and the two are saved as .mat files in `data_tumor` folder (resp. in `data_organs`).  

### 2. Clustering
In `segmentation` folder, `main_clustering.m` script:
- loads a ROI and its associated ground truth segmentation, 
- computes an image clustering,
- at least one cluster region covers the tumor (resp. organ) and a Dice score is computed between the cluster region(s) of interest and the ground truth segment.

The unsupervised clustering algorithm is based on mixtures-of-experts models with spatial constraints calculated upon functional observations as those of the energy levels. A dedicated expectation-maximization (EM) algorithm is written to estimate the joint maximum likelihood of the model parameters.

The clustering algorithm is an adaptation from Faicel Chamroukhi repository:  
Curve clustering with the MixFRHLP model and the EM (or a CEM-like) algorithm  
https://github.com/fchamroukhi/mixRHLP_m  

*Currently in fine-tuning stage*


### 3. Classification
In `segmentation` folder, `main_classif.m` script:
- loads cluster regions and ground truths,
- trains a model to classify voxels belonging to a tumor or not (resp. organ),
- determine the tumor region.

*Currently under development*.



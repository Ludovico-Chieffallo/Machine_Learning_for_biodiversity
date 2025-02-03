# Machine Learning for Biodiversity: UAV Images to Infer Bee Abundance

## Overview
This repository contains the R code used to analyze UAV (Unoccupied Aerial Vehicle) imagery and estimate flower abundance using machine learning models (**GBM, SVM, RF, and NNET**). The objective of this study is to evaluate the relationship between flower abundance and bee abundance/diversity using remote sensing techniques and artificial intelligence.

## Dataset
The dataset includes:
- **UAV Images:** High-resolution **JPG** and **GeoTIFF** orthomosaics.
- **Vector Data:** Shapefiles (`.shp`) containing buffered sampling areas.
- **Tabular Data:** **CSV** files with annotated field data on flower and bee abundance.
- **Metadata:** **ISO 19115-compliant metadata** in **XML** format.

## Dependencies
This analysis is performed using **R** and requires the following libraries:
```r
library(terra)
library(pROC)
library(caret)
library(sf)
library(purrr)
library(tools)
library(tidyverse)
library(rsample)
library(furrr)
library(pbapply)
library(ranger)
library(tidymodels)
library(e1071)
```

 ## Data Processing Steps

### Load and Preprocess Data
- Import **UAV orthomosaic** images.
- Load **buffer shapefiles** and overlay them with UAV imagery.
- Mask and crop the images based on the buffer regions.

### Feature Extraction
- Extract **spectral information** from UAV imagery at sampled locations.
- Label extracted data as `flower` or `grass` for classification.

### Dataset Preparation
- Merge extracted features and labels.
- Balance the dataset for training and testing.
- Split the data into **70% training** and **30% testing** sets.

### Model Training
- Train **machine learning models** (**GBM, SVM, RF, NNET**) using cross-validation.
- Optimize hyperparameters using a **grid search**.
- Evaluate models using **accuracy, precision, recall, F1-score, and AUC-ROC**.

### Prediction and Validation
- Apply trained models to **UAV imagery** to classify flower abundance.
- Compute **statistical correlations** between UAV-derived metrics and field observations.
- Generate **confusion matrices** and performance metrics.

### Visualization and Interpretation
- Plot **evaluation metrics**.
- Compare **model predictions** with ground-truth field data.
- Assess the relationship between **flower abundance and bee species richness, abundance, and Shannon diversity index**.

## Results
The results demonstrate the feasibility of using **UAV RGB imagery and machine learning** to predict flower abundance, which is highly correlated with **bee abundance and diversity**. The **best-performing model was GBM**, achieving the highest accuracy and AUC-ROC scores.


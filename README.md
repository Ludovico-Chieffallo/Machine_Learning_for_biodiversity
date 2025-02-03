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

###TO DO###
####Riorganizzare data preparation, prima creare unico dataset con fiori ed erba
####quindi EDA e solo dopo Rifare splitting con funzioni tidymodels e/o ID -> polygon_ID e creare ID vero (di riga)

#Library----
library(tidyverse)
library(randomForest)
library(sf)
library(terra)
library(caret)
library(pROC)
library(rsample)
library(GGally)


#Load data orthomosaic----
EY1_ortho<- rast("data/Ortho EY1.tif")
EY1_ortho

#Crop orthomosaic on transfect buffer
area_buffer<-st_read("data/Buffer EY1.shp")

EY1_ortho<-crop(EY1_ortho, ext(area_buffer))
EY1_ortho<-mask(EY1_ortho, area_buffer)
#EY1_ortho<-EY1_ortho[[1:3]] # uncomment when dealing with > 3 bands

#Load Flowers and Grass----
poligoni_fiori <- st_read("data/fiori EY1.shp")
poligoni_erba <- st_read("data/erba EY1.shp")

nrow(poligoni_fiori)
nrow(poligoni_erba)

#Remove empty geometries----
poligoni_fiori <- poligoni_fiori[!st_is_empty(poligoni_fiori), ]
nrow(poligoni_fiori)

poligoni_erba<- poligoni_erba[!st_is_empty(poligoni_erba), ]
nrow(poligoni_erba)

#Extraction for pixel values---- 
df_fiori <- extract(EY1_ortho, poligoni_fiori)
df_erba <- extract(EY1_ortho, poligoni_erba)

colnames(df_fiori) <- c("ID" ,paste0("banda_", 1:(ncol(df_fiori)-1)))
colnames(df_erba) <- c("ID" ,paste0("banda_", 1:(ncol(df_erba)-1)))

#Add label (fiori e erba)----
df_fiori$label <- "fiore"
df_erba$label <- "erba"

df_tot <- rbind(df_fiori, df_erba) %>% 
  drop_na()

#EDA----
summary(df_tot)
str(df_tot)
ggpairs(df_tot[, 2:5] %>% 
          mutate(banda_4 = sqrt(abs(banda_3 - banda_2))) %>% 
          mutate(banda_5=banda_3- banda_2) %>% 
          mutate(banda_6=(banda_3- banda_2)- banda_1) ,
        mapping = aes(color = label, alpha = 0.7))




#Splitting into train and test----
set.seed(123)
data_split <- initial_split(df_tot, prop = 0.9, strata = "label")
training_data <- training(data_split)
test_data <- testing(data_split)

####Riorganizzare data preparation, prima creare unico dataset con fiori ed erba
####quindi EDA e solo dopo Rifare splitting con funzioni tidymodels e/o ID -> polygon_ID e creare ID vero (di riga)



#Creation of train dataset----
set.seed(123)
df_fiori_train <- df_fiori %>% sample_n(., nrow(.) * 0.9)
df_erba_train <- df_erba %>% sample_n(., nrow(.) * 0.9)

df_train <- rbind(df_fiori_train, df_erba_train)
df_train$label <- as.factor(df_train$label)

#Creation of test dataset----
df_fiori_test <- df_fiori %>% anti_join(df_fiori_train)
df_erba_test <- df_erba %>% anti_join(df_erba_train)

df_test <- rbind(df_fiori_test, df_erba_test)
df_test$label <- as.factor(df_test$label)


# Addestramento del modello sul subset
RF<- randomForest(label ~ ., data = df_train)

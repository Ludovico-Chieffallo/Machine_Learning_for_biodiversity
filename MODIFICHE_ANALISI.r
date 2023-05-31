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
library(ggplot2)
library(patchwork)
library(raster)

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

colnames(df_fiori) <- c("ID_poly" ,paste0("banda_", 1:(ncol(df_fiori)-1)))
colnames(df_erba) <- c("ID_poly" ,paste0("banda_", 1:(ncol(df_erba)-1)))

#Add label (fiori e erba)----
df_fiori$label <- "fiore"
df_erba$label <- "erba"

df_tot <- rbind(df_fiori, df_erba) %>% 
  rowid_to_column(.) %>%
  rename(ID_pixel = rowid)

#Outliers and NA----

# Crea un dataframe 'lungo' 
df_long <- df_tot %>%
  tidyr::pivot_longer(cols =c("banda_1", "banda_2", "banda_3"),
                      names_to = "banda",
                      values_to = "valore")

#vediamo i valori NA
df_long %>% 
  filter(is.na(valore)) %>% 
  view()
#abbiamo solo NA per i fiori e solo per 2 bande(1,2), quindi conviene
#estrarre i valori medi della banda 1 e 2 e assegnarli agli NA

df_fiori_NA<-df_fiori %>% 
  na.omit(df_fiori)

# Calcolare la mediana per banda_1 e banda_2
mediana_banda_1 <- median(df_fiori_NA$banda_1)
mediana_banda_2 <- median(df_fiori_NA$banda_2)

print(paste("Mediana banda 1: ", mediana_banda_1))
print(paste("Mediana banda 2: ", mediana_banda_2))

df_tot <- df_tot %>%
  mutate(banda_1 = replace_na(banda_1, mediana_banda_1),
         banda_2 = replace_na(banda_2, mediana_banda_2))

df_long <- df_tot %>%
  tidyr::pivot_longer(cols =c("banda_1", "banda_2", "banda_3"),
                      names_to = "banda",
                      values_to = "valore")

#Density plot

#OVERLAPPED
#1
gg_lab<-ggplot(df_long, aes(x = valore, fill = label)) +
  geom_density(alpha = 0.5) +
  theme_minimal() +
  labs(x = "Valore", y = "Densità", fill = "Banda", 
       title = "Densità dei valori per label")
#2
gg_band<-ggplot(df_long, aes(x = valore, fill = banda)) +
  geom_density(alpha = 0.5) +
  theme_minimal() +
  labs(x = "Valore", y = "Densità", fill = "Banda", 
       title = "Densità dei valori per banda")

gg_lab/gg_band

#NOT OVERLAPPED
ggplot(df_long, aes(x = valore, fill = label)) +
  geom_density(alpha = 0.5) +
  theme_minimal() +
  facet_wrap(~banda, scales= "free")+
  labs(x = "Valore", y = "Densità", fill = "Banda", 
       title = "Densità dei valori per label")

#Histograms
hist_band<-ggplot(df_long, aes(x = valore)) +
  geom_histogram(bins = 30, alpha = 0.5, color = "black", fill = "skyblue") +
  facet_wrap(~banda, scales = "fixed") +
  theme_minimal() +
  labs(x = "Valore", y = "Frequenza", 
       title = "Istogrammi dei valori per banda")

hist_lab<-ggplot(df_long, aes(x = valore)) +
  geom_histogram(bins = 30, alpha = 0.5, color = "black", fill = "skyblue") +
  facet_wrap(~ label, scales = "fixed") +
  theme_minimal() +
  labs(x = "Valore", y = "Frequenza", 
       title = "Istogrammi dei valori per label")

hist_band+hist_lab

hist_lab_free<-ggplot(df_long, aes(x = valore)) +
  geom_histogram(bins = 30, alpha = 0.5, color = "black", fill = "skyblue") +
  facet_wrap(~ label, scales = "free") +
  theme_minimal() +
  labs(x = "Valore", y = "Frequenza", 
       title = "Istogrammi dei valori per label(scala libera)")

hist_lab +hist_lab_free


#EDA----
ggpairs(df_tot[, 3:6] %>% 
          mutate(banda_4 = sqrt(abs(banda_3 - banda_2))) %>% 
          mutate(banda_6=(banda_3- banda_2)- banda_1),
        mapping = aes(color = label, alpha = 0.6))































#DA QUI IL PANICO###################################################
#Splitting into train and test----
#Seleziono le bande 

df_tot$banda_4 <- sqrt(abs(df_tot$banda_3 - df_tot$banda_2))
df_tot$banda_5 <- (df_tot$banda_3 - df_tot$banda_2) - df_tot$banda_1

set.seed(123)

# Crea un nuovo dataframe con solo le colonne che vuoi usare
df_selected <- df_tot %>%
  select(banda_1,banda_2 , banda_3, label)

# Dividi il dataframe selezionato in set di addestramento e di test
data_split <- initial_split(df_selected, prop = 0.9, strata = "label")
training_data <- training(data_split)
test_data <- testing(data_split)

# Convertiamo in un fattore
training_data$label <- as.factor(training_data$label)
test_data$label <- as.factor(test_data$label)

#Model----
RF<- randomForest(label ~ ., data = training_data)


predizioni <- predict(RF, test_data)
matrice_confusione <- table(test_data$label, predizioni)
accuracy <- sum(diag(matrice_confusione)) / sum(matrice_confusione)
print(accuracy)

head(df_selected)
min(df_selected$banda_4)
max(df_selected$banda_4)
#-65 valore più alto?


#Rename columns for matches with dataset----
names(EY1_ortho) <- colnames(df_selected)[1:3]

#Application of model----
classificazione_RF <- predict(EY1_ortho, RF)
plot(classificazione_RF)


# Salva la classificazione come un file GeoTIFF
output_filename <- "classificazione_RF.tif"
writeRaster(classificazione_RF, output_filename, overwrite = TRUE)










# Definizione della griglia di parametri per la cross-validation
param_grid <- expand.grid(mtry = 1)

# Configurazione del controllo della cross-validation
control <- trainControl(method = "cv", number = 20)

# Addestramento del modello con cross-validation
model_cv <- train(label ~ ., data = training_data, method = "rf", trControl = control, tuneGrid = param_grid)

# Stampa dei risultati della cross-validation
print(model_cv)

# Utilizzo del modello ottimizzato per la valutazione
predizioni_cv <- predict(model_cv, test_data)

print(length(test_data$label))
print(length(predizioni_cv))


matrice_confusione_cv <- table(test_data$label, predizioni_cv)
accuracy_cv <- sum(diag(matrice_confusione_cv)) / sum(matrice_confusione_cv)
print(accuracy_cv)

###################################

# Calcola le probabilità di previsione
prob_pred <- predict(model_cv, newdata = test_data, type = "prob")

# Calcola la previsione delle classi
class_pred <- predict(model_cv, newdata = test_data)

# Calcola la confusione Matrix
cm <- confusionMatrix(class_pred, test_data$label)

# Ottieni precisione, recall e F1-score
precision <- posPredValue(class_pred, test_data$label)
recall <- sensitivity(class_pred, test_data$label)
F1 <- (2 * precision * recall) / (precision + recall)

# Calcola l'AUC-ROC
roc_obj <- roc(test_data$label, prob_pred[,2]) # assumendo che la seconda colonna corrisponda alla classe "positiva"
auc <- auc(roc_obj)

print(paste("Precision:", precision))
print(paste("Recall:", recall))
print(paste("F1-score:", F1))
print(paste("AUC-ROC:", auc))

















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

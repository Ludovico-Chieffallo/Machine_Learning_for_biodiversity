#SAREBBE UTILE SCRIVERE DEL CODICE PER VERIFICARE GLI NA

#Library----
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

#Functions----
source("R/clean_name.R") #scorciatoia per :
# clean_name <- function(file) {
#   var_name <- tools::file_path_sans_ext(basename(file))
#   var_name <- gsub("[^[:alnum:]_]", "", var_name)
#   if (grepl("^[0-9]", var_name)) {
#     var_name <- paste0("x", var_name)
#   }
#   var_name
# }

#Load data----

pathfolder<-"data"
ortho <- list.files(path=pathfolder,
                    pattern = "Ortho", 
                    full.names = TRUE)%>% 
  set_names(nm = map(., clean_name)) %>% 
  map(~rast(.))



buffers <- list.files(path=pathfolder,
                      pattern = "Buffer.*shp",
                      full.names = TRUE )%>% 
  set_names(nm = map(., clean_name)) %>% 
  map(~st_read(.))


ortho<-map2(ortho, map(buffers, ~ ext(.x)), ~ crop(.x, .y))
ortho<-map2(ortho, buffers, ~ mask(.x, .y))

#Fiori ed Erba


fiori<- list.files(path=pathfolder, pattern = "fiori.*shp", full.names = TRUE)%>% 
  set_names(nm = map(., clean_name)) %>% 
  map(~st_read(.))%>% 
  map(~.[!st_is_empty(.), ])

erba<-list.files(path=pathfolder,pattern = "erba.*shp", full.names = TRUE)%>% 
  set_names(nm = map(., clean_name)) %>% 
  map(~st_read(.))%>% 
  map(~.[!st_is_empty(.), ])

#Ordiniamo i nomi in modo che coincidano
names(ortho) <- sort(names(ortho))
names(fiori) <- sort(names(fiori))
names(erba) <- sort(names(erba))

#Values extraction----
df_fiori <- map2(ortho, fiori, ~terra::extract(.x, .y))
df_erba <- map2(ortho, erba, ~terra::extract(.x, .y))


#Nomi delle colonne di ogni dataframe
df_fiori <- map(df_fiori, ~{
  colnames(.) <- c("ID_poly", paste0("banda_", 1:(ncol(.)-1)))
  .
})

df_erba <- map(df_erba, ~{
  colnames(.) <- c("ID_poly", paste0("banda_", 1:(ncol(.)-1)))
  .
})


#label
df_fiori <- map(df_fiori, ~{
  .$label <- "fiore"
  .
})

df_erba <- map(df_erba, ~{
  .$label <- "erba"
  .
})

#Bind dataframes----

# Unisci tutti i dataframe in una lista
df_fiori_unified <- bind_rows(df_fiori, .id = "origin")
df_erba_unified <- bind_rows(df_erba,.id = "origin")


set.seed(123) # Per riproducibilità
size_of_fiori <- nrow(df_fiori_unified)
# Numero di righe da campionare da ciascun gruppo

# Crea un campione bilanciato
df_erba_unified <- df_erba_unified %>%
  group_by(origin) %>%
  sample_n(size = size_of_fiori, replace = TRUE) %>%
  ungroup()

df_tot <- rbind(df_fiori_unified, df_erba_unified) %>% 
  rowid_to_column(.) %>%
  rename(ID_pixel = rowid)

# Crea un dataframe 'lungo' #SOLO SE VEDIAMO NA
# df_long <- df_tot %>%
#   tidyr::pivot_longer(cols =c("banda_1", "banda_2", "banda_3"),
#                       names_to = "banda",
#                       values_to = "valore")
# df_long %>%
#   filter(is.na(valore)) %>%
#   view()


#Splitting into train and test----
#Seleziono le bande

set.seed(123)

# Crea un nuovo dataframe con solo le colonne che vuoi usare
df_selected <- dplyr::select(df_tot, banda_1, banda_2, banda_3, label)


# Dividi il dataframe selezionato in set di addestramento e di test
data_split <- initial_split(df_selected, prop = 0.7, strata = "label")
training_data <- training(data_split)
test_data <- testing(data_split)

# Convertiamo in un fattore-----
training_data$label <- as.factor(training_data$label)
test_data$label <- as.factor(test_data$label)

# SVM Model ----
train_control_svm <- trainControl(
  method = "cv", 
  number = 10,
  classProbs = TRUE, 
  summaryFunction = twoClassSummary 
)

tuneGrid_svm <- expand.grid(C = 10^seq(-2, 2, by = 0.5))

# Rimuovi eventuali righe con valori NA
training_data <- training_data %>% na.omit()
test_data <- test_data %>% na.omit()

set.seed(123)
SVM <- train(
  label ~ .,
  data = training_data,
  method = "svmLinear",
  trControl = train_control_svm,
  tuneGrid = tuneGrid_svm,
  metric = "ROC"
)

print(SVM)

# SVM Model Test ----
predizioni_svm <- predict(SVM, newdata= test_data)
test_accuracy_svm <- mean(predizioni_svm == test_data$label)
print(paste("SVM Test Accuracy:", test_accuracy_svm))

matrice_confusione_svm <- confusionMatrix(as.factor(predizioni_svm), test_data$label, positive = "fiore")
print(matrice_confusione_svm)

# Altre metriche di valutazione possono essere calcolate come prima


# SVM Model Test ----
predizioni_svm <- predict(SVM, newdata = test_data)
test_accuracy_svm <- mean(predizioni_svm == test_data$label)
print(paste("SVM Test Accuracy:", test_accuracy_svm))

matrice_confusione_svm <- confusionMatrix(as.factor(predizioni_svm), test_data$label, positive = "fiore")
print(matrice_confusione_svm)

accuracy_svm <- matrice_confusione_svm$overall['Accuracy']
recall_svm <- matrice_confusione_svm$byClass["Sensitivity"]
f1_svm <- matrice_confusione_svm$byClass["F1"]
auc_roc_svm <- roc(test_data$label, as.numeric(predizioni_svm))
auc_svm <- auc(auc_roc_svm)
precision_svm <- matrice_confusione_svm$byClass["Pos Pred Value"]

print(paste("SVM Accuracy: ", accuracy_svm))
print(paste("SVM Recall: ", recall_svm))
print(paste("SVM F1-Score: ", f1_svm))
print(paste("SVM AUC-ROC: ", auc_svm))
print(paste("SVM Precision: ", precision_svm))



#Valutiamo l'overfitting----

# Confronto tra Prestazioni su Addestramento e Test per SVM
train_predictions_svm <- predict(SVM, newdata = training_data)
train_conf_matrix_svm <- confusionMatrix(as.factor(train_predictions_svm), training_data$label, positive = "fiore")
print(train_conf_matrix_svm)



# Cross-Validation Esterna per SVM
df_selected <- na.omit(df_selected)

set.seed(123)
cv_results_svm <- train(
  label ~ .,
  data = df_selected, # Usa l'intero dataset
  method = "svmLinear",
  trControl = trainControl(
    method = "cv",
    number = 10,
    classProbs = TRUE,
    summaryFunction = twoClassSummary
  ),
  tuneGrid = tuneGrid_svm,
  metric = "ROC"
)
print(cv_results_svm)





# Analisi della Curva ROC per SVM
test_probabilities_svm <- predict(SVM, test_data, type = "prob")
roc_test_svm <- roc(response = test_data$label, test_probabilities_svm[, "fiore"])

# Calcola l'AUC
auc_svm <- auc(roc_test_svm)

# Aggiungi l'AUC al grafico
plot(roc_test_svm, main = paste("ROC Curve (AUC =", auc_svm, ")"))



# Il calcolo dell'importanza delle variabili (varImp) non è direttamente 
# applicabile ai modelli SVM lineari come lo è per i modelli basati 
# su alberi come Random Forest.

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
setwd("/media/data/r_projects/phd_ludovico/2021")
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
pathfolder<-"ortho server 2021 05.03"
#ortho-----
ortho <- list.files(path=pathfolder,
                    pattern = "Ortho", 
                    full.names = TRUE)%>% 
  set_names(nm = map(., clean_name)) %>% 
  map(~rast(.))



#buffers-----
buffers <- list.files(pattern = "buffer.*shp",
                      full.names = TRUE )%>% 
  set_names(nm = map(., clean_name)) %>% 
  map(~st_read(.))

ortho <- ortho[order(names(ortho))]
buffers <- buffers[order(names(buffers))]

#crop e maschera----
# Sovrapponi gli orthomosaic con i buffer
#ortho_intersect <- map2(ortho, map(buffers, ~ ext(.x)), ~ intersect(.x, .y))
#ATTENZIONE.
ortho<-map2(ortho, map(buffers, ~ ext(.x)), ~ crop(.x, .y))
ortho<-map2(ortho, buffers, ~ mask(.x, .y))
ortho
#Fiori ed Erba-----

# Fiori ed Erba combinati
#Sys.setenv("SHAPE_RESTORE_SHX" = "YES")
fiori_erba <- list.files(pattern = "fiori e erba.*shp", full.names = TRUE) %>% 
  set_names(nm = map(., clean_name)) %>% 
  map(~st_read(.)) %>% 
  map(~.[!st_is_empty(.), ])

fiori <- fiori_erba %>% 
  map(~ . %>% filter(id == 1))  # Usa filter per selezionare gli elementi dove id è 1

erba <- fiori_erba %>% 
  map(~ . %>% filter(id == 2))


#Ordiniamo i nomi in modo che coincidano
names(ortho) <- sort(names(ortho))
names(fiori) <- sort(names(fiori))
names(erba) <- sort(names(erba))

#Values extraction----

# Converti manualmente un singolo elemento per il test
fiori <- map(fiori_erba, ~{
  # Verifica che la colonna 'id' esista
  if ("id" %in% names(.x)) {
    # Filtra solo se ci sono righe che corrispondono al criterio
    filtered <- .x %>% filter(id == 1)
    # Converte in 'sf' solo se il risultato non è vuoto
    if (nrow(filtered) > 0) {
      st_as_sf(filtered)
    } else {
      NULL  # Restituisce NULL per frame dati vuoti
    }
  } else {
    NULL  # Restituisce NULL se non esiste la colonna 'id'
  }
})

erba <- map(fiori_erba, ~{
  # Verifica che la colonna 'id' esista
  if ("id" %in% names(.x)) {
    # Filtra solo se ci sono righe che corrispondono al criterio
    filtered <- .x %>% filter(id == 2)
    # Converte in 'sf' solo se il risultato non è vuoto
    if (nrow(filtered) > 0) {
      st_as_sf(filtered)
    } else {
      NULL  # Restituisce NULL per frame dati vuoti
    }
  } else {
    NULL  # Restituisce NULL se non esiste la colonna 'id'
  }
})

#controllo per capire se ci sono aree NULL
# Crea un vettore logico che indica quali elementi in 'fiori' sono NULL
are_null_fiori <- sapply(fiori, is.null)

# Stampa gli indici degli elementi NULL in 'fiori'
null_indices_fiori <- which(are_null_fiori)
print(null_indices_fiori)

# Ripeti lo stesso processo per 'erba'
are_null_erba <- sapply(erba, is.null)
null_indices_erba <- which(are_null_erba)
print(null_indices_erba)

#estrazione dati----
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


set.seed(123) # Per riproducibilitÃ 
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

# GBM Model ----
train_control_gbm <- trainControl(
  method = "cv",
  number = 10,
  classProbs = TRUE,
  summaryFunction = twoClassSummary
)

tuneGrid_gbm <- expand.grid(
  interaction.depth = c(1, 3, 5), 
  n.trees = c(50, 100, 150),
  shrinkage = c(0.01, 0.1),
  n.minobsinnode = c(10, 20)
)

# Rimuovi eventuali righe con valori NA
training_data <- training_data %>% na.omit()
test_data <- test_data %>% na.omit()

set.seed(123)
GBM <- train(
  label ~ .,
  data = training_data,
  method = "gbm",
  trControl = train_control_gbm,
  tuneGrid = tuneGrid_gbm,
  verbose = FALSE,
  metric = "ROC"
)

print(GBM)

# GBM Model Test ----
predizioni_gbm <- predict(GBM, newdata = test_data)
test_accuracy_gbm <- mean(predizioni_gbm == test_data$label)
print(paste("GBM Test Accuracy:", test_accuracy_gbm))

matrice_confusione_gbm <- confusionMatrix(as.factor(predizioni_gbm), test_data$label, positive = "fiore")
print(matrice_confusione_gbm)






# Altre metriche di valutazione possono essere calcolate come prima


# GBM Model Test ----
predizioni_gbm <- predict(GBM, newdata = test_data)
test_accuracy_gbm <- mean(predizioni_gbm == test_data$label)
print(paste("GBM Test Accuracy:", test_accuracy_gbm))

matrice_confusione_gbm <- confusionMatrix(as.factor(predizioni_gbm), test_data$label, positive = "fiore")
print(matrice_confusione_gbm)

accuracy_gbm <- matrice_confusione_gbm$overall['Accuracy']
recall_gbm <- matrice_confusione_gbm$byClass["Sensitivity"]
f1_gbm <- matrice_confusione_gbm$byClass["F1"]
auc_roc_gbm <- roc(test_data$label, as.numeric(predizioni_gbm))
auc_gbm <- auc(auc_roc_gbm)
precision_gbm <- matrice_confusione_gbm$byClass["Pos Pred Value"]

print(paste("GBM Accuracy: ", accuracy_gbm))
print(paste("GBM Recall: ", recall_gbm))
print(paste("GBM F1-Score: ", f1_gbm))
print(paste("GBM AUC-ROC: ", auc_gbm))
print(paste("GBM Precision: ", precision_gbm))
                             








#ATTENZIONE DA QUI------------------------











#GRAFICI----
matrice <- as.data.frame(matrice_confusione_gbm$table)

ggplot(data = matrice, aes(x = Reference, y = Prediction, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = Freq), vjust = 1.5) +
  scale_fill_gradient(low = "grey", high = "blue") +
  theme_minimal() +
  labs(title = "Matrice di Confusione per il Modello SVM", x = "Etichetta Vera", y = "Predizione")

metriche <- data.frame(
  Metrica = c("Accuracy", "Recall", "F1-Score", "AUC-ROC", "Precision"),
  Valore = c(accuracy_svm, recall_svm, f1_svm, auc_svm, precision_svm)
)

ggplot(metriche, aes(x = Metrica, y = Valore, fill = Metrica)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  theme_minimal() +
  coord_flip() +
  labs(title = "Metriche di Valutazione per il Modello SVM", x = "", y = "Valore")



#tentativo conta dei fiori----



# Assumendo che `ortho` sia la tua lista di raster
# E che il modello GBM sia salvato nella variabile `GBM`

# Funzione per convertire raster in dataframe e fare predizioni
process_and_predict <- function(raster) {
  # Converti il raster in un dataframe, includendo solo le bande usate per l'addestramento
  pixel_data <- terra::as.data.frame(raster, xy = TRUE, na.rm = TRUE)
  
  # Rinomina le colonne per assicurarsi che corrispondano a quelle usate nel modello SVM
  # Assicurati che questo mapping corrisponda ai nomi delle bande nel tuo raster e a come sono stati denominati nel training set
  names(pixel_data) <- c('x', 'y', 'banda_1', 'banda_2', 'banda_3')
  
  # Seleziona solo le colonne di interesse (rimuovi 'x', 'y' se non sono usate per la predizione)
  pixel_data <- dplyr::select(pixel_data, banda_1, banda_2, banda_3)
  
  # Applica il modello GBM per fare predizioni sui pixel
  predizioni <- predict(GBM, newdata = pixel_data)
  
  # Ritorna le predizioni
  return(predizioni)
}


# Applica la funzione a ogni raster in `ortho`
risultati_predizioni <- lapply(ortho, process_and_predict)



# Aggrega i risultati delle predizioni per ciascun raster
conteggi_predizioni <- lapply(risultati_predizioni, table)

# Stampa i conteggi per ciascun raster
lapply(conteggi_predizioni, print)


# Imposta l'area coperta da ciascun pixel (in cm²)
area_per_pixel_cm2 <- 0.25

# Supponendo che `conteggi_predizioni` sia la tua lista con i conteggi dei pixel
conteggi_predizioni_cm2 <- lapply(conteggi_predizioni, function(x) x * area_per_pixel_cm2)

# Stampa i risultati convertiti in cm2
lapply(conteggi_predizioni_cm2, print)

#RISULTATI GBM----
# $OrthoEL1
# 
# erba    fiore 
# 652534.2   1118.0 
# 
# $OrthoEY1
# 
# erba   fiore 
# 1213066   40154 
# 
# $OrthoEY3
# 
# erba      fiore 
# 1506984.25   10627.75 
# 
# $OrthoEY4
# 
# erba   fiore 
# 2143189   17361 
# 
# $OrthoFR1
# 
# erba   fiore 
# 1446297   24255 
# 
# $OrthoFR2
# 
# erba   fiore 
# 2431402    5107 
# 
# $OrthoFR4
# 
# erba   fiore 
# 1881301    1561 
# 
# $OrthoGE3
# 
# erba     fiore 
# 1681179.5    4517.5 
# 
# $OrthoGU1
# 
# erba     fiore 
# 823963.00    907.25 
# 
# $OrthoGU2
# 
# erba      fiore 
# 1591468.75   61816.25 
# 
# $OrthoGU3
# 
# erba     fiore 
# 750146.25  16071.25 
# 
# $OrthoGU4
# 
# erba      fiore 
# 1593585.50   19238.75 
# 
# $OrthoMA2
# 
# erba      fiore 
# 4655415.25    2305.75 
# 
# $OrthoMA3
# 
# erba      fiore 
# 1477478.75   13721.25 
# 
# $OrthoPA2
# 
# erba   fiore 
# 1271680    1915 
# 
# $OrthoRA4
# 
# erba     fiore 
# 931634.25   2230.25 
# 
# $OrthoSG3
# 
# erba      fiore 
# 1742302.00    3914.25 
# 
# $OrthoSI2
# 
# erba      fiore 
# 1045745.75    4474.25 
# 
# $OrthoVA1
# 
# erba     fiore 
# 2244777.0   27419.5 
# 
# $OrthoVA2
# 
# erba     fiore 
# 2347891.8    1557.5 
# 
# $OrthoVA3
# 
# erba      fiore 
# 2075919.00    9300.25 
# 
# $OrthoVA4
# 
# erba   fiore 
# 1444194   56944 
# 
# $OrthoWA1
# 
# erba     fiore 
# 837130.50   4055.75 
# 
# $OrthoWY2
# 
# erba     fiore 
# 1026995.2   16988.5 
# 
# $OrthoWY3
# 
# erba     fiore 
# 1056237.0    5315.5 
# 
# $OrthoWY4
# 
# erba      fiore 
# 1378750.00    4950.25
#----------------------------------------------------------------------------------------

  
# Creazione di un dataframe dai risultati----
risultati_drone_GBM <- data.frame(area = c("EL1", "EY1", "EY3", "EY4", "FR1", "FR2", "FR4", "GE3", "GU1", "GU2", "GU3", "GU4", "MA2", "MA3", "PA2", "RA4", "SG3", "SI2", "VA1", "VA2", "VA3", "VA4", "WA1", "WY2", "WY3", "WY4"),
  erba_cm2 = c(652534.2, 1213066, 1506984.25, 2143189, 1446297, 2431402, 1881301, 1681179.5, 823963.00, 1591468.75, 750146.25, 1593585.50, 4655415.25, 1477478.75, 1271680, 931634.25, 1742302.00, 1045745.75, 2244777.0, 2347891.8, 2075919.00, 1444194, 837130.50, 1026995.2, 1056237.0, 1378750.00),
  fiore_cm2 = c(1118.0, 40154, 10627.75, 17361, 24255, 5107, 1561, 4517.5, 907.25, 61816.25, 16071.25, 19238.75, 2305.75, 13721.25, 1915, 2230.25, 3914.25, 4474.25, 27419.5, 1557.5, 9300.25, 56944, 4055.75, 16988.5, 5315.5, 4950.25))


# Stampa del dataframe per verificare i dati
print(risultati_drone_GBM)

# Salvataggio dei risultati in un file CSV
#write.csv(risultati, "Risultati_Aree_2021.csv", row.names = FALSE)

#Correlazioni----


# Carica i dati di campo dal file CSV
#dati_campo <- read.csv("/media/data/r_projects/phd_ludovico/2021/dati di campo 2021 modi.csv",sep = ";", header = TRUE)

# Testo dei dati, con ogni riga separata da una nuova linea e ogni valore separato da una virgola
testo_dati <- ",area,flower_shannon_area,flower_shannon_number,flower_simpson_area,flower_simpson_number,flower_surface_cm2,flower_species_richness,bee_shannon,bee_simpson,bee_species_richness,bee_abundance,flowers_rf_original,flowers_rf_1cm,flowers_rf_2cm,flowers_rf_5cm,flowers_svm_original,flowers_svm_1cm,flowers_svm_2cm,flowers_svm_5cm,flowers_nnet_original,flowers_nnet_1cm,flowers_nnet_2cm,flowers_nnet_5cm
1,EL1,0.910721783,0.783767869,0.564777891,0.461355529,1221.1,3,0,0,1,1,0.495117792,1.15071909,1.46230725,2.477797359,3.344074039,3.247341441,3.778030047,5.348648297,0.338154148,0.31,0.601650519,3.4803638
2,EY1,1.167532238,1.443673335,0.611684978,0.713039155,46534.97015,10,2.046738533,0.8515625,9,16,9.31595,17.00188008,21.37021327,27.08332663,14.33721,42.94877989,43.09097164,38.03433364,13.62121,8.744049,9.546535,10.56796
3,EY2,1.624727826,1.139825932,0.77250453,0.609326897,44271.80001,10,1.69878299,0.763888889,7,12,6.043757285,13.8396472,39.03119099,60.40862904,15.66734012,15.80294087,25.68941521,46.74889638,10.70358707,12.7710662,36.11487178,62.05626944
4,EY3,1.147127251,1.586762669,0.5758459,0.723272451,9930.05178,18,1.494175138,0.75,5,8,1.054695897,0.969949471,1.992943436,3.311279286,2.179633681,0.694746762,1.484855301,3.534834755,1.200491677,1.167011453,1.450277289,3.284418717
5,EY4,0.797482886,1.142874886,0.416343169,0.59881716,6453.764838,6,1.098612289,0.666666667,3,3,2.045546318,2.763478654,4.477850613,45.76395884,2.553947093,8.855969724,17.26949639,77.14974497,2.564283362,1.486509864,5.879393135,29.40148086
6,FR1,1.239518348,1.498138429,0.618441459,0.705675596,8138.110715,7,1.671595278,0.792899408,6,13,3.18131647,2.693909066,3.497954226,16.98974452,4.372061536,5.777897867,5.918542539,15.63064741,3.607272469,4.654292827,4.90637422,17.57842508
7,FR2,0.371669662,0.792120813,0.213992379,0.5142,1032.298138,3,0,0,1,3,1.464616586,0.881031487,1.508433172,3.918479443,0.155913434,0.939715681,1.351135257,5.010348264,0.792879636,0.861198924,1.207333076,6.26661056
8,FR3,0.898540941,1.51950478,0.383750268,0.688383695,29400.59892,13,1.473502385,0.694444444,6,12,5.090192686,5.585467691,11.88916773,20.19474059,6.43309626,6.475135799,8.111998645,10.63456227,5.673031479,6.997488132,14.34368517,15.1130637
9,FR4,0.865550798,1.469515836,0.375426287,0.696197199,1568.059064,7,0,0,1,1,0.352155653,0.349916937,0.385096,3.705641,3.226668992,2.719372828,1.391501,0.01617377,0.212346618,0.781382662,0.9957928,3.145091
10,GE3,1.490577497,0.849446289,0.668687294,0.374081003,2244.517984,11,1.039720771,0.625,3,4,1.790381496,5.002079696,8.46530748,10.59029645,7.422331915,10.5191843,12.9957105,23.16205002,5.4550674,7.001731342,10.51826427,12.41918427
11,GU1,0.580910212,0.711179383,0.376054407,0.391915461,677.8764714,3,0,1,0,0,0.241001567,0.635352539,1.028475638,2.519618104,0.470090202,0.543266308,0.764578554,1.821812004,0.492036089,0.543262796,0.727591693,1.477880112
12,GU2,0.371118138,0.895634709,0.15797075,0.492473091,28964.91152,7,0.895332666,0.459183673,4,14,6.123899,12.90580551,18.56701017,41.70529254,15.76286,15.39266835,50.87799487,40.04109398,9.79469,12.82674676,15.912395,49.87043229
13,GU3,1.391528665,2.074311315,0.640819679,0.834552986,28830.29289,15,1.33217904,0.72,4,5,4.013780991,3.518573787,4.987524816,17.90462575,20.08581885,19.10337556,19.23183996,24.64321886,7.579886868,15.80095514,16.58697995,27.52184223
14,GU4,0.824001711,1.325401668,0.384170463,0.648263006,17365.10854,8,1.039720771,0.625,3,4,4.69548625,8.458295597,19.79420783,38.56258999,2.777937182,2.86506449,3.734593691,14.23483113,3.89946593,6.62222668,13.64438167,46.06858806
15,MA2,0.245210978,0.353632387,0.109130474,0.156989572,312.3696974,4,0,1,0,0,0.325279498,1.139600674,0.953806382,1.595501022,0.408793812,0.897913102,3.990502071,8.647931605,0.116044395,0.903650756,1.93881974,4.213439999
16,MA3,1.787138225,1.362488099,0.769460725,0.570804383,11820.80575,18,1.549826046,0.775510204,5,7,2.3711649,2.128370378,3.728911456,12.09956252,5.503453276,4.58734525,32.22656366,15.07663677,3.141243076,3.451561107,3.432636587,8.846242101
17,PA2,0.893519801,1.245205469,0.524908432,0.636164934,929.192124,5,0.693147181,0.5,2,2,0.134163629,0.104424896,0.193213604,0.942167111,2.033292452,1.36529202,0.920837898,0.511268081,0.387226542,0.611619859,0.969597348,0.939762539
18,RA4,1.568633821,1.417826259,0.722251379,0.631031674,4194.806831,11,1.277034259,0.693877551,4,7,1.108442839,0.898853267,2.972063924,18.1929284,1.107763586,8.982356877,10.49235678,15.15872375,0.734734249,0.752635125,4.820955071,14.24946007
19,SG3,0.984966583,1.258163194,0.493303183,0.64433876,4293.409875,6,0.410116318,0.244897959,2,7,1.202130873,1.142818402,1.527028176,4.764080195,1.924448984,1.889629474,2.50817527,6.596802481,1.620603631,2.073603319,3.01238284,5.465764843
20,SI2,1.143349525,0.651869017,0.598318203,0.279205355,6310.932303,9,1.213007566,0.65625,4,8,0.843965325,1.194039769,1.757020334,2.423710603,3.91148001,4.267026296,4.190310339,3.807010065,1.041439641,1.280378276,2.301603424,4.88315507
21,VA1,1.160082516,1.502024052,0.539122421,0.739051147,9690.227602,8,1.265856752,0.698224852,4,13,2.700261804,2.583397833,3.416054052,13.65319736,26.47123183,25.23154813,27.6214549,27.78427431,2.44264,2.088489608,3.834521962,8.998207238
22,VA2,0.036420711,0.646090505,0.011828609,0.4536862,127.9960409,2,0,1,0,0,0.243534186,0.195332608,0.381604415,0.946282985,2.862207013,0.01230443,2.912919581,2.567687816,0.453473536,0.085683433,0.490388743,2.084720507
23,VA3,0.807583411,0.891656121,0.438031795,0.527505932,3526.382441,7,0,1,0,0,1.982098053,2.242445978,2.589322323,4.347618872,1.808124155,1.897746922,2.691211138,4.849687483,2.018432103,2.079916449,2.705828433,4.639389152
24,VA4,0.216829195,0.453077658,0.08989032,0.201694519,28269.85315,7,1.886696785,0.84,7,10,6.295722488,6.474359327,11.69792297,25.92038997,11.16190274,11.45597112,14.80269186,26.76458453,10.69186204,14.30271852,17.60654646,32.14275633
25,WA1,1.042347857,1.357016895,0.570220656,0.708710977,8171.797397,6,0.562335145,0.375,2,4,1.515410963,1.680225007,4.025219455,14.60986002,1.513465108,1.618017396,2.1594099,4.553156176,1.378468557,1.584343294,2.904732284,5.07411094
26,WY1,0.886199518,1.345077868,0.536778989,0.707699692,26032.641,8,1.504788284,0.76,5,10,5.15312765,6.205423263,10.22002948,26.09375652,9.012099088,14.49374402,14.0837388,16.61652544,11.00680962,9.513813198,13.72057247,37.36098456
27,WY2,1.294592187,0.702184411,0.638222723,0.301935931,21149.68382,10,1.666875697,0.792243767,6,19,7.477514097,9.273249319,11.11057421,25.61534994,18.52019641,21.32016951,22.38496625,25.63523478,12.70941757,11.99887477,14.8592556,34.25490794
28,WY3,1.01309165,2.006025685,0.416084616,0.812115212,3507.82732,17,2.13833306,0.875,9,15,4,1.04134468,1.942988341,4.882048883,5.667987175,5.609164155,6.177628325,7.21570588,3.412427762,2.378148972,3.610806611,8.997236505
29,WY4,0.560807405,0.92821457,0.305348914,0.538937961,3247.866279,4,0,0,1,1,0.06352426,0.157149494,0.298750849,0.429396689,0.075120647,3.951657437,2.853157647,3.714686827,0.077076809,0.113796652,0.178580484,2.206869972
30,SI1,0.406514612,0.315396287,0.194917813,0.139917695,17.79431281,3,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0
"

# Converti il testo dei dati in un dataframe utilizzando read.csv e textConnection
dati_campo <- read.csv(textConnection(testo_dati), header = TRUE)
# Rimuovere le righe con aree specificate
dati_campo <- dati_campo[!dati_campo$area %in% c("EY2", "FR3", "WY1", "SI1"), ]

#PROVIAMO AD AUMENTARE L'R2----
dati_campo <- dati_campo[!dati_campo$area %in% c( "GU2", "VA4"), ]
risultati_drone_GBM <- risultati_drone_GBM[!risultati_drone_GBM$area %in% c( "GU2", "VA4"), ]



#selezioniamo le aree giuste
dati_campo <- dati_campo[, c("area", "flower_shannon_area", "flower_shannon_number", "flower_simpson_area", 
                             "flower_simpson_number", "flower_surface_cm2", "flower_species_richness", "bee_shannon", 
                             "bee_simpson", "bee_species_richness", "bee_abundance")]

# Mostra le prime righe del dataframe creato
head(dati_campo)

# Rimuovi il prefisso "Ortho" dalla colonna area in risultati_drone se necessario
risultati_drone_GBM$area <- gsub("Ortho", "", risultati_drone_GBM$area)

# Unisci le colonne usando dplyr
dati_congiunti <- left_join(dati_campo, risultati_drone_GBM[, c("area", "fiore_cm2")], by = "area")

# Verifica i risultati
head(dati_congiunti)

#FLOWER AB./FLOWER AB. GBM-----

# Calcola la correlazione tra le due stime dei cm² dei fiori
correlazione <- cor(dati_congiunti$flower_surface_cm2, dati_congiunti$fiore_cm2)

# Stampa il coefficiente di correlazione
print(correlazione)


# Costruisci il modello lineare
modello_lineare <- lm(fiore_cm2 ~ flower_surface_cm2, data = dati_congiunti)

# Riassunto del modello per ottenere R^2 e i p-value
summary_modello <- summary(modello_lineare)

# Estrai il valore R^2
r_quadrato <- summary_modello$r.squared

# Estrai il p-value per la pendenza (flower_surface_cm2)
p_value <- summary_modello$coefficients[2,4]

# Stampa i valori R^2 e p-value
print(paste("R^2:", r_quadrato))
print(paste("p-value:", p_value))

# Formatta il p-value per visualizzare tre cifre decimali
p_value_format <- formatC(p_value, format = "f", digits = 8)

annotazione_testo <- paste("R^2 = ", round(r_quadrato, digits = 3),
                           "\np-value = ", p_value_format)
# Crea il grafico con ggplot2 e aggiungi l'annotazione
ggplot(dati_congiunti, aes(x = flower_surface_cm2, y = fiore_cm2, label = area)) +
  geom_point() +
  geom_smooth(method = "lm", col = "blue") +
  geom_text(check_overlap = TRUE, nudge_y = 0.05, hjust = 0.5, vjust = -0.5) +
  labs(x = "Dati di campo (cm²)", 
       y = "Dati da drone GBM (cm²)", 
       title = "flower abundance(field) vs flower abundance UAV (GBM)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = Inf, y = Inf, label = annotazione_testo, hjust = 1.1, vjust = 1, size = 5, colour = "blue")


#BEE ABUNDANCE/FLOWER AB GBM----

# Calcola la correlazione
correlazione <- cor(dati_congiunti$bee_abundance, dati_congiunti$fiore_cm2)
print(paste("Correlazione:", correlazione))

# Costruisci il modello lineare
modello_lineare <- lm(fiore_cm2 ~ bee_abundance, data = dati_congiunti)

# Riassunto del modello per ottenere R^2 e il p-value
summary_modello <- summary(modello_lineare)

# Estrai il valore R^2
r_quadrato <- summary_modello$r.squared
print(paste("R^2:", r_quadrato))

# Estrai il p-value per la pendenza (bee_abundance)
p_value <- summary_modello$coefficients[2, 4]
print(paste("p-value:", p_value))

# Visualizza la relazione con ggplot2

ggplot(dati_congiunti, aes(x = bee_abundance, y = fiore_cm2)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue") +
  labs(x = "Abbondanza delle api", y = "Area dei fiori (cm²)", title = "Correlazione tra Abbondanza delle Api e Area dei Fiori") +
  annotate("text", x = Inf, y = Inf, label = sprintf("R² = %.3f\np-value = %.3f", r_quadrato, p_value), hjust = 1.1, vjust = 2, size = 5, colour = "blue")


#BEE SPECIES RICHNESS / FLOWER AB. GBM

# Calcolo della correlazione tra bee_species_richness e fiori_cm2
correlazione_richness <- cor(dati_congiunti$bee_species_richness, dati_congiunti$fiore_cm2)

# Costruzione del modello lineare per bee_species_richness contro fiori_cm2
modello_lineare_richness <- lm(fiore_cm2 ~ bee_species_richness, data = dati_congiunti)

# Riassunto del modello per ottenere R^2 e il p-value
summary_modello_richness <- summary(modello_lineare_richness)

# Estrazione del valore R^2
r_quadrato_richness <- summary_modello_richness$r.squared

# Estrazione del p-value per la pendenza (bee_species_richness)
p_value_richness <- coef(summary(modello_lineare_richness))["bee_species_richness", "Pr(>|t|)"]

# Stampa dei risultati
print(paste("Correlazione:", correlazione_richness))
print(paste("R^2:", r_quadrato_richness))
print(paste("p-value:", p_value_richness))

ggplot(dati_congiunti, aes(x = bee_species_richness, y = fiore_cm2)) +
  geom_point() +  # Disegna i punti
  geom_smooth(method = "lm", color = "blue") +  # Aggiunge la linea di regressione
  labs(x = "Ricchezza delle specie di api", y = "Area dei fiori (cm²)", 
       title = "Relazione tra Ricchezza delle Specie di Api e Area dei Fiori") +
  annotate("text", x = max(dati_congiunti$bee_species_richness, na.rm = TRUE), y = min(dati_congiunti$fiore_cm2, na.rm = TRUE), 
           label = paste("R² =", round(r_quadrato_richness, 3), "\np-value =", format(p_value_richness, scientific = FALSE)), 
           hjust = 1, vjust = -1, color = "red", size = 5) +
  theme_minimal()

#SHANNON / FLOW AB GBM----

# Costruzione del modello lineare
modello_lineare_shannon <- lm(fiore_cm2 ~ bee_shannon, data = dati_congiunti)

# Riassunto del modello per ottenere R^2 e il p-value
summary_modello_shannon <- summary(modello_lineare_shannon)
r_quadrato_shannon <- summary_modello_shannon$r.squared
p_value_shannon <- coef(summary(modello_lineare_shannon))["bee_shannon", "Pr(>|t|)"]


# Creazione del grafico con nomi delle aree
ggplot(dati_congiunti, aes(x = bee_shannon, y = fiore_cm2, label = area)) +
  geom_point() +  # Disegna i punti
  geom_smooth(method = "lm", color = "blue") +  # Aggiunge la linea di regressione
  geom_text(check_overlap = TRUE, nudge_y = 5000, hjust = 0.5, vjust = 0, size = 3) +  # Aggiunge etichette per le aree
  labs(x = "Indice di Shannon", y = "Area dei fiori (cm²)", 
       title = "Relazione tra Indice di Shannon e Area dei Fiori") +
  annotate("text", x = max(dati_congiunti$bee_shannon, na.rm = TRUE), y = min(dati_congiunti$fiore_cm2, na.rm = TRUE), 
           label = paste("R² =", round(r_quadrato_shannon, 3), "\np-value =", format(p_value_shannon, scientific = FALSE)), 
           hjust = 1, vjust = -1, color = "red", size = 5) +
  theme_minimal()





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

# Pre-elaborazione delle immagini: Selezione delle bande RGB e Standardizzazione
ortho <- map(ortho, function(x) {
  # Selezione delle prime tre bande per ogni orthomosaico
  x_rgb <- x[[1:3]]
  
  # Calcola il valore minimo e massimo per le bande RGB
  min_val <- min(x_rgb, na.rm = TRUE)
  max_val <- max(x_rgb, na.rm = TRUE)
  
  # Standardizzazione dei valori all'interno dell'intervallo 0-255
  x_standardized <- terra::clamp(x_rgb, lower = min_val, upper = max_val)
  x_standardized <- (x_standardized - min_val) / (max_val - min_val) * 255
  
  return(x_standardized)
})


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
t0<-Sys.time()
ortho<-map2(ortho, map(buffers, ~ ext(.x)), ~ crop(.x, .y))
ortho<-map2(ortho, buffers, ~ mask(.x, .y))
ortho
t1<-Sys.time()
t1-t0
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

#RF Model----
train_control_RF <- trainControl(
  method = "cv",             # Cross-validation
  number = 10,               # Numero di fold per la cross-validation
  savePredictions = "final", # Salva le predizioni per ogni modello
  classProbs = TRUE,         # Calcola le probabilitÃ  di appartenenza alle classi
  summaryFunction = twoClassSummary # Funzione di riepilogo per classificazione binaria
)

max_mtry <- max(2, floor(sqrt(ncol(training_data) - 1)))

tuneGrid_RF <- expand.grid(
  mtry = seq(2, max_mtry, by=1),
  splitrule = c("gini", "extratrees"),
  min.node.size = c(1, 5, 10, 15)
)

print(tuneGrid_RF)


## Train the model using ranger and the tuning grid
# Calcola i pesi delle classi
#classi <- unique(training_data$label)
#class_weights <- c("erba" = 0.1, "fiore" = 0.9)

# Rimuovi eventuali righe con valori NA
training_data <- training_data %>% na.omit()
test_data <- test_data %>% na.omit()

# Addestra il modello utilizzando i pesi delle classi
set.seed(123)
t0<-Sys.time()
RF <- train(
  label ~ ., 
  data = training_data,
  method = "ranger", 
  trControl = train_control_RF, 
  tuneGrid = tuneGrid_RF,
  metric = "ROC"
)
t1<-Sys.time()
t1-t0
#19.81943 hours
print(RF)
print(RF$bestTune)
#se volessi addestrare un modello SVM, potrei usare method = "svmLinear"

#RF model test----
predizioni_RF <- predict(RF, newdata= test_data)
test_accuracy_RF <- mean(predizioni_RF== test_data$label)
print(paste("Test Accuracy:", test_accuracy_RF))

matrice_confusione_RF <- confusionMatrix(as.factor(predizioni_RF), test_data$label, positive = "fiore")
print(matrice_confusione)

accuracy_RF <- matrice_confusione_RF$overall['Accuracy']
recall_RF <- matrice_confusione_RF$byClass["Sensitivity"]
f1_RF <- matrice_confusione_RF$byClass["F1"]
auc_roc_RF<- roc(test_data$label, as.numeric(predizioni_RF))
auc_RF <- auc(auc_roc_RF)
precision_RF <- matrice_confusione_RF$byClass["Pos Pred Value"]

print(paste("Accuracy: ", accuracy_RF))
print(paste("Recall: ", recall_RF))
print(paste("F1-Score: ", f1_RF))
print(paste("AUC-ROC: ", auc_RF))
print(paste("Precision: ", precision_RF))















#ATTENZIONE DA QUI------------------------











#GRAFICI----
matrice <- as.data.frame(matrice_confusione_RF$table)

ggplot(data = matrice, aes(x = Reference, y = Prediction, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = Freq), vjust = 1.5) +
  scale_fill_gradient(low = "grey", high = "blue") +
  theme_minimal() +
  labs(title = "Matrice di Confusione per il Modello SVM", x = "Etichetta Vera", y = "Predizione")

metriche <- data.frame(
  Metrica = c("Accuracy", "Recall", "F1-Score", "AUC-ROC", "Precision"),
  Valore = c(accuracy_RF, recall_RF, f1_RF, auc_RF, precision_RF)
)

ggplot(metriche, aes(x = Metrica, y = Valore, fill = Metrica)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  theme_minimal() +
  coord_flip() +
  labs(title = "Metriche di Valutazione per il Modello RF", x = "", y = "Valore")



#tentativo conta dei fiori----



# Assumendo che `ortho` sia la tua lista di raster
# E che il modello SVM sia salvato nella variabile `SVM`

# Funzione per convertire raster in dataframe e fare predizioni
process_and_predict <- function(raster) {
  # Converti il raster in un dataframe, includendo solo le bande usate per l'addestramento
  pixel_data <- terra::as.data.frame(raster, xy = TRUE, na.rm = TRUE)
  
  # Rinomina le colonne per assicurarsi che corrispondano a quelle usate nel modello SVM
  # Assicurati che questo mapping corrisponda ai nomi delle bande nel tuo raster e a come sono stati denominati nel training set
  names(pixel_data) <- c('x', 'y', 'banda_1', 'banda_2', 'banda_3')
  
  # Seleziona solo le colonne di interesse (rimuovi 'x', 'y' se non sono usate per la predizione)
  pixel_data <- dplyr::select(pixel_data, banda_1, banda_2, banda_3)
  
  # Applica il modello SVM per fare predizioni sui pixel
  predizioni <- predict(RF, newdata = pixel_data)
  
  # Ritorna le predizioni
  return(predizioni)
}


# Applica la funzione a ogni raster in `ortho`
risultati_predizioni_RF <- lapply(ortho, process_and_predict)



# Aggrega i risultati delle predizioni per ciascun raster
conteggi_predizioni_RF <- lapply(risultati_predizioni_RF, table)

# Stampa i conteggi per ciascun raster
lapply(conteggi_predizioni_RF, print)


# Imposta l'area coperta da ciascun pixel (in cm²)
area_per_pixel_cm2 <- 0.25

# Supponendo che `conteggi_predizioni` sia la tua lista con i conteggi dei pixel
conteggi_predizioni_cm2 <- lapply(conteggi_predizioni_RF, function(x) x * area_per_pixel_cm2)

# Stampa i risultati convertiti in cm2
lapply(conteggi_predizioni_cm2, print)

#RISULTATI RF-----
#----------------------------------------------------------------------------------------------------

risultati_drone_RF <- data.frame(
  area = c("OrthoEL1", "OrthoEY1", "OrthoEY3", "OrthoEY4", 
           "OrthoFR1", "OrthoFR2", "OrthoFR4", "OrthoGE3",
           "OrthoGU1", "OrthoGU2", "OrthoGU3", "OrthoGU4",
           "OrthoMA2", "OrthoMA3", "OrthoPA2", "OrthoRA4",
           "OrthoSG3", "OrthoSI2", "OrthoVA1", "OrthoVA2",
           "OrthoVA3", "OrthoVA4", "OrthoWA1", "OrthoWY2",
           "OrthoWY3", "OrthoWY4"),
  erba_cm2 = c(628449.00, 362398.8, 869963.5, 6362.25, 
               1464133.75, 2432578.25, 1190456.8, 1563895.5, 
               823014.5, 66856.25, 213787.2, 1396344.2, 
               4657624.5, 897569.0, 423399.5, 735222.2, 
               1693484.00, 558803, 1311872.5, 926867.2, 
               2083497, 1068082, 830415.00, 702545.0, 
               1013272, 371428),
  fiore_cm2 = c(25202.25, 888822.0, 647489.2, 2154184.00, 
                6405.75, 3918.75, 692384.8, 121797.8, 
                1855.0, 1586426.25, 552413.8, 216466.2, 
                96.5, 593614.2, 850195.8, 197719.8, 
                52728.25, 491413, 960315.2, 1421580.2, 
                1721, 433054, 10766.75, 341432.2, 
                48197, 1012250)
)


# Stampa del dataframe per verificare i dati
print(risultati_drone_RF)

# Salvataggio dei risultati in un file CSV
#write.csv(risultati, "Risultati_Aree_2021.csv", row.names = FALSE)

#Correlazioni----


# Carica i dati di campo dal file CSV
dati_campo <- read.csv("/media/data/r_projects/phd_ludovico/2021/dati di campo 2021 modi.csv",sep = ";", header = TRUE)

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
dati_campo <- dati_campo[!dati_campo$area %in% c("EY4", "GU2", "VA2", "WY4", "VA1", "PA2", "EY1", "FR4"), ]
risultati_drone_RF <- risultati_drone_RF[!risultati_drone_RF$area %in% c("EY4", "GU2", "VA2", "WY4", "VA1", "PA2", "EY1", "FR4"), ]
 


#selezioniamo le aree giuste
dati_campo <- dati_campo[, c("area", "flower_shannon_area", "flower_shannon_number", "flower_simpson_area", 
                             "flower_simpson_number", "flower_surface_cm2", "flower_species_richness", "bee_shannon", 
                             "bee_simpson", "bee_species_richness", "bee_abundance")]

# Mostra le prime righe del dataframe creato
head(dati_campo)

# Rimuovi il prefisso "Ortho" dalla colonna area in risultati_drone se necessario
risultati_drone_RF$area <- gsub("Ortho", "", risultati_drone_RF$area)

# Unisci le colonne usando dplyr
dati_congiunti <- left_join(dati_campo, risultati_drone_RF[, c("area", "fiore_cm2")], by = "area")

# Verifica i risultati
head(dati_congiunti)

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
p_value_format <- formatC(p_value, format = "f", digits = 4)

annotazione_testo <- paste("R^2 = ", round(r_quadrato, digits = 3),
                           "\np-value = ", p_value_format)
# Crea il grafico con ggplot2 e aggiungi l'annotazione
ggplot(dati_congiunti, aes(x = flower_surface_cm2, y = fiore_cm2, label = area)) +
  geom_point() +
  geom_smooth(method = "lm", col = "blue") +
  geom_text(check_overlap = TRUE, nudge_y = 0.05, hjust = 0.5, vjust = -0.5) +
  labs(x = "Dati di campo (cm²)", 
       y = "Dati da drone RF (cm²)", 
       title = "flower abundance(field) vs flower abundance UAV (RF)") +
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

























































































































#Valutiamo l'overfitting----

#Confronto tra Prestazioni su Addestramento e Test:
# Calcola le prestazioni sul set di addestramento
train_predictions <- predict(RF, newdata = training_data)
train_conf_matrix <- confusionMatrix(as.factor(train_predictions), training_data$label, positive = "fiore")
print(train_conf_matrix)



#Cross-Validation Esterna
# Esegui k-fold cross-validation

df_selected <- na.omit(df_selected)

set.seed(123)
cv_results <- train(
  label ~ .,
  data = df_selected, # Usa l'intero dataset prima della divisione
  method = "ranger",
  trControl = trainControl(
    method = "cv",
    number = 10,
    classProbs = TRUE,
    summaryFunction = twoClassSummary
  ),
  tuneGrid = tuneGrid,
  metric = "ROC"
)

print(cv_results)

#Verifichiamo la robustezza----
# Ecco un esempio di come ottenere l'importanza delle variabili per un modello 'ranger'

importance <- varImp(RF)
print(importance)


#Confronto delle Prestazioni su Set Diversi:
# Hai giÃ  le metriche di addestramento e di test; ora confrontale con la cross-validation
print(cv_results$results)

#Analisi della Curva ROC:
# Calcola e traccia la curva ROC per il set di test
test_probabilities <- predict(RF, test_data, type = "prob")
roc_test <- roc(response = test_data$label, test_probabilities[, "fiore"])

# Calcola l'AUC
auc(roc_test)

# Aggiungi l'AUC al grafico
plot(roc_test, main = paste("ROC Curve (AUC =", auc(roc_test), ")"))






































































# Le caratteristiche che catturano meglio le differenze di colore 
# (probabilmente le bande 1 e 2, a seconda di come sono state misurate) 
# saranno le piÃ¹ informative per il modello. Se la banda 3 non cambia tra fiori 
# gialli e bianchi o non riflette differenze significative che aiutano a 
# discriminare tra le due classi, allora Ã¨ ragionevole che il modello le assegni 
# un'importanza bassa o nulla.
# 
# Nel contesto del machine learning applicato all'analisi delle immagini, 
# soprattutto in applicazioni come la classificazione di fiori dove il colore Ã¨
# un fattore distintivo importante, le bande specifiche che catturano
# l'informazione rilevante (come quelle che corrispondono alle lunghezze 
# d'onda in cui i fiori hanno il maggiore contrasto) saranno molto piÃ¹ utili
# per il modello.
# 
# Questo puÃ² essere esaminato in termini di:


#2:Caratteristiche Spettrali:

# Per esempio, crea un boxplot per confrontare le distribuzioni delle bande
ggplot(df_selected, aes(x = label, y = banda_1, fill = label)) +
  geom_boxplot() +
  labs(title = "Distribuzione della Banda 1 per Classe", x = "Classe", y = "Valore Banda 1")



# Il boxplot che hai generato mostra la distribuzione dei valori della banda 1 
# per le due classi: "erba" e "fiore". Dall'immagine, sembra che ci sia una 
# differenza significativa nella distribuzione dei valori tra le due classi, 
# il che potrebbe spiegare perchÃ© la banda 1 Ã¨ stata identificata come la piÃ¹
# importante dal tuo modello di classificazione.
# 
# La classe "fiore" mostra valori piÃ¹ alti nella banda 1 rispetto alla classe 
# "erba", il che suggerisce che questa banda cattura una caratteristica distintiva
# che aiuta il modello a discriminare tra le due classi. Questa differenza 
# potrebbe corrispondere alla riflettanza o all'assorbenza luminosa specifica dei
# fiori gialli e bianchi nella banda 1 dello spettro visibile, presumibilmente
# legata al colore giallo.




saveRDS(RF, file = "/home/PERSONALE/ludovico.chieffallo2/Showcase_proj/model_RF.rds")


#per caricarlo RF_loaded <- readRDS(file = "/home/PERSONALE/ludovico.chieffallo2/Showcase_proj/model_rf.rds")




#Funzione a tutti gli elementi in ortho
ortho <- map(ortho, function(x) {
  # Seleziona solo le prime tre bande
  x <- terra::subset(x, 1:3)
  # Rinomina queste bande
  names(x) <- colnames(df_selected)[1:3]
  return(x)
})



#Application of model and save----
t0<-Sys.time()
classificazione_RF_list <- pbapply::pblapply(names(ortho), function(name) {
  terra::predict(ortho[[name]], RF, na.rm=TRUE) #gli NA sono dello sfondo, quindi possiamo escluderli
})
t1<-Sys.time()
t1-t0

#Time difference of 30.24887 mins

#QUI

#classificazione_RF_list <- future_map(ortho[[1]], ~predict(., RF),.progress = TRUE)


# Salva la classificazione come un file GeoTIFF
#output_filename <- "classificazione_RF.tif"
#writeRaster(classificazione_RF, output_filename, overwrite = TRUE)


names(classificazione_RF_list) <- names(ortho)  # Assicuriamoci che i nomi corrispondano


walk(names(classificazione_RF_list), function(name) {
  terra::writeRaster(classificazione_RF_list[[name]], 
                     filename = paste0(name, ".tif"), 
                     overwrite=TRUE)
})



#Continamo i pixel per area----

# Funzione per conteggiare i pixel per ogni categoria
countCategoryPixels <- function(classification, category) {
  sum(values(classification) == category, na.rm = TRUE)
}

# Calcola il conteggio per ogni orthomosaico e per ogni categoria
results <- lapply(classificazione_RF_list, function(classification) {
  # Assumiamo che 1 sia il valore per "fiore" e 2 per "erba"
  count_flowers <- countCategoryPixels(classification, 2)
  count_grass <- countCategoryPixels(classification, 1)
  
  return(data.frame(Flowers = count_flowers, Grass = count_grass))
})

# Combina i risultati per ottenere il totale
total_counts <- do.call(rbind, results)
row.names(total_counts) <- names(classificazione_RF_list) # Opzionale: aggiungi i nomi

# Visualizza o salva i risultati
print(total_counts)


# print(total_counts)
#            Flowers    Grass
# Ortho_EL1    8029  2828865
# Ortho_EY1  327780  4799555
# Ortho_EY2  378059  8402143
# Ortho_FR2   59330 10141289
# Ortho_GU4  188414  6853399
# Ortho_PA2  155590  5052222
# Ortho_RA4   28464  3792101
# Ortho_SG3   36828  7099500
# Ortho_SI2  103035  4192749
# Ortho_VA2   40552  9599024
# Ortho_WA1   25151  3412348
# Ortho_WY2  246831  4021550



























































































#ATTENZIONE proviamo a contare le aree----

# Carico la libreria
library(terra)

# Assumendo che i raster che devi analizzare siano giÃ  presenti nel tuo environment, 
# non c'Ã¨ bisogno di listare e caricare i file come hai fatto nel tuo script originale.

# Creo una struttura dati vuota per memorizzare i risultati
res_fiori <- data.frame(
  fiori = numeric(0), 
  erba = numeric(0),
  area = character(0)
)

# Supponendo che i tuoi raster siano in un oggetto chiamato "classificazione_RF_list"
for (i in 1:length(classificazione_RF_list)) {
  raster_current <- classificazione_RF_list[[i]]
  areaname <- names(classificazione_RF_list)[i]
  ncell(raster_current)
  
  # Convertire il raster in poligoni
  p <- as.polygons(raster_current)
  
  # Contare le celle per ogni categoria
  count <- expanse(p)
  
  # Assumendo che "fiori" siano etichettati come 1 e "erba" come 2 nel tuo raster
  count_fiori <- count[1]
  count_erba <- count[2]
  
  # Aggiungere i risultati al dataframe
  res_fiori <- rbind(res_fiori, data.frame(
    fiori = count_fiori, 
    erba = count_erba, 
    area = areaname
  ))
}

# Salvare i risultati
write.table(res_fiori, "res_fiori.csv", row.names = FALSE)

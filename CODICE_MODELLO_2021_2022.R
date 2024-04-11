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

GBM <- readRDS(file = "/home/PERSONALE/ludovico.chieffallo2/Showcase_proj/models/GBM_model.rds")

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

pathfolder<-"data1"
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


# Funzione per preparare i dati e fare la predizione
predici_ortho <- function(ortho) {
  # Seleziona solo le prime tre bande per ogni immagine
  ortho_rgb <- ortho[[1:3]]
  
  # Converti il raster in dataframe per la predizione
  dati_predizione <- as.data.frame(ortho_rgb, xy = FALSE, na.rm = TRUE)
  
  # Controlla se il dataframe è vuoto
  if (nrow(dati_predizione) == 0) {
    return(NULL)
  }
  
  # Assicurati che i nomi delle colonne corrispondano a quelli attesi dal modello
  names(dati_predizione) <- c("banda_1", "banda_2", "banda_3")
  
  # Applica il modello per la predizione
  predizioni <- predict(GBM, newdata = dati_predizione)
  
  return(predizioni)
}

# Applica la funzione a ciascuna immagine orto
risultati_predizioni <- lapply(ortho, predici_ortho)

# Visualizza i risultati delle predizioni
print(risultati_predizioni)


conteggi_predizioni <- lapply(risultati_predizioni, table)

# Calcola l'area in cm2
area_per_pixel_cm2 <- 0.25
conteggi_predizioni_cm2 <- lapply(conteggi_predizioni, function(x) x * area_per_pixel_cm2)

# Prepara i dati di campo
dati_campo <- data.frame(
    area = c("Ortho_EL1", "Ortho_EY1", "Ortho_EY2", "Ortho_FR2", "Ortho_FR3", 
             "Ortho_GU4", "Ortho_PA2", "Ortho_RA4", "Ortho_SG3", "Ortho_SI2", 
             "Ortho_VA2", "Ortho_WA1", "Ortho_WY1", "Ortho_WY2"),
    erba = c(136.8554864, 172.8088145, 185.3238983, 177.5062205, 179.6626651,
             192.5570632, 197.4387587, 171.6467693, 201.3373421, 164.845817,
             204.4391318, 195.6002123, 178.8172751, 191.2315377),
    fiori = c(1.272015184, 1.54877876, 3.126231016, 1.711925487, 0.129040724,
              10.20119567, 1.42807106, 2.131151049, 2.93460803, 4.724712543,
              0.560481177, 1.862174039, 1.42491097, 10.87105428),
    fiori_rescaled = c(12720.15184, 15487.7876, 31262.31016, 17119.25487, 1290.407237,
                       102011.9567, 14280.7106, 21311.51049, 29346.0803, 47247.12543,
                       5604.811768, 18621.74039, 14249.1097, 108710.5428),
    FC_per_site = c(1293.922476, 20215.86551, 31188.99036, 14410.38212, 925.2031281,
                    81587.75116, 3926.135726, 24150.63874, 113560.394, 48382.92472,
                    4300.175342, 17965.69816, 8897.594522, 32428.54389),
    ABBONDANZA = c(2, 3, 11, 13, 3, 2, 0, 8, 11, 21, 0, 1, 5, 16),
    RICHNESS = c(2, 3, 5, 4, 3, 2, 0, 7, 8, 13, 0, 1, 5, 8),
    SHANNON = c(0.693147181, 1.098612289, 1.159588814, 1.031767111, 1.098612289,
                0.693147181, NA, 1.906154747, 1.972246979, 2.441660565,
                NA, 0, 1.609437912, 1.981332515)
  )
dati_campo$area
names(ortho)
dati_campo <- dati_campo[-c(5,7, 11, 13), ]
dati_campo
# Stampa il dataframe per verificarne il contenuto
print(dati_campo)


risultati_drone_GBM <- data.frame(
  area = c("Ortho_EL1", "Ortho_EY1", "Ortho_EY2", "Ortho_FR2",
           "Ortho_GU4", "Ortho_PA2", "Ortho_RA4", "Ortho_SG3",
           "Ortho_SI2", "Ortho_VA2", "Ortho_WA1", "Ortho_WY2"),
  erba_cm2 = c(708035.5, 1240966.50, 2142361.25, 2544872,
           1739698.8, 1300004.75, 952865.2, 1780104.5,
           1069406, 2407908.75, 855280.2, 1049888.50),
  fiore_cm2 = c(1188.0, 40867.25, 52689.25, 5283,
            20754.5, 1948.25, 2276.0, 3977.5,
            4540, 1985.25, 4094.5, 17206.75)
)

risultati_drone_GBM <- risultati_drone_GBM[-c(6,10), ]


print(risultati_drone_GBM)
risultati_drone_GBM$area <- gsub("Ortho_", "", risultati_drone_GBM$area)
risultati_drone_GBM$area <- gsub("^_", "", risultati_drone_GBM$area)
risultati_drone_GBM
dati_campo$area<- gsub("Ortho_", "", risultati_drone_GBM$area)
dati_campo$area <- gsub("^_", "", risultati_drone_GBM$area)
dati_campo

dati_congiunti_2021_2022 <- left_join(dati_campo, risultati_drone_GBM[, c("area", "fiore_cm2")], by = "area")
dati_congiunti_2021_2022


#------------------------------------------------------------------------------------
#FLOWER AB./FLOWER AB. GBM-----

# Calcola la correlazione tra le due stime dei cm² dei fiori
correlazione <- cor(dati_congiunti_2021_2022$fiori_rescaled, dati_congiunti_2021_2022$fiore_cm2)

# Stampa il coefficiente di correlazione
print(correlazione)


# Costruisci il modello lineare
modello_lineare <- lm(fiore_cm2 ~ flower_surface_cm2, data = dati_congiunti_2021_2022)

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
  annotate("text", x = max(dati_congiunti$bee_shannon, na.rm = TRUE), y = min(dati_congiunti_2021_2022$fiore_cm2, na.rm = TRUE), 
           label = paste("R² =", round(r_quadrato_shannon, 3), "\np-value =", format(p_value_shannon, scientific = FALSE)), 
           hjust = 1, vjust = -1, color = "red", size = 5) +
  theme_minimal()

library(caret)
library(raster)
install.packages("randomForest") 
library(randomForest)


brick_for_prediction_norm<-raster("D:/Desktop/DOTTORATO/1.SHOWCASE DATI/DATI 2022/orthofoto/DAG5/27.EY1/gis/file da testare machine learning ey1_1.tif")
brick_for_prediction_norm
brick_for_prediction_norm<-stack("D:/Desktop/ey1_buffer.tif")
brick_for_prediction_norm <- dropLayer(brick_for_prediction_norm, c(4))
plot(brick_for_prediction_norm)
names(brick_for_prediction_norm)<- c("EY1_rgb_buffer1m.1","EY1_rgb_buffer1m.2","EY1_rgb_buffer1m.3")
#modello----

# Imposta il percorso del file
file_path <- "D:/Desktop/model_rf_EY1.1.rds"

# Carica il file .RDS
model_rf_ey1 <- readRDS(file_path)

model_rf_ey1$coefnames


#----

#crop----
area_buffer<-shapefile("D:/Desktop/DOTTORATO/1.SHOWCASE DATI/DATI 2022/orthofoto/DAG5/27.EY1/gis/buffer.shp")
area_buffer

ey1_cropped<-crop(brick_for_prediction_norm, extent(area_buffer))
ey1_cropped_2<-mask(ey1_cropped, area_buffer)
#-----
writeRaster(ey1_cropped_2,"ey1_buffer.tif")
getwd()


predict_rf <- raster::predict(object = brick_for_prediction_norm,
                              model = model_rf_ey1, type = 'raw')
setwd("D:/Desktop")
writeRaster(predict_rf, "r_ey1_rf.tif")
brick_for_prediction_norm

predict_rf








#system.time({
#  predict_rf <- raster::predict(object = brick_for_prediction_norm,
#                                model = model_rf, type = 'raw')
#  predict_svm <- raster::predict(object = brick_for_prediction_norm,
#                                 model = model_svm, type = 'raw')
#  predict_nnet <- raster::predict(object = brick_for_prediction_norm,
#                                  model = model_nnet, type = 'raw')
#})


writeRaster(predict_rf, paste0("G:/SHOWCASE_dati/result/machine_learning/1_2_fixed/rf/original/",nome,"_new_rf_1500.tiff"),overwrite=T )
writeRaster(predict_svm, paste0("G:/SHOWCASE_dati/result/machine_learning/1_2_fixed/svm/original//",nome,"_new_svm_1500.tiff"),overwrite=T )
writeRaster(predict_nnet, paste0("G:/SHOWCASE_dati/result/machine_learning/1_2_fixed/nnet/original/",nome,"_new_nnet_1500.tiff"),overwrite=T )


#bick_for_prediction_norm<-inserire il brick del raster rgb clippato al buffer. vedere se ha 3 o 4 bande (togli la quarta xon brick_for_prediction_norm <- dropLayer(brick_for_prediction_norm, c(4)))


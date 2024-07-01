rm(list = ls())
gc(reset=TRUE)
source('snow_functions.R')
library(matrixStats)
library(landsat)
library(imager)
library(raster)
library(R.matlab)
library(gsignal)
library(pracma)



maindir <- '~/dev/TopoVegSnow/CleanCode/'

#Location of each image: snow depth, ground surface elevation, 
#distance to shrubs with height = 1m+, and vegetation height 
file_snow <- paste(maindir,"data/snow.txt",sep = '')
file_topo <- paste(maindir,"data/topo.txt",sep = '')
file_d100 <- paste(maindir,"data/d100.txt",sep = '')
file_vegheight <- paste(maindir,"data/veg_height.txt",sep = '')

#Definition of the lengths used for the stacked directional filtering
lengths <- c(3, 5, 7, 9, 11, 15, seq(from = 21, to = 231, by = 10))
#Number of pixels to remove along boundary after filtering
cut_idx <- 56

#------------------------------------------------------------------------------#
#Perform stacked directional filtering of ground surface elevation and snow
#depth, or load files if they has already been processed
#Inputs: ground surface elevation, snow depth

filename = paste(maindir,"outputs/FilteredImages.rds",sep = '')
if (file.exists(filename)) {
  fi <- readRDS(filename)
  topo_filtered<-fi$topo_filtered
  snow_filtered<-fi$snow_filtered
  image <- readMat(file_snow)
  image <- image$snow
  snow  <- cut_edges(image, cut_idx)
  rm(fi) 
  
} else {
  image <- unname(as.matrix(read.table(file_topo)))
  topo_filtered <- stacked_directional_filter(image, lengths,TRUE)
  names(topo_filtered$filtered_images_1) <- paste('TF1_', lengths, sep = "")
  names(topo_filtered$filtered_images_2) <- paste('TF2_', lengths, sep = "")
  names(topo_filtered$filtered_slopes) <- c(paste('SL_', lengths, sep = ""),'SL_BASE')
  names(topo_filtered$filtered_asps) <- c(paste('ASP_', lengths, sep = ""),'ASP_BASE')
  
  image <-  unname(as.matrix(read.table(file_snow)))
  snow  <- cut_edges(image, cut_idx)
  snow_filtered <- stacked_directional_filter(image, lengths,FALSE)
  names(snow_filtered$filtered_images_1) <- paste('SF1_', lengths, sep = "")
  names(snow_filtered$filtered_images_2) <- paste('SF2_', lengths, sep = "")
  saveRDS(list('topo_filtered'=topo_filtered,'snow_filtered'=snow_filtered),
          file=filename)
}

#------------------------------------------------------------------------------#
#Fit simple model, or load file if this has already been processed
#Inputs: filtered topography, filtered snow depth, distance from 1m vegetation
filename = paste(maindir,"outputs/SimpleModel.rds",sep = '')
if (file.exists(filename)) {
  sm <- readRDS(filename)
  sm_outs <- sm$sm_outs
  sm_eval <- sm$sm_eval
  rm(sm)

} else {
  d100 <-  unname(as.matrix(read.table(file_d100)))
  sm_outs <- calculate_simple_model(topo_filtered$filtered_images_1,
                                    topo_filtered$filtered_images_2,
                                    snow_filtered$filtered_images_1,
                                    snow_filtered$filtered_images_2,
                                    d100, cut_idx)
  
  sm_eval <- evaluate_simple_model(sm_outs$SM,image)
  saveRDS(list('sm_outs'=sm_outs,'sm_eval'=sm_eval),file=filename)
}

#------------------------------------------------------------------------------#
#Transform vegetation height into canopy trapping fields, or load file if this
#has already been processed
#Inputs: vegetation height

filename = paste(maindir,"outputs/VegFields.rds",sep = '')
if (file.exists(filename)) {
  veg <- readRDS(filename)
  ctfield1 <- veg$ctfield1
  ctfield2 <- veg$ctfield2
  ctgrad1 <- veg$ctgrad1
  ctgrad2 <- veg$ctgrad2
  rm(veg)
  
} else {
  veg_height <- unname(as.matrix(read.table(file_vegheight)))
  veg_height[veg_height<0] <- 0
  ctfield1 <- calculate_canopy_trapping_field(veg_height,1)
  ctgrad1 <- gradient(ctfield1,1,1)
  ctfield2 <-  calculate_canopy_trapping_field(veg_height,2)
  ctgrad2 <- gradient(ctfield2,1,1)
  
  saveRDS(list('ctfield1'=ctfield1,'ctfield2'=ctfield2,'ctgrad1'=ctgrad1,
               'ctgrad2'=ctgrad2),file=filename)
}

#------------------------------------------------------------------------------#
#Train machine learning ensemble and predict snow depths in the presence and 
#in the absence of tall shrubs


#Combine data needed for ML model into single dataframe
data <- data.frame('SM'= as.vector(cut_edges(sm_outs$SM,cut_idx)))
data <- cbind(data,data.frame(lapply(lapply(topo_filtered$filtered_images_1,cut_edges,cut_idx =cut_idx),as.vector)))
data <- cbind(data,data.frame(lapply(lapply(topo_filtered$filtered_images_2,cut_edges,cut_idx =cut_idx),as.vector)))
data <- cbind(data,data.frame(as.vector(cut_edges(topo_filtered$filtered_base,cut_idx))))
data <- cbind(data,data.frame(lapply(lapply(topo_filtered$filtered_slopes,cut_edges,cut_idx =cut_idx),as.vector)))
data <- cbind(data,data.frame(lapply(lapply(topo_filtered$filtered_asps,cut_edges,cut_idx =cut_idx),as.vector)))
data <-cbind(data, data.frame('ctfield1'=as.vector(cut_edges(ctfield1,cut_idx)),
                              'ctgrad1.X'=as.vector(cut_edges(ctgrad1$X,cut_idx)),
                              'ctgrad1.Y'=as.vector(cut_edges(ctgrad1$Y,cut_idx)),
                              'ctfield2'=as.vector(cut_edges(ctfield2,cut_idx)),
                              'ctgrad2.X'=as.vector(cut_edges(ctgrad2$X,cut_idx)),
                              'ctgrad2.Y'=as.vector(cut_edges(ctgrad2$Y,cut_idx))),
             data.frame('Snow' = as.vector(snow)))


ML_output_list <- list()
gc()
# Number of iterations
num_iterations <- 10

# Run the function iteratively and save the output
for (i in 1:num_iterations) {
  
  output <- snow_ML_model(data,15000)
  ML_output_list[[i]] <- output
}

saveRDS(ML_output_list,file='../outputs/ML_outputs.rds')

ML_output_list <- readRDS('../outputs/ML_outputs.rds')


#Set canopy snow trapping potential to low, constant value (i.e. predict snow 
#depth in the absence of shrubs)
data$ctfield1 <- quantile(data$ctfield1,0.05)
data$ctfield2 <- quantile(data$ctfield2,0.05)
data$ctgradmag1 <-0
data$ctgradmag2 <-0
#data$ctgrad1.X <- 0
#data$ctgrad1.Y <- 0
#data$ctgrad2.X <- 0
#data$ctgrad2.Y <- 0
preds_noveg <- list()
for (i in 1:num_iterations) {
  brt.snow <- ML_output_list[[i]]$brt.snow
  gbm.predict.grids(brt.snow, data, want.grids = F, sp.name =
                      "preds")
  preds_noveg[[i]] <- preds
}
saveRDS(preds_noveg,file='../outputs/preds_noveg.rds')

#ML predictions of shrub canopy snow trapping can be determined by subtracting
#the ML ensemble predictions in the absence of shrubs from the predictions in
#the presence of shrubs

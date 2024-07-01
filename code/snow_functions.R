library(matrixStats)
library(landsat)
library(imager)
library(raster)
library(gsignal)
source('brt.functions.R')

mirror_horizontal <- function(array) {
  return(array[, ncol(array):1])
}

mirror_vertical <- function(array) {
  return(array[nrow(array):1, ])
}


expand_array <- function(array){
  out<- rbind(mirror_vertical(cbind(mirror_horizontal(array),array,mirror_horizontal(array))),
              cbind(mirror_horizontal(array),array,mirror_horizontal(array)),
              mirror_vertical(cbind(mirror_horizontal(array),array,mirror_horizontal(array))))
  return(out)
}

reduce_array <- function(array,shape){
  out<- array[(shape[1]+1):(shape[1]*2),(shape[2]+1):(shape[2]*2)]  
  return(out)
}

cut_edges <- function(array, cut_idx){
  return(array[cut_idx:(nrow(array)-cut_idx),cut_idx:(ncol(array)-cut_idx)])
}

#------------------------------------------------------------------------------#
stacked_directional_filter <- function(image,lengths,slope_asp){
  #image: input array that will be filtered
  #lengths: array of progressively increasing filter lengths
  #lats: array of lats (same shape as image)
  #lons: array of lons (same shape as image)
  #slope_asp: if True, return slope and aspect of filtered image at each scale
  temp <- image
  shape <- dim(image)
  #tile image to reduce boundary effects
  temp <- expand_array(temp)
  imf1 <- list()
  imf2 <- list()
  imf_sum <- matrix(0, nrow = nrow(temp), ncol = ncol(temp))
  if (slope_asp) {
    slopes <- list()
    aspects <- list()
  }
  
  for (l in 1:length(lengths)) {
    if (slope_asp){
      aspect_slope_grad <- slopeasp(reduce_array(temp,shape),1,1)
      aspect <- aspect_slope_grad$aspect
      slope <- aspect_slope_grad$slope
      aspects[[l]] <- aspect
      slopes[[l]] <- slope
    } 
    
    # Create the linear filter kernel in direction 1
    filter_kernel <- matrix(1, nrow = lengths[l], ncol = 1) / lengths[l]
    # Apply the filter using convolution
    mapfilt_l <- conv2(temp, filter_kernel,shape='same')
    imf1 [[l]] <- reduce_array(temp - mapfilt_l,shape)
    imf_sum <- imf_sum + temp- mapfilt_l
    temp<- mapfilt_l
    # Create the linear filter kernel in direction 2
    filter_kernel <- matrix(1, nrow = 1, ncol = lengths[l]) / lengths[l]
    # Apply the filter using convolution
    mapfilt_l <- conv2(temp, filter_kernel,shape='same')
    imf2 [[l]] <- reduce_array(temp- mapfilt_l,shape)
    imf_sum <- imf_sum + temp- mapfilt_l
    temp<- mapfilt_l
  }
  imf_base <- reduce_array(temp,shape)
  out <- list('filtered_images_1'=imf1,'filtered_images_2' = imf2,
              'filtered_base'=imf_base)
  if (slope_asp){
    aspect_slope_grad <- slopeasp(imf_base,1,1)
    aspects[[length(lengths) + 1]] <- aspect_slope_grad$aspect
    slopes[[length(lengths) + 1]] <- aspect_slope_grad$slope
    out <- list('filtered_images_1'=imf1,'filtered_images_2' = imf2,
                'filtered_base'=imf_base,
                'filtered_slopes'=slopes,'filtered_asps' = aspects)
  }
  return(out)
}

#------------------------------------------------------------------------------#
calculate_simple_model <- function(tf1, tf2, sf1, sf2,d100,cut_idx){
  # Initialize lists to store results
  lms <- list()
  ms <- numeric()
  bs <- numeric()
  r2s <- numeric()
  rmses <- numeric()
  sds <- numeric()
  sdms <- numeric()
  
  # Initialize sm matrix
  shape <- dim(tf1[[1]])
  sm <- matrix(0, nrow = shape[1], ncol = shape[2])
  
  # Iterate over j values
  for (j in 1:length(tf1)) {
    X <- cut_edges(tf1[[j]] + tf2[[j]],cut_idx)
    
    Y <- cut_edges(sf1[[j]] + sf2[[j]],cut_idx)
    
    X2 <- tf1[[j]] + tf2[[j]] 
    lm <- lm(Y[d100 > 25] ~ 0+X[d100 > 25])
    #lms[[j]] <- lm
    ms[j] <- coef(lm)
    r2s[j] <- summary(lm)$r.squared
    rmses[j] <- summary(lm)$sigma / sd(Y)
    
    sm <- sm + X2 * coef(lm)
  }
  
  out <- list('SM'= sm,  'ms'=ms,'r2s'=r2s,'rmses'=rmses)
  return(out)
}

#------------------------------------------------------------------------------#
evaluate_simple_model <- function(sm, snow){
  # Initialize arrays for ms, bs, r2s, rmses
  shape <- dim(sm)
  idxs1  <- seq(from = 51, to = shape[1]-50, by = 10)
  idxs2 <- seq(from = 51, to = shape[2]-50, by = 10)
  ms <- matrix(0, nrow = length(idxs1), ncol = length(idxs2))
  bs <- matrix(0, nrow = length(idxs1), ncol = length(idxs2))
  r2s <- matrix(0, nrow = length(idxs1), ncol = length(idxs2))
  rmses <- matrix(0, nrow = length(idxs1), ncol = length(idxs2))
  
  # Iterate over m and n values
  for (m in 1:length(idxs1)) {
    for (n in 1:length(idxs2)) {
      idx1 <- idxs1[m]
      idx2 <- idxs2[n]
      X <- sm[(idx1 - 50):(idx1 + 49), (idx2 - 50):(idx2 + 49)]
      Y <- snow[(idx1 - 50):(idx1 + 49), (idx2 - 50):(idx2 + 49)]
      
      lm <- lm(as.vector(Y) ~ as.vector( X))
      ms[m, n] <- coef(lm)[2]
      bs[m, n] <- coef(lm)[1]
      r2s[m, n] <- summary(lm)$r.squared
      rmses[m, n] <- summary(lm)$sigma / sd(Y)
    }
  }
  out <- list('ms'=ms, 'bs'=bs, 'r2s'=r2s, 'rmses'=rmses)
  return(out)
}

#------------------------------------------------------------------------------#
calculate_canopy_trapping_field <- function(veg_height,p){
  # Initialize ctfield with zeros
  ctfield <- matrix(0, nrow = nrow(veg_height), ncol = ncol(veg_height))
  
  # Get the shape of the full array
  shape_full <- dim(veg_height)
  q_all <- matrix( pmin(1.5, veg_height), nrow = shape_full[1], ncol = shape_full[2])
  # Loop through rows and columns
  for (i in 1:shape_full[1]) {
    for (j in 1:shape_full[2]) {
      q <- q_all[i, j]
#      for (k in 0:49) {
#        for (l in 0:49) {
      for (k in 0:49) {
        for (l in 0:49) {
          if ((k==0)&(l==0)){
            val <- 2 *q
          }else{
          val <- min(2 * q, q / (k^2 + l^2)^(0.5*p))
          }
          if ((i + k) <= shape_full[1] && (j + l) <= shape_full[2]) {
            ctfield[i + k, j + l] <- ctfield[i + k, j + l] + val
          }
          if ((i - k) >= 1 && (j + l) <= shape_full[2]) {
            ctfield[i - k, j + l] <- ctfield[i - k, j + l] + val
          }
          if ((i - k) >= 1 && (j - l) >= 1) {
            ctfield[i - k, j - l] <- ctfield[i - k, j - l] + val
          }
          if ((i + k) <= shape_full[1] && (j - l) >= 1) {
            ctfield[i + k, j - l] <- ctfield[i + k, j - l] + val
          }
        }
      }
    }
  }
  
  ctfield <- ctfield / max(ctfield, na.rm = TRUE)
  return(ctfield)
  
}

#------------------------------------------------------------------------------#
snow_ML_model <- function(data,numpoints){
  
  dt = sort(sample(nrow(na.omit(data)),numpoints))
  train<-na.omit(data)[dt,]
  
  brt.snow <- gbm.step(data=train,
                            gbm.x = c(1:(ncol(train)-1)),
                            gbm.y = ncol(train),
                            family = "gaussian",
                            tree.complexity = 3,
                            learning.rate = 0.5,
                            bag.fraction = 0.5)
  
  gbm.predict.grids(brt.snow, data, want.grids = F, sp.name =
                      "preds")
  
  return(list('brt.snow'=brt.snow,'preds'=preds))
}

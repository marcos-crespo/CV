library(fda.usc)
library(fsemipar)
library(dplyr)
library(tidyr)
library(ggplot2)


# Auxiliary functions -----------------------------------------------------


# Auxiliary function to compute integral from a fda object
compute_integral <- function(fd_obj, curve_index) {
  # Extract the curve from the fd object
  curve_fun <- function(x) {
    eval.fd(x, fd_obj[curve_index])
  }
  
  # Compute the integral over the range
  integral_value <- integrate(curve_fun, lower = 1, upper = 100)$value
  return(integral_value)
}

# Define the MAR mechanism function
apply_mar_mechanism <- function(data, gamma_value) {
  
  # Extract the predictors and the integral values
  X1 <- data$Protein
  X2 <- data$Water
  spectral_integrals <- data$Spectral_Integral
  
  # Calculate the logistic probability of missingness
  linear_predictor <- gamma_value * (X1 + X2 + spectral_integrals)
  prob_missing <- 1 / (1 + exp(linear_predictor))
  
  # Introduce missingness in the 'Fat' column
  set.seed(123)
  missing_indices <- which(runif(nrow(data)) < prob_missing)
  data$Fat[missing_indices] <- NA
  
  return(data)
}


# Tecator fitting ---------------------------------------------------------

fit_tecator<- function() {
  
  data('tecator', package = "fda.usc")
  data <- cbind(tecator$absorp.fdata$data, tecator$y)
  
  # Train - test split of the data
  trainData <- data[1:160,]
  testData <- data[161:215,]
  
  print('Data loaded and split')
  
  ### Plotting
  # Convert the data into a long format for ggplot2
  train_long <- as.data.frame(trainData[1:160,]) %>%
    mutate(id = row_number()) %>%
    pivot_longer(cols = 1:100, names_to = "Wavelength", values_to = "Absorbance") %>%
    mutate(Wavelength = as.numeric(Wavelength))
  # Plot
  ggplot(train_long, aes(x = Wavelength, y = Absorbance, group = id)) +
    geom_line(color = "darkslategray4", alpha = 0.5) +
    labs(
      x = "Wavelength discretization point",
      y = "Absorbance"
    ) +
    theme_minimal(base_size = 15) +
    scale_y_continuous(limits = c(2, 5.5))
  
  ### Outlier detection
  # Basis expansion
  basis <- create.bspline.basis(rangeval = c(1, 100), nbasis = 10)
  # Smooth the data to create a functional data object
  fd_object <- smooth.basis(1:100, t(trainData[, 1:100]), basis)$fd
  # Detect outliers using depth trimming
  outlier_detection <- outliers.depth.trim(fd_object)
  # Extract the indices of detected outliers
  outlier_indices <- outlier_detection$outliers
  non_outlier_indices <- setdiff(1:nrow(trainData), outlier_indices)
  
  # Check the detected outliers
  print('Outlier analysis done')
  print(outlier_indices)
  # Convert the data into a long format for ggplot2
  train_long <- as.data.frame(trainData) %>%
    mutate(id = row_number()) %>%
    pivot_longer(cols = 1:100, names_to = "Wavelength", values_to = "Absorbance") %>%
    mutate(Wavelength = as.numeric(Wavelength),
           Outlier = ifelse(id %in% outlier_indices, "Outlier", "Non-Outlier"))
  # Plot
  ggplot(train_long, aes(x = Wavelength, y = Absorbance, group = id, color = Outlier)) +
    geom_line(alpha = 0.8) +
    labs(
      x = "Wavelength discretization point",
      y = "Absorbance"
    ) +
    scale_color_manual(values = c("Non-Outlier" = "darkslategray4", "Outlier" = "red")) +
    theme_minimal(base_size = 15) +
    scale_y_continuous(limits = c(2, 5.5))
  
  ### Introduce some MAR responses.
  # Compute integrals for each curve
  train_Outlier <- trainData[non_outlier_indices,]
  # Smooth the data to create a functional data object
  fd_object2 <- smooth.basis(1:100, t(train_Outlier[, 1:100]), basis)$fd
  integrals <- sapply(1:length(non_outlier_indices), function(i) compute_integral(fd_object2, i))
  
  train_Outlier$Spectral_Integral <- integrals
  # Apply the MAR mechanism to the training data. 0.0035 -> 30 missings. 0.0008 -> 70
  train_MAR <- apply_mar_mechanism(train_Outlier,0.0008)
  
  print(paste0('NAs introduced: ', sum(is.na(train_MAR$Fat))))
  
  train_filtered <- na.omit(train_MAR)
  
  ###### Kernel part #######
  MAR_model_kernel <- sfplsim.kernel.fit(
    train_filtered[, 1:100],
    train_filtered[, 102:103],
    train_filtered$Fat,
    order.Bspline = 4,
    nknot.theta = 4,
    nknot = 20,
    range.grid = c(850, 1050)
  )
  print('kernel mar fit')

  pred <- predict(MAR_model_kernel,
                  newdata.x = as.matrix(train_MAR[is.na(train_MAR$Fat), 1:100]),
                  newdata.z = as.matrix(train_MAR[is.na(train_MAR$Fat), 102:103]))
  na_indices <- which(is.na(train_MAR$Fat))
  train_imputed <- train_MAR
  train_imputed$Fat[na_indices] <- pred
  
  Fitted_kernel_Model <- sfplsim.kernel.fit(
    train_imputed[, 1:100],
    train_imputed[, 102:103],
    train_imputed$Fat,
    order.Bspline = 4,
    nknot.theta = 4,
    nknot = 20,
    range.grid = c(850, 1050)
  )
  
  print('kernel imputed fit')

  
  # SE between real and estimated MAR responses
  impValues_kernel_mse <- sum((pred - train_Outlier$Fat[na_indices])^2)/length(na_indices)
  
  ### Test predictions with MAR and Fitted models
  pred_kernel_mar <- predict(
    MAR_model_kernel,
    newdata.x = as.matrix(testData[, 1:100]),
    newdata.z = as.matrix(testData[, 102:103]),
    y.test = testData$Fat
  )
  pred_kernel_imp <- predict(
    Fitted_kernel_Model,
    newdata.x = as.matrix(testData[, 1:100]),
    newdata.z = as.matrix(testData[, 102:103]),
    y.test = testData$Fat
  )
  
  ###### KNN part ######
  
  MAR_model_knn <- sfplsim.kNN.fit(
    train_filtered[, 1:100],
    train_filtered[, 102:103],
    train_filtered$Fat,
    order.Bspline = 4,
    nknot.theta = 4,
    nknot = 20,
    range.grid = c(850, 1050)
  )
  print('knn mar fit')

  pred <- predict(MAR_model_knn,
                  newdata.x = as.matrix(train_MAR[is.na(train_MAR$Fat), 1:100]),
                  newdata.z = as.matrix(train_MAR[is.na(train_MAR$Fat), 102:103]))
  
  na_indices <- which(is.na(train_MAR$Fat))
  train_imputed <- train_MAR
  train_imputed$Fat[na_indices] <- pred
  
  Fitted_knn_Model <- sfplsim.kNN.fit(
    train_imputed[, 1:100],
    train_imputed[, 102:103],
    train_imputed$Fat,
    order.Bspline = 4,
    nknot.theta = 4,
    nknot = 20,
    range.grid = c(850, 1050)
  )
  
  print('knn imputed fit')

  # se between real and estimated MAR responses
  impValues_knn_mse <- sum((pred - train_Outlier$Fat[na_indices])^2)/length(na_indices)
  
  ### Test predictions with MAR and Fitted models
  pred_knn_mar <- predict(
    MAR_model_knn,
    newdata.x = as.matrix(testData[, 1:100]),
    newdata.z = as.matrix(testData[, 102:103]),
    y.test = testData$Fat
  )
  pred_knn_imp <- predict(
    Fitted_knn_Model,
    newdata.x = as.matrix(testData[, 1:100]),
    newdata.z = as.matrix(testData[, 102:103]),
    y.test = testData$Fat
  )
  
  ### Full model Benchmark ###
  Full_kernel_Model <- sfplsim.kernel.fit(
    trainData[, 1:100],
    trainData[, 102:103],
    trainData$Fat,
    order.Bspline = 4,
    nknot.theta = 4,
    nknot = 20,
    range.grid = c(850, 1050)
  )
  print('full kernel fitted')

  
  Full_knn_Model <- sfplsim.kNN.fit(
    trainData[, 1:100],
    trainData[, 102:103],
    trainData$Fat,
    order.Bspline = 4,
    nknot.theta = 4,
    nknot = 20,
    range.grid = c(850, 1050)
  )
  
  print('full knn fitted')

  
  pred_kernel_full <- predict(
    Full_kernel_Model,
    newdata.x = as.matrix(testData[, 1:100]),
    newdata.z = as.matrix(testData[, 102:103]),
    y.test = testData$Fat
  )
  pred_knn_full <- predict(
    Full_knn_Model,
    newdata.x = as.matrix(testData[, 1:100]),
    newdata.z = as.matrix(testData[, 102:103]),
    y.test = testData$Fat
  )
  
  
  list(
    imputation_kernel_error = impValues_kernel_mse,
    imputation_knn_error = impValues_knn_mse,
    kernelMAR_test_MSEP = pred_kernel_mar$MSEP.1,
    kernelIMP_test_MSEP = pred_kernel_imp$MSEP.1,
    kernelFULL_MSEP = pred_kernel_full$MSEP.1,
    knnMAR_test_MSEP = pred_knn_mar$MSEP.1,
    knnIMP_test_MSEP = pred_knn_imp$MSEP.1,
    knnFULL_MSEP = pred_knn_full$MSEP.1
  )
  
}


fit_tecator()


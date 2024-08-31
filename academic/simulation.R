# Load necessary libraries
library(fda)
library(fsemipar)   
library(dplyr)
library(MASS)


# Set wd for saving intermediate steps
setwd('C:/Users/marco/Documents/Universidad/Master/02. TFM/scripts')


# Auxiliary functions -----------------------------------------------------

# Create the grid of n observations in t time stamps for function fun
create_grid <- function(fun, a, b, t) {
  outer(a, t, function(ai, ti)
    fun(ti, ai, b))
}

# Our functional random variable
f <- function(x, a, b) {
  a * (x - 0.5)^2 + b
}

# Define theta as a function formed by a B-spline basis
basis <- create.bspline.basis(rangeval = c(0, 1), nbasis = 6)
alpha_0 <- c(0, 1.741539, 0, 1.741539, -1.741539, -1.741539)

theta_0 <- function(t) {
  eval.basis(t, basis) %*% alpha_0
}

# Calculates r(<chi,theta>) as int^3.
r <- function(f, theta0, a, b) {
  integrand <- function(x)
    f(x, a, b) * theta0(x)
  result <- integrate(integrand, lower = 0, upper = 1)$value
  result ^ 3
}

# SFPLSIM formula to one row.
model_index <- function(i, x_lin, betas, r, f, theta0, a, b) {
  sum(x_lin[i, ] * betas) + r(f, theta0, a, b)
}

# Expit formula to probability
expit <- function(u) {
  exp(u) / (1 + exp(u))
}

# Function to create one dataset. The true betas must be adjusted if k is changed.
# generates n+25 observations. n train and 20 test observations.
generate_df <- function(f, n, t, j, cor, c) {
  dim <- n + 25
  # function parameters
  a <- runif(dim, min = -2, max = 2)
  b <- rnorm(dim)
  # function discretization
  grid <- create_grid(f, a, b, t)
  
  # linear covariates
  covariance_matrix <- diag(1, j) + matrix(cor, j, j) - diag(cor, j)
  x_lin <- mvrnorm(dim, mu = rep(0, j), Sigma = covariance_matrix)
  betas <- c(-1, .7, 2)
  
  # generate response
  y <- sapply(1:dim, function(i)
    model_index(i, x_lin, betas, r, f, theta_0, a[i], b[i]))
  err <- rnorm(dim, mean = 0, sd = c * sd(y)) # empirical variance
  y <- y + err
  
  # integral value to use in the MAR introduction
  funs2_value <- c()
  for (i in 1:dim) {
    temp<- function(x) f(x,a[i],b[i])^2
    funs2_value[i]= integrate(temp, 0, 1)$value
  }
  
  temp <- data.frame(y, x_lin, grid, funs2_value)
  
  # create train test partition
  test_indices <- sample(seq_len(nrow(temp)), size = 25)
  list(train = temp[-test_indices, ], test = temp[test_indices, ])
}


# Function to create MAR in a copy of the train dataset df
generate_MAR <- function(df, alpha) {
  train.MAR <- df
  combined_value <- df$X1 + df$X2 + df$X3 + df$funs2_value
  probs <- expit(2 * alpha / pi * combined_value)
  miss_indices <- which(probs < 1 - mean(probs))
  train.MAR$y[miss_indices] <- NA
  train.MAR
}



# Data generation ---------------------------------------------------------

# Define the Data-Generating Process
set.seed(123)

# Parameters
n <- c(50, 200)
t <- seq(0, 1, length.out = 100)
j <- 3
p <- c(.01,2)
cor <- c(0, .5)
c <- c(.01, .05)

# Generate Data Sets
num_simulations <- 100
data_list <- list()

# Generate a list of lists. In each element we have a list containing train-test partition in data, train with MAR in mar,
# and the parameters in dim, prob, cor and c.
index <- 1
for (size in n) {
  for (prob in p) {
    for (corr in cor) {
      for (sig_noise in c) {
        for (i in 1:num_simulations) {
          data <- generate_df(f, size, t, j, corr, sig_noise)
          print(paste0('data gen ', index))
          mar <- generate_MAR(data$train,prob)
          data_list[[index]] <- list(data = data,
                                     mar = mar,
                                     dim = size,
                                     probability = prob,
                                     cor = corr,
                                     c = sig_noise)
          index <- index +1
        }
      }
    }
  }
}

# Save the data
save(data_list, file='data_list.Rda')


# Models fit --------------------------------------------------------------

# Function to impute the MAR and then fit the model.
fit_MAR<- function(data_list_elem) {
  # Drop the integral values from the df and extract train and test from the list of data.
  last_index <- ncol(data_list_elem$data$train)
  train <- data_list_elem$data$train[,-last_index]
  test <- data_list_elem$data$test[,-last_index]
  train_MAR <- data_list_elem$mar[,-last_index]
  
  # Train with ommited MAR responses
  train_filtered <- na.omit(train_MAR)
  
  
  ###### Kernel part ########
  MAR_model_kernel <- sfplsim.kernel.fit(
    as.matrix(train_filtered[, 5:ncol(train_filtered)]),
    as.matrix(train_filtered[, 2:4]),
    train_filtered$y,
    nknot.theta = 3,
    nknot = 20,
    range.grid = c(0, 1)
  )
  print('kernel mar fit')
  pred <- predict(MAR_model_kernel,
                  newdata.x = as.matrix(train_MAR[is.na(train_MAR$y), 5:ncol(train_MAR)]),
                  newdata.z = as.matrix(train_MAR[is.na(train_MAR$y), 2:4]))
  
  na_indices <- which(is.na(train_MAR$y))
  train_imputed <- train_MAR
  train_imputed$y[na_indices] <- pred
  
  Fitted_kernel_Model <- sfplsim.kernel.fit(
    as.matrix(train_imputed[, 5:ncol(train_imputed)]),
    as.matrix(train_imputed[, 2:4]),
    train_imputed$y,
    nknot.theta = 3,
    nknot = 20,
    range.grid = c(0, 1)
  )
  
  print('imputed fit')
  
  # se between real and (full) estimated betas
  beta_kernel_se <- sum((Fitted_kernel_Model$beta.est - c(-1, .7, 2)) ^ 2)
  
  # mse between real and estimated MAR responses
  impValues_kernel_mse <- sum((pred - train$y[na_indices])^2)/length(na_indices)
  
  ### Test predictions with MAR and Fitted models
  pred_kernel_mar <- predict(
    MAR_model_kernel,
    newdata.x = as.matrix(test[, 5:ncol(test)]),
    newdata.z = as.matrix(test[, 2:4]),
    y.test = test$y
  )
  pred_kernel_imp <- predict(
    Fitted_kernel_Model,
    newdata.x = as.matrix(test[, 5:ncol(test)]),
    newdata.z = as.matrix(test[, 2:4]),
    y.test = test$y
  )
  
  ###### KNN part ########
  
  MAR_model_knn <- sfplsim.kNN.fit(
    as.matrix(train_filtered[, 5:ncol(train_filtered)]),
    as.matrix(train_filtered[, 2:4]),
    train_filtered$y,
    nknot.theta = 3,
    nknot = 20,
    range.grid = c(0, 1)
  )
  print('knn mar fit')
  pred <- predict(MAR_model_knn,
                  newdata.x = as.matrix(train_MAR[is.na(train_MAR$y), 5:ncol(train_MAR)]),
                  newdata.z = as.matrix(train_MAR[is.na(train_MAR$y), 2:4]))
  
  na_indices <- which(is.na(train_MAR$y))
  train_imputed <- train_MAR
  train_imputed$y[na_indices] <- pred
  
  Fitted_knn_Model <- sfplsim.kNN.fit(
    as.matrix(train_imputed[, 5:ncol(train_imputed)]),
    as.matrix(train_imputed[, 2:4]),
    train_imputed$y,
    nknot.theta = 3,
    nknot = 20,
    range.grid = c(0, 1)
  )
  
  print('knn imputed fit')
  
  # se between real and (full) estimated betas
  beta_knn_se <- sum((Fitted_knn_Model$beta.est - c(-1, .7, 2)) ^ 2)
  
  # mse between real and estimated MAR responses
  impValues_knn_mse <- sum((pred - train$y[na_indices])^2)/length(na_indices)
  
  ### Test predictions with MAR and Fitted models
  pred_knn_mar <- predict(
    MAR_model_knn,
    newdata.x = as.matrix(test[, 5:ncol(test)]),
    newdata.z = as.matrix(test[, 2:4]),
    y.test = test$y
  )
  pred_knn_imp <- predict(
    Fitted_knn_Model,
    newdata.x = as.matrix(test[, 5:ncol(test)]),
    newdata.z = as.matrix(test[, 2:4]),
    y.test = test$y
  )
  
  list(
    beta_kernel_se = beta_kernel_se,
    theta_kernel_mise = theta_kernel_mise,
    imputation_kernel_error = impValues_kernel_mse,
    kernelMAR_test_MSEP = pred_kernel_mar$MSEP.1,
    kernelIMP_test_MSEP = pred_kernel_imp$MSEP.1,
    beta_knn_se = beta_knn_se,
    imputation_knn_error = impValues_knn_mse,
    knnMAR_test_MSEP = pred_knn_mar$MSEP.1,
    knnIMP_test_MSEP = pred_knn_imp$MSEP.1
  )
}



# Metric generation -------------------------------------------------------

load('data_list.Rda')

metrics_df <- data.frame()

load('metrics_df.Rda')
### Change index for running in batches. Very computational expensive. Batches of 100 are recomended from 1:1600
for (i in 1401:1500){
  print(i)
  info <- fit_MAR(data_list[[i]])
  metrics_df <- rbind(
    metrics_df,
    data.frame(
      beta_kernel_se = info$beta_kernel_se,
      imputation_kernel_error = info$imputation_kernel_error,
      kernelMAR_test_MSEP = info$kernelMAR_test_MSEP,
      kernelIMP_test_MSEP = info$kernelIMP_test_MSEP,
      beta_knn_se = info$beta_knn_se,
      imputation_knn_error = info$imputation_knn_error,
      knnMAR_test_MSEP = info$knnMAR_test_MSEP,
      knnIMP_test_MSEP = info$knnIMP_test_MSEP,
      size = data_list[[i]]$dim,
      probability = data_list[[i]]$probability,
      cor = data_list[[i]]$cor,
      c = data_list[[i]]$c
    )
  )
}

save(metrics_df, file='metrics_df.Rda')


##################### errors #################################

# data sets (98,114) has no MAR. We fill with NA
null_vec <- rep(NA,14)
metrics_df <- rbind(metrics_df, null_vec)


# metric means by groups (100 iterations) ------------------------------------------------------------

# Convert size and probability to factors for better plotting
metrics_df$size <- factor(metrics_df$size)
metrics_df$probability <- factor(metrics_df$probability)
metrics_df$cor <- factor(metrics_df$cor)
metrics_df$c <- factor(metrics_df$c)

# Aggregate the data to calculate mean metrics for each combination of size and probability
aggregated_metrics <- metrics_df %>%
  group_by(size, probability, cor, c) %>%
  summarise(
    beta_kernel_se = mean(beta_kernel_se),
    imputation_kernel_error = mean(imputation_kernel_error),
    kernelMAR_test_MSEP = mean(kernelMAR_test_MSEP),
    kernelIMP_test_MSEP = mean(kernelIMP_test_MSEP),
    beta_knn_se = mean(beta_knn_se),
    imputation_knn_error = mean(imputation_knn_error),
    knnMAR_test_MSEP = mean(knnMAR_test_MSEP),
    knnIMP_test_MSEP = mean(knnIMP_test_MSEP)
  ) %>%
  ungroup()

save(aggregated_metrics, file='agg_metrics.Rda')


# Some visualizations -----------------------------------------------------

load('agg_metrics.Rda')
# Load necessary libraries
# install.packages("xtable")
library(xtable)


######### Tables #################
# Create a LaTeX table from the data frame
latex_table <- xtable(aggregated_metrics, digits = c(0, 0, 2, 1, 2, 6, 6, 6, 6, 6,6,6,6))

# Print the table to the console in LaTeX format
print(latex_table)


############### Graphs ####################

library(ggplot2)
library(patchwork)
library(tidyr)
library(stringr)

# Gather the metrics columns into key-value pairs
metrics_long <- aggregated_metrics %>%
  pivot_longer(cols = beta_kernel_se:knnIMP_test_MSEP, 
               names_to = "metric", 
               values_to = "value")

# Create separate columns for method and metric
metrics_long <- metrics_long %>%
  mutate(method = ifelse(grepl("kernel", metric), "kernel", "knn"),
         metric = gsub("kernel", "", metric),
         metric = gsub("knn", "", metric))

aggregated_metrics <- drop_na(aggregated_metrics)

##################### n comparison #########################

# Custom palette and labels
palette <- c('#c7522a',
             '#e5c185',
             '#fbf2c4',
             'bisque4',
             '#74a892',
             '#008585',
             '#004a4a',
             'cyan3')
custom_labels <- c(
  "cor =  0, alpha =  0.01, c = 0.01",
  "cor =  0.5, alpha =  0.01, c = 0.05",
  "cor =  0, alpha =  2, c = 0.01",
  "cor =  0.5, alpha =  2, c = 0.05",
  "cor =  0, alpha =  0.01, c = 0.01",
  "cor =  0.5, alpha =  0.01, c = 0.05",
  "cor =  0, alpha =  2, c = 0.01",
  "cor =  0.5, alpha =  2, c = 0.05"
)

my_labeller <- as_labeller(
  c(
    beta__se = "beta[0[MSE]]",
    imputation__error = "Imputation[MSE]",
    MAR_test_MSEP = "MAR_Model[MSEP]",
    IMP_test_MSEP = "Fitted_Model[MSEP]"
  ),
  default = label_parsed
)

# Plot
ggplot(
  metrics_long %>% filter(
    metric %in% c(
      "beta__se",
      "imputation__error",
      "MAR_test_MSEP",
      "IMP_test_MSEP"
    )
  ),
  aes(
    x = interaction(size, method),
    y = value,
    color = interaction(probability, cor, c)
  )
) +  labs(x = "Size.Method interaction", y = "Metric Value") +
  geom_point(size = 2) +
  facet_wrap( ~ metric, scales = "free", labeller = my_labeller) +
  theme_minimal() +
  scale_color_manual(labels = custom_labels, values = palette)

############################ Kernel Knn comparison ##################

library(dplyr)

# Compute the rate of change between Kernel and KNN estimations
comparison_table <- aggregated_metrics %>%
  mutate(
    rate_change_kernelMAR_knnMAR = (knnMAR_test_MSEP - kernelMAR_test_MSEP) / kernelMAR_test_MSEP,
    rate_change_kernelIMP_knnIMP = (knnIMP_test_MSEP - kernelIMP_test_MSEP) / kernelIMP_test_MSEP
  ) %>%
  select(
    size, probability, cor, c,
    kernelMAR_test_MSEP, knnMAR_test_MSEP, rate_change_kernelMAR_knnMAR,
    kernelIMP_test_MSEP, knnIMP_test_MSEP, rate_change_kernelIMP_knnIMP
  ) %>%
  mutate(Group = row_number()) %>%
  select(
    Group,kernelIMP_test_MSEP, knnIMP_test_MSEP, rate_change_kernelIMP_knnIMP
  )

# Print the table
print(comparison_table)


#################### alpha comparison ##########################


# Create the comparison table for KernelIMP_test_MSEP and KnnIMP_test_MSEP
comparison_table <- aggregated_metrics %>%
  select(size, probability, cor, c, kernelIMP_test_MSEP, knnIMP_test_MSEP) %>%
  pivot_wider(names_from = probability, values_from = c(kernelIMP_test_MSEP, knnIMP_test_MSEP)) %>%
  mutate(
    kernelIMP_rate_change = (`kernelIMP_test_MSEP_2` - `kernelIMP_test_MSEP_0.01`) / `kernelIMP_test_MSEP_0.01`,
    knnIMP_rate_change = (`knnIMP_test_MSEP_2` - `knnIMP_test_MSEP_0.01`) / `knnIMP_test_MSEP_0.01`
  ) %>%
  select(size, cor, c, kernelIMP_test_MSEP_0.01, kernelIMP_test_MSEP_2, kernelIMP_rate_change, knnIMP_test_MSEP_0.01, knnIMP_test_MSEP_2, knnIMP_rate_change)


# Create a LaTeX table from the data frame
latex_table <- xtable(comparison_table, digits = c(0, 2, 2, 6, 6, 6, 6, 6, 6,6))

# Print the table to the console in LaTeX format
print(latex_table)


######################### gamma comparation##########################

custom_labels <- c(
  "c =  0.01, alpha =  0.01, n = 50",
  "c =  0.05, alpha =  0.01, n = 50",
  "c =  0.01, alpha =  2, n = 50",
  "c =  0.05, alpha =  2, n = 50",
  "c =  0.01, alpha =  0.01, n = 200",
  "c =  0.05, alpha =  0.01, n = 200",
  "c =  0.01, alpha =  2, n = 200",
  "c =  0.05, alpha =  2, n = 200"
)
ggplot(metrics_long %>% filter(metric %in% c('beta__se', 'IMP_test_MSEP', 'imputation__error')), aes(x = interaction(cor,method), y = value, color =interaction(c,probability,size))) +
  geom_point(size=2) +
  facet_wrap(~metric, scales = "free", labeller = my_labeller) +
  labs(x = "Gamma.Method interaction", y = "Metric Value") +
  theme_minimal() + 
  scale_color_manual(labels=custom_labels, values = palette) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


######################### signal to noise comparation##########################
custom_labels <- c(
  "cor =  0, alpha =  0.01, n = 50",
  "cor =  0.5, alpha =  0.01, n = 50",
  "cor =  0, alpha =  2, n = 50",
  "cor =  0.5, alpha =  2, n = 50",
  "cor =  0, alpha =  0.01, n = 200",
  "cor =  0.5, alpha =  0.01, n = 200",
  "cor =  0, alpha =  2, n = 200",
  "cor =  0.5, alpha =  2, n = 200"
)

ggplot(metrics_long %>% filter(metric %in% c('beta__se', 'IMP_test_MSEP', 'imputation__error')), 
       aes(x = interaction(c, method), y = value, color = interaction(cor, probability, size))) +
  geom_point(size = 2) +
  facet_wrap(~metric, scales = "free", labeller = my_labeller) +
  labs(x = "c.method interaction", y = "Metric Value") +
  theme_minimal() +
  scale_color_manual(labels = custom_labels, values = palette) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



#################### MSEPs comparison ####################################

ggplot(metrics_long%>% filter(metric %in% c('MAR_test_MSEP', 'IMP_test_MSEP')), aes(x = metric, y = value, color = interaction(probability, cor, c))) +
  geom_line() +
  geom_point(size=2) +
  facet_wrap(~metric, scales = "free") +
  scale_color_manual(values = green_palette, labels = custom_labels) +
  labs(x = "Size", y = "Metric Value", color = "Groups") +
  theme_minimal()


# Compute the rate of change between Kernel and KNN estimations
comparison_table <- aggregated_metrics %>%
  mutate(
    rate_change_kernel = (kernelMAR_test_MSEP - kernelIMP_test_MSEP) / kernelIMP_test_MSEP,
    rate_change_knn = (knnMAR_test_MSEP - knnIMP_test_MSEP) / knnIMP_test_MSEP
  ) %>%
  select(
    size, probability, cor, c,
    kernelMAR_test_MSEP,kernelIMP_test_MSEP,rate_change_kernel,
    knnMAR_test_MSEP, knnIMP_test_MSEP, rate_change_knn
  )

# Create a LaTeX table from the data frame
latex_table <- xtable(comparison_table, digits = c(0, 2, 2, 6, 6, 6, 6, 6, 6,6,6))

# Print the table to the console in LaTeX format
print(latex_table)

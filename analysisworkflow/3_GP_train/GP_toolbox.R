#############################
# script GP_toolbox.R
#
# Contains helper functions for training a GP regression model based on a 
# database of OpenMalaria simulations, and for performing a sobol sensitivity 
# analysis
#
# Created 17.06.2019
# monica.golumbeanu@unibas.ch
#
# Adapted October 2021
# lydia.braunack-mayer@swisstph.ch
#
#############################

##############
### HEADER ###
##############

# Load required pacakges
library(tgp)
library(hetGP)
library(sensitivity)


###############################
### ACTIVE HELPER FUNCTIONS ###
###############################

rep.row <- function(x, n){
  return(matrix(rep(x, each = n), nrow = n))
}


rep.col <- function(x, n){
  matrix(rep(x, each = n), ncol = n, byrow = TRUE)
}


test_GP_plot <- function(cv_result, save) {
  
  # ----------------------------------------------------------
  # This function visualises accuracy of a fitted gaussian process regression model
  #
  # Inputs
  # cv_result: outputs of the function cv_train_matern
  # save: file path to save function outputs
  #
  # Required packages:
  # hetGP
  #
  # Outputs
  # Predicted vs. true plot for a fitted gaussian process regression model
  #
  # ----------------------------------------------------------
  
  # Set up
  test_data <- cv_result[["test_data"]]
  train_data <- cv_result[["train_data"]]
  
  # Generate model predictions
  test_predict <- predict(x = as.matrix(test_data[, 1:(ncol(test_data) - 1)]), object = cv_result[["GP_model"]])
  train_predict <- predict(x = as.matrix(train_data[, 1:(ncol(train_data) - 1)]), object = cv_result[["GP_model"]])
  
  # Calculate out of sample Mean Squared Error
  # mse <- sum((test_data[, ncol(test_data)] - test_predict$mean)^2) / nrow(test_data)
  # print(paste0("Out-of-Sample Mean Square Error: ", mse))
  
  # Generate predicted vs. true plot
  out_test <- data.frame("predicted" = test_predict$mean, "true" = test_data[, ncol(test_data)])
  out_train <- data.frame("predicted" = train_predict$mean, "true" = train_data[, ncol(train_data)])
  
  jpeg(save, width = 400, height = 400)
  plot(out_train,
       col = rgb(red = 70, green = 138, blue = 178, alpha = 100, maxColorValue = 255), #densCols(out_train),
       xlab = "Predicted value", ylab = "True value",
       main = paste0("Within-Sample R2 = ", round(cor(out_train[, "predicted"], out_train[, "true"]), 2),
                     ", Out-of-Sample R2 = ", round(cor(out_test[, "predicted"], out_test[, "true"]), 2)), 
       pch = 20)
  points(out_test,
         col = rgb(red = 191, green = 50, blue = 39, alpha = 200, maxColorValue = 255),
         pch = 20)
  legend("bottomright", 
         legend = c("Training data", "Testing data"), 
         col = c(rgb(red = 70, green = 138, blue = 178, maxColorValue = 255),
                 rgb(red = 191, green = 50, blue = 39, maxColorValue = 255)),
         pch = 20)
  dev.off()
    
}


train_GP_matern <- function(train_data, lower, upper) {
  
  # ----------------------------------------------------------
  # This function uses gaussian process regression to train an OpenMalaria emulator
  #
  # Inputs
  # train_data: data used to fit the gaussian process regression model
  # lower: lower input parameter bounds passed to mleHetGP() for mle optimisation
  # upper: upper input parameter bounds passed to mleHetGP() for mle optimisation
  #
  # Required packages:
  # hetGP
  #
  # Outputs
  # GP_model, a list containing all training data and the trained gaussian process regression model
  #
  # ----------------------------------------------------------
  
  # Set up function
  D <- ncol(train_data) - 1    
  param_col <- c(1:D)
  response_col <- D + 1
  
  # Debug
  if (is.null(lower) | is.null(upper)) {
      print("Error in train_GP_matern(): object 'lower' or 'upper' not found.")
  } else if (length(lower) != D | length(upper) != D) {
    print("Error in train_GP_matern(): the numer of bounds specified in object 'lower' or 'upper' is not the same as the number of parameters in the input data.")
  }
  
  # Format data as required for mleHetGP
  prdata <- find_reps(X = as.matrix(train_data[, param_col]), 
                     Z = as.matrix(train_data[, response_col]), 
                     rescale = FALSE, normalize = FALSE)
  
  # Fit GP regression model
   
  GP_model <- mleHetGP(X = list(X0 = as.matrix(prdata$X0), Z0 = as.matrix(prdata$Z0), mult = prdata$mult),
                       Z = prdata$Z, lower = lower, upper = upper, covtype = "Matern5_2")
  
  # Return function outputs
  return(GP_model)
  
}


cv_train_matern <- function(input_data, lower, upper, scale, test_prop = 0.2) {

  # ----------------------------------------------------------
  # This function splits data into test and train sets, and performs gaussian process regression on the train set
  #
  # Inputs
  # input_data: a data frame containing data on which the gaussian process regression will be trained
  # lower: argument passed to train_GP_matern, lower input parameter bounds passed to mleHetGP() for mle optimisation
  # upper: argument passed to train_GP_matern, upper input parameter bounds passed to mleHetGP() for mle optimisation
  # test_prop: the proportion of data set aside as test data, set by default to 10%
  #
  # Required packages:
  # hetGP
  #
  # Outputs
  # list containing 'train_data', the dataset used to perform gaussian process regression,
  # 'test_data', the dataset set aside as test data, and 
  # 'GP_model', the gaussian process regression model fitted to train_data
  #
  # ----------------------------------------------------------
  
  # Set up function
  D <- ncol(input_data) - 1    
  param_col <- c(1:D)
  response_col <- D + 1

  if (is.null(scale)) {
    print("Fitting gaussian process regression. Inputs used on their original scale.")
  } else {
    print("Fitting gaussian process regression. Inputs have been scaled.")
  }
  
  # Scale data according to bounds defined in 'scale'
  if (!is.null(scale)) {
    for (i in 1:D) {
      input_data[, i] <- (input_data[, i] - scale[i, 1]) / (scale[i, 2] - scale[i, 1])
    }
  }
  
  # Split input data into test and train sets
  test_indices <- sample(1:nrow(input_data), size = test_prop*nrow(input_data))
  test_data <- input_data[test_indices, ]
  
  train_indices <- setdiff(1:nrow(input_data), test_indices)
  train_data <- input_data[train_indices, ]
  
  # Train model on train set
  trained_model <- train_GP_matern(train_data, lower, upper)
  print(summary(trained_model))
  
  # Return function outputs
  return(list(train_data = train_data, test_data = test_data, GP_model = trained_model))
  
}


calc_sobol_idx <- function(GP_model, param_spec, num_points){
  
  # ----------------------------------------------------------
  # This function performs a sobol global sensitivity analysis based on an OpenMalaria emulator
  #
  # Inputs
  # GP_model: a model object, e.g. outputs of fitting gaussian process regression
  # param_spec: a matrix specifying ranges for input parameters for GP_model
  # num_points: number of samples passed to lhs()
  #
  # Required packages:
  # hetGP
  # sensitivity
  # tgp
  #
  # Outputs
  # list containing 'S_eff', a vector containing sobol first-order effects for each input parameter for GP_model, and 'T_eff', sobol total-order effects
  #
  # ----------------------------------------------------------
  
    # Set up function
    S_mat <- T_mat <- NULL
    print(param_spec)
    
    # Define wrapper for the GP_model prediction to be called by soboljansen()
    GP_f <- function(X){
        out <- predict(x = as.matrix(X), object = GP_model)
        return(out$mean)
    }
    
    # Construct the two random parameter samples
    X1 <- lhs(num_points, as.matrix(param_spec))
    X2 <- lhs(num_points, as.matrix(param_spec))
    
    # Compute the Sobol indices
    SA <- sensitivity::soboljansen(model = GP_f, as.data.frame(X1), as.data.frame(X2), nboot = 1000)
    S_eff <- SA$S$original
    S_eff[S_eff < 0] <- 0
    T_eff<- SA$T$original
    
    # Return function outputs
    #return(list(S_eff = S_eff, T_eff = T_eff, X1=as.data.frame(X1), X2 = as.data.frame(X2)))
    return(list(SA = SA))
    
}


#################################
### ARCHIVED HELPER FUNCTIONS ###
#################################

# plot_gradient = function(data_tab, plot_file = NULL, labels_names) {
#   a = interp.loess(data_tab[,1],data_tab[,2],data_tab[,3])
#   #image(a, col=terrain.colors(10))
#   d = as.data.frame(cbind(expand.grid(a$x, a$y), c(a$z)))
#   colnames(d) = c("X", "Y", "Z")
#   colfunc = colorRampPalette(c("#2c7fb8", "#7fcdbb", "#fec44f"))
#   ggplot() + theme_bw() + geom_tile(data = d, aes(x = X, y = Y, fill = Z)) +
#     xlab(labels_names[1]) + ylab(labels_names[2]) +
#     #scale_y_continuous(limits=c(0,1)) + scale_x_continuous(limits=c(0,10)) +
#     #scale_fill_gradientn(colours = viridis(50), name = labels_names[3]) + 
#     scale_fill_gradient(low = "darkgrey", high = "white", name = labels_names[3])+
#     theme(panel.grid.minor = element_blank(),
#           panel.grid.major = element_blank(), panel.background = element_blank()) 
# }


# # Function that trains a GP regression given a training set
# train_GP = function(train_data) {
#   D = ncol(train_data) - 1    
#   param_col = c(1:D)
#   response_col = D + 1
#   
#   print("Training GP regression model ...")
#   
#   prdata = find_reps(X = as.matrix(train_data[, param_col]), 
#                      Z = as.matrix(train_data[, response_col]), 
#                      rescale = FALSE, normalize = FALSE)
#   GP_model = mleHetGP(X = list(X0 = as.matrix(prdata$X0), 
#                                Z0 = as.matrix(prdata$Z0), mult = prdata$mult), 
#                       Z = prdata$Z, lower = rep(0.0001, D), upper = rep(10, D), 
#                       covtype = "Gaussian")
#   return(GP_model)
# }


# cv_train = function(input_data, K) {
#   # retrieve all data points indices
#   indices = c(1:nrow(input_data))
#   
#   # split the data in K groups
#   cv_indices = sample(rep(1:K, length.out = nrow(input_data)))
#   train_data = test_data = NULL
#   for (i in 1:K) {
#     print(paste("Performing CV run",i,"..."))
#     test_points = which(cv_indices == i)
#     train_points = setdiff(indices, test_points)
#     trained_model = train_GP(input_data[train_points,]) 
#     print(paste("Computing errors run",i,"..."))
#     data_tabs = test_GP(trained_model, input_data[train_points,], 
#                         input_data[test_points,])
#     train_data = rbind.data.frame(train_data, cbind.data.frame(i, data_tabs$train_data))
#     test_data = rbind.data.frame(test_data, cbind.data.frame(i, data_tabs$test_data))
#   }
#   # Train the final GP model on the entire training set
#   trained_model = train_GP(input_data) 
#   return(list(train_data = train_data, test_data = test_data, GP_model = trained_model))
# }



# test_GP = function(GP_model, train_data, test_data) {
#   D = ncol(train_data) - 1
#   param_col = c(1:D)
#   response_col = D + 1
#   # apply the model on the train data and calculate the error
#   prediction_train = predict(x = as.matrix(train_data[, param_col]), object = GP_model)
#   predicted = prediction_train$mean
#   train_data = cbind.data.frame(train_data, predicted)
#   
#   # apply the model on the train data and calculate the error
#   if(!is.null(test_data)) {
#     prediction_test = predict(x = as.matrix(test_data[, param_col]), object = GP_model)
#     predicted = prediction_test$mean
#     test_data = cbind.data.frame(test_data, predicted)
#   } else {
#     test_data = NULL
#   }
#   
#   return(list(train_data = train_data, test_data = test_data))
# }
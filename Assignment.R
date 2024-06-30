
# Load required libraries
library(matlib)
library(ggplot2)
library(rsample)

# Read the CSV data
data <- read.csv("C:/Users/Master Fluffy/Documents/gene_data.csv")


# Convert data to long format
data_long <- pivot_longer(data, cols = c(x1, x2, x3, x4, x5), 
                          names_to = "Signal", values_to = "Value")

# Create the time series plot
ggplot(data_long, aes(x = Time..min., y = Value, color = Signal)) +
  geom_line() +
  labs(title = "Time Series Plot of X Signals",
       x = "Time (min)",
       y = "Value") +
  theme_minimal() +
  scale_color_brewer(palette = "Set1")

# If you prefer separate plots for each signal:
ggplot(data_long, aes(x = Time..min., y = Value, color = Signal)) +
  geom_line() +
  facet_wrap(~Signal, scales = "free_y", ncol = 1) +
  labs(title = "Time Series Plots of X Signals",
       x = "Time (min)",
       y = "Value") +
  theme_minimal() +
  theme(legend.position = "none")

# Save the plots (optional)
ggsave("combined_time_series_plot.png", width = 10, height = 6)
ggsave("faceted_time_series_plots.png", width = 10, height = 12)


# Create the time series plot for X1 (treating it as Y series)
ggplot(data, aes(x = Time..min., y = x1)) +
  geom_line(color = "blue") +
  labs(title = "Time Series Plot of Y Signal",
       x = "Time (min)",
       y = "Y Value") +
  theme_minimal()

# Save the plot (optional)
ggsave("y_series_time_plot.png", width = 10, height = 6)


# Combine all X signals into a single vector
X_CSV_DATA <- c(data$x1, data$x2, data$x3, data$x4, data$x5)

# Create density plot of X signal
density_of_X <- density(X_CSV_DATA)

# Plot the density
plot(density_of_X, main = "Density plot of whole input signal", 
     xlab = "Value", ylab = "Density")

#=================================================================================
# 1. Density plot of whole input signal
density_of_X <- density(X_CSV_DATA)
plot(density_of_X, main = "Density plot of whole input signal", 
     xlab = "X signal", ylab = "Density")

# 2. Histogram with density overlay
hist(X_CSV_DATA, freq = FALSE, main = "Density and Histogram of X signal",
     xlab = "X signal", ylab = "Density")
lines(density_of_X, lwd = 2, col = "brown")
rug(jitter(X_CSV_DATA))
 
#===================================================================================  

# Create plots for each variable

  var_name <- names(data)[2]
  plot_hist_density(data[[2]], var_name)

  var_name <- names(data)[3]
  plot_hist_density(data[[3]], var_name)
  
  var_name <- names(data)[4]
  plot_hist_density(data[[4]], var_name)
  
  var_name <- names(data)[5]
  plot_hist_density(data[[5]], var_name)
  
  var_name <- names(data)[6]
  plot_hist_density(data[[6]], var_name)

#===================================================================================  
  
  # Plot the density Y
  # Extract time data
  time <- data[, 1]
  
  # Set up a 1x2 plot layout
  par(mfrow = c(1, 2))
  
  # Density plot for time
  plot_density <- function(x, var_name) {
    d <- density(x)
    plot(d, main = paste("Density Plot of", var_name),
         xlab = var_name, ylab = "Density",
         col = "blue", lwd = 2)
    polygon(d, col = rgb(0, 0, 1, 0.2), border = "blue")
    rug(jitter(x), col = rgb(0, 0, 0, 0.1))
  }
  
  plot_density(time, "Time (min)")


var_name <- names(data)[1]
plot_hist_density(data[[1]], var_name)

#===================================================================================  

# Task 1.3 Correlation and scatter plots
# to arrange plot in a single screen


# Load the data
gene_data <- read.csv("C:/Users/Master Fluffy/Documents/gene_data.csv")



# Fix the column name to remove the space
colnames(gene_data)[1] <- "Time_min"

# Create the scatter plot using ggplot2
library(ggplot2)

ggplot(gene_data, aes(x = x1, y = x2)) +
  geom_point() +
  labs(title = "Scatter Plot of x1 vs x2",
       x = "Gene x1 Expression", 
       y = "Gene x2 Expression")

ggplot(gene_data, aes(x = x1, y = x3)) +
  geom_point() +
  labs(title = "Scatter Plot of x1 vs x3",
       x = "Gene x1 Expression", 
       y = "Gene x3 Expression")

ggplot(gene_data, aes(x = x1, y = x4)) +
  geom_point() +
  labs(title = "Scatter Plot of x1 vs x4",
       x = "Gene x1 Expression", 
       y = "Gene x4 Expression")

ggplot(gene_data, aes(x = x1, y = x5)) +
  geom_point() +
  labs(title = "Scatter Plot of x1 vs x5",
       x = "Gene x1 Expression", 
       y = "Gene x5 Expression")

#==2

ggplot(gene_data, aes(x = x2, y = x3)) +
  geom_point() +
  labs(title = "Scatter Plot of x2 vs x3",
       x = "Gene x1 Expression", 
       y = "Gene x3 Expression")

ggplot(gene_data, aes(x = x2, y = x4)) +
  geom_point() +
  labs(title = "Scatter Plot of x2 vs x4",
       x = "Gene x2 Expression", 
       y = "Gene x4 Expression")


ggplot(gene_data, aes(x = x2, y = x5)) +
  geom_point() +
  labs(title = "Scatter Plot of x2 vs x5",
       x = "Gene x2 Expression", 
       y = "Gene x5 Expression")
#==3

ggplot(gene_data, aes(x = x3, y = x4)) +
  geom_point() +
  labs(title = "Scatter Plot of x3 vs x4",
       x = "Gene x3 Expression", 
       y = "Gene x4 Expression")

ggplot(gene_data, aes(x = x3, y = x5)) +
  geom_point() +
  labs(title = "Scatter Plot of x3 vs x5",
       x = "Gene x3 Expression", 
       y = "Gene x5 Expression")

#==4
ggplot(gene_data, aes(x = x4, y = x5)) +
  geom_point() +
  labs(title = "Scatter Plot of x4 vs x5",
       x = "Gene x4 Expression", 
       y = "Gene x5 Expression")


#===================================================================================  

#Task 2 begins

# Load necessary libraries
library(readxl)
library(ggplot2)
library(GGally)
library(cowplot)

# Load the Data

file_path <- "C:/Users/Master Fluffy/Documents/cleaned_dataset 1.xlsx"
data <- read_excel(file_path)

# Extracting time and gene expression data from the tibbles
time <- data$`Time (min)`
genes <- data[, -1] 
GenesData <- as.matrix(genes)


#Task 2.1

# Define input variables according to each model
x1 <- GenesData[, "x1"]
x3 <- GenesData[, "x3"]
x4 <- GenesData[, "x4"]
x5 <- GenesData[, "x5"]
y <- GenesData[, "x2"]

# Model 1: y = θ₁x₄ + θ₂x₃² + θbias
X_model1 <- cbind(x4, x3^2, rep(1, length(x4)))  # Design matrix including bias term
theta_model1 <- solve(t(X_model1) %*% X_model1) %*% t(X_model1) %*% y  # Least squares estimation

# Model 2: y = θ₁x₄ + θ₂x₃² + θ₃x₅ + θbias
X_model2 <- cbind(x4, x3^2, x5, rep(1, length(x4)))
theta_model2 <- solve(t(X_model2) %*% X_model2) %*% t(X_model2) %*% y

# Model 3: y = θ₁x₃ + θ₂x₄ + θ₃x₅³
X_model3 <- cbind(x3, x4, x5^3)
theta_model3 <- solve(t(X_model3) %*% X_model3) %*% t(X_model3) %*% y

# Model 4: y = θ₁x₄ + θ₂x₃² + θ₃x₅³ + θ bias
X_model4 <- cbind(x4, x3^2, x5^3, rep(1, length(x4)))
theta_model4 <- solve(t(X_model4) %*% X_model4) %*% t(X_model4) %*% y

# Model 5: y = θ₁x₄ + θ₂x₁² + θ₃x₃² + θ bias
X_model5 <- cbind(x4, x1^2, x3^2, rep(1, length(x4)))
theta_model5 <- solve(t(X_model5) %*% X_model5) %*% t(X_model5) %*% y

# Display estimated parameters
print("Model 1 Theta:")
print(theta_model1)
print("Model 2 Theta:")
print(theta_model2)
print("Model 3 Theta:")
print(theta_model3)
print("Model 4 Theta:")
print(theta_model4)
print("Model 5 Theta:")
print(theta_model5)


# Model 1: y = θ₁x₄ + θ₂x₃² + θbias
X_model1 <- cbind(1, x4, x3^2)  # Design matrix including bias term
theta_model1 <- solve(t(X_model1) %*% X_model1) %*% t(X_model1) %*% y  # Least squares estimation

# Model 2: y = θ₁x₄ + θ₂x₃² + θ₃x₅ + θbias
X_model2 <- cbind(1, x4, x3^2, x5)
theta_model2 <- solve(t(X_model2) %*% X_model2) %*% t(X_model2) %*% y

# Model 3: y = θ₁x₃ + θ₂x₄ + θ₃x₅³
X_model3 <- cbind(x3, x4, x5^3)
theta_model3 <- solve(t(X_model3) %*% X_model3) %*% t(X_model3) %*% y

# Model 4: y = θ₁x₄ + θ₂x₃² + θ₃x₅³ + θ bias
X_model4 <- cbind(1, x4, x3^2, x5^3)
theta_model4 <- solve(t(X_model4) %*% X_model4) %*% t(X_model4) %*% y

# Model 5: y = θ₁x₄ + θ₂x₁² + θ₃x₃² + θ bias
X_model5 <- cbind(1, x4, x1^2, x3^2)
theta_model5 <- solve(t(X_model5) %*% X_model5) %*% t(X_model5) %*% y

# Display estimated parameters
print("Model 1 Theta:")
print(theta_model1)
print("Model 2 Theta:")
print(theta_model2)
print("Model 3 Theta:")
print(theta_model3)
print("Model 4 Theta:")
print(theta_model4)
print("Model 5 Theta:")
print(theta_model5)


#### Task 2.1

#Define the matrix of ones (bias term)
onesMatrix <- matrix(1, nrow = length(y), ncol = 1)

# Model 1: y = θ₁x₄ + θ₂x₃² + θbias
X_model1 <- cbind(x4, x3^2, onesMatrix)
theta_model1 <- solve(t(X_model1) %*% X_model1) %*% t(X_model1) %*% y

# Model 2: y = θ₁x₄ + θ₂x₃² + θ₃x₅ + θbias
X_model2 <- cbind(x4, x3^2, x5, onesMatrix)
theta_model2 <- solve(t(X_model2) %*% X_model2) %*% t(X_model2) %*% y

# Model 3: y = θ₁x₃ + θ₂x₄ + θ₃x₅³
X_model3 <- cbind(x3, x4, x5^3)
theta_model3 <- solve(t(X_model3) %*% X_model3) %*% t(X_model3) %*% y

# Model 4: y = θ₁x₄ + θ₂x₃² + θ₃x₅³ + θ bias
X_model4 <- cbind(x4, x3^2, x5^3, onesMatrix)
theta_model4 <- solve(t(X_model4) %*% X_model4) %*% t(X_model4) %*% y

# Model 5: y = θ₁x₄ + θ₂x₁² + θ₃x₃² + θ bias
X_model5 <- cbind(onesMatrix,x4, x1^2, x3^2)
theta_model5 <- solve(t(X_model5) %*% X_model5) %*% t(X_model5) %*% y

# Display estimated parameters
print("Model 1 Theta:")
print(theta_model1)
print("Model 2 Theta:")
print(theta_model2)
print("Model 3 Theta:")
print(theta_model3)
print("Model 4 Theta:")
print(theta_model4)
print("Model 5 Theta:")
print(theta_model5)


### Task 2.2 RSS

# Compute Predicted values for each model
y_hat_model1 <- X_model1 %*% theta_model1
y_hat_model2 <- X_model2 %*% theta_model2
y_hat_model3 <- X_model3 %*% theta_model3
y_hat_model4 <- X_model4 %*% theta_model4
y_hat_model5 <- X_model5 %*% theta_model5

# Residuals for each model
residuals_model1 <- y - y_hat_model1
residuals_model2 <- y - y_hat_model2
residuals_model3 <- y - y_hat_model3
residuals_model4 <- y - y_hat_model4
residuals_model5 <- y - y_hat_model5

# Compute RSS for each model
RSS_model1 <- sum(residuals_model1^2)
RSS_model2 <- sum(residuals_model2^2)
RSS_model3 <- sum(residuals_model3^2)
RSS_model4 <- sum(residuals_model4^2)
RSS_model5 <- sum(residuals_model5^2)

# Display RSS for each model
print("RSS for Model 1:")
print(RSS_model1)
print("RSS for Model 2:")
print(RSS_model2)
print("RSS for Model 3:")
print(RSS_model3)
print("RSS for Model 4:")
print(RSS_model4)
print("RSS for Model 5:")
print(RSS_model5)



## Task 2.3
n <- length(y)
# Calculating variance
# Model 1
variance_model1 = RSS_model1/(n-1)
# Model 2
variance_model2 = RSS_model2/(n-1)
# Model 3
variance_model3 = RSS_model3/(n-1)
# Model 4
variance_model4 = RSS_model4/(n-1)
# Model 5
variance_model5 = RSS_model5/(n-1)

# Calculation log-likelihood for each model
log_likelihood_model1 <- -n/2 * log(2 * pi) - n/2 * log(variance_model1) - 1/(2 * variance_model1) * RSS_model1
log_likelihood_model2 <- -n/2 * log(2 * pi) - n/2 * log(variance_model2) - 1/(2 * variance_model2) * RSS_model2
log_likelihood_model3 <- -n/2 * log(2 * pi) - n/2 * log(variance_model3) - 1/(2 * variance_model3) * RSS_model3
log_likelihood_model4 <- -n/2 * log(2 * pi) - n/2 * log(variance_model4) - 1/(2 * variance_model4) * RSS_model4
log_likelihood_model5 <- -n/2 * log(2 * pi) - n/2 * log(variance_model5) - 1/(2 * variance_model5) * RSS_model5

# Display log-likelihood for each model
print("Log-Likelihood for Model 1:")
print(log_likelihood_model1)
print("Log-Likelihood for Model 2:")
print(log_likelihood_model2)
print("Log-Likelihood for Model 3:")
print(log_likelihood_model3)
print("Log-Likelihood for Model 4:")
print(log_likelihood_model4)
print("Log-Likelihood for Model 5:")
print(log_likelihood_model5)




## Task 2.4

# Number of parameters (k) for each model
k_model1 <- length(theta_model1)
k_model2 <- length(theta_model2)
k_model3 <- length(theta_model3)
k_model4 <- length(theta_model4)
k_model5 <- length(theta_model5)

# Compute AIC for each model
AIC_model1 <- 2 * k_model1 - 2 * log_likelihood_model1
AIC_model2 <- 2 * k_model2 - 2 * log_likelihood_model2
AIC_model3 <- 2 * k_model3 - 2 * log_likelihood_model3
AIC_model4 <- 2 * k_model4 - 2 * log_likelihood_model4
AIC_model5 <- 2 * k_model5 - 2 * log_likelihood_model5

# Compute BIC for each model
BIC_model1 <- k_model1 * log(n) - 2 * log_likelihood_model1
BIC_model2 <- k_model2 * log(n) - 2 * log_likelihood_model2
BIC_model3 <- k_model3 * log(n) - 2 * log_likelihood_model3
BIC_model4 <- k_model4 * log(n) - 2 * log_likelihood_model4
BIC_model5 <- k_model5 * log(n) - 2 * log_likelihood_model5

# Display AIC and BIC for each model
print("AIC and BIC for Model 1:")
print(paste("AIC:", AIC_model1))
print(paste("BIC:", BIC_model1))
print("AIC and BIC for Model 2:")
print(paste("AIC:", AIC_model2))
print(paste("BIC:", BIC_model2))
print("AIC and BIC for Model 3:")
print(paste("AIC:", AIC_model3))
print(paste("BIC:", BIC_model3))
print("AIC and BIC for Model 4:")
print(paste("AIC:", AIC_model4))
print(paste("BIC:", BIC_model4))
print("AIC and BIC for Model 5:")
print(paste("AIC:", AIC_model5))
print(paste("BIC:", BIC_model5))


## Task 2.5
#Plot Histograms of Residuals
par(mfrow = c(3, 2))  # Set up a 3x2 grid for plots
hist(residuals_model1, main = "Model 1 Residuals", xlab = "Residuals")
hist(residuals_model2, main = "Model 2 Residuals", xlab = "Residuals")
hist(residuals_model3, main = "Model 3 Residuals", xlab = "Residuals")
hist(residuals_model4, main = "Model 4 Residuals", xlab = "Residuals")
hist(residuals_model5, main = "Model 5 Residuals", xlab = "Residuals")
hist(residuals_model6, main = "Model 6 Residuals", xlab = "Residuals")

par(mfrow = c(1,1))  # Reset plotting layout

# Q-Q Plot for Residuals
par(mfrow = c(3, 2))
qqnorm(residuals_model1, col = "#336600",main = "Q-Q Plot: Model 1 Residuals")
qqline(residuals_model1)

qqnorm(residuals_model2, col = "#336600",main = "Q-Q Plot: Model 2 Residuals")
qqline(residuals_model2)

qqnorm(residuals_model3, col = "#336600",main = "Q-Q Plot: Model 3 Residuals")
qqline(residuals_model3)

qqnorm(residuals_model4, col = "#336600",main = "Q-Q Plot: Model 4 Residuals")
qqline(residuals_model4)

qqnorm(residuals_model5,col = "#336600", main = "Q-Q Plot: Model 5 Residuals")
qqline(residuals_model5)

qqnorm(residuals_model6, col = "#336600",main = "Q-Q Plot: Model 6 Residuals")
qqline(residuals_model6)
par(mfrow = c(1,1))


### Task 2.7

# Set seed for reproducibility
set.seed(123)

# Split data into training (70%) and testing (30%)
train_indices <- sample(1:nrow(GenesData), 0.7 * nrow(GenesData))
train_data <- GenesData[train_indices, ]
test_data <- GenesData[-train_indices, ]

# Define input variables for training and testing
x1_train <- train_data[, "x1"]
x3_train <- train_data[, "x3"]
x4_train <- train_data[, "x4"]
x5_train <- train_data[, "x5"]
y_train <- train_data[, "x2"]

x1_test <- test_data[, "x1"]
x3_test <- test_data[, "x3"]
x4_test <- test_data[, "x4"]
x5_test <- test_data[, "x5"]
y_test <- test_data[, "x2"]

# Define design matrix for Model 5 using training data
X_train <- cbind(x4_train, x1_train^2, x3_train^2, rep(1, length(x4_train)))

# Estimate model parameters using least squares method
theta_model5 <- solve(t(X_train) %*% X_train) %*% t(X_train) %*% y_train
print(theta_model5)
# Define design matrix for testing data
X_test <- cbind(x4_test, x1_test^2, x3_test^2, rep(1, length(x4_test)))

# Compute model predictions on testing data
y_pred <- X_test %*% theta_model5

# Compute residuals
residuals_test <- y_test - y_pred

# Compute RSS (Residual Sum of Squares)
RSS_Testing_Data <- sum(residuals_test^2)

# Compute standard error of residuals
se <- sqrt(sum(residuals_test^2) / (length(y_test) - length(theta_model5)))
alpha <- 0.05  # Significance level for 95% confidence interval
t_value <- qt(1 - alpha/2, df = length(y_test) - length(theta_model5))  # t-value for alpha/2 and degrees of freedom

# Compute lower and upper bounds of confidence interval
lower_bound <- y_pred - t_value * se
upper_bound <- y_pred + t_value * se

# Combine data into a data frame
plot_data <- data.frame(x4_test = x4_test, y_test = y_test, y_pred = y_pred, lower_bound = lower_bound, upper_bound = upper_bound)

# Plot with ggplot2
ggplot(plot_data, aes(x = x4_test)) +
  geom_point(aes(y = y_test, color = "blue"), size = 3, alpha = 0.6, show.legend = TRUE) +
  geom_point(aes(y = y_pred, color = "red"), size = 3, show.legend = TRUE) +
  geom_errorbar(aes(ymin = lower_bound, ymax = upper_bound, color = "green"), width = 0.1, size = 0.3, show.legend = TRUE) +
  labs(x = "x4 (Testing Data)", y = "y (Model Prediction)", 
       title = "Model Prediction with 95% Confidence Intervals") +
  scale_color_manual(values = c("blue" = "blue", "red" = "red", "green" = "green"),
                     labels = c("Testing Data", "Model Prediction", "95% CI")) +
  theme_minimal()

# Display estimated parameters and RSS
print("Estimated Parameters for Model 5:")
print(theta_model5)
print("RSS for Model 5:")
print(RSS_Testing_Data)

# Additional Plots

# Distribution of output signal (training data) with density
ggplot(train_data, aes(x = x2)) +
  geom_histogram(aes(y = ..density..), bins = 20, fill = "blue", alpha = 0.5) +
  geom_density(col = "red", lwd = 1.5) +
  labs(title = "Distribution of Output Signal (Training Data)", x = "Gene x2", y = "Density") +
  theme_minimal()

# Error distribution (residuals) with density
ggplot(data.frame(residuals_test), aes(x = residuals_test)) +
  geom_histogram(aes(y = ..density..), bins = 20, fill = "blue", alpha = 0.5) +
  geom_density(col = "red", lwd = 1.5) +
  labs(title = "Error Distribution (Residuals) - Testing Data", x = "Residuals", y = "Density") +
  theme_minimal()

# Q-Q plot for residuals
qqnorm(residuals_test)
qqline(residuals_test, col = "red", lwd = 2)
title("Q-Q Plot of Residuals - Testing Data")
par(mfrow = c(1,1))

# Prepare data for plotting
plot_data <- data.frame(x4_test = x4_test, y_test = y_test, y_pred = y_pred, lower_bound = lower_bound, upper_bound = upper_bound)

# Plot using ggplot2
library(ggplot2)

ggplot(plot_data, aes(x = x4_test)) +
  geom_point(aes(y = y_test), color = "blue", size = 3) +  # Actual values
  geom_point(aes(y = y_pred), color = "red", size = 3) +   # Predicted values
  geom_errorbar(aes(ymin = lower_bound, ymax = upper_bound), width = 0.2, color = "green") +  # Error bars for CI
  labs(x = "x4 (Testing Data)", y = "y (Model Prediction)", 
       title = "Model 5 Prediction with 95% Confidence Intervals") +
  theme_minimal()

## Task 3
RSS <- 8.62241676218405  
# Step 1: Identify the two parameters with the largest absolute values
abs_theta <- abs(theta_model5)
largest_params_indices <- order(abs_theta, decreasing = TRUE)[1:2]
fixed_params_indices <- setdiff(1:length(theta_model5), largest_params_indices)

# Extract the indices of the largest parameters
index1 <- largest_params_indices[1]
index2 <- largest_params_indices[2]

# Estimated values of the largest parameters
theta1_hat <- theta_model5[index1]
theta2_hat <- theta_model5[index2]

# Step 2: Define Uniform priors around the estimated values
range_factor <- 0.5  # Define the range factor for the Uniform distribution
theta1_range <- c(theta1_hat - range_factor * abs(theta1_hat), theta1_hat + range_factor * abs(theta1_hat))
theta2_range <- c(theta2_hat - range_factor * abs(theta2_hat), theta2_hat + range_factor * abs(theta2_hat))

# Step 3: Perform rejection ABC
n_samples <- 1000  # Number of samples to draw from the prior
epsilon <- RSS* 2  # Tolerance level for acceptance

# Draw samples from the prior distributions
theta1_samples <- runif(n_samples, min = theta1_range[1], max = theta1_range[2])
theta2_samples <- runif(n_samples, min = theta2_range[1], max = theta2_range[2])

# Function to compute RSS given parameter values
compute_rss <- function(theta1, theta2) {
  theta <- theta_model5
  theta[index1] <- theta1
  theta[index2] <- theta2
  y_pred <- X_train %*% theta
  rss <- sum((y_train - y_pred)^2)
  return(rss)
}

# Perform rejection sampling
accepted_samples <- data.frame(theta1 = numeric(), theta2 = numeric())
rejected_samples <- data.frame(theta1 = numeric(), theta2 = numeric())

for (i in 1:n_samples) {
  theta1 <- theta1_samples[i]
  theta2 <- theta2_samples[i]
  rss <- compute_rss(theta1, theta2)
  if (rss < RSS + epsilon) {
    accepted_samples <- rbind(accepted_samples, data.frame(theta1 = theta1, theta2 = theta2))
  } else {
    rejected_samples <- rbind(rejected_samples, data.frame(theta1 = theta1, theta2 = theta2))
  }
}

# Joint posterior distribution plot
plot_joint <- ggplot() +
  geom_point(data = accepted_samples, aes(x = theta1, y = theta2), color = "#3366ff", alpha = 0.5, size = 2) +
  labs(title = "Joint Posterior Distribution of theta1 and theta2",
       x = "Theta1", y = "Theta2") +
  theme_minimal()

# Marginal posterior distribution of theta1
plot_theta1 <- ggplot(accepted_samples, aes(x = theta1)) +
  geom_histogram(aes(y = ..density..), bins = 30, color = "black", alpha = 0.5) +
  geom_density(size = 1) +
  labs(title = "Marginal Posterior Distribution of theta1", x = "Theta1", y = "Density") +
  theme_minimal()

# Marginal posterior distribution of theta2
plot_theta2 <- ggplot(accepted_samples, aes(x = theta2)) +
  geom_histogram(aes(y = ..density..), bins = 20, color = "black", alpha = 0.5) +
  geom_density(size = 1) +
  labs(title = "Marginal Posterior Distribution of theta2", x = "Theta2", y = "Density") +
  theme_minimal()

# Combine plots using cowplot for a cleaner presentation (optional)
library(cowplot)
plot_grid(plot_joint, plot_theta1, plot_theta2, ncol = 1)

library(cowplot)

# Step 4: Plot the joint and marginal posterior distributions
# Joint posterior distribution plot
ggplot(accepted_samples, aes(x = theta1, y = theta2)) +
  geom_point(alpha = 0.5) +
  labs(title = "Joint Posterior Distribution of theta1 and theta2",
       x = "theta1", y = "theta2") +
  theme_minimal()

# Marginal posterior distributions
ggplot(accepted_samples, aes(x = theta1)) +
  geom_histogram(aes(y = ..density..), bins = 30, color = "black", alpha = 0.5) +
  geom_density( size = 1) +
  labs(title = "Marginal Posterior Distribution of theta1", x = "theta1", y = "Density") +
  theme_minimal()

ggplot(accepted_samples, aes(x = theta2)) +
  geom_histogram(aes(y = ..density..), bins = 20,color = "black", alpha = 0.5) +
  geom_density( size = 1) +
  labs(title = "Marginal Posterior Distribution of theta2", x = "theta2", y = "Density") +
  theme_minimal()
par(mfrow = c(1,1))

#second way
thetabias = 1.36054805
thetaA = 0.89403436
thetaB = 0.52925488
thetaC = 1.36054805

array1 = 0
array2 = 0
f_value = 0
s_value = 0

number = 1000
counter = 0
for (i in 1:number){
  range11 = runif(1,-1.36054805,1.36054805)
  range22 = runif(1, -0.89403436, 0.89403436)
  newThetahat = matrix(c(range11,range22, thetaB,thetaC))
  newYHat <- X_model5 %*% newThetahat
  newRSS = sum((y- newYHat)^2)
  if(newRSS > epsilon){
    array1[i] = range11
    array2[i] = range22
    counter = counter +1
    Fvalue = matrix(array1)
    Svalue = matrix(array2)
  }
}
plot(Fvalue, Svalue, col = c("#ff1a1a", "#3366ff"), main = "Joint and Marginal Posterior Distribution Model 5")
par(mfrow = c(1,1))



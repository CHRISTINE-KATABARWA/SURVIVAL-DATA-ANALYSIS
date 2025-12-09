data<-read.csv("C:\\Users\\chris\\OneDrive\\Documents\\BIOSTAT 2\\SEM 1\\sda\\Homework\\HW1\\hodgkin2b.csv")
View(data)
## exploring the data##
#Age variable
avg_age<-mean(data$age)
avg_age
min_age<-min(data$age)
max_age<-max(data$age)

# age categories
# Calculating quantiles for 4 groups (quartiles)
quantiles <- quantile(data$age, probs = seq(0, 1, by = 0.25))  # Adjust `by` for more/less groups
print(quantiles)  # View the boundaries of each quantile

# Using the calculated quantiles to group the ages
data$age_group <- cut(data$age, 
                        breaks = quantiles,  # calculated quantile boundaries
                        include.lowest = TRUE,  # to include the lowest value
                        labels = c("Group 1", "Group 2", "Group 3", "Group 4"))  # Adjusting labels as needed
table(data$age_group)  # Checking how many people are in each age group

#age distribution in the groups
age_range <- tapply(data$age, data$age_group, range)
print(age_range)

#Kaplan Meier survival curves for the different age groups
library(survival)
library(survminer)
# Create age groups
data$age_group1 <- cut(data$age, breaks=c(0, 40, 60, 100), labels=c("Young", "Middle", "Old"))

# Fit Kaplan-Meier model based on age group
km_fit <- survfit(Surv(firsttime, event) ~ age_group1, data = data)

# Print summary of the Kaplan-Meier fit
summary(km_fit)


# Plot KM curves by age group
plot(km_fit, xlab = "Time", ylab = "Survival Probability", col = c("red", "blue", "green"))
legend("topright", legend = c("Young", "Middle-aged", "Old"), col = c("red", "blue", "green"), lty = 1)

############################
km_fit <- survfit(Surv(firsttime, event) ~ age_group, data = data)

# Print summary of the Kaplan-Meier fit
summary(km_fit)


# Plot KM curves by age group
plot(km_fit, xlab = "Time", ylab = "Survival Probability", col = c("red", "blue", "green","purple"))
legend("bottomleft", legend = c("g1", "g2", "g3","g4"), col = c("red", "blue", "green","purple"), lty = 1)

##################################################################
# Load required libraries
library(dplyr)

# Sort the data by time (assuming firsttime = survival time, event = event status)
data <- data %>% arrange(firsttime)

# Initialize variables
n <- nrow(data) # Total number of individuals
km_estimate <- numeric(n) # To store survival estimates
km_estimate[1] <- 1 # The first survival estimate starts at 1 (100%)

# Loop through time points to calculate KM estimates
for (i in 2:n) {
  # Calculate the survival probability at each time point
  risk_set_size <- n - (i - 1) # Number of people at risk
  num_events <- sum(data$event[data$firsttime == data$firsttime[i]]) # Number of events at that time
  km_estimate[i] <- km_estimate[i - 1] * (1 - (num_events / risk_set_size))
}

# Add the KM estimates to the dataset
data$km_estimate <- km_estimate

# Base R plotting for Kaplan-Meier
plot(data$firsttime, data$km_estimate, type = "s", 
     xlab = "Time", ylab = "Survival Probability",
     main = "Manually Computed Kaplan-Meier Survival Curve",
     lwd = 2, col = "blue")
#####################################
# Sort the data by time
data <- data %>% arrange(firsttime)

# Initialize Kaplan-Meier estimates for both groups (no metastasis and metastasis)
km_estimate_no_met <- 1
km_estimate_met <- 1

# Loop through survival times for each group and calculate the KM estimate manually
for (i in 2:nrow(data)) {
  if (data$cmt[i] == 0) {  # No metastasis group
    risk_set_size <- sum(data$cmt == 0) - (i - 1)
    num_events <- sum(data$event[data$firsttime == data$firsttime[i] & data$cmt == 0])
    km_estimate_no_met[i] <- km_estimate_no_met[i - 1] * (1 - (num_events / risk_set_size))
  } else {  # Metastasis group
    risk_set_size <- sum(data$cmt == 1) - (i - 1)
    num_events <- sum(data$event[data$firsttime == data$firsttime[i] & data$cmt == 1])
    km_estimate_met[i] <- km_estimate_met[i - 1] * (1 - (num_events / risk_set_size))
  }
}

# Add KM estimates to the dataset
data$km_no_met <- ifelse(data$cmt == 0, km_estimate_no_met, NA)
data$km_met <- ifelse(data$cmt == 1, km_estimate_met, NA)


# Plot manually computed Kaplan-Meier curves for both groups
ggplot() +
  geom_step(data = data[!is.na(data$km_no_met), ], aes(x = firsttime, y = km_no_met), color = "blue") +
  geom_step(data = data[!is.na(data$km_met), ], aes(x = firsttime, y = km_met), color = "red") +
  labs(title = "Kaplan-Meier Survival Curves (Manual) by Metastasis Status", 
       x = "Time", y = "Survival Probability") +
  theme_minimal() +
  scale_color_manual(values = c("blue", "red"), labels = c("No Metastasis", "Metastasis"))



#######################################################################################
# Install and load flexsurv package if you haven't already
# install.packages("flexsurv")
library(flexsurv)

# Load your data (assuming it is named 'hodgkin' in R)
# Convert categorical variables to factors if necessary
data$male <- factor(data$male, levels = c(0, 1))
data$cmt <- factor(data$cmt, levels = c(0, 1))
data$mediast <- factor(data$mediast, levels = c(0, 1, 2))
data$nodes <- factor(data$nodes, levels = c(0, 1))
data$clinstg <- factor(data$clinstg, levels = c(1, 2))  # Adjust levels as appropriate

# Fit the AFT model with a generalized gamma distribution
aft_model <- flexsurvreg(Surv(firsttime, event) ~ age + male + cmt + nodes + clinstg + mediast,
                         data = data,
                         dist = "gengamma",y=TRUE)  # "gengamma" is generalized gamma
aft_model

######Cox residuals
res.GG1 <-resid(aft_model,type="coxsnell")
# Standardizing the Cox-Snell residuals
coxsnell_mean <- mean(res.GG1)
coxsnell_sd <- sd(res.GG1)

# Standardized residuals
std_res <- (res.GG1 - coxsnell_mean) / coxsnell_sd


# Step 3: Get predicted survival probabilities using predict() for flexsurvreg model
surv_probs <- predict(aft_model, type = "response")

# Step 4: Plot standardized residuals vs survival probabilities
plot(surv_probs, std_res, 
     xlab = "Survival Probability", 
     ylab = "Standardized Residuals", 
     main = "Standardized Residuals vs Survival Probability", 
     pch = 19, col = "blue")


survplot(npsurv(res.GG1 ~1),conf ="none",ylab="Survival probability", xlab="Residual")
cs <- coxsnell_flexsurvreg(aft_model)

head(cs)
surv <- survfit(Surv(cs$est, data$event) ~ 1)
plot(surv, fun="cumhaz")
abline(0, 1, col="red")

length(cs$est)
length(data$event)

################################################################################"""
  surv_probs <- predict(aft_model, type = "survival")

# The 'surv_probs' might be a data frame, so extract the survival probabilities (column)
surv_probs_values <- surv_probs$surv  # Assuming 'surv' contains the survival probabilities

# Compute the Cox-Snell residuals: -log(survival probabilities)
cox_snell <- -log(surv_probs_values)

# Check the structure of the prediction output
surv_probs <- predict(aft_model, type = "survival")

# Check the structure of the prediction result
str(surv_probs)
head(surv_probs)
# Unnest the .pred column to get the survival probabilities
library(tidyr)

# Unnest the predictions
surv_probs_unnested <- unnest(surv_probs, .pred)

# Inspect the unnested data
head(surv_probs_unnested)
# Extract survival probabilities from the unnested data
surv_probs_values <- surv_probs_unnested$.pred_survival  # Assuming 'surv' is the column name for survival probabilities

# Compute Cox-Snell residuals
cox_snell <- -log(surv_probs_values)


# Fit Kaplan-Meier survival curve for comparison

# Step 1: Recalculate Cox-Snell residuals for individual observations
cox_snell <- residuals(aft_model, type = "coxsnell")
length(cox_snell)  # Length should be 865

# Step 2: Extract survival probabilities for individual observations
surv_probs_individual <- predict(aft_model, type = "survival", newdata = data)
length(surv_probs_individual)  # Length should be 865

# Step 3: Plot Cox-Snell residuals vs -ln(K-M Survival)
if (length(cox_snell) == length(surv_probs_individual)) {
  plot(-log(surv_probs_individual), cox_snell,
       xlab = "-ln(K-M Survival)",
       ylab = "Cox-Snell Residuals",
       main = "Cox-Snell Residuals vs Kaplan-Meier Survival",
       pch = 19)
} else {
  cat("Lengths do not match, check your data extraction.")
}

# Inspect the structure of the survival probabilities
str(surv_probs_individual)
library(dplyr)
library(tidyr)

# Unnest the .pred list-column
unnested_surv_probs <- surv_probs_individual %>%
  unnest(.pred) # This will unnest the list column to a long format

# Add the Cox-Snell residuals
unnested_surv_probs <- unnested_surv_probs %>%
  mutate(cox_snell = cox_snell) # Add the Cox-Snell residuals

# Now you can plot the predicted survival probabilities against the Cox-Snell residuals
library(ggplot2)

ggplot(unnested_surv_probs, aes(x = cox_snell, y = .pred_survival)) +
  geom_point() +
  labs(x = "Cox-Snell Residuals", y = "Predicted Survival Probability") +
  theme_minimal()


##################""
# Assuming cox_snell has a corresponding ID or identifier to join on
# If cox_snell is a separate vector, make sure it aligns by adding an ID to both dataframes

# Add an ID to the unnested data for merging
unnested_surv_probs <- surv_probs_individual %>%
  unnest(.pred) %>%
  mutate(id = row_number()) # Create an ID for each row

# Assuming cox_snell is a vector with the same length as the original data (or can be mapped to the data)
cox_snell_data <- tibble(id = 1:length(cox_snell), cox_snell = cox_snell)

# Join the data on 'id'
unnested_surv_probs <- unnested_surv_probs %>%
  left_join(cox_snell_data, by = "id")

# Plot the predicted survival probabilities against the Cox-Snell residuals
ggplot(-log(unnested_surv_probs), aes(x = cox_snell, y = .pred_survival)) +
  geom_point() +
  labs(x = "Cox-Snell Residuals", y = "Predicted Survival Probability") +
  theme_minimal()

# Adjust margins to reduce whitespace
par(mar=c(4, 4, 2, 2))  # Set margins: bottom, left, top, right

# Kaplan-Meier model
km_fit <- survfit(Surv(firsttime, event) ~ 1, data = data)



# Get standardized residuals from the AFT model
std_residuals <- residuals(aft_model, type = "response")
length(std_residuals)
# Get the fitted survival probabilities from the AFT model
surv_probs <- summary(aft_model, type = "survival", t = data$firsttime)

# Extract the survival probabilities from the summary as a vector
fitted_probs <- sapply(surv_probs, function(x) x$est)
# Check lengths
length(std_residuals)  # Should match the length of fitted_probs
length(fitted_probs)


# Step 4: Plot standardized residuals vs survival probabilities
plot(std_residuals, fitted_probs, 
     xlab = "Standardized Residuals", 
     ylab = "Predicted Survival Probability", 
     main = "Standardized Residuals vs Predicted Survival Probability", 
     col = "blue", pch = 19)
# Step 3: Generate predicted survival probabilities at each unique `firsttime`
unique_times <- sort(unique(data$firsttime))
surv_probs_summary <- summary(aft_model, type = "survival", t = unique_times)

# Extract survival probabilities from the summary as a vector
fitted_probs <- sapply(surv_probs_summary, function(x) x$est)

# Step 4: Map each `firsttime` in data to the corresponding probability
fitted_probs_by_time <- setNames(fitted_probs, unique_times)
data$fitted_prob <- fitted_probs_by_time[as.character(data$firsttime)]

# Now we should have `std_residuals` and `data$fitted_prob` with matching lengths

# Step 5: Plot standardized residuals vs survival probabilities
plot(std_residuals, data$fitted_prob, 
     xlab = "Standardized Residuals", 
     ylab = "Predicted Survival Probability", 
     main = "Standardized Residuals vs Predicted Survival Probability", 
     col = "blue", pch = 19)

# Step 6: Overlay Kaplan-Meier curve and AFT model's fitted survival curve
plot(km_fit, xlab = "Time", ylab = "Survival Probability", 
     main = "Kaplan-Meier vs Fitted Survival (AFT Model)", 
     col = "black", lwd = 2)

# Adding the AFT model's predicted survival as a line
lines(unique_times, fitted_probs, col = "red", lwd = 2)

# Add legend to distinguish curves
legend("bottomleft", legend = c("Kaplan-Meier", "AFT Model"), col = c("black", "red"), lty = 1, lwd = 2)

# Step 5: Plot standardized residuals vs survival probabilities
plot(std_residuals, data$fitted_prob, 
     xlab = "Standardized Residuals", 
     ylab = "Predicted Survival Probability", 
     main = "Standardized Residuals vs Predicted Survival Probability", 
     col = "blue", pch = 19)

# Add a loess smoothed line
loess_fit <- loess(data$fitted_prob ~ std_residuals)
lines(std_residuals, predict(loess_fit), col = "red", lwd = 2)


library(ggplot2)

# Create a data frame with the residuals and survival probabilities
plot_data <- data.frame(std_residuals = std_residuals, fitted_prob = data$fitted_prob)

# Generate the plot with a smooth line
ggplot(plot_data, aes(x = std_residuals, y = fitted_prob)) +
  geom_point(color = "blue") +
  geom_smooth(method = "loess", color = "red", se = FALSE) +
  labs(x = "Standardized Residuals", 
       y = "Predicted Survival Probability", 
       title = "Standardized Residuals vs Predicted Survival Probability")

library(survival)

# Fit the Kaplan-Meier survival curve
km_fit <- survfit(Surv(firsttime, event) ~ 1, data = data)

# Create a data frame with time and survival probabilities from Kaplan-Meier
km_data <- data.frame(time = km_fit$time, km_surv_prob = km_fit$surv)

# Assuming `data$fitted_prob` contains the predicted survival probabilities
plot_data <- data.frame(std_residuals = std_residuals, fitted_prob = data$fitted_prob)

# Base plot with standardized residuals vs predicted survival probabilities
plot(std_residuals, data$fitted_prob, 
     xlab = "Standardized Residuals", 
     ylab = "Survival Probability", 
     main = "Standardized Residuals vs Survival Probability",
     col = "blue", pch = 19)

# Add a smoothed line to represent the fitted AFT model probabilities
loess_fit <- loess(fitted_prob ~ std_residuals, data = plot_data)
lines(plot_data$std_residuals, predict(loess_fit), col = "red", lwd = 2)

# Add Kaplan-Meier survival probability line
lines(km_data$time, km_data$km_surv_prob, col = "black", lwd = 2, lty = 2)
legend("topright", legend = c("Predicted AFT Model", "Kaplan-Meier"), 
       col = c("red", "black"), lty = c(1, 2), lwd = 2)

library(ggplot2)

# Combine all data into one data frame for ggplot2
plot_data <- data.frame(std_residuals = std_residuals, 
                        fitted_prob = data$fitted_prob,
                        km_time = km_data$time[1:length(std_residuals)], 
                        km_surv_prob = km_data$km_surv_prob[1:length(std_residuals)])
plot_data
# Generate the overlay plot
ggplot(plot_data) +
  geom_point(aes(x = std_residuals, y = fitted_prob), color = "blue") +
  geom_smooth(aes(x = std_residuals, y = fitted_prob), method = "loess", color = "red", se = FALSE) +
  geom_line(aes(x = std_residuals, y = km_surv_prob), color = "black", linetype = "dashed", size = 1) +
  labs(x = "Standardized Residuals", 
       y = "Survival Probability", 
       title = "Standardized Residuals vs Survival Probability") +
  scale_color_manual(name = "Legend", values = c("Predicted AFT Model" = "red", "Kaplan-Meier" = "black")) +
  theme_minimal()

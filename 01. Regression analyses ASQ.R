### ASQ REGRESSION ANALYSIS ###

# Load libraries
library(readxl)
library(tidyr)
library(dplyr)
library(Hmisc)
library(ggplot2)
library(broom)
library(ggpubr)
library(quantreg)
library(mgcv)
library(lmtest)
library(car)
library(AER)
library(mice)
library(kableExtra)


# Load data
data <- read_excel("PATHWAY/ASQ data 6 month.xlsx")

str(data) # check variables 
data <- data %>% 
  mutate(across(c(comm_total, gm_total, fm_total, ps_total, perso_total), as.numeric)) # total score variables need to be numeric

# Summary of SARS-CoV-2 infections in pregnancy
data$preg_pos[is.na(data$preg_pos)] <- "No"
preg_pos <- data %>%
  group_by(preg_pos) %>%
  summarise(Count = n()) %>%
  mutate(Percentage = (Count / sum(Count)) * 100)
print(preg_pos)


# Summary of ASQ scores 
asq <- data %>%
  group_by(preg_pos) %>%
  summarise(across(c(comm_total, gm_total, fm_total, ps_total, perso_total),  #### PLEASE NOTE: One participant (ID=4223) is missing Fine Motor Skill score
                   list(
                     mean = ~ sprintf("%.1f", mean(.x, na.rm = TRUE)),
                     sd = ~ sprintf("%.1f", sd(.x, na.rm = TRUE))   
  )
  ),
  .groups = 'drop' 
  )
#print(asq)

# Create a table with the means and SDs
formatted_table_mean <- asq %>%
  kable("html",  col.names = c("SARS-CoV-2 positive", "Mean", "SD", "Mean", "SD","Mean", "SD","Mean", "SD","Mean", "SD"), align = "rccc") %>%
  kable_styling("hover", full_width = F, font_size = 13, html_font = 'helvetica') %>%
  column_spec(1, bold = T) %>%
  row_spec(0, bold = TRUE) %>%
  add_header_above(c(" ", "Communication" = 2, "Gross Motor" = 2, "Fine Motor" = 2, "Problem Solving" = 2, "Personal Social" = 2))
formatted_table_mean

# Format table
formatted_table_mean <- asq %>%
  kable("html", col.names = c("SARS-CoV-2 positive", "Mean", "SD", "Mean", "SD", "Mean", "SD", "Mean", "SD", "Mean", "SD"), align = "rccc") %>%
  kable_styling("hover", full_width = F, font_size = 13, html_font = 'helvetica') %>%
  column_spec(1, bold = T) %>%
  row_spec(0, bold = TRUE) %>%
  add_header_above(c(" ", "Communication" = 2, "Gross Motor" = 2, "Fine Motor" = 2, "Problem Solving" = 2, "Personal Social" = 2))

# Make sure appropriate variables are factors/dates
data$comm <-  factor(data$comm, order=T, levels= c("Typical", "Clinical"))
data$gm <-  factor(data$gm, order=T, levels= c("Typical", "Clinical"))
data$fm <-  factor(data$fm, order=T, levels= c("Typical", "Clinical"))
data$ps <-  factor(data$ps, order=T, levels= c("Typical", "Clinical"))
data$perso <-  factor(data$perso, order=T, levels= c("Typical", "Clinical"))
data$raceethnicitycombined = factor(data$raceethnicitycombined)
data$edu = factor(data$edu)
data$preg_pos = factor(data$preg_pos, levels = c("No", "Yes"))
data$min_pos_timing_cat = factor(data$min_pos_timing_cat)
data$anxdep = factor(data$anxdep)
data$sex = factor(data$sex)

### Summary of ASQ cut offs
# Numbers of individuals who are typical or clinical for the 5 subscales by SARS-CoV-2 group
# Get counts and percentages for each variable
# For your own information to quickly view
counts <- data %>%
  group_by(preg_pos) %>%
  summarise(
    count=n(), #count total occurences per group
    across(c(comm, gm, fm, ps, perso), ~paste0(table(.x), collapse = ", ")),  # Count for each level 
  )
#print(counts)

percentages <- asq_cat_per <- data %>%
  group_by(preg_pos) %>%
  summarise(
    count=n(), #count total occurrences per group
    across(c(comm, gm, fm, ps, perso), ~paste0(round(prop.table(table(.x)) * 100, 1), collapse = ", "))  # Percentage for each level
  )
#print(percentages)

# Prepare the data to create a table
counts_1 <- data %>% 
   group_by(preg_pos) %>%
  count(comm, name = "Communication") %>%
  mutate('Communication %'= round((Communication / sum(Communication)) * 100, 1))%>%   # Percentage rounded to 1 decimal
  rename(asq = comm)
counts_2 <- data %>% 
   group_by(preg_pos) %>%
  count(gm, name = "GrossMotor") %>% 
  mutate('Gross Motor %' = round((GrossMotor / sum(GrossMotor)) * 100, 1)) %>% 
  rename(asq = gm)
counts_3 <- data %>% 
  filter(!is.na(fm))%>%
  group_by(preg_pos) %>%
  count(fm, name = "FineMotor")  %>%
  mutate('Fine Motor %' = round((FineMotor / sum(FineMotor)) * 100, 1)) %>% 
  rename(asq = fm)
counts_4 <- data %>% 
  group_by(preg_pos) %>%
  count(ps, name = "ProblemSolving")  %>%
  mutate('Problem Solving %' = round((ProblemSolving / sum(ProblemSolving)) * 100, 1)) %>% 
  rename(asq = ps)
counts_5 <- data %>% 
  group_by(preg_pos) %>%
  count(perso, name = "PersonalSocial")  %>%
  mutate('Personal Social %' = round((PersonalSocial / sum(PersonalSocial)) * 100, 1))%>% 
  rename(asq = perso)

# Merge all of the counts and percentages together
counts <-  counts_1 %>%                                 
  inner_join(counts_2, by = c("preg_pos", "asq")) %>%
  inner_join(counts_3, by = c("preg_pos", "asq")) %>%
  inner_join(counts_4, by = c("preg_pos", "asq")) %>%
  inner_join(counts_5, by = c("preg_pos", "asq")) 
  
# Renaming is unneccessary and tedious but avoids picking the wrong column by accident
counts_name <- counts%>%
    rename('SARS-CoV-2 positive'= preg_pos, 'ASQ category' = asq,'Gross Motor' = GrossMotor, 'Fine Motor' = FineMotor, 'Problem Solving' = ProblemSolving, 'Personal Social' = PersonalSocial) 
 
# Pivot the table so that ASQ category becomes a column rather than being presented long
wide <- counts_name %>%
    pivot_wider(names_from = 'ASQ category' , values_from = c(Communication, 'Communication %', 'Gross Motor', 'Gross Motor %', 'Fine Motor', 'Fine Motor %', 'Problem Solving', 'Problem Solving %', 'Personal Social', 'Personal Social %'))   %>%  
  dplyr::select(`SARS-CoV-2 positive`, Communication_Typical, `Communication %_Typical`, Communication_Clinical, `Communication %_Clinical`, 
                `Gross Motor_Typical`, `Gross Motor %_Typical`, `Gross Motor_Clinical`, `Gross Motor %_Clinical`, 
                 `Fine Motor_Typical`, `Fine Motor %_Typical`, `Fine Motor_Clinical`, `Fine Motor %_Clinical`, 
                 `Problem Solving_Typical`, `Problem Solving %_Typical`, `Problem Solving_Clinical`, `Problem Solving %_Clinical`, 
                  `Personal Social_Typical`, `Personal Social %_Typical`, `Personal Social_Clinical`, `Personal Social %_Clinical`) 
 
# Create a table with the counts and percentages
formatted_table_counts <- wide %>%
  kable("html",  col.names = c("SARS-CoV-2 positive", "n", "%", "n", "%","n", "%","n", "%","n", "%", "n", "%", "n", "%","n", "%","n", "%","n", "%"), align = "rccc") %>%
  column_spec(1, bold = T) %>%
  row_spec(0, bold = TRUE) %>%
  add_header_above(c(" ", "Typical" = 2, "Clinical" = 2, "Typical" = 2, "Clinical" = 2, "Typical" = 2, "Clinical" = 2, "Typical" = 2, "Clinical" = 2, "Typical" = 2, "Clinical" = 2), align = "l")%>%
add_header_above(c(" ", "Communication" = 4, "Gross Motor" = 4, "Fine Motor" = 4, "Problem Solving" = 4, "Personal Social" = 4)) %>%
#  collapse_rows(columns = 1, align = "middle")%>%
kable_styling("condensed","hover", full_width = F, font_size = 13, html_font = 'helvetica') 
formatted_table_counts


#Check the data for missing values and to see what it looks like (different methods to check the same thing)
#describe(data) 
#md.pattern(data)
sapply(data, function(x) sum(is.na(x))) # Missing ASQ items were imputed using the method described in their manual before being loaded here. 
#If more than two ASQ items are missing for a subscale, the subscale cannot be imputed, which is why there is one person who doesn't have fm scores (this participant had 3 items missing)
#Importantly, we have 35 participants with missing education data and 34 participants with missing EPDS data, which we will impute using MICE (see below). 
# Please note: The 34 participants with missing ASQ age confirmed that their child was between 5 months and 6 months and 30 days on the questionnaire. 

# aggr_plot <- aggr(data, col=c('navyblue','red'),
#                   numbers=TRUE,
#                  sortVars=TRUE,
#                  labels=names(data),
#                  cex.axis=.7,
#                  gap=3,
#                   ylab=c("Histogram of Missing dat","Pattern"))

# Impute missing data
# Exclude Infection Timing variable as imputation will not be reliable and we aren't interested in it here by setting its method to "none"
no_timing <- data[, !(names(data) %in% "min_pos_timing_cat")]

# Run imputations
imputed_data <- mice(no_timing, m = 50, method = 'rf', seed = 123, printFlag = TRUE) # Imputing data using Random Forest because it is effective for both categorical (education) and continuous (EPDS) variables, especially when dealing with complex relationships and interactions.
# If you were interested in only working with one imputation: Combine the imputed data with the excluded variable. But really we aren't going to use this below, this is just FYI
#completed_data <- complete(imputed_data)
#completed_data <- cbind(completed_data, data[["min_pos_timing_cat"]])

# Test normality of the ASQ scores
hist(data$comm_total)
hist(data$gm_total)
hist(data$fm_total)
hist(data$ps_total)
hist(data$perso_total)

min(data$perso_total, na.rm = TRUE)

# LOGISTIC REGRESSIONS for categorical ASQ scores
# Communication
a <- with(imputed_data, glm(comm ~ preg_pos + raceethnicitycombined + edu + gestationalagedays + ageatdelivery + epds + sex, data = data, family = "binomial"))
pool_1 <- pool(a)
summary_1 <- summary(pool_1, conf.int = TRUE)
# Use broom's tidy function to get a clean summary of coefficients
tidy_1 <- tidy(pool_1)
# Extract p-values
p_values_a <- tidy_1$p.value
#p_values_a
# Apply FDR correction
fdr_1 <- p.adjust(p_values_a, method = "fdr")
# Add FDR adjusted p-values to the summary
summary_1$coefficients <- cbind(summary_1$coefficients, fdr_1)
print(summary_1)

# GM 
b <- with(imputed_data, glm(gm ~ preg_pos + raceethnicitycombined + edu + gestationalagedays + ageatdelivery + epds + sex, data = data , family = "binomial"))
pool_2 <- pool(b)
summary_2 <- summary(pool_2, conf.int = TRUE)
# Use broom's tidy function to get a clean summary of coefficients
tidy_2 <- tidy(pool_2)
# Extract p-values
p_values_b <- tidy_2$p.value
p_values_b
# Apply FDR correction
fdr_2 <- p.adjust(p_values_b, method = "fdr")
# Add FDR adjusted p-values to the summary
summary_2$coefficients <- cbind(summary_2$coefficients, fdr_2)
print(summary_2)

# FM
c <- with(imputed_data, glm(fm ~ preg_pos + raceethnicitycombined + edu + gestationalagedays + ageatdelivery + epds + sex,  data = data , family = "binomial"))
pool_3 <- pool(c)
summary_3 <- summary(pool_3, conf.int = TRUE)
# Use broom's tidy function to get a clean summary of coefficients
tidy_3 <- tidy(pool_3)
# Extract p-values
p_values_c <- tidy_3$p.value
p_values_c
# Apply FDR correction
fdr_3 <- p.adjust(p_values_c, method = "fdr")
# Add FDR adjusted p-values to the summary
summary_3$coefficients <- cbind(summary_3$coefficients, fdr_3)
print(summary_3)

# PS
d <- with(imputed_data, glm(ps ~ preg_pos + raceethnicitycombined + edu + gestationalagedays + ageatdelivery + epds + sex, data = data , family = "binomial"))
pool_4 <- pool(d)
summary_4 <- summary(pool_4, conf.int = TRUE)
# Use broom's tidy function to get a clean summary of coefficients
tidy_4 <- tidy(pool_4)
# Extract p-values
p_values_d <- tidy_4$p.value
p_values_d
# Apply FDR correction
fdr_4 <- p.adjust(p_values_d, method = "fdr")
# Add FDR adjusted p-values to the summary
summary_4$coefficients <- cbind(summary_4$coefficients, fdr_4)
print(summary_4)

# Perso
e <- with(imputed_data, glm(perso ~ preg_pos + raceethnicitycombined + edu + gestationalagedays + ageatdelivery + epds + sex, data = data , family = "binomial"))
pool_5 <- pool(e)
summary_5 <- summary(pool_5, conf.int = TRUE)
# Use broom's tidy function to get a clean summary of coefficients
tidy_5 <- tidy(pool_5)
# Extract p-values
p_values_e <- tidy_5$p.value
p_values_e
# Apply FDR correction
fdr_5 <- p.adjust(p_values_e, method = "fdr")
# Add FDR adjusted p-values to the summary
summary_5$coefficients <- cbind(summary_5$coefficients, fdr_5)
print(summary_5)

# Combine the outputs
ms <- list(Communication = summary_1, `Gross Motor` = summary_2, `Fine Motor` = summary_3,`Problem Solving` = summary_4, `Personal Social` = summary_5)
head(ms)
combined_summaries <- bind_rows(ms, .id = "model")

# Rename columns
renamed_columns <- combined_summaries %>%
  mutate(
    estimate = round(estimate, 3),
    p.value = round(p.value, 3),
    `2.5 %` = round(`2.5 %`, 2) ,
    `97.5 %` = round(`97.5 %`, 2),
    coefficients  = round(coefficients, 3)) 

# Select important rows and columns 
filtered <-  renamed_columns  %>%
  filter(term == "preg_posYes") %>%
  dplyr::select(c(model, estimate, p.value, `2.5 %`, `97.5 %`, coefficients))

# Regression table for the categorical outcomes 
formatted_table_cat <- filtered %>%
  kable("html", col.names = c("ASQ Subscale", "Estimate", "p-value", "Lower", "Upper", "p-value" ), align = "lccc") %>%
  kable_styling("condensed", "hover", full_width = F, font_size = 13, html_font = 'helvetica') %>%
  add_header_above(c(" ", " ", " ",  "95% CI" = 2, "FDR-corrected" = 1)) %>%
  column_spec(1, bold = T) %>%
  row_spec(0, bold = TRUE) 
formatted_table_cat

## LINEAR REGRESSIONS for total scores
# Communication
com_rg <- with(imputed_data, glm(comm_total ~ preg_pos + raceethnicitycombined + edu + gestationalagedays + ageatdelivery + epds + sex, data = data, family = Gamma(link = "log")))
pool_1 <- pool(com_rg)               
summary_1 <- summary(pool_1, conf.int = TRUE)
# Use broom's tidy function to get a clean summary of coefficients
tidy_1 <- tidy(pool_1)
# Extract p-values
p_values_a <- tidy_1$p.value
#p_values_a
# Apply FDR correction
fdr_1 <- p.adjust(p_values_a, method = "fdr")
# Add FDR adjusted p-values to the summary
summary_1$coefficients <- cbind(summary_1$coefficients, fdr_1)
print(summary_1)

# Gross Motor
data$gm_adjusted <- data$gm_total + 0.001  # Adding a small constant because GM contains 0s
gm_rg <- with(imputed_data, glm(gm_adjusted ~ preg_pos + raceethnicitycombined + edu + gestationalagedays + ageatdelivery + epds + sex, data = data, family = Gamma(link = "log")))
pool_2 <- pool(gm_rg)               
summary_2 <- summary(pool_2, conf.int = TRUE)
# Use broom's tidy function to get a clean summary of coefficients
tidy_2 <- tidy(pool_2)
# Extract p-values
p_values_b <- tidy_2$p.value
p_values_b
# Apply FDR correction
fdr_2 <- p.adjust(p_values_b, method = "fdr")
# Add FDR adjusted p-values to the summary
summary_2$coefficients <- cbind(summary_2$coefficients, fdr_2)
print(summary_2)

# Fine Motor
data$fm_adjusted <- data$fm_total + 0.001  # Adding a small constant because GM contains 0s
fm_rg <- with(imputed_data, glm(fm_adjusted ~ preg_pos + raceethnicitycombined + edu + gestationalagedays + ageatdelivery + epds + sex, data = data, family = Gamma(link = "log")))
pool_3 <- pool(fm_rg)               
summary_3 <- summary(pool_3, conf.int = TRUE)
# Use broom's tidy function to get a clean summary of coefficients
tidy_3 <- tidy(pool_3)
# Extract p-values
p_values_c <- tidy_3$p.value
p_values_c
# Apply FDR correction
fdr_3 <- p.adjust(p_values_c, method = "fdr")
# Add FDR adjusted p-values to the summary
summary_3$coefficients <- cbind(summary_3$coefficients, fdr_3)
print(summary_3)

# Problem Solving
data$ps_adjusted <- data$ps_total + 0.001  # Adding a small constant because GM contains 0s
ps_rg <- with(imputed_data, glm(ps_adjusted ~ preg_pos + raceethnicitycombined + edu + gestationalagedays + ageatdelivery + epds + sex, data = data, family = Gamma(link = "log")))
pool_4 <- pool(ps_rg) 
summary_4 <- summary(pool_4, conf.int = TRUE)
# Use broom's tidy function to get a clean summary of coefficients
tidy_4 <- tidy(pool_4)
# Extract p-values
p_values_d <- tidy_4$p.value
p_values_d
# Apply FDR correction
fdr_4 <- p.adjust(p_values_d, method = "fdr")
# Add FDR adjusted p-values to the summary
summary_4$coefficients <- cbind(summary_4$coefficients, fdr_4)
print(summary_4)

# Personal Social
data$perso_adjusted <- data$perso_total + 0.001  # Adding a small constant because GM contains 0s
perso_rg <- with(imputed_data, glm(perso_adjusted ~ preg_pos + raceethnicitycombined + edu + gestationalagedays + ageatdelivery + epds + sex, data = data, family = Gamma(link = "log")))
pool_5 <- pool(perso_rg) 
summary_5 <- summary(pool_5, conf.int = TRUE)
# Use broom's tidy function to get a clean summary of coefficients
tidy_5 <- tidy(pool_5)
# Extract p-values
p_values_e <- tidy_5$p.value
p_values_e
# Apply FDR correction
fdr_5 <- p.adjust(p_values_e, method = "fdr")
# Add FDR adjusted p-values to the summary
summary_5$coefficients <- cbind(summary_5$coefficients, fdr_5)
print(summary_5)

ms1 <- list(Communication = summary_1, `Gross Motor` = summary_2, `Fine Motor` = summary_3,`Problem Solving` = summary_4, `Personal Social` = summary_5)
head(ms1)

combined_summaries1 <- bind_rows(ms1, .id = "model")

# Rename columns
renamed_columns1 <- combined_summaries1 %>%
  mutate(
    estimate = round(estimate, 3),
    p.value = round(p.value, 3),
    `2.5 %` = round(`2.5 %`, 2) ,
    `97.5 %` = round(`97.5 %`, 2),
    coefficients  = round(coefficients, 3)) 

# Select only the important rows and columns 
filtered1 <-  renamed_columns1  %>%
  filter(term == "preg_posYes") %>%
  dplyr::select(c(model, estimate, p.value, `2.5 %`, `97.5 %`, coefficients))
colnames(filtered1)

# Create a regression table for the continuous scores
formatted_table_total <- filtered1 %>%
  kable("html", col.names = c("ASQ Subscale", "Estimate", "p-value", "Lower", "Upper", "p-value" ), align = "lccc") %>%
  kable_styling("condensed", "hover", full_width = F, font_size = 13, html_font = 'helvetica') %>%
  add_header_above(c(" ", " ", " ",  "95% CI" = 2, "FDR-corrected" = 1)) %>%
  column_spec(1, bold = T) %>%
  row_spec(0, bold = TRUE) 
formatted_table_total


#How to Test Model Assumptions
###Linearity
#Residuals vs. Fitted Values Plot: You can create a scatter plot of the residuals versus fitted values to check for linearity.
# Check the residual plot for a random scatter around zero.

# Residuals vs. Fitted Plot
#plot(com_rg, which = 1)  # This will show a plot of residuals vs fitted values

###Independence
#Durbin-Watson Test: This test checks for autocorrelation in the residuals.
#The Durbin-Watson statistic should be close to 2 (values between 1.5 and 2.5 are generally acceptable).
dw_test <- dwtest(com_rg)
print(dw_test)

###Homoscedasticity
# Breusch-Pagan Test: This test can formally test for homoscedasticity.
#Look for a non-significant p-value in the Breusch-Pagan test.
bptest(com_rg)

###Normality of Residuals (for continuous outcomes)
# Shapiro-Wilk Test
# The Shapiro-Wilk test should have a non-significant p-value (p > 0.05).
shapiro.test(residuals(com_rg))

###No Multicollinearity (only for multiple regression!)
#Variance Inflation Factor (VIF): This quantifies how much the variance is inflated due to multicollinearity.
# A VIF value greater than 5-10 indicates potential multicollinearity.
# Calculate VIF
#vif_values <- vif(com_rg)
#print(vif_values)
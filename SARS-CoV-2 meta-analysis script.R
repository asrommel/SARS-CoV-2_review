## Load packages
library(dplyr)
library(mice)
library(tidyverse)
library(Hmisc)
library(VIM)
require(foreign)
require(ggplot2)
require(MASS)
require(Hmisc)
require(reshape2)
require(metafor)

library(metafor)

##########################################################################################################
################################### Meta-analysis of hearing data ########################################
##########################################################################################################

dft <- read.csv("PATHWAY/Auditory.csv", header=T)


### Add Q-test, I^2, and tau^2 estimate info
mlabfun <- function(text, x) {
  list(bquote(paste(.(text),
                    " (Q = ", .(fmtx(x$QE, digits=2)),
                    ", df = ", .(x$k - x$p), ", ",
                    .(fmtp(x$QEp, digits=3, pname="p", add0=TRUE, sep=TRUE, equal=TRUE)), "; ",
                    I^2, " = ", .(fmtx(x$I2, digits=1)), "%, ",
                    tau^2, " = ", .(fmtx(x$tau2, digits=2)), ")")))}

##############################################################################################################################################################################################################
## FIRST SCREENING
### Calculate log risk ratios and corresponding sample variances for all hearing studies - FIRST SCREENING
dft_1 <- escalc(measure="RR", ai=ex_f_1, bi=ex_p_1,
                ci=un_f_1, un_p_1,
                slab=paste(Author,", ", Year, sep=""), data=dft)
dft_1
# Random-effects model (using log risk ratios and variances as input)
res_1 <- rma(yi, vi, data=dft_1)
res_1
# Predict pooled risk ratio and corresponding CI/PI
predict(res_1, transf=exp, digits=2)

# Egger's test
regtest(res_1)

#word
#reporter(res_1, format="word_document")

##############################################################################################################################################################################################################
## SECOND SCREENING
### Calculate log risk ratios and corresponding sample variances for all hearing studies - SECOND SCREENING
dft_2 <- escalc(measure="RR", ai=ex_f_2, bi=ex_p_2,
                ci=un_f_2, un_p_2,
                slab=paste(Author,", ", Year, sep=""), data=dft)
dft_2
# Random-effects model (using log risk ratios and variances as input)
res_2 <- rma(yi, vi, data=dft_2)
res_2
# Predict pooled risk ratio and corresponding CI/PI
predict(res_2, transf=exp, digits=2)

# Egger's test
regtest(res_2)

#word
#reporter(res_2, format="word_document")
##############################################################################################################################################################################################################

#COMBINED FIGURE
##############################################################################################################################################################################################################
folder <- "PATHWAY/Figures"
file_name <- "forestplot_auditory_combined.png"

# Save high resolution 300dpi
png(file.path(folder, file_name), width = 20, height = 24, units = "in", res = 700)

# Set up a 3-row, 2-column layout
par(mfrow=c(2, 1))
par(mar = c(7, 5, 4, 2))
par(mgp = c(5.5, 3, 0))

forest(res_1, addpred=T, header="Authors and Years", atransf=exp, at=log(c(0.5, 1, 2, 4, 8)), 
       ilab=cbind(dft_1$ex_f_1, dft_1$ex_p_1, dft_1$un_f_1, dft_1$un_p_1),
       ilab.xpos=c(-2.6, -2.2, -1.5, -1), cex=2.2, psize=1.2, mlab=mlabfun("RE Model", res_1),
       xlim=c(-6, 4),xlab="Risk Ratio for First Hearing Screening (log scale)" )


### add additional column headings to the plot
text(c(-2.65, -2.1, -1.45, -0.9), 12.5, font=2, cex=2,  c("Refer", "Pass", "Refer", "Pass"))
text(c(-2.4,-1.2), 11.5, font=2, cex=2, c("Exposed", "Unexposed"))
text(-5.9, 13, "a", cex=2.5, font=2)

# switch to bold italic font
par(font=4)

forest(res_2, addpred=T, header="Authors and Years", atransf=exp, at=log(c(0.5, 1, 2, 4, 8)), 
       ilab=cbind(dft_2$ex_f_2, dft_2$ex_p_2, dft_2$un_f_2, dft_2$un_p_2),
       ilab.xpos=c(-2.6, -2.2, -1.5, -1), cex=2.2, psize=1.2, mlab=mlabfun("RE Model", res_2),
       xlim=c(-6, 4),xlab="Risk Ratio for Second Hearing Screening (log scale)" )

### add additional column headings to the plot
text(c(-2.65, -2.1, -1.45, -0.9), 12.5, font=2, cex=2,  c("Refer", "Pass", "Refer", "Pass"))
text(c(-2.4,-1.2), 11.5, font=2, cex=2, c("Exposed", "Unexposed"))
text(-5.9, 13, "b", cex=2.5, font=2)
# switch to bold italic font
par(font=4)

dev.off()


##########################################################################################################
######################### Meta-analysis of categorical ASQ scores ########################################
##########################################################################################################

## Read categorical ASQ-3 data
data <- read.csv("PATHWAY/ASQ cat.csv", header=T)

## Read categorical ASQ-SE data
data_se <- read.csv("PATHWAY/ASQSE cat.csv", header=T)


### Add Q-test, I^2, and tau^2 estimate info
mlabfun <- function(text, x) {
  list(bquote(paste(.(text),
                    " (Q = ", .(fmtx(x$QE, digits=2)),
                    ", df = ", .(x$k - x$p), ", ",
                    .(fmtp(x$QEp, digits=3, pname="p", add0=TRUE, sep=TRUE, equal=TRUE)), "; ",
                    I^2, " = ", .(fmtx(x$I2, digits=1)), "%, ",
                    tau^2, " = ", .(fmtx(x$tau2, digits=2)), ")")))}

##############################################################################################################################################################################################################
### Calculate log risk ratios and corresponding sample variances for communication skills
# 'slab' refers to the study labels of the dataset
df_com <- escalc(measure="RR", ai=com_ex_d, bi=com_ex_t,
                 ci=com_un_d, di=com_un_t,
                 slab=paste(Author,", ", Year, sep=""), data=data)

df_com

## Random-effects model (using log risk ratios and variances as input)
res_com <- rma(yi, vi, data=df_com)
res_com

# Predict pooled risk ratio and corresponding CI/PI
predict(res_com, transf=exp, digits=2)

# Egger
regtest(res_com)

#word
#reporter(res_com, format="word_document")

##############################################################################################################################################################################################################
### Calculate log risk ratios and corresponding sample variances for gross motor skills
# 'slab' refers to the study labels of the dataset
df_gm <- escalc(measure="RR", ai=gm_ex_d, bi=gm_ex_t,
                ci=gm_un_d, di=gm_un_t,
                slab=paste(Author,", ", Year, sep=""), data=data)
df_gm

# Random-effects model (using log risk ratios and variances as input)
res_gm <- rma(yi, vi, data=df_gm)
res_gm

# Predict pooled risk ratio and corresponding CI/PI
predict(res_gm, transf=exp, digits=2)

# Egger
regtest(res_gm)

#word
#reporter(res_gm, format="word_document")

##############################################################################################################################################################################################################
### Calculate log risk ratios and corresponding sample variances for fine motor skills
# 'slab' refers to the study labels of the dataset
df_fm <- escalc(measure="RR", ai=fm_ex_d, bi=fm_ex_t,
                ci=fm_un_d, di=fm_un_t,
                slab=paste(Author,", ", Year, sep=""), data=data)

df_fm

# Random-effects model (using log risk ratios and variances as input)
res_fm <- rma(yi, vi, data=df_fm)
res_fm

# Predict pooled risk ratio and corresponding CI/PI
predict(res_fm, transf=exp, digits=2)

# Egger
regtest(res_fm)

#word
#reporter(res_fm, format="word_document")

##############################################################################################################################################################################################################
### Calculate log risk ratios and corresponding sample variances for problem solving skills
df_ps <- escalc(measure="RR", ai=ps_ex_d, bi=ps_ex_t,
                ci=ps_un_d, di=ps_un_t,
                slab=paste(Author,", ", Year, sep=""), data=data)

df_ps

# Random-effects model (using log risk ratios and variances as input)
res_ps<- rma(yi, vi, data=df_ps)
res_ps


# Predict pooled risk ratio and corresponding CI/PI
predict(res_ps, transf=exp, digits=2)

# Egger
regtest(res_ps)

#word
#reporter(res_ps, format="word_document")

##############################################################################################################################################################################################################
### Calculate log risk ratios and corresponding sample variances for personal-social skills
df_perso <- escalc(measure="RR", ai=perso_ex_d, bi=perso_ex_t,
                   ci=perso_un_d, di=perso_un_t,
                   slab=paste(Author,", ", Year, sep=""), data=data)

df_perso

# Random-effects model (using log risk ratios and variances as input)
res_perso<- rma(yi, vi, data=df_perso)
res_perso


# Predict pooled risk ratio and corresponding CI/PI
predict(res_perso, transf=exp, digits=2)

# Egger
regtest(res_perso)

#word
#reporter(res_perso, format="word_document")

##############################################################################################################################################################################################################
### Calculate log risk ratios and corresponding sample variances for ASQ social-emotional
df_se <- escalc(measure="RR", ai=asqse_ex_d, bi=asqse_ex_t,
                ci=asqse_un_d, di=asqse_un_t,
                slab=paste(Author,", ", Year, sep=""), data=data_se)

df_se

## Random-effects model (using log risk ratios and variances as input)
res_se <- rma(yi, vi, data=df_se)
res_se


# Predict pooled risk ratio and corresponding CI/PI
predict(res_se, transf=exp, digits=2)

# Egger
regtest(res_se)

#word
#reporter(res_se, format="word_document")

##############################################################################################################################################################################################################
## FIGURES
##############################################################################################################################################################################################################

# Set up the folder and file name for the combined plot
folder <- "PATHWAY/Figures"
file_name <- "forestplot_asq_cat_combined.png"

# Save high resolution 700dpi
png(file.path(folder, file_name), width = 27, height = 23, units = "in", res = 700)

# Set up a 3-row, 2-column layout
par(mfrow=c(3, 2))
par(mgp = c(6, 3, 0))
par(mar = c(8, 4, 4, 3))


# Function to create a forest plot with a label
create_forest_plot <- function(res, data, ex_d, ex_t, un_d, un_t, xlab, label) {
  forest(res, addpred=T, header="Authors and Year", atransf=exp, at=log(c(0.05, 0.25, 1, 4)), 
         ilab=cbind(data[[ex_d]], data[[ex_t]], data[[un_d]], data[[un_t]]),
         ilab.xpos=c(-6.1, -5.5, -4.6, -4), cex=1.6, psize=1.1, mlab=mlabfun("RE Model", res),
         xlim=c(-10, 4), xlab=xlab)
  text(-9.5, 9, label, cex=1.5, font=2)
}

# Plot 1: Communication Skills (df_com)
res_com <- rma(yi, vi, data=df_com)
create_forest_plot(res_com, data, "com_ex_d", "com_ex_t", "com_un_d", "com_un_t", 
                   "Risk Ratio for Communication Skills (log scale)","")
### add additional column headings to the plot
text(c(-6.3, -5.6, -4.7, -4), 10.5, font=2, cex=1.5,  c("Clinical", "Typical", "Clinical", "Typical"))
text(c(-6,-4.4), 11, font=2, cex=1.5, c("Exposed", "Unexposed"))
text(-9.5, 12, "a", cex=2.5, font=2)

#######################################
# Plot 2: Gross Motor Skills (df_gm)
res_gm <- rma(yi, vi, data=df_gm)
create_forest_plot(res_gm, data, "gm_ex_d", "gm_ex_t", "gm_un_d", "gm_un_t", 
                   "Risk Ratio for Gross Motor Skills (log scale)", "")
### add additional column headings to the plot
text(c(-6.3, -5.6, -4.7, -4), 10.5, font=2, cex=1.5,  c("Clinical", "Typical", "Clinical", "Typical"))
text(c(-6,-4.4), 11, font=2, cex=1.5, c("Exposed", "Unexposed"))
text(-9.5, 12, "b", cex=2.5, font=2)

######################################
# Plot 3: Fine Motor Skills (df_fm)
res_fm <- rma(yi, vi, data=df_fm)
create_forest_plot(res_fm, data, "fm_ex_d", "fm_ex_t", "fm_un_d", "fm_un_t", 
                   "Risk Ratio for Fine Motor Skills (log scale)", "")
### add additional column headings to the plot
text(c(-6.3, -5.6, -4.7, -4), 10.5, font=2, cex=1.5,  c("Clinical", "Typical", "Clinical", "Typical"))
text(c(-6,-4.4), 11, font=2, cex=1.5, c("Exposed", "Unexposed"))
text(-9.5, 12, "c", cex=2.5, font=2)

########################################
# Plot 4: Problem-Solving Skills (df_ps)
res_ps <- rma(yi, vi, data=df_ps)
create_forest_plot(res_ps, data, "ps_ex_d", "ps_ex_t", "ps_un_d", "ps_un_t", 
                   "Risk Ratio for Problem-Solving Skills (log scale)", "")
### add additional column headings to the plot
text(c(-6.3, -5.6, -4.7, -4), 10.5, font=2, cex=1.5,  c("Clinical", "Typical", "Clinical", "Typical"))
text(c(-6,-4.4), 11, font=2, cex=1.5, c("Exposed", "Unexposed"))
text(-9.5, 12, "d", cex=2.5, font=2)

#########################################
# Plot 5: Personal-Social Skills (df_perso)
res_perso <- rma(yi, vi, data=df_perso)
create_forest_plot(res_perso, data, "perso_ex_d", "perso_ex_t", "perso_un_d", "perso_un_t", 
                   "Risk Ratio for Personal-Social Skills (log scale)", "")
### add additional column headings to the plot
text(c(-6.3, -5.6, -4.7, -4), 10.5, font=2, cex=1.5,  c("Clinical", "Typical", "Clinical", "Typical"))
text(c(-6,-4.4), 11, font=2, cex=1.5, c("Exposed", "Unexposed"))
text(-9.5, 12, "e", cex=2.5, font=2)

#########################################
# Plot 6: Social-Emotional Skills (df_se)
res_se <- rma(yi, vi, data=df_se)
create_forest_plot(res_se, data_se, "asqse_ex_d", "asqse_ex_t", "asqse_un_d", "asqse_un_t", 
                   "Risk Ratio for ASQ Social-Emotional (log scale)", "")
### add additional column headings to the plot
text(c(-6.3, -5.6, -4.7, -4), 4.5, font=2, cex=1.5,  c("Clinical", "Typical", "Clinical", "Typical"))
text(c(-6,-4.4), 5, font=2, cex=1.5, c("Exposed", "Unexposed"))
text(-9.5, 9, "f", cex=2.5, font=2)

dev.off()

##########################################################################################################
################################ Meta-analysis of mean ASQ scores ########################################
##########################################################################################################


## Read data
dat_asq <- read.csv("PATHWAY/ASQ scores.csv", header=T)

### a little helper function to add Q-test, I^2, and tau^2 estimate info
mlabfun <- function(text, x) {
  list(bquote(paste(.(text),
                    " (Q = ", .(fmtx(x$QE, digits=2)),
                    ", df = ", .(x$k - x$p), ", ",
                    .(fmtp(x$QEp, digits=3, pname="p", add0=TRUE, sep=TRUE, equal=TRUE)), "; ",
                    I^2, " = ", .(fmtx(x$I2, digits=1)), "%, ",
                    tau^2, " = ", .(fmtx(x$tau2, digits=2)), ")")))}

##############################################################################################################################################################################################################
## Communication skills
da_com <- escalc(measure="SMD", m1i=com_ex_m, sd1i=com_ex_sd, n1i=ex_n,
                 m2i=com_un_m, sd2i=com_un_sd, n2i=un_n, slab=paste(Author,", ", Year, sep=""), data=dat_asq)
da_com

# Random-effects model
res_co <- rma(yi, vi, data=da_com)
res_co

# Egger's test
regtest(res_co)

#word
#reporter(res_co, format="word_document")

##############################################################################################################################################################################################################
# Gross Motor
da_gm <- escalc(measure="SMD", m1i=gm_ex_m, sd1i=gm_ex_sd, n1i=ex_n,
                m2i=gm_un_m, sd2i=gm_un_sd, n2i=un_n, slab=paste(Author,", ", Year, sep=""), data=dat_asq)
da_gm

res_g <- rma(yi, vi, data=da_gm)
res_g

# Egger
regtest(res_g)

#word
#reporter(res_g, format="word_document")

##############################################################################################################################################################################################################
# Fine Motor
da_fm <- escalc(measure="SMD", m1i=fm_ex_m, sd1i=fm_ex_sd, n1i=ex_n,
                m2i=fm_un_m, sd2i=fm_un_sd, n2i=un_n, slab=paste(Author,", ", Year, sep=""), data=dat_asq)
da_fm

res_f <- rma(yi, vi, data=da_fm)
res_f

# Egger's test
regtest(res_f)

#word
#reporter(res_f, format="word_document")

##############################################################################################################################################################################################################
# Problem solving
da_ps <- escalc(measure="SMD", m1i=ps_ex_m, sd1i=ps_ex_sd, n1i=ex_n,
                m2i=ps_un_m, sd2i=ps_un_sd, n2i=un_n, slab=paste(Author,", ", Year, sep=""), data=dat_asq)
da_ps

res_p <- rma(yi, vi, data=da_ps)
res_p

# Egger's test
regtest(res_p)

#word
#reporter(res_p, format="word_document")

##############################################################################################################################################################################################################
# Personal-social
da_perso <- escalc(measure="SMD", m1i=perso_ex_m, sd1i=perso_ex_sd, n1i=ex_n,
                   m2i=perso_un_m, sd2i=perso_un_sd, n2i=un_n, slab=paste(Author,", ", Year, sep=""), data=dat_asq)
da_perso

res_per <- rma(yi, vi, data=da_perso)
res_per

# Egger's test
regtest(res_per)

#word
#reporter(res_per, format="word_document")

##############################################################################################################################################################################################################
#COMBINED PLOT 5x2
##############################################################################################################################################################################################################
folder <- "PATHWAY/Figures"
file_name <- "forestplot_asq_msd_combined1.png"

# Save high resolution 700dpi
png(file.path(folder, file_name), width = 15, height = 25, units = "in", res = 700)

# Set up a 3-row, 2-column layout
#par(mfrow=c(5, 2))

#layout
layout(matrix(c(1,2,3,4,5,6,7,8,9,10), nrow = 5, ncol = 2, byrow = TRUE), widths = c(2,1))

#Communication
#Forest
forest(res_co, addpred=T, header="Authors and Years", atransf=exp, at=log(c(0.25, 1, 4)), 
       ilab=cbind(dat_asq$com_ex_m, dat_asq$com_ex_sd, dat_asq$com_un_m, dat_asq$com_un_sd),
       ilab.xpos=c(-5.3, -4.6, -3.7, -3), cex=1.5, psize=1.1, mlab=mlabfun("RE Model", res_co),
       xlim=c(-11, 4), xlab="Standardized Mean Difference for Communication Skills" )

### add additional column headings to the plot
text(c(-5.3, -4.6, -3.7, -3), 8.5, font=2, cex=1.5,  c("Mean", "SD", "Mean", "SD"))
text(c(-5,-3.4), 9.5, font=2, cex=1.5, c("Exposed", "Unexposed"))
text(-9.5, 10, "a", cex=2.5, font=2)
# switch to bold italic font
par(font=4)

#Funnel
funnel(res_co, ylim=c(0, 1.0), las=1, digits=list(1L,1))
text(-2, 0.1, "f", cex=2.5, font=2)

#Gross motor
#Forest
forest(res_g, addpred=T, header="Authors and Years", atransf=exp, at=log(c(0.25, 1, 4)), 
       ilab=cbind(dat_asq$gm_ex_m, dat_asq$gm_ex_sd, dat_asq$gm_un_m, dat_asq$gm_un_sd),
       ilab.xpos=c(-5.3, -4.6, -3.7, -3), cex=1.5, psize=1.1, mlab=mlabfun("RE Model", res_g),
       xlim=c(-10.5, 4), xlab="Standardized Mean Difference for Gross Motor Skills" )

### add additional column headings to the plot
text(c(-5.3, -4.6, -3.7, -3), 8.5, font=2, cex=1.5,  c("Mean", "SD", "Mean", "SD"))
text(c(-5,-3.4), 9.5, font=2, cex=1.5, c("Exposed", "Unexposed"))
text(-9.5, 10, "b", cex=2.5, font=2)
# switch to bold italic font
par(font=4)

#Funnel
funnel(res_g, ylim=c(0, 1.0), las=1, digits=list(1L,1))
text(-2, 0.1, "g", cex=2.5, font=2)

#Fine motor
#Forest
forest(res_f, addpred=T, header="Authors and Years", atransf=exp, at=log(c(0.25, 1, 4)), 
       ilab=cbind(dat_asq$fm_ex_m, dat_asq$fm_ex_sd, dat_asq$fm_un_m, dat_asq$fm_un_sd),
       ilab.xpos=c(-5.3, -4.6, -3.7, -3), cex=1.5, psize=1.1, mlab=mlabfun("RE Model", res_f),
       xlim=c(-10.5, 4), xlab="Standardized Mean Difference for Fine Motor Skills" )

### add additional column headings to the plot
text(c(-5.3, -4.6, -3.7, -3), 8.5, font=2, cex=1.5,  c("Mean", "SD", "Mean", "SD"))
text(c(-5,-3.4), 9.5, font=2, cex=1.5, c("Exposed", "Unexposed"))
text(-9.5, 10, "c", cex=2.5, font=2)
# switch to bold italic font
par(font=4)

#Funnel
funnel(res_f, ylim=c(0, 1.0), las=1, digits=list(1L,1))
text(-2, 0.1, "h", cex=2.5, font=2)

#Problem solving
#Forest
forest(res_p, addpred=T, header="Authors and Years", atransf=exp, at=log(c(0.25, 1, 4)), 
       ilab=cbind(dat_asq$ps_ex_m, dat_asq$ps_ex_sd, dat_asq$ps_un_m, dat_asq$ps_un_sd),
       ilab.xpos=c(-5.3, -4.6, -3.7, -3), cex=1.5, psize=1.1, mlab=mlabfun("RE Model", res_p),
       xlim=c(-10.5, 4), xlab="Standardized Mean Difference for Problem Solving Skills" )
### add additional column headings to the plot
text(c(-5.3, -4.6, -3.7, -3), 8.5, font=2, cex=1.5,  c("Mean", "SD", "Mean", "SD"))
text(c(-5,-3.4), 9.5, font=2, cex=1.5, c("Exposed", "Unexposed"))
text(-9.5, 10, "d", cex=2.5, font=2)
# switch to bold italic font
par(font=4)

#Funnel
funnel(res_p, ylim=c(0, 1.0), las=1, digits=list(1L,1))
text(-2, 0.1, "i", cex=2.5, font=2)

#Personal social
#Forest
forest(res_per, addpred=T, header="Authors and Years", atransf=exp, at=log(c(0.25, 1, 4)), 
       ilab=cbind(dat_asq$perso_ex_m, dat_asq$perso_ex_sd, dat_asq$perso_un_m, dat_asq$perso_un_sd),
       ilab.xpos=c(-5.3, -4.6, -3.7, -3), cex=1.5, psize=1.1, mlab=mlabfun("RE Model", res_per),
       xlim=c(-10.5, 4), xlab="Standardized Mean Difference for Communication Skills" )
### add additional column headings to the plot
text(c(-5.3, -4.6, -3.7, -3), 8.5, font=2, cex=1.5,  c("Mean", "SD", "Mean", "SD"))
text(c(-5,-3.4), 9.5, font=2, cex=1.5, c("Exposed", "Unexposed"))
text(-9.5, 10, "e", cex=2.5, font=2)
# switch to bold italic font
par(font=4)

#Funnel
funnel(res_per, ylim=c(0, 1.0), las=1, digits=list(1L,1))
text(-2, 0.1, "j", cex=2.5, font=2)

dev.off()


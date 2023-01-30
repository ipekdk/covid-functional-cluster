# This script has the codes used in the article:
# Functional Cluster and Canonical Correlation Analysis of EU Countries by Number of Daily Deaths and Stringency Index During Covid-19 Pandemic
# Please cite as:
#
# Keser, I. K., & Deveci Kocakoc, İ. (2021). Functional cluster and canonical correlation analysis of EU countries by number of daily deaths and stringency index during Covid-19 pandemic. Electronic Journal of Applied Statistical Analysis, 14(1), 197-216.

# Necessary libraries:

library(dplyr)
library(readxl)
library(lubridate, warn.conflicts = FALSE)
library(tidyselect)
library(magrittr)
library(rlang)
library(tidyr)
library(rio)
library(readr)
library(fda.usc)
library(tidyverse) # programming tools
library(fda) # functional data analysis
library(functional) # functional programming tools
library(MASS)
library(data.table)


# You should download the files below into your computer first.

# Main data file: https://github.com/ipekdk/covid-functional-cluster/blob/main/owid-covid-data_23.06.2020.xlsx
# Country selection list: https://github.com/ipekdk/covid-functional-cluster/blob/main/country_list-eu.xlsx
# Variable selection list: https://github.com/ipekdk/covid-functional-cluster/blob/main/variable_list.xlsx
# stringency index data: https://github.com/ipekdk/covid-functional-cluster/blob/main/covid-stringency-index.xlsx 
# All of the data above as an Rdata file: https://github.com/ipekdk/covid-functional-cluster/blob/main/data_all_github.RData
# Modified functional data plot function: https://github.com/ipekdk/covid-functional-cluster/blob/main/fd_plot.R

# You can select the countries and variables you like by writing 1 next to them in the selection list files.

# Importing data from your computer:

mdata <- read_excel("~/owid-covid-data_23.06.2020.xlsx")
country_list <- read_excel("~/country_list-eu.xlsx")
variable_list <- read_excel("~/variable_list.xlsx")
str_ind<-read_excel("~/covid-stringency-index.xlsx")

# or, you can load the data file from Rdata:

load(~/data_all_github.RData)

# cleaning variable names if necessary
library("janitor")
mdata <- clean_names(mdata)
mdata$date<-ymd(mdata$date)

# Looking at the number of data for each country
countries<-mdata %>% count(location)

# Take counrty and variable choices from the selection lists:
selected_c<-country_list[ which(country_list$select==1),"country" ]
selected_v<-t(variable_list[ which(variable_list$select==1),"v_name" ])

# Make a subset of selected countries 
small<- mdata[mdata$location %in% selected_c$country, ]
ind<-str_ind[str_ind$Entity %in% selected_c$country, ]

# Determine date range
b_date<-"2020-03-01"
e_date<-"2020-06-23" #max(mdata$date)-1 if all data has the same end

# Filter data in the selected date range
small<-filter(small, (date>b_date&date<e_date) )
ind<-filter(ind, (Date>b_date&Date<e_date) )

# Make a subset of selected variables 
small<-subset(small, select=as.character(selected_v))

# Make tables of date vs countries for the selected variables 

mdata_sets <- list()
for (i in 3:length(selected_v))
{c<-selected_v[i]
mdata_set<-subset(small, select=c("location","date",as.character(c)))
mdata_set<-spread(mdata_set,location,as.character(c))
mdata_sets[[paste0("v_", as.character(c))]] <- mdata_set
}

# Make a table of stringency index and add it to the main data set
stringency_index <- list()
c2<-"stringency_index"
mdata_set2<-subset(ind, select=c("Entity","Date",as.character(c2)))
mdata_set2<-spread(mdata_set2,Entity,as.character(c2))
mdata_sets[["v_stringency_index"]] <- mdata_set2

###----------------------------------------------------------###

# Preparing data for functional conversion

source("~/fd_plot.R", encoding = 'UTF-8')

# Get the data for the selected variable as "data1"

s<-4 # We selected new deaths per million variable, which is the 4th variable in mdata_sets
data1<-as.data.frame(mdata_sets[s])
rownames(data1)<-as.character(t(data1[1]))
data1[1]<-NULL
colnames(data1)<-c(t(selected_c))
data1<-na.omit(data1) # This is the ready-to-analyze data for the selected variable

# Get the data for stringency index as "data2" 

s<-5
data2<-as.data.frame(mdata_sets[s])
rownames(data2)<-as.character(t(data2[1]))
data2[1]<-NULL
dates<-rownames(data2)
data2 <- as.data.frame(sapply(data2, as.numeric)) 
rownames(data2)<-dates
data2<-na.omit(data2)

# data1 and data2 are combined based on "date" 
data2<-data2[rownames(data2) %in% rownames(data), ] 

###----------------------------------------------------------###

# preparing functional basis 

nb=10 #number of basis,selected subjectively by trial and error 
l=0.01 #smoothing parameter, selected subjectively by trial and error. You can change it to see its effect
dates<-as.Date(rownames(data2))
days <- 1:dim(data2)[1] # days 
basis <- create.bspline.basis(rangeval = range(days), nbasis = nb) # B-spline basis with 10 basis functions

# Functional object for the selected variable (Ours is the number of new deaths per million)
fdobj <- as.fd(smooth.basis(argvals = days,y=as.matrix(data),basis))
# Functional object for the stringency index
fdobj2 <- as.fd(smooth.basis(argvals = days,y=as.matrix(data2),basis))

###----------------------------------------------------------###

# plotting functional curves for the selected variable
fd_plot(fdobj,1)
title(main="Functional Curves")

# plotting functional curves for the stringency index
fd_plot(fdobj2,fl=3)
title(main="Stringency Index Curves")

# Taking the first derivatives and plotting them
f_deriv1<-deriv.fd(fdobj,1)
fd_plot(f_deriv1,2)
title(main="First Derivatives")

# Taking the second derivatives and plotting them
f_deriv2<-deriv.fd(fdobj,2)
fd_plot(f_deriv2,2)
title(main="Second Derivatives")

###----------------------------------------------------------###

# Curve registration (This is not a must, you can skip it if you think your data doesn't need registering)
# For more information on curve registration, please refer to: Ramsay et al.(2009) (https://link.springer.com/content/pdf/10.1007/978-0-387-98185-7.pdf, Section 8.)
# Also: https://www.rdocumentation.org/packages/fda/versions/5.1.9/topics/register.fd

# Curve registration for the selected variable

palette("default")

smBv <- deriv.fd(fdobj, 1)
ncountry <- 26
#  Define the target function as the mean of the first ncountry records
smBv0 = mean.fd(smBv[1:ncountry])
# Register these curves. The default choice for the functional
# parameter object WfdParObj is used.
smB.reg.0 <- register.fd(smBv0, smBv[1:ncountry])
# plot each curve. Click on the R Graphics window to show each plot.
# The left panel contains:
# -- the unregistered curve (dashed blue line)
# -- the target function (dashed red line)
# -- the registered curve (solid blue line)
# The right panel contains:
# -- the warping function h(t)
# -- the linear function corresponding to no warping

plotreg.fd(smB.reg.0)
# Notice that all the warping functions all have simple shapes
# due to the use of the simplest possible basis

# Now, define a more flexible basis for the warping functions
Wnbasis   <- nb
Wbasis    <- create.bspline.basis(rangeval = range(days), Wnbasis)
Wfd0      <- fd(matrix(0,Wnbasis,1),Wbasis)
#  set up the functional parameter object using only
#      a light amount smoothing
WfdParobj <- fdPar(Wfd0, Lfdobj=2, lambda=0.01)
#  register the curves
smB.reg.1 <- register.fd(smBv0, smBv[1:ncountry], WfdParobj)
plotreg.fd(smB.reg.1)

# Notice that now the warping functions can have more complex shapes

# Change the target to the mean of the registered functions ...
# this should provide a better target for registration

smBv1 <- mean.fd(smB.reg.1$regfd)
#  plot the old and the new targets
par(mfrow=c(1,1),ask=FALSE)
plot(smBv1)
lines(smBv0, lty=2)

#  Notice how the new target (solid line) has sharper features and
#  a stronger peek relative to the old target
#  (dashed line).  Now register to the new target
smB.reg.2 <- register.fd(smBv1, smBv[1:ncountry], WfdParobj) 
plotreg.fd(smB.reg.2)
#  Plot the mean of these curves as well as the first and second targets
par(mfrow=c(1,1),ask=FALSE)
plot(mean.fd(smB.reg.2$regfd))
lines(smBv0, lty=2)
lines(smBv1, lty=3)

#  Now register the previously registered functions to the new target.
smB.reg.3 <- register.fd(smBv1, smB.reg.1$regfd, WfdParobj)
plotreg.fd(smB.reg.3)

reg3<-smB.reg.3$regfd #registered curves

fd_plot(reg3,1)
title(main="Registered Functional Curves")

###----------------------------------------------------------###

# Functional cluster analysis

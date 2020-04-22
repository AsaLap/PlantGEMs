setwd("~/INRAE/Work/Plant-GEMs")
data_tom = read.csv("Data/tomato_values.csv", sep=',', header=TRUE)
data_kiw = read.csv("Data/kiwi_values.csv", sep=',', header=TRUE)
data_cuc = read.csv("Data/cucumber_values.csv", sep=',', header=TRUE)
data_che = read.csv("Data/cherry_values.csv", sep=',', header=TRUE)
data_cam = read.csv("Data/camelina_values.csv", sep=',', header=TRUE)

#Evaluating the linearity
library(carData)
scatterplot(Identity~E_Value, data = data_tom)

#Checking the parameters


#Linear regression

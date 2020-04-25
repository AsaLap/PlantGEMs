setwd("~/INRAE/Work/Plant-GEMs")
data_tom = read.csv("AllValues/Data/tomato_values.csv", sep=',', header=TRUE)
data_kiw = read.csv("AllValues/Data/kiwi_values.csv", sep=',', header=TRUE)
data_cuc = read.csv("AllValues/Data/cucumber_values.csv", sep=',', header=TRUE)
data_che = read.csv("AllValues/Data/cherry_values.csv", sep=',', header=TRUE)
data_cam = read.csv("AllValues/Data/camelina_values.csv", sep=',', header=TRUE)

library(ggplot2)

#Evaluating the linearity
#Between identity and E_Value:
plot(x = data_tom$Identity, y = data_tom$E_Value, xlab = 'Percentage of identity', ylab = 'E_Value', ylim=c(0, 1e-100))
plot(x = data_tom$Identity, y = data_tom$E_Value, ylim=c(0, 1e-100))
plot(x = data_kiw$Identity, y = data_kiw$E_Value)
plot(x = data_cuc$Identity, y = data_cuc$E_Value)
plot(x = data_che$Identity, y = data_che$E_Value)
plot(x = data_cam$Identity, y = data_cam$E_Value)
plot(x = data_cam$Identity, y = data_cam$E_Value, ylim=c(0, 1e-10))

#Between identity and score:
plot(x = data_tom$Identity, y = data_tom$Score, xlab = 'Percentage of identity', ylab = 'Score', ylim = c(0,6000))


#Linear regression
#Tomato
#Between identity and E_Value:
reg_tom_Evalue = lm(data_tom$Identity~data_tom$E_Value)
summary(reg_tom_Evalue)
#Between identity and Score:
reg_tom_Score = lm(data_tom$Identity~data_tom$Score)
summary(reg_tom_Score)

#Kiwi
#Between identity and E_Value:
reg_kiw_Evalue = lm(data_kiw$Identity~data_kiw$E_Value)
summary(reg_kiw_Evalue)
#Between identity and Score:
reg_kiw_Score = lm(data_kiw$Identity~data_kiw$Score)
summary(reg_kiw_Score)

#Cucumber
#Between identity and E_Value:
reg_cuc_Evalue = lm(data_cuc$Identity~data_cuc$E_Value)
summary(reg_cuc_Evalue)
#Between identity and Score:
reg_cuc_Score = lm(data_cuc$Identity~data_cuc$Score)
summary(reg_cuc_Score)

#Cherry
#Between identity and E_Value:
reg_che_Evalue = lm(data_che$Identity~data_che$E_Value)
summary(reg_kiw_Evalue)
#Between identity and Score:
reg_che_Score = lm(data_che$Identity~data_che$Score)
summary(reg_che_Score)

#Camelina
#Between identity and E_Value:
reg_cam_Evalue = lm(data_cam$Identity~data_cam$E_Value)
summary(reg_cam_Evalue)
#Between identity and Score:
reg_cam_Score = lm(data_cam$Identity~data_cam$Score)
summary(reg_cam_Score)
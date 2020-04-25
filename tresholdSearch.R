setwd("~/INRAE/Work/Plant-GEMs")
data_kiw = read.csv("Data/Treshold_50_30_1e-100_20_Score/KiwiTreshold.csv", sep=',', header=TRUE)
data_tom = read.csv("Data/Treshold_50_30_1e-100_20_Score/TomatoTreshold.csv", sep=',', header=TRUE)
data_cuc = read.csv("Data/Treshold_50_30_1e-100_20_Score/CucumberTreshold.csv", sep=',', header=TRUE)
data_che = read.csv("Data/Treshold_50_30_1e-100_20_Score/CherryTreshold.csv", sep=',', header=TRUE)
data_cam = read.csv("Data/Treshold_50_30_1e-100_20_Score/CamelinaTreshold.csv", sep=',', header=TRUE)

#Tomato double y axis
dataFrame = data_tom
par(mar = c(5, 5, 3, 5))
plot(dataFrame$Score, dataFrame$Nb.genes, type = "l", pch = 19,
     col = "blue", xlab = "Score", ylab = "Nombre de gènes du draft")
par(new = TRUE)
plot(dataFrame$Nb.reactions, type = "l", pch = 19,
     col = "red", lty = 2, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
axis(side = 4)
mtext("Nombre de réactions du draft", side = 4, line = 3)
legend("bottomleft", c("Gènes", "Réactions"), col = c("blue", "red"), lty = c(1, 2), cex = 1.2)
title("Draft Tomate 50/30/1e-100/20")

#Kiwi double y axis
dataFrame = data_kiw
par(mar = c(5, 5, 3, 5))
plot(dataFrame$Score, dataFrame$Nb.genes, type = "l", pch = 19,
     col = "blue", xlab = "Score", ylab = "Nombre de gènes du draft")
par(new = TRUE)
plot(dataFrame$Nb.reactions, type = "l", pch = 19,
     col = "red", lty = 2, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
axis(side = 4)
mtext("Nombre de réactions du draft", side = 4, line = 3)
legend("bottomleft", c("Gènes", "Réactions"), col = c("blue", "red"), lty = c(1, 2), cex = 1.2)
title("Draft Kiwi 50/30/1e-100/20")

#Cucumber double y axis
dataFrame = data_cuc
par(mar = c(5, 5, 3, 5))
plot(dataFrame$Score, dataFrame$Nb.genes, type = "l", pch = 19,
     col = "blue", xlab = "Score", ylab = "Nombre de gènes du draft")
par(new = TRUE)
plot(dataFrame$Nb.reactions, type = "l", pch = 19,
     col = "red", lty = 2, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
axis(side = 4)
mtext("Nombre de réactions du draft", side = 4, line = 3)
legend("bottomleft", c("Gènes", "Réactions"), col = c("blue", "red"), lty = c(1, 2), cex = 1.2)
title("Draft Concombre 50/30/1e-100/20")

#Cherry double y axis
dataFrame = data_che
par(mar = c(5, 5, 3, 5))
plot(dataFrame$Score, dataFrame$Nb.genes, type = "l", pch = 19,
     col = "blue", xlab = "Score", ylab = "Nombre de gènes du draft")
par(new = TRUE)
plot(dataFrame$Nb.reactions, type = "l", pch = 19,
     col = "red", lty = 2, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
axis(side = 4)
mtext("Nombre de réactions du draft", side = 4, line = 3)
legend("bottomleft", c("Gènes", "Réactions"), col = c("blue", "red"), lty = c(1, 2), cex = 1.2)
title("Draft Cerise 50/30/1e-100/20")

#Camelina double y axis
dataFrame = data_cam
par(mar = c(5, 5, 3, 5))
plot(dataFrame$Score, dataFrame$Nb.genes, type = "l", pch = 19,
     col = "blue", xlab = "Score", ylab = "Nombre de gènes du draft")
par(new = TRUE)
plot(dataFrame$Nb.reactions, type = "l", pch = 19,
     col = "red", lty = 2, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
axis(side = 4)
mtext("Nombre de réactions du draft", side = 4, line = 3)
legend("bottomleft", c("Gènes", "Réactions"), col = c("blue", "red"), lty = c(1, 2), cex = 1.2)
title("Draft Cameline 50/30/1e-100/20")


# #Tomato
# dataFrame = data_tom
# min(dataFrame[,2:3])
# plot(dataFrame$EValue..1e.x., dataFrame$Nb.genes, type = "l", pch = 19,
#      col = "red", xlab = "E-Value (1e-X)", ylab = "Quantité de genes/reactions conservés",
#      ylim = c(min(dataFrame[,2:3]), max(dataFrame[,2:3])))
# lines(dataFrame$EValue..1e.x., dataFrame$Nb.reactions, type = "l", pch = 19,
#       col = "blue", lty = 2)
# legend("bottomleft", legend=c("Gènes", "Réactions"),
#        col=c("red", "blue"), lty = 1:2, cex=1.2)
# title("Tomato 80/20")
# 
# #Kiwi
# dataFrame = data_kiw
# plot(dataFrame$EValue..1e.x., dataFrame$Nb.genes, type = "l", pch = 19,
#      col = "red", xlab = "E-Value (1e-X)", ylab = "Quantité de genes/reactions conservés",
#      ylim = c(min(dataFrame[,2:3]), max(dataFrame[,2:3])))
# lines(dataFrame$EValue..1e.x., dataFrame$Nb.reactions, type = "l", pch = 19,
#       col = "blue", lty = 2)
# legend("bottomleft", legend=c("Gènes", "Réactions"),
#        col=c("red", "blue"), lty = 1:2, cex=1.2)
# title("Kiwi 80/20")
# 
# #Cucumber
# dataFrame = data_cuc
# plot(dataFrame$EValue..1e.x., dataFrame$Nb.genes, type = "l", pch = 19,
#      col = "red", xlab = "E-Value (1e-X)", ylab = "Quantité de genes/reactions conservés",
#      ylim = c(min(dataFrame[,2:3]), max(dataFrame[,2:3])))
# lines(dataFrame$EValue..1e.x., dataFrame$Nb.reactions, type = "l", pch = 19,
#       col = "blue", lty = 2)
# legend("bottomleft", legend=c("Gènes", "Réactions"),
#        col=c("red", "blue"), lty = 1:2, cex=1.2)
# title("Cucumber 80/20")
# 
# #Cherry
# dataFrame = data_che
# plot(dataFrame$EValue..1e.x., dataFrame$Nb.genes, type = "l", pch = 19,
#      col = "red", xlab = "E-Value (1e-X)", ylab = "Quantité de genes/reactions conservés",
#      ylim = c(min(dataFrame[,2:3]), max(dataFrame[,2:3])))
# lines(dataFrame$EValue..1e.x., dataFrame$Nb.reactions, type = "l", pch = 19,
#       col = "blue", lty = 2)
# legend("bottomleft", legend=c("Gènes", "Réactions"),
#        col=c("red", "blue"), lty = 1:2, cex=1.2)
# title("Cherry 80/20")
# 
# #Cameline
# dataFrame = data_cam
# plot(dataFrame$EValue..1e.x., dataFrame$Nb.genes, type = "l", pch = 19,
#      col = "red", xlab = "E-Value (1e-X)", ylab = "Quantité de genes/reactions conservés",
#      ylim = c(min(dataFrame[,2:3]), max(dataFrame[,2:3])))
# lines(dataFrame$EValue..1e.x., dataFrame$Nb.reactions, type = "l", pch = 19,
#       col = "blue", lty = 2)
# legend("bottomleft", legend=c("Gènes", "Réactions"),
#        col=c("red", "blue"), lty = 1:2, cex=1.2)
# title("Camelina 80/20")

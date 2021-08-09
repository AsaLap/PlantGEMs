setwd("~/INRAE/Work/Plant-GEMs")
data_kiw = read.csv("Data/Treshold_0_100_1_0_BitScore/KiwiTreshold.csv", sep=',', header=TRUE)
data_tom = read.csv("Data/Treshold_0_100_1_0_BitScore/TomatoTreshold.csv", sep=',', header=TRUE)
data_cuc = read.csv("Data/Treshold_0_100_1_0_BitScore/CucumberTreshold.csv", sep=',', header=TRUE)
data_che = read.csv("Data/Treshold_0_100_1_0_BitScore/CherryTreshold.csv", sep=',', header=TRUE)
data_cam = read.csv("Data/Treshold_0_100_1_0_BitScore/CamelinaTreshold.csv", sep=',', header=TRUE)

#Tomato double y axis
dataFrame = data_tom
par(mar = c(5, 5, 3, 5))
plot(dataFrame$Bit_Score, dataFrame$Nb.genes, type = "l", pch = 19,
     col = "blue", xlab = "Score", ylab = "Nombre de gènes du draft", ylim = c(0, max(dataFrame$Nb.genes)))
par(new = TRUE)
plot(dataFrame$Nb.reactions, type = "l", pch = 19,
     col = "red", lty = 2, xaxt = "n", yaxt = "n", xlab = "", ylab = "", ylim = c(0, max(dataFrame$Nb.reactions)))
axis(side = 4)
mtext("Nombre de réactions du draft", side = 4, line = 3)
legend("bottomleft", c("Gènes", "Réactions"), col = c("blue", "red"), lty = c(1, 2), cex = 1.2)
title("Draft Tomate 0/100/1/0 BitScore")

#Kiwi double y axis
dataFrame = data_kiw
par(mar = c(5, 5, 3, 5))
plot(dataFrame$Bit_Score, dataFrame$Nb.genes, type = "l", pch = 19,
     col = "blue", xlab = "Score", ylab = "Nombre de gènes du draft", ylim = c(0, max(dataFrame$Nb.genes)))
par(new = TRUE)
plot(dataFrame$Nb.reactions, type = "l", pch = 19,
     col = "red", lty = 2, xaxt = "n", yaxt = "n", xlab = "", ylab = "", ylim = c(0, max(dataFrame$Nb.reactions)))
axis(side = 4)
mtext("Nombre de réactions du draft", side = 4, line = 3)
legend("bottomleft", c("Gènes", "Réactions"), col = c("blue", "red"), lty = c(1, 2), cex = 1.2)
title("Draft Kiwi 0/100/1/0 BitScore")

#Cucumber double y axis
dataFrame = data_cuc
par(mar = c(5, 5, 3, 5))
plot(dataFrame$Bit_Score, dataFrame$Nb.genes, type = "l", pch = 19,
     col = "blue", xlab = "Score", ylab = "Nombre de gènes du draft", ylim = c(0, max(dataFrame$Nb.genes)))
par(new = TRUE)
plot(dataFrame$Nb.reactions, type = "l", pch = 19,
     col = "red", lty = 2, xaxt = "n", yaxt = "n", xlab = "", ylab = "", ylim = c(0, max(dataFrame$Nb.reactions)))
axis(side = 4)
mtext("Nombre de réactions du draft", side = 4, line = 3)
legend("bottomleft", c("Gènes", "Réactions"), col = c("blue", "red"), lty = c(1, 2), cex = 1.2)
title("Draft Concombre 0/100/1/0 BitScore")

#Cherry double y axis
dataFrame = data_che
par(mar = c(5, 5, 3, 5))
plot(dataFrame$Bit_Score, dataFrame$Nb.genes, type = "l", pch = 19,
     col = "blue", xlab = "Score", ylab = "Nombre de gènes du draft", ylim = c(0, max(dataFrame$Nb.genes)))
par(new = TRUE)
plot(dataFrame$Nb.reactions, type = "l", pch = 19,
     col = "red", lty = 2, xaxt = "n", yaxt = "n", xlab = "", ylab = "", ylim = c(0, max(dataFrame$Nb.reactions)))
axis(side = 4)
mtext("Nombre de réactions du draft", side = 4, line = 3)
legend("bottomleft", c("Gènes", "Réactions"), col = c("blue", "red"), lty = c(1, 2), cex = 1.2)
title("Draft Cerise 0/100/1/0 BitScore")

#Camelina double y axis
dataFrame = data_cam
par(mar = c(5, 5, 3, 5))
plot(dataFrame$Bit_Score, dataFrame$Nb.genes, type = "l", pch = 19,
     col = "blue", xlab = "Score", ylab = "Nombre de gènes du draft", ylim = c(0, max(dataFrame$Nb.genes)))
par(new = TRUE)
plot(dataFrame$Nb.reactions, type = "l", pch = 19,
     col = "red", lty = 2, xaxt = "n", yaxt = "n", xlab = "", ylab = "", ylim = c(0, max(dataFrame$Nb.reactions)))
axis(side = 4)
mtext("Nombre de réactions du draft", side = 4, line = 3)
legend("bottomleft", c("Gènes", "Réactions"), col = c("blue", "red"), lty = c(1, 2), cex = 1.2)
title("Draft Cameline 0/100/1/0 BitScore")


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

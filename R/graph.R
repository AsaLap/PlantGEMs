setwd("~/INRAE/Work/Plant-GEMs/")
df_tot = read.csv("Data/AllValues/values.csv", sep=',', header=TRUE)
#rownames(e_values) = e_values[,1]
#e_values = e_values[,-1]
#e_values = t(e_values)

library(ggplot2)
library(plyr)

#For all organisms
muID <- ddply(df_tot, "Organism", summarise, grp.mean=mean(Identity))
muSCORE <- ddply(df_tot, "Organism", summarise, grp.mean=mean(Score))
muEVAL <- ddply(df_tot, "Organism", summarise, grp.mean=mean(E_Value))
muBIT <- ddply(df_tot, "Organism", summarise, grp.mean=mean(Bit_Score))


percID = ggplot(df_tot, aes(x=Identity, color=Organism, fill=Organism)) +
  geom_histogram(position="identity", alpha=0.1, binwidth = 5)+
  geom_vline(data=muID, aes(xintercept=grp.mean, color=Organism),
             linetype="dashed")+
  labs(title="Percentage of identity histogram plot",x="Identity(%)", y = "Number of counts")+
  theme_classic()
percID

score = ggplot(df_tot, aes(x=Score, color=Organism, fill=Organism)) +
  geom_histogram(position="identity", alpha=0.1, binwidth = 5)+
  geom_vline(data=muSCORE, aes(xintercept=grp.mean, color=Organism),
             linetype="dashed")+
  labs(title="Score histogram plot",x="Score", y = "Number of counts")+
  theme_classic()
score

e_value = ggplot(df_tot, aes(x=E_Value, color=Organism, fill=Organism)) +
  geom_histogram(position="identity", alpha=0.1, binwidth = 1)+
  geom_vline(data=muEVAL, aes(xintercept=grp.mean, color=Organism),
             linetype="dashed")+
  labs(title="E_Value histogram plot",x="E_Value", y = "Number of counts")+
  theme_classic()
e_value

bit_score = ggplot(df_tot, aes(x=Bit_Score, color=Organism, fill=Organism)) +
  geom_histogram(position="identity", alpha=0.1, binwidth = 5)+
  geom_vline(data=muBIT, aes(xintercept=grp.mean, color=Organism),
             linetype="dashed")+
  labs(title="Bit Score histogram plot",x="Bit Score", y = "Number of counts")+
  theme_classic()
bit_score

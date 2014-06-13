library(plyr)
library(ggplot2)

In_put <- read.table("qpcr_test//qpcr.csv", sep = ",", header = T)
In_put$Ct <- as.numeric(as.character( In_put$Ct))
Genes <- levels(In_put$Detector)
Ref_Gen <- "HPRT"
Test_Gen <- Genes[-which(Genes == Ref_Gen)]
Start_Data <- data.frame(In_put$Position, In_put$Detector, In_put$Sample, In_put$RNA, In_put$Ct)
Temp <- colnames(Start_Data)
colnames(Start_Data) <- substring(Temp, first = 8)
Data_Sum <- ddply(Start_Data, summarise, .variables = .(Detector, Sample, RNA), mean = mean(Ct), Err = sd(Ct))
Gen_Data <- subset(Data_Sum, Detector == Test_Gen, drop = T)
Ref_Data <- subset(Data_Sum, Detector == Ref_Gen, drop = T)
# test of RT- values
if (sum((subset(Gen_Data, RNA == "RT-")$mean), na.rm = T) != 0) 
	print(summary(subset(Gen_Data, RNA == "RT-")$mean)) else print("No detection in no RT samples")
# test of NTC samples
if (sum((subset(Data_Sum, Sample == "NTC")$mean), na.rm = T) != 0) 
	print(subset(Gen_Data, RNA == "RT-")) else print("No detection in NTC samples")

Ref_Data <- subset(Ref_Data, Sample != "NTC")
Gen_Data <- subset(Gen_Data, Sample != "NTC")
Gen_Data <- subset(Gen_Data, RNA == "RT+")
Comb_Data <- data.frame(Gen_Data$Sample, Gen_Data$mean, Gen_Data$Err, Ref_Data$mean, Ref_Data$Err)

Ct_Diff <- Comb_Data[2] - Comb_Data[4]
Ct_Error <- Ct_Diff + (sqrt(Comb_Data[3]^2 + Comb_Data[5]^2))


Results <- 2^-Ct_Diff
ErBar <- 2^-Ct_Error

Final_Gen_Data <- data.frame(Comb_Data[1], Results, ErBar)
Table_ColNames <- c("Sample", "Values", "Errors")
colnames(Final_Gen_Data) <- Table_ColNames

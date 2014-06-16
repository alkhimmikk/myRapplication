### - Libraries section - ####
library(ggplot2)
library(plyr)
library(reshape2)
library(plotrix)
library(grid)
library(extrafont)

### - Presets section - ####
rm(list=ls())
loadfonts(device="postscript", quiet=T)

# input data - temporarily manually changed in script
Samples<-c("UT","ZAP","CT\nDNA","poly\ndA","poly\ndT","poly\ndAdT")
SamOrder<-c("UT","ZAP","CT\nDNA","poly\ndA","poly\ndT","poly\ndAdT")
Label<-c("29-01-14 MTT Assay on HD-11 cells 1h")
Label2<-c("29-01-14 MTT Assay on HD-11 cells 3h")
GenLabel<-c("29-01-14 MTT Assay: polydA/dT/dA:dT")
File<-"14-01-30_HD-11_polydAdT"
Hours<-c(1,3) # time points

#defining plot options
TxtSize<-element_text(family = "Ubuntu Condensed", size = rel(1.7), face=2)
LegTxtSize<-element_text(family = "Ubuntu Condensed", size = rel(23), angle=0, face=2)
GridType<-element_blank()
AxLines<-element_line(size=2)
AR=2.75/length(SamOrder)

# creating derivative object by input data
RawSamples<-rep(Samples, each=2)
RealFile<-paste(File,"_1h.pda.txt",sep="")
RealFile2<-paste(File,"_3h.pda.txt",sep="")
Labels<-c(Label,Label2)
Files<-c(RealFile, RealFile2)

### - Functions section - ####

# Object - data.frame, 
#       a - # of column, 
#       No - new name of column (character)
RenameColumn<- function (Object, a, Name) 
{
RCtemp<-colnames(Object)
RCtemp[a]<-Name
colnames(Object)<-RCtemp
return(Object)
}

# File - name of file for analysis, 
#       Samples - list of treatment,
#       Order - list of treatments in desirable order,
      
TimePointFun<-function (File, Samples) 
{
        TPFV1<-head(read.table(file=File, sep="\t", fill=T), -2)
        InDeX<-which(TPFV1=="Wells", arr.ind=T)
        TPFV1<-head(read.table(file=File, sep="\t", fill=T, skip=InDeX[1]), -2)
        TPFV1[,1]<-factor(TPFV1[,1])
        TPFtemp<-c(1:ncol(TPFV1))
        for (i in 1:ncol(TPFV1)) TPFtemp[i]<-as.character(TPFV1[1,i])
        TPFtemp[3]<-"Sample#"
        TPFtemp[4]<-"Values"
        colnames(TPFV1)<-TPFtemp
        TPFV2<-tail(TPFV1,-1)
        TPFV2[,1]<-factor(TPFV2[,1])
        for (i in 2:nrow(TPFV2)) 
        {
                if (TPFV2[i,1]=="") TPFV2[i,1]<-TPFV2[i-1,1]
        }
        TPFV3<-ddply(TPFV2, .(Sample), summarize, DuplMean=mean(Values), 
                     DuplMax=max(Values), DuplMin=min(Values))
        TPFRawSamples<-rep(Samples, each=2)
        Cuvette<-rep(c("1st","2nd"), times=length(Samples))
        TPFV3$Sample<-RawSamples
        TPFV4<-cbind(TPFV3[,1], Cuvette, TPFV3[,2:4])
        TPFV4<-RenameColumn(TPFV4, 1, "Samples")
        TPFV5<-ddply(TPFV4, .(Samples), summarize, Mean=mean(DuplMean), 
                     Max=max(DuplMean), Min=min(DuplMean))
        
       TPFV4EBars<-aes(ymin=TPFV4$DuplMin, ymax=TPFV4$DuplMax)
       TPFV5EBars<-aes(ymin=TPFV5$Min, ymax=TPFV5$Max)
       TPFV4Ylim<-(ceiling(max(TPFV4$DuplMax)*10)+1)/10
       TPFV5Ylim<-(ceiling(max(TPFV5$Max)*10)+1)/10
       
       list(TPFV4,TPFV4Ylim, TPFV5,TPFV5Ylim)
}

# Datum - data frame of results,
#       ErrBars - aestetics with error bars description
#       Title - name of plot
#       YLim - max of Y-axis
DuplPlot<-function (Datum, ErrBars, YLim, Title, a)
{
        g<-ggplot(Datum, aes(x=Samples, y=DuplMean, fill=Cuvette, width=0.7 ))
        g1<-g+geom_bar(position = position_dodge(width = 0.75), stat="identity")+
                scale_fill_manual(values=c("gray90","gray20"))
        g2<-g1+scale_x_discrete(limits=SamOrder)
        g3<-g2+geom_errorbar(ErrBars, position=position_dodge(width=0.9), width=0.2, size=1.2)
        g4<-g3+scale_y_continuous(breaks=c(seq(0, YLim, by=0.2)),limits=c(0,YLim),expand=c(0,0))
        g5<-g4+xlab("")+ylab("Relative viability")+ggtitle(Title)
        g6<-g5+theme_classic()+theme(axis.title = TxtSize, axis.text = TxtSize,aspect.ratio=AR,
                                     plot.margin=unit(c(1,1,1,3), "line"), plot.title=TxtSize,
                                     axis.title.y=element_text(vjust=0), 
                                     axis.line = AxLines, axis.ticks = AxLines,
                                     panel.grid.major = GridType, panel.grid.minor = GridType)
        g7<-g6+ guides(fill = guide_legend(title=a, title.position="left", 
                                             keywidth=1.5, keyheight=1, 
                                             title.hjust=0.5, title.theme=LegTxtSize, 
                                             label.theme=LegTxtSize, label.position = "right",
                                             ncol=2))+
                theme(legend.position="top", legend.key=element_rect(colour="black", size=1.5))
        g8<-g7+geom_bar(position = position_dodge(width = 0.75),
        		colour="black", stat="identity", size=1.5, show_guide=FALSE)+
                scale_fill_manual(values=c("gray90","gray20"))+
        	geom_errorbar(ErrBars, position=position_dodge(width=0.9), width=0.2, size=1.2)
        return(g8)
}

# Datum - data frame of results,
#       ErrBars - aestetics with error bars description
#       Title - name of plot
#       YLim - max of Y-axis
PointPlot<-function (Datum, ErrBars, YLim, Title)
{
        g<-ggplot(Datum, aes(x=Samples, y=Mean, width=0.75))
        g1<-g+geom_bar(position="dodge", colour="black", stat="identity", fill="gray50", size=1.2)
        g2<-g1+scale_x_discrete(limits=SamOrder)
        g3<-g2+geom_errorbar(ErrBars, position=position_dodge(width=0.9), width=0.2, size=1.5)
        g4<-g3+scale_y_continuous(breaks=c(seq(0, YLim, by=0.2)),limits=c(0,YLim),expand=c(0,0))
        g5<-g4+xlab("")+ylab("Relative viability")+ggtitle(Title)
        g6<-g5+theme_classic()+theme(axis.title = TxtSize, axis.text = TxtSize,aspect.ratio=AR*1.4,
                                     plot.margin=unit(c(3,3,3,3), "mm"), plot.title=TxtSize,
                                     axis.title.y=element_text(vjust=0), 
                                     axis.line = AxLines, axis.ticks = AxLines,
                                     panel.grid.major = GridType, panel.grid.minor = GridType)
        return(g6)
}

#  
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y, gp=gpar(mar=c(1,1,1,1)))

# create lists-conteiners of plots and data tables
Plots<-list() 
Datum<-list()

# filing containers
for (i in c(1:length(Files)))
{
        Data<-TimePointFun(Files[i], Samples)
        for ( j in c(1:4)) Datum[[4*(i-1)+j]]<-Data[[j]]
}

###############
# PROBLEM: "i" is changed in cycles and this prohibits set up list
#	of plots with "i" as parameter
###############
# This is temporarily solution, should be changed as soon as 
# I find how to set up ggplot2-s without dependencies to "i"

for (i in c(1:(length(Datum)/4)))
{
        ifelse (i==1, 
{
        Plots[[1]]<-DuplPlot(Datum=Datum[[1]],
                             ErrBars=aes(ymin=Datum[[1]]$DuplMin,
                                         ymax=Datum[[1]]$DuplMax),
                             Datum[[2]],  Labels[1], "Cuvette")
        Plots[[2]]<-PointPlot(Datum=Datum[[3]],
                              ErrBars=aes(ymin=Datum[[3]]$Min,
                                          ymax=Datum[[3]]$Max),
                              Datum[[4]],  Labels[1])
}, 
ifelse (i==2,
{
        Plots[[3]]<-DuplPlot(Datum=Datum[[5]],
                             ErrBars=aes(ymin=Datum[[5]]$DuplMin,
                                         ymax=Datum[[5]]$DuplMax),
                             Datum[[6]],  Labels[2], "Cuvette")
        Plots[[4]]<-PointPlot(Datum=Datum[[7]],
                              ErrBars=aes(ymin=Datum[[7]]$Min,
                                          ymax=Datum[[7]]$Max),
                              Datum[[8]],  Labels[2])
},
ifelse (i==3,
{
        Plots[[5]]<-DuplPlot(Datum=Datum[[9]],
                             ErrBars=aes(ymin=Datum[[9]]$DuplMin,
                                         ymax=Datum[[9]]$DuplMax),
                             Datum[[10]],  Labels[3], "Cuvette")
        Plots[[6]]<-PointPlot(Datum=Datum[[11]],
                              ErrBars=aes(ymin=Datum[[11]]$Min,
                                          ymax=Datum[[11]]$Max),
                              Datum[[12]],  Labels[3])
},
print("Too much data, only 3 timepoint can be analysed in a row"))))
}

###############

# creating of lists contained data from duplicates & average per timepoints
DuplData<-as.list(c(1:(length(Datum)/4)))
AvData<-as.list(c(1:(length(Datum)/4)))

for (i in c(1:(length(Datum)/4))) 
{
	DuplData[[i]]<-Datum[[4*i-3]]
	AvData[[i]]<-Datum[[4*i-1]]
}

Samples<-AvData[[1]]$Samples # copy order of samples

# creating a vector with timepoint in form Nh

Times<-character() 
for (i in c(1:length(Hours)))  Times<-c(Times, paste(Hours[i],"h", sep=""))

#naming of data.frames with data
names(AvData)<-Times
names(DuplData)<-Times

### writing files with data for duplicated data - DataNh.csv
for (k in c(1:length(DuplData))) 
{
        Data1H<-as.data.frame(DuplData[[k]])
        CuvetteSplit<-split(Data1H, levels(Data1H$Cuvette))
        Exps<-levels(Data1H$Cuvette)
        for (i in c(1:length(CuvetteSplit))) CuvetteSplit[[i]]$Cuvette<-NULL
        for (i in c(1:length(CuvetteSplit))) CuvetteSplit[[i]]$Samples<-NULL
        
        for (i in c(1:length(CuvetteSplit)))  
        {
                a<-colnames(CuvetteSplit[[i]])
                for (j in c(1:length(a))) a[j]<-substring((paste(a[j],Exps[i],sep="_")),5)
                colnames(CuvetteSplit[[i]])<-a
        }
        FinData1<-as.data.frame(Samples)
        for (i in c(1:length(CuvetteSplit))) FinData1<-cbind(FinData1, CuvetteSplit[[i]])
        FinData1<-FinData1[match(SamOrder, FinData1$Samples),]
        FinData1$Samples<-gsub("\n"," ", FinData1$Samples)
        write.table(FinData1, file=paste(File,"_Values_", Times[k], ".csv", sep=""), 
                    sep=",", col.names=T, row.names=F)
}

#transforming average data for final plot

for (i in c(1:length(AvData))) AvData[[i]]$Samples<-NULL
for (i in c(1:length(AvData))) 
{
        TP<-rep(Times[i], times=length(Samples))
        AvData[[i]]<-cbind(Samples,TP, AvData[[i]])
}

PlotTable<-data.frame()
for (i in c(1:length(Samples)))
{
        Names<-Samples[i]
        for (j in c(1:length(AvData)))
        {
                PlotTable<-rbind(PlotTable, AvData[[j]][i,])
        }
}
colnames(PlotTable)<-colnames(Datum[[1]])
TotErrBar<-aes(ymin=PlotTable$DuplMin, ymax=PlotTable$DuplMax)
TotYlim<-(ceiling(max(PlotTable$DuplMax)*10)+1)/10

FinPlot<-DuplPlot(PlotTable,ErrBars=TotErrBar, YLim=TotYlim, GenLabel, a="Time :")

#transforming average data for final table
for (i in c(1:length(AvData))) AvData[[i]]$Samples<-NULL
Range<-vector()
for (i in c(1:(length(AvData)))) 
{
	Range<-(AvData[[i]]$Max-AvData[[i]]$Min)
	AvData[[i]]<-cbind(AvData[[i]],Range)
 	AvData[[i]]<-AvData[[i]][,-(3:4)]
	AvData[[i]]$TP<-NULL
}

for (i in c(1:length(AvData)))  
{
	a<-colnames(AvData[[i]])
	for (j in c(1:length(a))) a[j]<-(paste(Times[i],a[j],sep="_"))
	colnames(AvData[[i]])<-a
}
FinalData<-as.data.frame(Samples)
for (i in c(1:length(AvData)))
{
	FinalData<-cbind(FinalData, AvData[[i]])
}
FinalData<-FinalData[match(SamOrder, FinalData$Samples),]
FinalData$Samples<-gsub("\n"," ", FinalData$Samples)
write.table(FinalData, file=paste(File, "_Values.csv", sep=""), 
	    sep=",", col.names=T, row.names=F)


############################## Before plotting:
# 	object FinPlot contains final plot of both (?) tp
# 	object Plots contains 
# 		plots of 2 cuvettes in odd position 
# 		and average data in even position

########################
postscript(fonts="Ubuntu Condensed")
unlink("Rplots.ps")
for (j in c(1:(length(Plots)/2)))
{
	dev.new()
grid.newpage()
pushViewport(viewport(layout = grid.layout(2,1)))

print(Plots[[2*j-1]], vp=vplayout(1,1))
	print(Plots[[2*j]], vp=vplayout(2,1))

postscript(file=paste(File,"_", Times[j],".ps",sep=""), family="Ubuntu Condensed")
pushViewport(viewport(layout = grid.layout(2,1)))
	print(Plots[[2*j-1]], vp=vplayout(1,1))
	print(Plots[[2*j]], vp=vplayout(2,1))
dev.off()
}	     


dev.new()
print(FinPlot)
	
postscript(file=paste(File, "_Summary.ps",sep=""), family="Ubuntu Condensed")
print(FinPlot)
dev.off()
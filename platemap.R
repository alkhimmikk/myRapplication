library(dplyr)
library(ggplot2)

TxtSize<-element_text( size = rel(1.0), face=2, colour="black")
LegTxtSize<-element_text(size = rel(23), angle=0, face=2)
GridType<-element_blank()
AxLines<-element_line(size=.1)

platemap <- read.csv("SamplesForTest//Plate.csv")

platemap <- mutate(platemap,
		   Row=as.numeric(match(toupper(substr(Position, 1, 1)), LETTERS)),
		   Column=as.numeric(substr(Position, 2, 5)))

platemap$Samples<-factor(platemap$Sample, levels=unique(as.character(platemap$Sample)))

	
# ggplot(data=platemap, aes(x=Column, y=Row)) +
# 	geom_point(size=10) +
# 	labs(title="Test platemap")

g<-ggplot(data=platemap, aes(x=Column, y=Row, fill=Target) )+
	geom_point(data=expand.grid(seq(1, 24), seq(1, 16)), aes(x=Var1, y=Var2),
		   color="black", fill="white", shape=22, size=11) +
	geom_point(shape=22, size=10)+
	scale_fill_manual(name="Primers",
			  values=c("gray90", "gray75", "gray60", "gray45"),
			  labels=levels(platemap$Target))+
	geom_point(aes(colour=cDNA, shape=Samples), size=4, fill="grey15")+
	scale_shape_manual(name="Treatment",
			   values=LETTERS[1:19],
			   labels=levels(platemap$Samples))+
	scale_colour_manual(name="cDNA", values=c("white","black"), labels=c("RT+","RT-"))+
	scale_y_reverse(breaks=seq(1, 16), labels=LETTERS[1:16]) +
	scale_x_continuous(breaks=seq(1, 24)) +
	coord_fixed(ratio=(25/24)/(17/16), xlim = c(0.5, 24.5), ylim=c(0.5, 16.5))+
	xlab("Columns")+ylab("Rows")+ggtitle("Platemap")+
	theme(axis.title=TxtSize, axis.text=TxtSize, plot.title=TxtSize,
	      axis.title.y=element_text(vjust=0),axis.line = AxLines, axis.ticks = AxLines)+
	guides(fill = guide_legend(title.position="left", 
				   keywidth=1.5, keyheight=1, 
				   title.hjust=0.5, label.position = "right",
				   nrow=1))+
	guides(colour = guide_legend(legend.positon="bottom",title.position="left", 
				     keywidth=1.5, keyheight=1, 
				     title.hjust=0.5, label.position = "right",
				     ncol=2))+
	guides(shape = guide_legend(legend.positon="bottom",title.position="left", 
				    keywidth=1.5, keyheight=1, 
				    title.hjust=0.5, label.position = "right",
				    ncol=10, byrow=T))+
	theme(legend.position="bottom")
dev.new()
print(g)

postscript("14-05-30_Platemap.ps")
print(g)
dev.off()

system("gv 14-05-30_Platemap.ps")

# print(g)
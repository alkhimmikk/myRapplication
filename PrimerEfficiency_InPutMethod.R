# x<-c(1,.1,.01,.001,0.0001)
x<-c(1,.1,.01)
y1<-c(32.028,	35.164,	37.310)
y2<-c(31.897,	34.127,	36.916)	
Ymean<-apply(rbind(y1,y2),2,mean)
	
a<-data.frame(x,y1,y2, Ymean)


b<-lm(formula=Ymean~log10(x), data=a)
d<-coefficients(b)
Slope<-d[2]
AmplEff<-(10^(-1*1/Slope)-1)*100
PowOfAmpl<-2*AmplEff/100


dev.new()
plot(Ymean~x,a, type="b", log="x", col="grey20", lwd=4 )
abline(b, col="black", lwd=5)
text(x=x[3]*2, y=mean(a$Ymean), paste("Slope=", round(d[2],3), sep=""))
text(x=x[3]*2, y=mean(a$Ymean)-.5, paste("Power of\n Amplification=", round(PowOfAmpl,3), sep=""))

svg("STING-efficiency.svg")
plot(Ymean~x,a, type="b", log="x", col="grey20", lwd=4 )
abline(b, col="black", lwd=5)
text(x=x[3]*2, y=mean(a$Ymean), paste("Slope=", round(d[2],3), sep=""))
text(x=x[3]*2, y=mean(a$Ymean)-.5, paste("Power of\n Amplification=", round(PowOfAmpl,3), sep=""))
dev.off()
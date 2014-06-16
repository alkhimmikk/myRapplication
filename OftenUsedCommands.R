options(device="X11")
if (length(dev.list()>0))(for (j in (1:(length(dev.list())))) dev.off())
rm(list=ls())


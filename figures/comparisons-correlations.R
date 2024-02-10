# Clinico-pathological analysis

# Function to plot relation
two_groups_plot = function(data_min, data_plus, file, ymax, ystep, ylab, yline, sign, ylabels=NA, mar=c(6,7,4,3)) {
  
  df = data.frame(rbind(cbind(data_min, rep("min", length(data_min))), cbind(data_plus, rep("plus", length(data_plus)))))
  colnames(df)= c("data", "group")
  
  data = ggplot(df, aes(x=factor(group), y=as.numeric(data), fill=factor(group))) +
    geom_boxplot() +
    geom_dotplot(binaxis='y', stackdir='center', dotsize=1)
  pos = ggplot_build(data)$data[[2]]
  data

  png(file, res=600, 3300, 2300)  
  par(mar=mar)
  xlim = c(0,4)
  ylim = c(0,ymax)
  
  plot(xlim, 
       ylim, 
       type = "n", 
       axes = F, 
       xlab = "",
       ylab = "",
       xaxs = "i", 
       yaxs = "i")
  xlim[2] = 3
  xat = seq(xlim[1], xlim[2], by=1)
  yat = seq(ylim[1], ylim[2], by=ystep)
  #ylim[2]=55
  #axis(side = 1, at = xat, labels = xat, col = "#b1b1b1", lwd = 1.4, col.axis="#333333")
  if (is.na(ylabels)) {
    ylabels = paste0(yat,"")
  }
  axis(side = 2, at = yat, labels = ylabels, las=2, col = "#b1b1b1", lwd = 1.4, col.axis="#333333",cex.axis=1.1)
  #mtext(side=2,"Concentration",line=4.55,xpd=T, col="#333333")
  mtext(side=2,ylab,line=yline,xpd=T, col="#333333", cex=1.1)
  
  #segments(xat, ylim[1], xat, ylim[2], col="#eeeeee", lwd=1.4, xpd=T)
  segments(xlim[1], yat, xlim[2], col="#eeeeee", lwd=1.4, xpd=T)
  segments(xlim[1], ylim[1], xlim[1], ylim[2], col="#b1b1b1", lwd=1.4, xpd=T)
  
  #dist = table(vitreous_data$GEP_class[which(vitreous_data$Classification %in% c("vfDNA-/UM-","vfDNA+/UM-"))])
  #rect(0.25, 0, 0.75, bruch_min*100, border=NA, col=RColorBrewer::brewer.pal(10,"Paired")[10])
  #rect(0.25, dist[2]/sum(dist)*100, 0.75, 100, border=NA, col=RColorBrewer::brewer.pal(10,"Paired")[9])
  text(0.8, -ylim[2]/10, col="#333333", labels=paste0("vfDNA+/UM\u2212"), xpd=T, cex=1.1)
  text(0.8, -ylim[2]/10*2.1, col="#333333", labels=paste0("(n=",length(data_min),")"), xpd=T, cex=1.1)
  
  #dist = table(vitreous_data$GEP_class[which(vitreous_data$Classification == "vfDNA+/UM+")])
  #rect(1, 0, 1.5, bruch_plus*100, border=NA, col=RColorBrewer::brewer.pal(10,"Paired")[10])
  #rect(1, dist[2]/sum(dist)*100, 1.5, 100, border=NA, col=RColorBrewer::brewer.pal(10,"Paired")[9])
  text(2.2, -ylim[2]/10, col="#333333", labels="vfDNA+/UM+", xpd=T, cex=1.1)
  text(2.2, -ylim[2]/10*2.1, col="#333333", labels=paste0("(n=",length(data_plus),")"), xpd=T, cex=1.1)
  
  color = c("#F44336","#2196F3") #RColorBrewer::brewer.pal(9,"Paired")[c(6,2)]
  pos$group[which(pos$group == 1)] = 0.8
  pos$group[which(pos$group == 2)] = 2.2
  points(pos$group+pos$stackpos/10, pos$y, pch=16, col=color[round(pos$group)])
  m = median(data_min, na.rm=T)
  segments(0.8-0.4,m,0.8+0.4,m, lwd=2.8)
  
  m = median(data_plus, na.rm = T)
  segments(2.2-0.4,m,2.2+0.4,m, lwd=2.8)
  
  max_1 = max(data_min, na.rm=T)
  max_2 = max(data_plus, na.rm=T)
  max = max(c(max_1,max_2))+(ylim[2]-ylim[1])/10*1.5
  segments(0.8,max,2.2,lwd=1.4,col="#333333",xpd=T)
  segments(0.8,max_1+(ylim[2]-ylim[1])/10,0.8,max,lwd=1.4,col="#333333",xpd=T)
  segments(2.2,max_2+(ylim[2]-ylim[1])/10,2.2,max,lwd=1.4,col="#333333",xpd=T)
  text(1.5, max+ylim[2]/10, col="#333333", labels=sign, xpd=T, cex=1.1)
  print(sign)
  segments(xlim[1], ylim[1], xlim[2], ylim[1], col="#b1b1b1", lwd=1.4, xpd=T)
  
  
  system("taskkill /F /IM Microsoft.Photos.exe /T")
  dev.off()
  system(paste("open",file))
}

# Load libraries
library(ggplot2)
library(stringr)
library(readxl)

# Load data
vitreous_data = readxl::read_xlsx("data/data - public.xlsx", skip=0)[1:65,]

# Set groups
um_min = which(vitreous_data$vfDNA_classification == "vfDNA+/UM-")
um_plus = which(vitreous_data$vfDNA_classification == "vfDNA+/UM+")

# Prominence
data_min = as.numeric(vitreous_data$prominence[um_min])
data_min = data_min[!is.na(data_min)]
data_plus = as.numeric(vitreous_data$prominence[um_plus])
data_plus = data_plus[!is.na(data_plus)]
median(data_min)
median(data_plus)
wilcox.test(data_min, data_plus)
max(c(data_min,data_plus))
two_groups_plot(data_min, data_plus, "figures/prominence.png", ymax=15, ystep=5, ylab="Prominence (mm)", yline=2.75, sign = expression(italic(p)~"= 0.018"))

# Diameter
data_min = as.numeric(vitreous_data$diameter[um_min])
data_min = data_min[!is.na(data_min)]
data_plus = as.numeric(vitreous_data$diameter[um_plus])
data_plus = data_plus[!is.na(data_plus)]
median(data_min)
median(data_plus)
wilcox.test(data_min, data_plus)
max(c(data_min,data_plus))
two_groups_plot(data_min, data_plus, "figures/diameter.png", ymax=25, ystep=5, ylab="Diameter (mm)", yline=2.75, sign = "n/s")

# Total concentration vfDNA in vfDNA+/UM- vs. vfDNA+/UM+
data_min = log10(vitreous_data$total_copies_per_uL_vf[um_min])+1
data_plus = log10(vitreous_data$total_copies_per_uL_vf[um_plus])+1
wilcox.test(data_min, data_plus)
max(c(data_min,data_plus))
two_groups_plot(data_min, data_plus, "figures/vfDNA.png", ymax=6, ystep=1, ylab="Total vfDNA\nconcentration\n(copies/uL)", yline=3.75, sign = expression(italic(p)~"= 0.002"), ylabels = c(0.01,0.1,1,10,100,1000,10000))

# Total concentration non-melanoma cell derived vfDNA in vfDNA+/UM- vs. vfDNA+/UM+
data_min = log10(vitreous_data$total_copies_per_uL_vf[um_min]*(1-as.numeric(vitreous_data$Gaq_melanoma_cell_derived_fraction_vfDNA[um_min])))+1
data_plus = log10(vitreous_data$total_copies_per_uL_vf[um_plus]*(1-as.numeric(vitreous_data$Gaq_melanoma_cell_derived_fraction_vfDNA[um_plus])))+1
wilcox.test(data_min, data_plus)
max(c(data_min,data_plus))
two_groups_plot(data_min, data_plus, "figures/vfDNA_melanoma.png", ymax=6, ystep=1, ylab="Non-melanoma\ncell-derived vfDNA\nconcentration\n(copies/uL)", yline=3.75, sign = expression(italic(p)~"= 0.011"), ylabels = c(0.01,0.1,1,10,100,1000,10000), mar = c(6,9,4,1))

# Correlation between melanoma-cell derived vfDNA and tumour prominence
y = log10(as.numeric(vitreous_data$mutant_copies_per_ul_vf[um_plus])*10)
x = as.numeric(vitreous_data$prominence[um_plus])/15
cor.test(x,y, method = "spearman")
file = "figures/correlation.png"
png(file, res=600, 3300, 2300)  
par(mar=c(5,8,5,8))
xlim = c(0,1)
ylim = c(0,5)
plot(xlim, 
     ylim, 
     type = "n", 
     axes = F, 
     xlab = "",
     ylab = "",
     xaxs = "i", 
     yaxs = "i")
xat = seq(xlim[1], xlim[2], by=1/3)
yat = seq(ylim[1], ylim[2], by=1)
axis(side = 1, at = xat[c(1,3,5)], labels = paste0(xat*3*5,"")[c(1,3,5)], col = "#b1b1b1", lwd = 1.4, col.axis="#333333", cex.axis=1.1)
axis(side = 1, at = xat[c(2,4)], labels = paste0(xat*3*5,"")[c(2,4)], col = "#b1b1b1", lwd = 1.4, col.axis="#333333", cex.axis=1.1)
axis(side = 2, at = yat, labels = c(0.1,1,10,100,1000,10000), las=2, col = "#b1b1b1", lwd = 1.4, col.axis="#333333", cex.axis=1.1)
mtext(side=2,"Melanoma-cell\nderived vfDNA\nconcentration\n(copies/uL)",line=3.75,xpd=T, col="#333333", cex=1.1)
mtext(side=1,"Prominence (mm)",line=2.75,xpd=T, col="#333333", cex=1.1)
segments(xat, ylim[1], xat, ylim[2], col="#eeeeee", lwd=1.4, xpd=T)
segments(xlim[1], yat, xlim[2], col="#eeeeee", lwd=1.4, xpd=T)
segments(xlim[1], ylim[1], xlim[1], ylim[2], col="#b1b1b1", lwd=1.4, xpd=T)
segments(xlim[1], ylim[1], xlim[2], ylim[1], col="#b1b1b1", lwd=1.4, xpd=T)
points(x,y, pch=16, cex=1, col="#2196F3",xpd=T)
corr = lm(y~x)
segments(min(x),corr$coefficients[1]+corr$coefficients[2]*min(x),max(x),corr$coefficients[1]+corr$coefficients[2]*max(x),lwd=1.4,col="#333333",xpd=T)
text(1.05,2.25,labels=substitute(italic("p")~"< 0.001"), xpd=T, col="#333333",pos=4, cex=1.1)
text(1.05,2.75,labels=substitute(italic("rho")~"= 0.52"), xpd=T, col="#333333",pos=4, cex=1.1)
dev.off()
system(paste("open",file))

# Function to plot relation
two_bars_plot = function(prop_min, prop_plus, length_min, length_plus, file, ylab, yline, sign, ylabels=NA) {
  
  png(file, res=600, 3300, 2300)  
  par(mar=c(6,7,4,3))
  xlim = c(0,4)
  ylim = c(0,100)
  
  plot(xlim, 
       ylim, 
       type = "n", 
       axes = F, 
       xlab = "",
       ylab = "",
       xaxs = "i", 
       yaxs = "i")
  xlim[2] = 3
  xat = seq(xlim[1], xlim[2], by=1)
  yat = seq(ylim[1], ylim[2], by=25)
  #ylim[2]=55
  #axis(side = 1, at = xat, labels = xat, col = "#b1b1b1", lwd = 1.4, col.axis="#333333")
  if (is.na(ylabels)) {
    ylabels = paste0(yat,"%")
  }
  axis(side = 2, at = yat, labels = ylabels, las=2, col = "#b1b1b1", lwd = 1.4, col.axis="#333333",cex.axis=1.1)
  mtext(side=2,ylab,line=yline,xpd=T, col="#333333",cex=1.1)
  segments(xlim[1], yat, xlim[2], col="#eeeeee", lwd=1.4, xpd=T)
  segments(xlim[1], ylim[1], xlim[1], ylim[2], col="#b1b1b1", lwd=1.4, xpd=T)
  color = c("#F44336","#2196F3")
  rect(0.8-0.4, 0, 0.8+0.4, prop_min*100, border=NA, col=color[1])
  text(0.8, -ylim[2]/10, col="#333333", labels=paste0("vfDNA+/UM\u2212"), xpd=T, cex=1.1)
  text(0.8, -ylim[2]/10*2.1, col="#333333", labels=paste0("(n=",length_min,")"), xpd=T, cex=1.1)
  rect(2.2-0.4, 0, 2.2+0.4, prop_plus*100, border=NA, col=color[2])
  text(2.2, -ylim[2]/10, col="#333333", labels="vfDNA+/UM+", xpd=T, cex=1.1)
  text(2.2, -ylim[2]/10*2.1, col="#333333", labels=paste0("(n=",length_plus,")"), xpd=T, cex=1.1)
  
  max_1 = max(prop_min*100, na.rm=T)
  max_2 = max(prop_plus*100, na.rm=T)
  max = max(c(max_1,max_2))+(ylim[2]-ylim[1])/10*1.5
  segments(0.8,max,2.2,lwd=1.4,col="#333333")
  segments(0.8,max_1+(ylim[2]-ylim[1])/10,0.8,max,lwd=1.4,col="#333333")
  segments(2.2,max_2+(ylim[2]-ylim[1])/10,2.2,max,lwd=1.4,col="#333333")
  text(1.5, max+ylim[2]/10, col="#333333", labels=sign, xpd=T, cex=1.1)
  
  segments(xlim[1], ylim[1], xlim[2], ylim[1], col="#b1b1b1", lwd=1.4, xpd=T)
  
  system("taskkill /F /IM Microsoft.Photos.exe /T")
  dev.off()
  system(paste("open",file))
}

# Bruch:
#-: 13x0 : 6x1
#+ 15x 0 : 17x 1
two_bars_plot(6/19, 17/32, 19, 32, "figures/bruchs_membrane.png", ylab="Broken through\nBruch's membrane", yline=3.75, sign = "n/s")

# Classes of primary tumours
#clI 9x : clII 15x
#clI 15x : clII 24x
prop_min = 15/24
prop_plus = 24/39
length_min = 24
length_plus = 39
file = "figures/class_II.png"
yline=3.75
sign = "n/s"
ylabels=NA
png(file, res=600, 3300, 2300)  
par(mar=c(6,7,4,3))
xlim = c(0,4)
ylim = c(0,100)
plot(xlim, 
     ylim, 
     type = "n", 
     axes = F, 
     xlab = "",
     ylab = "",
     xaxs = "i", 
     yaxs = "i")
xlim[2] = 3
xat = seq(xlim[1], xlim[2], by=1)
yat = seq(ylim[1], ylim[2], by=25)
if (is.na(ylabels)) {
  ylabels = paste0(yat,"%")
}
axis(side = 2, at = yat, labels = ylabels, las=2, col = "#b1b1b1", lwd = 1.4, col.axis="#333333",cex.axis=1.1)
mtext(side=2,substitute(italic("BAP1")~"mutation"),line=yline+2.2,xpd=T, col="#333333",cex=1.1)
mtext(side=2,"and/or monosomy 3p\nin primary tumour",line=yline,xpd=T, col="#333333",cex=1.1)
segments(xlim[1], yat, xlim[2], col="#eeeeee", lwd=1.4, xpd=T)
segments(xlim[1], ylim[1], xlim[1], ylim[2], col="#b1b1b1", lwd=1.4, xpd=T)
color = c("#F44336","#2196F3")
rect(0.8-0.4, 0, 0.8+0.4, prop_min*100, border=NA, col=color[1])
text(0.8, -ylim[2]/10, col="#333333", labels=paste0("vfDNA+/UM\u2212"), xpd=T, cex=1.1)
text(0.8, -ylim[2]/10*2.1, col="#333333", labels=paste0("(n=",length_min,")"), xpd=T, cex=1.1)
rect(2.2-0.4, 0, 2.2+0.4, prop_plus*100, border=NA, col=color[2])
text(2.2, -ylim[2]/10, col="#333333", labels="vfDNA+/UM+", xpd=T, cex=1.1)
text(2.2, -ylim[2]/10*2.1, col="#333333", labels=paste0("(n=",length_plus,")"), xpd=T, cex=1.1)
max_1 = max(prop_min*100, na.rm=T)
max_2 = max(prop_plus*100, na.rm=T)
max = max(c(max_1,max_2))+(ylim[2]-ylim[1])/10*1.5
segments(0.8,max,2.2,lwd=1.4,col="#333333")
segments(0.8,max_1+(ylim[2]-ylim[1])/10,0.8,max,lwd=1.4,col="#333333")
segments(2.2,max_2+(ylim[2]-ylim[1])/10,2.2,max,lwd=1.4,col="#333333")
text(1.5, max+ylim[2]/10, col="#333333", labels=sign, xpd=T, cex=1.1)
segments(xlim[1], ylim[1], xlim[2], ylim[1], col="#b1b1b1", lwd=1.4, xpd=T)
system("taskkill /F /IM Microsoft.Photos.exe /T")
dev.off()
system(paste("open",file))
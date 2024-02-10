# Close current images
system("taskkill /F /IM Microsoft.Photos.exe /T")
# Init library
library(stringr)
# Load data
vitreous_data = readxl::read_xlsx("data/data - public.xlsx")
# Order data
vitreous_data = vitreous_data[order(vitreous_data$order),]
# Init PNG
graphics.off()

png("figures/figure-2.png",res=600,width=8000,height=6000)
par(mar=c(6,8,6,4))
plot(c(0,75),c(-10,10),type="n",axes=F,ylab="",xlab="",xaxs = "i", yaxs = "i",)
selection = 1:65
segments(0,5+(0:5)/5*3,77,col="#eeeeee",xpd=T,lwd=1.4)
x = 1.5
for (i in selection) {
  x = x+1
  if (i %in% c(3,27)) {
    x=x+4
  }
  if (vitreous_data$total_copies_per_uL_vf[i]==0) {
    text(x,5.1, labels="no vfDNA",srt=90,cex=0.8,adj=0,col="#333333")
  }
  else {
    rect(x-0.4, 5, x+0.4, 5.6+log10(as.numeric(vitreous_data$total_copies_per_uL_vf[i]))/5*3, border=NA, col="#777777")
  }
}
axis(side=2,at=5+(0:5)/5*3,labels=c(0.1,1,10,100,1000,10000),las=2,col="#b1b1b1",cex.axis=1.1,lwd=1.4,col.axis="#333333")
mtext(side=2,"Total vfDNA",line=6.0,xpd=T,at=6.5,cex=1.1,cex.axis=1.1,col.axis="#333333")
mtext(side=2,"concentration",line=5.2,xpd=T,at=6.5,cex=1.1,cex.axis=1.1,col.axis="#333333")
mtext(side=2,"(copies/uL)",line=4.4,xpd=T,at=6.5,cex=1.1,cex.axis=1.1,col.axis="#333333")

# Labels on top of plot for three groups
text(3,12.4-1.5, labels="", cex=1.1,col="#333333",xpd=T)
text(3,11.9-1.5, labels="vfDNA\u2212/UM\u2212", cex=1.1,col="#333333",xpd=T)
text(3,11.4-1.5, labels="(n=2)", cex=1.1,col="#333333",xpd=T)
text(20,12.4-1.5, labels="", cex=1.1,col="#333333",xpd=T)
text(20,11.9-1.5, labels="vfDNA+/UM\u2212", cex=1.1,col="#333333",xpd=T)
text(20,11.4-1.5, labels="(n=24)", cex=1.1,col="#333333",xpd=T)
text(55.5,12.4-1.5, labels="", cex=1.1,col="#333333",xpd=T)
text(55.5,11.9-1.5, labels="vfDNA+/UM+", cex=1.1,col="#333333",xpd=T)
text(55.5,11.4-1.5, labels="(n=39)", cex=1.1,col="#333333",xpd=T)


acco = function(x1,x2, y1,y2) {
  
  segments(x1, y2, x1, mean(c(y1,y2)), col='#b1b1b1', lwd=1.4, xpd=T)
  segments(x2, y2, x2, mean(c(y1,y2)), col='#b1b1b1', lwd=1.4, xpd=T)
  segments(x1, mean(c(y1,y2)), x2, col='#b1b1b1', lwd=1.4, xpd=T)
  segments(mean(c(x1,x2)), mean(c(y1,y2)), mean(c(x1,x2)), y1, col='#b1b1b1', lwd=1.4, xpd=T)
}

acco (2.5-0.4, 3.5+0.4, 9.4,8.8)
acco (8.5-0.4, 31.5+0.4, 9.4,8.8)
acco (36.5-0.4, 74.5+0.4, 9.4,8.8)

# Fraction cancer-cell derived
ys = -2
segments(34,-2+(0:4)/4*3,77,col="#eeeeee",xpd=T,lwd=1.4)
x = 1.5
for (i in selection) {
  x = x+1
  if (i %in% c(3,27)) {
    x=x+4
  }
  rect(x-0.4, 0+ys, x+0.4, 0+ys+3*as.numeric(vitreous_data$Gaq_melanoma_cell_derived_fraction_vfDNA[i]), border=NA, col="#777777")
}
axis(side=2,at=-2+(0:4)/4*3,pos = 34,las=2,labels=c("0%","25%","50%","75%","100%"),col="#b1b1b1",cex.axis=1.1,lwd=1.4,col.axis="#333333")
xx=-24.5
mtext(side=2,"Fraction",line=5.2+xx,xpd=T,at = -.5, cex=1.1,cex.axis=1.1,col.axis="#333333")
mtext(side=2,"melanoma-cell",line=4.5+xx,xpd=T, at =-.5, cex=1.1,cex.axis=1.1,col.axis="#333333")
mtext(side=2,"derived",line=3.8+xx,xpd=T, at =-.5, cex=1.1,cex.axis=1.1,col.axis="#333333")
axis(side=2,at=5-0.2-1+c(0,-0.5,-1,-1.5),las=2,labels=c("GNAQ","GNA11","CYSLTR2","PLCB4"),font=3,col="#b1b1b1",cex.axis=1.1,lwd=1.4,col.axis="#333333")
axis(side=2,pos = 34,at=2.5-0.2-7+c(0,-0.5,-1)+1.5,las=2,labels=c("EIF1AX","SF3B1","BAP1"),font=3,col="#b1b1b1",cex.axis=1.1,lwd=1.4,col.axis="#333333")
axis(side=2,pos = 34,at=0.5-0.2-7+c(0,-0.5)+1.5-.5,las=2,labels=c("Chr. 3p","Chr. 8q"),font=1,col="#b1b1b1",cex.axis=1.1,lwd=1.4,col.axis="#333333")

color_check = function(val) {
  if (!is.na(val)) {
    if (str_detect(val,"not detected")) {
      col="#F44336" #col=RColorBrewer::brewer.pal(9,"Paired")[6]
    }
    else if (str_detect(val,"no CNA detected")) {
      col="#F44336" #RColorBrewer::brewer.pal(9,"Paired")[6]
    }
    else if (str_detect(val,"detected")) {
      col="#2196F3" #RColorBrewer::brewer.pal(9,"Paired")[2]
    }
    else if (str_detect(val,"no assay")) {
      col="#CCCCCC"
    }
    else if (str_detect(val,"not measured")) {
      col=RColorBrewer::brewer.pal(9,"Paired")[1]
    }
    else if (str_detect(val,"underpowered")) {
      col="#ffffff" #RColorBrewer::brewer.pal(9,"Paired")[5] #5
    }
  }
  else {
    col="#CCCCCC"
  }
  return(col)
}
# Tumour characteristics
# Start x
x = 1.5
# Iterate through samples
for (i in selection) {
  # Determine x position
  x = x+1
  if (i %in% c(3,27)) {
    x=x+4
  }
  # Start y
  y = 5-1
  # GNAQ
  col="#EEEEEE"
  if (str_detect(vitreous_data$Gaq_mutation_tDNA[i],"GNAQ")) {
    if (vitreous_data$mutant_copies_per_ul_vf[i] == "n/a") {
      #col = color_check("not detected")
    }
    else if (vitreous_data$mutant_copies_per_ul_vf[i] != "0") {
      col = color_check("detected")
    }
    else {
      col = color_check("not detected")
    }
  }
  rect(x-0.4, y, x+0.4, y-.4, border=NA, col=col)
  y=y-0.5
  
  # GNA11
  col="#EEEEEE"
  if (str_detect(vitreous_data$Gaq_mutation_tDNA[i],"GNA11")) {
    if (vitreous_data$mutant_copies_per_ul_vf[i] == "n/a") {
      #col = color_check("not detected")
    }
    else if (vitreous_data$mutant_copies_per_ul_vf[i] != "0") {
      col = color_check("detected")
    }
    else {
      col = color_check("not detected")
    }
  }
  rect(x-0.4, y, x+0.4, y-.4, border=NA, col=col)
  y=y-0.5
  
  # CYSLTR2
  col="#EEEEEE"
  if (str_detect(vitreous_data$Gaq_mutation_tDNA[i],"CYSLTR2")) {
    if (vitreous_data$mutant_copies_per_ul_vf[i] == "n/a") {
      #col = color_check("not detected")
    }
    else if (vitreous_data$mutant_copies_per_ul_vf[i] != "0") {
      col = color_check("detected")
    }
    else {
      col = color_check("not detected")
    }
  }
  rect(x-0.4, y, x+0.4, y-.4, border=NA, col=col)
  y=y-0.5
  
  # PLCB4
  col="#EEEEEE"
  if (str_detect(vitreous_data$Gaq_mutation_tDNA[i],"PLCB4")) {
    if (vitreous_data$mutant_copies_per_ul_vf[i] == "n/a") {
      #col = color_check("not detected")
    }
    else if (vitreous_data$mutant_copies_per_ul_vf[i] != "0") {
      col = color_check("detected")
    }
    else {
      col = color_check("not detected")
    }
  }
  rect(x-0.4, y, x+0.4, y-.4, border=NA, col=col)
  y=y-0.5
}
x = 1.5+34
for (i in 27:66) {
  x = x+1
  y = 5-8-1.5+1.5
  # EIF1AX
  col = "#EEEEEE"
  if (!is.na(vitreous_data$EIF1AX_mutation_tDNA[i])) {
    col = color_check(vitreous_data$EIF1AX_mutation_vfDNA[i])
  }
  rect(x-0.4, y, x+0.4, y-.4, border=NA, col=col)
  y=y-0.5
  col = "#EEEEEE"
  if (!is.na(vitreous_data$SF3B1_mutation_tDNA[i])) {
    col = color_check(vitreous_data$SF3B1_mutation_vfDNA[i])
  }
  rect(x-0.4, y, x+0.4, y-.4, border=NA, col=col)
  #text(0,y-0.2,labels="GNA11",pos=2,xpd=T,font=3,cex=0.8)
  y=y-0.5
  col = "#EEEEEE"
  if (!is.na(vitreous_data$BAP1_mutation_tDNA[i])) {
    col = color_check(vitreous_data$BAP1_mutation_vfDNA[i])
    if (vitreous_data$BAP1_mutation_tDNA[i] == "unknown") {
      col="#CCCCCC"
    }
  }
  rect(x-0.4, y, x+0.4, y-.4, border=NA, col=col)
  y=y-1.5
  col = "#EEEEEE"
  col = color_check(vitreous_data$chr3p_vfDNA[i])
  if (vitreous_data$chr3p_match[i] == "match" & vitreous_data$chr3p_vfDNA[i] == "no CNA detected") {
    col = "#EEEEEE"
  }
  rect(x-0.4, y, x+0.4, y-.4, border=NA, col=col)
  y=y-0.5
  col = "#EEEEEE"
  col = color_check(vitreous_data$chr8q_vfDNA[i])
  if (vitreous_data$chr8q_match[i] == "match" & vitreous_data$chr8q_vfDNA[i] == "no CNA detected") {
    col = "#EEEEEE"
  }
  rect(x-0.4, y, x+0.4, y-.4, border=NA, col=col)
  #text(x, y-1.5, srt=90, labels=vitreous_data$vfDNA_ID[i], cex=0.5)
}
x=0.5
y=0+0.5
text(x-1.25,y,labels="Mutations and CNAs",pos=4,xpd=T,font=2,cex=1.1, col='#333333')
y=y-0.5
rect(x-0.4, y, x+0.4, y-.4, border=NA, col=color_check("detected"))
text(x+0.3,y-0.2,labels="detected",pos=4,xpd=T,font=1,cex=1.1, col='#333333')
y=y-0.5
rect(x-0.4, y, x+0.4, y-.4, border=NA, col=color_check("not detected"))
text(x+0.3,y-0.2,labels="not detected",pos=4,xpd=T,font=1,cex=1.1, col='#333333')
y=y-0.5
rect(x-0.4, y, x+0.4, y-.4, border=NA, col=color_check("no assay"))
text(x+0.3,y-0.2,labels="no assay",pos=4,xpd=T,font=1,cex=1.1, col='#333333')
y=y-0.5
rect(x-0.4, y, x+0.4, y-.4, border=NA, col="#EEEEEE")
text(x+0.3,y-0.2,labels="no alteration",pos=4,xpd=T,font=1,cex=1.1, col='#333333')

y=y-1.5
#text(x-1.4,y,labels="CNAs",pos=4,xpd=T,font=2,cex=0.8)
y=y-0.5
#rect(x-0.4, y, x+0.4, y-.4, border=NA, col="#EEEEEE")
#text(x,y-0.2,labels="no CNA",pos=4,xpd=T,font=1,cex=0.8)
y=y-0.5
#rect(x-0.4, y, x+0.4, y-.4, border=NA, col=color_check("detected"))
#text(x,y-0.2,labels="detected",pos=4,xpd=T,font=1,cex=0.8)
y=y-0.5
#rect(x-0.4, y, x+0.4, y-.4, border=NA, col=color_check("not detected"))
#text(x,y-0.2,labels="not detected",pos=4,xpd=T,font=1,cex=0.8)
y=y-0.5
#rect(x-0.4, y, x+0.4, y-.4, border=NA, col=color_check("underpowered"))
#text(x,y-0.2,labels="underpowered",pos=4,xpd=T,font=1,cex=0.8)
y=y-0.5
#rect(x-0.4, y, x+0.4, y-.4, border=NA, col=color_check("not measured"))
#text(x,y-0.2,labels="not measured",pos=4,xpd=T,font=1,cex=0.8)


segments(0,4,0,2.1,col="#b1b1b1",lty=1,xpd=T,lwd=1.4)
segments(6,c(5,4),6,c(8,2.1),col="#b1b1b1",lty=1,lwd=1.4)
segments(34,c(5,4),34,c(8,2.1),col="#b1b1b1",lty=1,lwd=1.4)
segments(77,c(5,4,1,-3,-5-.5),77,c(8,2.1,-2,-4.4,-5.9-.5),col="#b1b1b1",lty=1,xpd=T,lwd=1.4)
segments(34,c(5,4,1,-3,-5-.5),34,c(8,2.1,-2,-4.4,-5.9-.5),col="#b1b1b1",lty=1,xpd=T,lwd=1.4)

segments(36.2,-5.2-.5,54.8-1,col="#b1b1b1",lty=1,xpd=T,lwd=1.4)
rect(49.8-.5,-5.4-.5,41.2-.5,-5.0-.5,col="#FFFFFF",border=NA,xpd=T)
text(mean(c(36.1,53.9)),-5.2-.5,labels="underpowered",cex=1.1,col="#333333",font=3)

segments(36.2,-5.7-.5,55.8-1,col="#b1b1b1",lty=1,xpd=T,lwd=1.4)
rect(49.8-.5,-5.6-.5,41.2-.5,-5.9-.5,col="#FFFFFF",border=NA,xpd=T)
text(mean(c(36.1,53.9)),-5.7-.5,labels="underpowered",cex=1.1,col="#333333",font=3)

dev.off()

system("open figures/figure-2.png")

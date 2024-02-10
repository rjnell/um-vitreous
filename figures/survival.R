# Load libraries
library(stringr)
library(readxl)
library(survival)

# Load data
vitreous_data = readxl::read_xlsx("data/data.xlsx", skip=0)[1:65,]
survival = vitreous_data$`follow-up_days`/365.25
status = vitreous_data$melanoma_death
status[which(status == "yes")] = 1
status[which(status == "no")] = 0
status = as.numeric(status)

# Function create figures
survival_plot = function(filename, r=F, width=3400, height=2500, xlab = "Years", p_label="", labels = rep("",10), ylab = "Melanoma-related survival", tick = 0.01, colors, legend_positions) {
  
  # Close current images
  system("taskkill /F /IM Microsoft.Photos.exe /T")
  
  # Init PNG and plot
  png(filename, res=600, width=3600, height=2300)
  par(mar=c(5,6,5,10),xpd=T)
  xlim = c(0,16)
  ylim = c(0,1)
  plot(xlim, ylim, type="n", axes=F, ylab="", xlab="", xaxs = "i", yaxs = "i", xlim=xlim, ylim=ylim)
  xlim = c(0,14)
  
  # Plot y grid
  yat = seq(0, 1, by=0.25)
  segments(xlim[1],yat,xlim[2],lwd=1.4,col="#EEEEEE", xpd=T)
  
  # Plot x grid
  xat = seq(xlim[1],xlim[2],by=2)
  segments(xat,ylim[1],xat,ylim[2],lwd=1.4,col="#EEEEEE", xpd=T)
  
  # Plot y axis
  axis(side = 2, at = yat, labels=paste0(yat*100,"%"),las=2, col = "#b1b1b1", lwd = 1.4, col.axis="#333333", cex.axis=1.1, line = 0)
  mtext(side = 2, text = "Melanoma-related", line = 4.9, at = mean(yat), cex=1.1, col="#333333")
  mtext(side = 2, text = "survival", line = 3.9, at = mean(yat), cex=1.1, col="#333333")
  mtext(side = 2, text = "", line = 4.5, at = mean(yat), cex=1.1)
  
  # Plot x axis
  axis(side = 1, at = xat[c(1,3,5,7)], col = "#b1b1b1", lwd = 1.4, col.axis="#333333", cex.axis=1.1, line = 0)
  axis(side = 1, at = xat[c(2,4,6,8)], col = "#b1b1b1", lwd = 1.4, col.axis="#333333", cex.axis=1.1, line = 0)
  mtext(side = 1, text = xlab, line=3,at = mean(xat), cex=1.1, col="#333333")
  
  # Plot survival curves
  survival = survfit(Surv(survival, status) ~ survival_groups)
  lines(survival, col=colors, lwd=1.4, xpd=T, lty=c(1,1))
  lines(survival, col=colors, lwd=1.4, xpd=T, lty=c(1,1))
  
  # Plot censor ticks
  tick = 0.01
  jj = 0
  for (k in 1:length(survival$strata)) {
    for (i in 1:survival$strata[k]) {
      j=i+jj
      if (survival$n.censor[j] > 0) {
        segments(survival$time[j],survival$surv[j]+tick,survival$time[j],survival$surv[j]-tick,lwd=1.4,col=colors[k],xpd=T)
      }
    }
    jj = j
  }
  
  # Plot legend
  for (k in 1:length(colors)) {
    y = legend_positions[k]
    x = 16
    segments(x, y, x+1, y,lwd=1.4,col=colors[k],xpd=T)
    segments(x+0.5, y+tick, x+0.5, y-tick,lwd=1.4,col=colors[k],xpd=T)
    text(x+1, y, labels=labels[k], cex=1.1, pos=4, col='#333333',xpd=T)  
  }
  
  rect(1.5,0.08,2.5,0.17,col="white",border=NA)
  if (r) {
    rect(1.5,0.08,4.5,0.17,col="white",border=NA)
    text(0.25, .105, labels=p_label, cex=1.1, pos=4, col='#333333',xpd=T) 
  }else {
    text(0.7, .125, labels=p_label, cex=1.1, pos=4, col='#333333',xpd=T) 
  }
  
  # Finalize plot and PNG
  dev.off()
  system(paste0("open ", filename))
}

# Define groups (1) and plot results
survival_groups = rep(NA, 65)
survival_groups[which(vitreous_data$Classification == "vfDNA+/UM-")] = 1
survival_groups[which(vitreous_data$Classification == "vfDNA+/UM+")] = 2
survdiff(Surv(survival, status) ~ survival_groups)
plot(survfit(Surv(survival, status) ~ survival_groups))
dev.off()
survival_plot("figures/survival_primary.png", 
              colors =c("#F44336","#2196F3"), 
              legend_positions = c(0.56,0.44),
              p_label = "n/s",
              labels = c("vfDNA+/UM\U2212 (n=24)","vfDNA+/UM+ (n=39)"),
              r=F)

# Define groups (2) and plot results
survival_groups = rep(NA, 65)
survival_groups[which(vitreous_data$vfDNA_survival_group == "1")] = 1
survival_groups[which(vitreous_data$vfDNA_survival_group == "2")] = 2
diff = survdiff(Surv(survival, status) ~ survival_groups)
pchisq(diff$chisq, length(diff$n)-1, lower.tail = FALSE)
survival_plot("figures/survival_secondary.png", 
              p_label=expression(italic(p)~"= 0.020"), 
              r = T,
              colors = RColorBrewer::brewer.pal(9,"Paired")[c(4,8)], 
              labels = c("",""),
              legend_positions = c(0.76,0.29))
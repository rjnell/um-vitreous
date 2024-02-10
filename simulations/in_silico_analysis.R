# Load libraries
library(digitalPCRsimulations)
library(readxl)

# Set seed
seed = 123

# Function to calculate stats
calculate_stats = function(cnv, n_simulations, n_droplets_range, alpha) {
  
  # Initialise results
  results = NULL
  
  # Classic CNV
  universe_cnv = universe(input_ng * cnv/2)
  universe_ref = universe(input_ng)
  
  # SNP-CNV
  cnv_snp1 = 1
  cnv_snp2 = cnv-cnv_snp1
  universe_snp1 = universe(input_ng * cnv_snp1/2)
  universe_snp2 = universe(input_ng * cnv_snp2/2)
  
  for (n_droplets in n_droplets_range) {   
    
    # Classic CNV
    classic_cnv = simulate_ratios(universe_cnv,
                                  universe_ref,
                                  n_droplets = n_droplets,
                                  n_simulations = n_simulations,
                                  alpha = alpha)
    classic_cnv = classic_cnv*2

    # SNP-CNV
    snp_cnv = simulate_ratios(universe_snp1,
                              universe_snp2,
                              n_droplets = n_droplets,
                              n_simulations = n_simulations,
                              alpha = alpha)
    for (sim in 1:nrow(snp_cnv)) {
      snp_cnv[sim,] = 1/snp_cnv[sim,c(1,3,2)]+1
    }
   
    # Bind results and print progress
    results = rbind(results, c(n_droplets, (1-stats(classic_cnv,2)$coverage)*100, (1-stats(snp_cnv,2)$coverage)*100))
    print(n_droplets)
  }
  
  return(results)
}

# Function to plot results
plot_results = function(results, main) {
  plot(c(min(n_droplets_range),max(n_droplets_range)),
       c(0,100),
       type="n",
       main=main,
       ylab="% sensitivity",
       xlab="droplets",
       xaxs = "i", yaxs = "i",
       axes=F)
  yat = seq(0,100,by=20)
  segments(n_droplets_range,0,n_droplets_range,100,col="#EEEEEE", xpd=T)
  segments(min(n_droplets_range),yat,max(n_droplets_range),col="#EEEEEE", xpd=T)
  axis(side=1,col = "#b1b1b1",at = n_droplets_range)
  axis(side=2,las=2,col = "#b1b1b1")
  points(results[,1], results[,2], pch=16, col="#2196F3", xpd=T)
  lines(results[,1],predict(loess(results[,2] ~ results[,1])), lty=2, col="#2196F3")
  points(results[,1], results[,3], pch=16, col="#F44336", xpd=T)
  lines(results[,1],predict(loess(results[,3] ~ results[,1])), lty=2, col="#F44336")
}

# Set droplets range
n_droplets_range = seq(15000, 45000, by=2500)

# Tumour ID
tumour = "vfDNA-24"

# Input values (from Gaq measurement)
measured_copiesperuL = 6.90
measured_ccf = 0.391
input_ng = measured_copiesperuL * 0.00085 / (10^-9) * (3.59 * 10^-12) * 20000
  
# How many simulations per condition should be performed?
n_simulations = 1000
  
# Which alpha should be used?
alpha = 0.01
  
# Perform simulations
# Clonal loss (CNV=1)
cnv = 1*measured_ccf+2*(1-measured_ccf)  
clonal_loss = calculate_stats(cnv, n_simulations, n_droplets_range, alpha)

# Clonal gain (CNV=3)
cnv = 3*measured_ccf+2*(1-measured_ccf)  
clonal_gain = calculate_stats(cnv, n_simulations, n_droplets_range, alpha)
  
# Combine and save results
combined_results = cbind(clonal_loss, clonal_gain)
colnames(combined_results) = c("clonal_loss_droplets","clonal_loss_classic_cnv","clonal_loss_snp_cnv",
                               "clonal_gain_droplets","clonal_gain_classic_cnv","clonal_gain_snp_cnv")
write.table(combined_results, paste0("simulations/",tumour,".tsv"), sep="\t", row.names = F)
  
# Plot PDF
graphics.off()
pdf(paste0("simulations/",tumour,".pdf"),width = 8, height=2.5)
par(mfrow=c(1,3))
  
# Plot info box
par(mar=c(0,0,0,0))
plot(c(0,120),c(0,100), type="n",axes=F, xlab="", ylab="")
text(0,92.5,labels=tumour, font=2, pos=4,cex=1.75)
text(0,85,labels="Concentration", pos=4)
text(50,85,labels=paste(measured_copiesperuL,"copies/uL"), pos=4)
text(0,80,labels="CCF", pos=4)
text(50,80,labels=paste0(format(measured_ccf*100,nsmall=1,digits=1),"%"), pos=4)
text(0,70,labels="n simulations", pos=4)
text(50,70,labels=n_simulations, pos=4)
text(0,65,labels="random seed", pos=4)
text(50,65,labels=seed, pos=4)
points(5+3, 10, pch=16, col="#2196F3", xpd=T)
lines(c(0,10)+3,c(10,10), lty=2, col="#2196F3")
text(10+3,10,labels="Classic CNV", pos=4)
points(65, 10, pch=16, col="#F44336", xpd=T)
lines(c(60,70),c(10,10), lty=2, col="#F44336")
text(70,10,labels="SNP-CNV", pos=4)
  
# Results
par(mar=c(5.1,4.1,4.1,2.1))
plot_results(clonal_loss, "Clonal loss (CNV=1)")
plot_results(clonal_gain, "Clonal gain (CNV=3)")

# Finalise PDF
dev.off()
system(paste0("open simulations/",tumour,".pdf"))
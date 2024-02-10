# Load library
library(digitalPCRsimulations)

# Shortcut function to generate copy number value (CNV) simulations
simulate_CNV = function(
    input_ng = NA, 
    cnv_snp1, 
    cnv_snp2 = 1, 
    n_droplets = 17500, 
    n_simulations = 100,
    alpha = 0.01, 
    input_copies = NA) {
  
  # Determine input (in copies or ng)
  if (is.na(input_copies)) {
    input = input_ng
  }
  else {
    input = input_copies * 0.00085 / (10^-9) * (3.59 * 10^-12) * 20000
  }
  
  # Calculate total number of copies target per haploid genome equivalent
  cnv = cnv_snp1/2 + cnv_snp2/2
  
  # Create universes
  universe_snp1 = universe(input * cnv_snp1/2)
  universe_snp2 = universe(input * cnv_snp2/2)
  universe_cnv = universe(input * cnv)
  universe_ref = universe(input)
  
  # Simulate classic CNV experiment
  classic_cnv = simulate_ratios(universe_cnv,
                                universe_ref,
                                n_droplets = n_droplets,
                                n_simulations = n_simulations,
                                alpha = alpha)
  classic_cnv = classic_cnv * 2
  
  # Simulate new SNP experiment
  snp_cnv = simulate_ratios(universe_snp1,
                            universe_snp2,
                            n_droplets = n_droplets,
                            n_simulations = n_simulations,
                            alpha = alpha)
  snp_cnv = snp_cnv+1
  
  # Return results
  res = list()
  res[["cnv_snp1"]] = cnv_snp1
  res[["classic"]] = classic_cnv
  res[["snp"]] = snp_cnv
  
  return(res)
}


# Function to plot CNV simulations and observations
plot_CNV = function(
    file, 
    cnv_snp1 = 0, 
    true = NULL,
    truetext = "",
    n_droplets, 
    subheader = "", 
    approach, 
    ylim, 
    ysteps = 0.1, 
    observed = NULL, 
    colors = TRUE,
    coverage = FALSE, 
    input_ng = NA, 
    input_copies = NA, 
    n_simulations = 100, 
    cnv_data = NULL) {
  
  # Initialise PNG image
  png(file, res = 600, width = 5000, height = 2200)
  
  # Initialise plot
  par(mar = c(4,7,6,7))
  xlim = c(0,6.5)
  plot(xlim, ylim, type = "n", axes = F, ylab = "", xlab = "", xaxs = "i", yaxs = "i")
  ys = 0.3
  
  #
  x1 = 00
  x2 = 3.575+(length(observed)-1)*0.125+0.4
  y2 = ylim[2]+ (ylim[2]-ylim[1])/5*1.25
  y1 = ylim[2]+ (ylim[2]-ylim[1])/5*2.
  
  segments(x1, y2, x1, mean(c(y1,y2)), col='#b1b1b1', lwd=1.4, xpd=T)
  segments(x2, y2, x2, mean(c(y1,y2)), col='#b1b1b1', lwd=1.4, xpd=T)
  segments(x1, mean(c(y1,y2)), x2, col='#b1b1b1', lwd=1.4, xpd=T)
  segments(mean(c(x1,x2)), mean(c(y1,y2)), mean(c(x1,x2)), y1, col='#b1b1b1', lwd=1.4, xpd=T)
  
  if (approach  ==  "classic") {
    lab = "Classic approach" }
  else {
    lab = "SNP-based approach"
  }
  text(labels = lab, x = mean(c(x1,x2)), y = ylim[2]+(ylim[2]-ylim[1])/5*3.1, cex = 1.1, col = "#333333",xpd=T)
  text(labels = subheader, x = mean(c(x1,x2)), y = ylim[2]+(ylim[2]-ylim[1])/5*2.5, cex = 1.1, col = "#333333",xpd=T)
  
  # Plot y axis
  yat = seq(ylim[1], ylim[2], by = ysteps)
  segments(0,yat,3.175,lwd = 1.4,col = "#EEEEEE", xpd = T)
  axis(side = 2, at = yat, labels = format(yat, nsmall = 2), las = 2, 
       col = "#b1b1b1", lwd = 1.4, col.axis = "#333333", cex.axis = 1.1, line = 0)
  
  
  # Get in silico data
  if (is.null(cnv_data)) {
    cnv = simulate_CNV(input_copies = input_copies, input_ng = input_ng, 
                       cnv_snp1 = cnv_snp1, n_droplets = n_droplets, 
                       n_simulations = n_simulations)
  }
  else {
    cnv = cnv_data
  }
  cnv_approach = cnv[[approach]][1:20,]
  cnv_approach = cnv_approach[order(cnv_approach[,1]),]
  
  if (!is.null(true)) {
    segments(0, true, 3.175, lwd = 1.4, lty = 3, col = "#b1b1b1")  
  }
  
  
  # Plot simulations
  for (i in 1:20) {
    x = 0.4+(i-1)*0.125
    col = "#808080"
    if (colors) {
      if (cnv_approach[i,2] < 2 & cnv_approach[i,3] > 2) {
        col = "#F88E86"
      }
    }
    arrows(x, cnv_approach[i,2], x, cnv_approach[i,3], length = 0.05, 
           angle = 90, code = 3, col = col, lwd = 1.4, xpd = T)
    if (approach  ==  "classic") {
      points(x, cnv_approach[i,1], pch = 19, cex = 0.8, col = "#333333")
    } else {
      points(x, cnv_approach[i,1], pch = 15, cex = 0.8, col = "#333333")
    }
  }
  
  # Plot subheaders
  mtext(text = paste0("Simulated"), at = 1.5875, line = 0.5, side = 3, xpd = T, 
        cex = 1.1, col = "#333333", font = 3)
  s = stats(cnv[[approach]], true_value = cnv_snp1+1)
  m = format(round(s$point_estimate_mean, 2), nsmall = 2)
  sd = format(round(s$point_estimate_sd, 2), nsmall = 2)
  cc = format(round(s$coverage*100, 1), nsmall = 1)
  ss = stats(cnv[[approach]], true_value = 2)
  sens = format(round((1-ss$coverage)*100, 1), nsmall = 1)
  mtext(text = paste0("Mean = ", m), at = 1.5875, line = 0.5, 
        side = 1, xpd = T, cex = 1.1, col = "#333333", font = 1)
  mtext(text = paste0("Standard deviation = ", sd), at = 1.5875, line = 1.5, 
        side = 1, xpd = T, cex = 1.1, col = "#333333", font = 1)
  if (coverage) {
    mtext(text = paste0("99%-CI coverage = ", cc, "%"), at = 1.5875, line = 2.5, 
          side = 1, xpd = T, cex = 1.1, col = "#333333", font = 1)
  } else {
    mtext(text = paste0("Sensitivity = ", sens, "%"), at = 1.5875, line = 2.5, 
          side = 1, xpd = T, cex = 1.1, col = "#333333", font = 1)
  }
  
  # Plot observations
  if (!is.null(observed)) {
    segments(3.175, yat, 3.575+(length(observed)-1)*0.125+0.4, lwd = 1.4,
             col = "#EEEEEE", xpd = T)
    if (!is.null(true)) {
      segments(3.175, true, 3.575+(length(observed)-1)*0.125+0.4, lwd = 1.4, lty = 3, col = "#b1b1b1")  
    }
    for (i in 1:length(observed)) {
      x = 3.575+(i-1)*0.125
      col = "#808080"
      if (colors) {
        if (observed[[i]][2] < 2 & observed[[i]][3] > 2) {
          col = "#F88E86"
        }
      }
      arrows(x, observed[[i]][2], x, observed[[i]][3], length = 0.05, 
             angle = 90, code = 3, col = col, lwd = 1.4, xpd = T)
      if (approach  ==  "classic") {
        points(x, observed[[i]][1], pch = 19, cex = 0.8, col = "#333333")
      } else {
        points(x, observed[[i]][1], pch = 15, cex = 0.8, col = "#333333")
      }
      
    }
    
    # Plot subheaders
    mtext(text = paste0("Observed"), at = mean(c(3.175,x+0.4)), line = 0.5, 
          side = 3, xpd = T, cex = 1.1, col = "#333333", font = 3)
    segments(3.175,ylim[1],3.175,ylim[2],lwd = 1.4,col = "#b1b1b1", xpd = T)
  }
  
  # Plot axis label
  mtext(side = 2, text = "Copy number value", line = 4, at = mean(yat), cex = 1.1, col="#333333", col.axis="#333333")
  
  
  dev.off()
  system(paste("open", file))
}


###
# Figure 5B: vfDNA-41 (loss)
###

# Perform simulations
set.seed(9045)
ccf = 36.1
cnv_snp1 = (2 * (1-ccf/100) + 1 * ccf/100) - 1 
cnv_data = simulate_CNV(input_copies = 153.3, cnv_snp1 = cnv_snp1, n_droplets = 15000, n_simulations = 1000)

# Save plot 1
plot_CNV(
  cnv_data = cnv_data,
  file = "figures/vfDNA_41_gain_classic_single.png",
  subheader = "(single experiment)",
  approach = "classic",
  ylim = c(1.3,2.0),
  ystep = 0.1,
  observed = list(
    c(1.64,1.50,1.80)
  )
)

# Save plot 2
plot_CNV(
  cnv_data = cnv_data,
  file = "figures/vfDNA_41_gain_snp_single.png",
  subheader = "(single experiment)",
  approach = "snp",
  ylim = c(1.3,2.0),
  ystep = 0.1,
  observed = list(
    c(1.63,1.55,1.73)
  )
)

###
# Figure 5C: vfDNA-24 (loss)
###

# Perform simulations
set.seed(8030)
ccf = 39.1
cnv_snp1 = (2 * (1-ccf/100) + 1 * ccf/100 ) - 1 
cnv_data_single = simulate_CNV(input_copies = 6.90, cnv_snp1 = cnv_snp1, n_droplets = 15000, n_simulations = 1000)
cnv_data_merged = simulate_CNV(input_copies = 6.90, cnv_snp1 = cnv_snp1, n_droplets = 45000, n_simulations = 1000)

# Save plot 0 (not in publication)
plot_CNV(
  cnv_data = cnv_data_single,
  file = "figures/vfDNA_24_loss_classic_single.png",
  subheader = "(single experiment)",
  approach = "classic",
  ylim = c(0.5,3.5),
  ystep = 0.5
)

# Save plot 1
plot_CNV(
  cnv_data = cnv_data_single,
  file = "figures/vfDNA_24_loss_snp_single.png",
  subheader = "(single experiment)",
  approach = "snp",
  ylim = c(0.5,3.5),
  ystep = 0.5,
  observed = list(
    c(1.51,1.23,1.96),
    c(1.61,1.28,2.18),
    c(1.64,1.29,2.25)
  )
)

# Save plot 2
plot_CNV(
  cnv_data = cnv_data_merged,
  file = "figures/vfDNA_24_loss_snp_merged.png",
  subheader = "(merged [3x] experiment)",
  approach = "snp",
  ylim = c(0.5,3.5),
  ystep = 0.5,
  observed = list(
    c(1.58,1.39,1.84)
  )
)

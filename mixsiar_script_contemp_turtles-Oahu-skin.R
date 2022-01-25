# Adriana Velasquez-Vacca
# February 1st, 2021
# Script for Contemporary Green turtle analysis including 3 isotopes in Big Island
# because we have different sources per island I had to make a different script for each island.

# Files needed:
# contemp_turtles_Oahu_consumer_skin.csv
# contemp_turtles_Oahu_sources_skin.csv
# contemp_turtles_Oahu_discrimination_skin.csv

rm(list=ls()) # clear variables
graphics.off() # close figures

start1 <- Sys.time() # To record starting time and then know how long it took to run the script.

setwd("~/Dropbox/PHD/Isotopes/results/Bones/MixSIAR")


set.seed(1) # que carajos es esto? no entiendo como funciona generar un solo n�mero random para que? The random numbers are the same, and they would continue to be the same no matter how far out in the sequence we went. Tip. Use the set.seed function when running simulations to ensure all results, figures, etc are reproducible.


library(MixSIAR)

# Load consumer
#mix.filename <- system.file("data", "contemp_turtles_Oahu_consumer_skin.csv", package = "MixSIAR")
mix <- load_mix_data(filename="data/skin/contemp_turtles_Oahu_consumer_skin.csv",
           iso_names=c("d13C","d15N","d34S"),
           factors=NULL,
           fac_random=NULL,
           fac_nested=NULL,
           cont_effects=NULL)

# Load source
source <- load_source_data(filename="data/skin/contemp_turtles_Oahu_sources_skin.csv",
               source_factors=NULL,
               conc_dep=TRUE,
               data_type="means",
               mix)

# Load trophic discrimination factor (TDF)
discr <- load_discr_data(filename="data/skin/contemp_turtles_Oahu_discrimination_skin.csv", mix)

# Plot isospace

plot_data(filename="results/CNS/skin/contemp_turtles_Oahu_isospace_plot_skin",
      plot_save_pdf=FALSE, # If I put yes to pdf it doesn't work and doesnt even save the png
      plot_save_png=TRUE,
      mix,source,discr)
          
# Calculate standardized convex hull area
if(mix$n.iso==2) calc_area(source=source,mix=mix,discr=discr)


################################################################################

# Output options for the plots:
output_options <- list(summary_save = TRUE,
                       summary_name = "results/CNS/skin/contemp_turtles_Oahu_summary_statistics_skin_uninf", # 
                       sup_post = FALSE,
                       plot_post_save_pdf = TRUE,
                       plot_post_name = "results/CNS/skin/contemp_turtles_Oahu_posterior_density_skin_uninf", # 
                       sup_pairs = FALSE,
                       plot_pairs_save_pdf = TRUE,
                       plot_pairs_name = "results/CNS/skin/contemp_turtles_Oahu_pairs_plot_skin_uninf", # 
                       sup_xy = TRUE,
                       plot_xy_save_pdf = TRUE,
                       plot_xy_name = "results/CNS/skin/contemp_turtles_Oahu_xy_plot_skin_uninf", # 
                       gelman = TRUE,
                       heidel = TRUE,
                       geweke = TRUE,
                       diag_save = TRUE,
                       diag_name = "results/CNS/skin/contemp_turtles_Oahu_diag_skin_uninf",
                       indiv_effect = FALSE,
                       plot_post_save_png = TRUE,
                       plot_pairs_save_png = TRUE,
                       plot_xy_save_png = TRUE,
                       diag_save_ggmcmc = TRUE)

# default "UNINFORMATIVE" / GENERALIST prior (alpha = 1)

par(mar=c(1,1,1,1)) # margins size
plot_prior(alpha.prior=1, source=source, plot_save_pdf=TRUE, plot_save_png=TRUE,filename = "results/CNS/skin/contemp_turtles_Oahu_prior_plot_kw_uninf_skin") # saves it as pdf

# Define model structure and write JAGS model file
model_filename <- "results/CNS/skin/contemp_turtles_Oahu_MixSIAR_model_kw_uninf_skin.txt"   # Name of the JAGS model file
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

# Run the JAGS model (test, normal, long, very long) Normal took
jags.uninf <- run_model(run="normal",mix,source,discr,model_filename,alpha.prior = 1, resid_err, process_err)
# Process diagnostics, summary stats, and posterior plots
output_JAGS(jags.uninf, mix, source, output_options)

# ################################################################################
# # INFORMATIVE prior (construct alpha with prior data)
# 
# # Output options for the plots:
# output_options <- list(summary_save = TRUE,
#                        summary_name = "results/CNS/skin/contemp_turtles_Oahu_summary_statistics_skin",
#                        sup_post = FALSE,
#                        plot_post_save_pdf = TRUE,
#                        plot_post_name = "results/CNS/skin/contemp_turtles_Oahu_posterior_density_skin",
#                        sup_pairs = FALSE,
#                        plot_pairs_save_pdf = TRUE,
#                        plot_pairs_name = "results/CNS/skin/contemp_turtles_Oahu_pairs_plot_skin",
#                        sup_xy = TRUE,
#                        plot_xy_save_pdf = TRUE,
#                        plot_xy_name = "results/CNS/skin/contemp_turtles_Oahu_xy_plot_skin",
#                        gelman = TRUE,
#                        heidel = TRUE,
#                        geweke = TRUE,
#                        diag_save = TRUE,
#                        diag_name = "results/CNS/skin/contemp_turtles_Oahu_diag_skin",
#                        indiv_effect = FALSE,
#                        plot_post_save_png = TRUE,
#                        plot_pairs_save_png = TRUE,
#                        plot_xy_save_png = TRUE,
#                        diag_save_ggmcmc = TRUE)
# 
# # Prior sample values # Halophila is 20%. Make sure is in alphabetical order, not the file order.
# kw.alpha <- c(1,1,1,1)
# 
# # Generate alpha hyperparameters scaling sum(alpha)=n.sources
#  kw.alpha <- kw.alpha*length(kw.alpha)/sum(kw.alpha)
# 
# # the Dirichlet hyperparameters for the alpha.prior cannot be 0 (but can set = .01)
#  kw.alpha[which(kw.alpha==0)] <- 0.01
# 
#  # Plot your informative prior
#  par(mar=c(1,1,1,1)) # margins size
#  plot_prior(alpha.prior=kw.alpha, source=source, plot_save_pdf=TRUE, plot_save_png=TRUE,filename = "results/CNS/skin/contemp_turtles_Oahu_prior_plot_kw_skin") # saves it as pdf
# 
# # Define model structure and write JAGS model file
#  model_filename <- "results/CNS/skin/contemp_turtles_Oahu_MixSIAR_model_kw_skin.txt"   # Name of the JAGS model file
#  resid_err <- TRUE
#  process_err <- TRUE
#  write_JAGS_model(model_filename, resid_err, process_err, mix, source)
# 
# # Run the JAGS model (test, normal, long, very long)
#  jags.inf <- run_model(run="test",mix,source,discr,model_filename,alpha.prior=kw.alpha)
# 
#  # Process diagnostics, summary stats, and posterior plots
#  output_JAGS(jags.inf, mix, source, output_options)

graphics.off() # To close all graph windows


end1 <- Sys.time() # Record ending time to run the script. 
time1 <- end1 - start1 # Calculate dif star and end and save it as time1.
time1 # Show how long it took

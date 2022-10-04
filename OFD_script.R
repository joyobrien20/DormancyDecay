# Joy O'Brien
# Lennon Lab 
# Attempting to do OFD analysis for the total DNA community and active community
# from Locey et al., 2020

# Notes from Jay & Nathan's thought document
# With these data, we would work with the site-by-species matrix. 

# Load libraries needed
library(vegan)
library(dplyr)
library(readxl)



# Convert to presence-absence matrix. 
# Presence absence for cDNA and DNA combined is the PondsPA dataset (INPond_Initia.RData), prettry sure its pondsPA

# We probably need to separate cDNA from DNA here
# use subset function?

# Read in the active community PA table 
PondsPA_active <- subset(PondsPA, rownames(PondsPA %in% "-cDNA"))

#PondsPA_total <- subset(PondsPA, rownames(PondsPA %in% "DNA")) #trying to get this line of code to work; subset must be logical error
PondsPA_total <- read_excel#(insert file path here)

# For each row (species), we would take the sum across sites (columns), which would result in a vector that 
#could be used to make frequency histograms: one for DNA and one for RNA. 

# Transpose the data
PondsPA_active_t <- t(PondsPA_active)
PondsPA_total_t <- t(PondsPA_total)

# Export the transposed data and re-import to get the right format (can't figure it out in R!)
write.csv(PondsPA_active_t,"~/Desktop/PondsPA_active_t.csv")
write.csv(PondsPA_total_t, "~/Desktop/PondsPA_total_t.csv")
#Delete the added row at the top of each data frame 
# PondsPA_active_t <- PondsPA_active_t[-1,] this line of code doesn't work, so I did it manually (files "PondsPA_active_OFD.csv" and "PondsPA_total_OFD.csv")

# Step 1: Sum across sites (columns)
# For active community
active <- colSums(PondsPA_active_OFD[,-1]) 
summary(active)
qqnorm(active)
histogram(active)

# For total community 
total <- colSums(PondsPA_total_OFD[,-1])
summary(total)
qqnorm(total)
histogram(total)


# We could then statistically compare by 
# fitting models or some other approach. I could imagine that we might also be able to do some resampling of the data. 
#We might want to control for differences in species richness. Perhaps this could be done by relativizing y-axis so it’s % species in the DNA or RNA pool. But we could also consider bootstrapping, either by species or sites. 
# For example pull out 1000 species randomly with replacement?


# Response from Nathan:
# Just compare colSums(DNA)/nrow(DNA) to colSums(RNA)/nrow(RNA) on your site by species matrix
#nrow = the number of rows 


# Performed with the "MOStest" function in the vegan package with a presence/absence OTU table 
test <- MOStest(active, total)
test <- MOStest(log(active), , family=quasipoisson())

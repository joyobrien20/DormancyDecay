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
library(writexl)
library(tidyverse)
library(ggplot2)

# Convert to presence-absence matrix. 
# Presence absence for cDNA and DNA combined is the PondsPA dataset (INPond_Initia.RData), prettry sure its pondsPA

# Subset the data based on active community and total community
# cDNA
active_PA <- PondsPA[grep('-cDNA', rownames(PondsPA)),]

# DNA
total_PA <- PondsPA[grep('-DNA', rownames(PondsPA)),]

# save these data matrices 
#saveRDS(active_PA,"~/Desktop/IU/Lennon Lab/DormancyDecay/active_PA.rds")
#saveRDS(total_PA, "~/Desktop/IU/Lennon Lab/DormancyDecay/total_PA.rds")

# For each row (species), we would take the sum across sites (columns), which would result in a vector that 
#could be used to make frequency histograms: one for DNA and one for RNA. 

# Transpose the data
PondsPA_active_t <- t(active_PA)
PondsPA_total_t <- t(total_PA)

# save these transposed matrices here
saveRDS(PondsPA_total_t, "~/Desktop/IU/Lennon Lab/DormancyDecay/PondsPA_total_t.rds")
saveRDS(PondsPA_active_t, "~/Desktop/IU/Lennon Lab/DormancyDecay/PondsPA_active_t.rds")

#*******************************************************************************************

# Step 1: Sum across sites (columns)
# For active community
active <- colSums(PondsPA_active_t[,-1]) 
summary(active)
print(active)

qqnorm(active)
hist(active)

# for total community 
total <- colSums(PondsPA_total_t[,-1])
print(total)

qqnorm(total)
hist(total)

# Step 2: Divide col sums RNA/n row RNA
nrow(PondsPA_active_t)
active_OFD <- active/34059 #34059 is the number of OTUs which is the number of rows in the site by species matrix 
print(active_OFD)

hist(active_OFD)

# For total community 
nrow(PondsPA_total_t)
total_OFD <- total/34059
summary(total_OFD)
print(total)

qqnorm(total)
hist(total_OFD)
plot(total)
barplot(active_OFD)
barplot(total_OFD)


# Now divide col sums DNA/nrow DNA 
total_OFD <- total/16383 #16383 is the number of OTUs
print(DNA_comp)
summary(DNA_comp)
summary(RNA_comp)
hist(DNA_comp)

# Write the DNA and RNA comps to merge in excel 
write.csv(total_OFD, "~/Desktop/total_OFD.csv")
write.csv(active_OFD, "~/Desktop/active_OFD.csv")


# Read in data for OFD


# https://rdrr.io/rforge/vegan/man/MOStest.html
test <- MOStest(DNA_comp, RNA_comp) #this isnt working because there are more samples for DNA than RNA #this is a whole problem where I need to start from sq.1
# Visualize the relationship 
# Merge the values?
# We could then statistically compare by 
# fitting models or some other approach. I could imagine that we might also be able to do some resampling of the data. 
#We might want to control for differences in species richness. Perhaps this could be done by relativizing y-axis so it’s % species in the DNA or RNA pool. But we could also consider bootstrapping, either by species or sites. 
# For example pull out 1000 species randomly with replacement?


# Response from Nathan:
# Just compare colSums(DNA)/nrow(DNA) to colSums(RNA)/nrow(RNA) on your site by species matrix
#nrow = the number of rows 

# Figure out how to merge this with Git- https://stackoverflow.com/questions/57637795/fatal-could-not-read-username-when-pushing-to-git-using-a-cron-job

# Performed with the "MOStest" function in the vegan package with a presence/absence OTU table 
test <- MOStest(total_OFD_amend_valonly, active_OFD_valonly)
plot(total_OFD_amend_valonly)
plot(active_OFD_valonly)

# trial


test <- MOStest(log(active), , family=quasipoisson())

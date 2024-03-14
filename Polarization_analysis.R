library(stringr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(tidyr)
library(ggh4x)
library(ggsignif)


# Read in data from ImageJ measurements, obtained with Measure_slices_with_predefined_ROIsets.ijm
# use both replicates (II: 240203, and II: 240212)
repII <- list.files(path = "./mesodiff_II_patterning_analysis/Measurements",full.names = TRUE)
repIII <- list.files(path = "./mesodiff_III_patterning_analysis/Measurements",full.names = TRUE)
files <- c(repII, repIII)

dat <- data.frame()
for (i in 1:length(files)){
  temp_dat <- read.csv(files[i], row.names = "X.1")
  temp_dat$Celltype <- str_extract(files[i], "(?<=crop_).*?(?=_cells)")
  
  dat <- rbind(dat, temp_dat)
}

# Distribute label information over individual meta data columns
dat$Channel <- factor(str_sub(dat$Label, start = 1, end = 2))
levels(dat$Channel) <- c("Nuclei", "Cer1", "T")
dat$Replicate <- factor(str_sub(dat$Label, start = 4, end = 9))
levels(dat$Replicate) <- c("Rep2", "Rep3")
dat$Position <- str_extract(dat$Label, "(?<=Position\\s)\\d+(?=_z)")
dat$Zslice <- str_extract(dat$Label, "(?<=_z)\\d{2}")
dat$Cell <- str_sub(dat$Label, start = -9, end = -1)

# Filter cells for size and circularity
dat <- subset(dat, Area < 200 & Area > 15) # throws out 1100/3 cells
dat <- subset(dat, Circ. > 0.8) # throws out another 1500/3 cells

# Simplify dataframe
dat <- dat[,c(6,7,13,20:25)]


# Add a normalized median (Median / max(Median) group for every BELA in every Channel in every slice in every replicate in every celltype)
# Calculate vector from mean position of all cells (center) to individual cells
normalized_dat <- dat %>%
  group_by(Channel, Replicate, Position, Zslice,Celltype) %>%
  mutate(normalized_median = Median / max(Median),
         Mean_x = mean(X), Mean_y = mean(Y)) %>% # for mean position grouping by channel does not matter
  ungroup()

# Plot calculated normalized intentisity for every BELA and zslice
normalized_dat$identity <- paste(normalized_dat$Channel, normalized_dat$Position, normalized_dat$Replicate, normalized_dat$Celltype)

# Calculate vector between each cell and center of corresponding BELA
normalized_dat$Dif_x <- normalized_dat$X-normalized_dat$Mean_x
normalized_dat$Dif_y <- normalized_dat$Y-normalized_dat$Mean_y


# weight individual cells by randomized brightness
# randomize medians within group
number_of_permutations <- 100
sampled_median <- data.frame(row.names = 1:length(normalized_dat$X))
for (i in 1:number_of_permutations){
  normalized_dat_temp <- normalized_dat %>%
    group_by(Channel, Replicate, Position, Zslice, Celltype) %>%
    mutate(random = sample(normalized_median))
  sampled_median = cbind(sampled_median, normalized_dat_temp$random)
  # add meaningful name for added column
  column_names <- colnames(sampled_median)
  column_names[length(column_names)] <- paste0("random_", i)
  colnames(sampled_median) <- column_names
}

# Merge randomized and actual data
normalized_dat <- cbind(normalized_dat, sampled_median)
# reformat combined data
normalized_dat_long <- melt(normalized_dat, id.vars = c("X", "Y","Median", "Celltype", "Channel", "Replicate", "Position",
                                                        "Zslice","Cell","Mean_x","Mean_y", "identity","Dif_x","Dif_y"), 
                            variable.name = "datatype", value.name = "normalized_median")

# weighting individual cell directions by brightness
normalized_dat_long$cell_dipol_vector_x <- normalized_dat_long$Dif_x * normalized_dat_long$normalized_median
normalized_dat_long$cell_dipol_vector_y <- normalized_dat_long$Dif_y * normalized_dat_long$normalized_median

# Sum up.cell individual vectors per BELA
# Calculate Rgyr per BELA and per celltype and per datatype
Rgyr_dipole_dat <- normalized_dat_long %>%
  group_by(Channel, Replicate, Position, Zslice, Celltype, datatype) %>%
  summarise(cell_dipol_vector_sum_x = sum(cell_dipol_vector_x), 
            cell_dipol_vector_sum_y = sum(cell_dipol_vector_y),
            sum_of_vectors = sum(Dif_x^2 + Dif_y^2), # for calculating rgyr (sum of (ri-r)^2)
            n_cells = n()) # for calculating rgyr (number of cells per BELA per slice per celll type)

Rgyr_dipole_dat$Rgyr <- sqrt(Rgyr_dipole_dat$sum_of_vectors/Rgyr_dipole_dat$n_cells)

# normalize dipole vector with size estimate for each BELA and cell type
Rgyr_dipole_dat$relative_cell_dipol_vector_sum_x <- Rgyr_dipole_dat$cell_dipol_vector_sum_x/(Rgyr_dipole_dat$Rgyr^2)
Rgyr_dipole_dat$relative_cell_dipol_vector_sum_y <- Rgyr_dipole_dat$cell_dipol_vector_sum_y/(Rgyr_dipole_dat$Rgyr^2)

# get length of vectors, both actual data and permutated
Rgyr_dipole_dat <- Rgyr_dipole_dat %>%
  group_by(Channel, Replicate, Position, Zslice, Celltype) %>%
  mutate(relative_dipole_magnitude = sqrt(relative_cell_dipol_vector_sum_x^2 + relative_cell_dipol_vector_sum_y^2)) %>%
  ungroup()

# Clean up dataframe and calculate mean of all permutated data for each BELA and celltype
Rgyr_dipole_dat <- Rgyr_dipole_dat[,c(1:6,12:14)]
Rgyr_dipole_dat$datatype_general <- ifelse(Rgyr_dipole_dat$datatype == "normalized_median", yes = "data", no = "permutated data")

Rgyr_dipole_dat <- Rgyr_dipole_dat %>%
  group_by(Channel, Replicate, Position, Zslice, Celltype, datatype_general) %>%
  summarise(summed_x = mean(relative_cell_dipol_vector_sum_x),
            summed_y = mean(relative_cell_dipol_vector_sum_y),
            dipole_magnitude = mean(relative_dipole_magnitude)) %>%
  ungroup()

# For plotting subset data to only include meaningful combinations of Cer1 in VE and T in Epi (Nuclei in both)
Rgyr_dipole_dat_sub <- subset(Rgyr_dipole_dat, (Channel == "Cer1" & Celltype == "VE") | (Channel == "T" & Celltype == "Epi"))
Rgyr_dipole_dat_sub$Celltype_rdclass <- paste(Rgyr_dipole_dat_sub$Celltype, Rgyr_dipole_dat_sub$datatype_general)
Rgyr_dipole_dat_sub$BELA <- paste(Rgyr_dipole_dat_sub$Replicate, Rgyr_dipole_dat_sub$Position)

# Plot dipole magnitude
ggplot(Rgyr_dipole_dat_sub, aes(x=datatype_general, y=dipole_magnitude, fill = Celltype_rdclass)) + 
  geom_violin(linewidth = 0.2) +
  scale_fill_manual(values = c("#ffa6ff", "lightgrey", "#ffff00", "lightgrey"))+
  geom_jitter(size = 0.1, shape = 16, width = 0.3) + 
  geom_hline(aes(yintercept=0.02779585), linewidth = 0.4, color = "red")+
  xlab("datatype") + 
  ylab("polarization vector") + 
  scale_x_discrete(labels=c("data" = "exp", "permutated data" = "ctrl")) +
  facet_grid(cols =vars(Channel, Celltype), scales = "free_x") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        text = element_text(color = 'black', size = 8),
        axis.text = element_text(colour = "black", size = 8),
        strip.text = element_text(colour = "black", size = 8),
        axis.line = element_line(colour = "black", linewidth = 0.48),
        axis.ticks = element_line(color = "black", linewidth = 0.48),
        legend.position = "none",
        panel.background = element_rect(fill = NA)) + 
  force_panelsizes(rows = unit(2.5, "cm"),
                   cols = unit(1, "cm"))

# restructure for calculating angle between vectors
Epi <- subset(Rgyr_dipole_dat_sub, Celltype == "Epi")
VE <- subset(Rgyr_dipole_dat_sub, Celltype == "VE")

Epi <- pivot_wider(Epi, 
                   id_cols = c("Replicate","Position","Zslice","Celltype","BELA","Celltype_rdclass", "datatype_general"),
                   names_from = "Channel",
                   values_from = c("summed_x", "summed_y", "dipole_magnitude"))
Epi$Celltype <- NULL
Epi$Celltype_rdclass <- NULL

VE <- pivot_wider(VE, 
                  id_cols = c("Replicate","Position","Zslice","Celltype","BELA","Celltype_rdclass", "datatype_general"),
                  names_from = "Channel",
                  values_from = c("summed_x", "summed_y", "dipole_magnitude"))
VE$Celltype <- NULL
VE$Celltype_rdclass <- NULL

EpiandVE <- merge(Epi, VE, by = c("Replicate","Position","Zslice","BELA", "datatype_general"), suffixes = c("_Epi", "_VE"))


# Threshold based on 95% percent quantile in permutated data in Cer1 polarisation
EpiandVE$Cer_pos <- ifelse(EpiandVE$dipole_magnitude_Cer1>quantile(subset(EpiandVE, datatype_general == "permutated data")$dipole_magnitude_Cer1, probs = 0.95), 
                           yes = "polarized", no = "non-polarized")
EpiandVE <- subset(EpiandVE, datatype_general == "data")


# Plot angle between channel-specific dipoles 
# Between Cer and T
EpiandVE$dotproduct_Cer1_T <- EpiandVE$summed_x_Cer1*EpiandVE$summed_x_T+EpiandVE$summed_y_Cer1*EpiandVE$summed_y_T
EpiandVE$angle_Cer1_T <- acos(EpiandVE$dotproduct_Cer1_T/(EpiandVE$dipole_magnitude_Cer1*EpiandVE$dipole_magnitude_T))* (180/pi)
ggplot(EpiandVE, aes(x=Cer_pos, y=angle_Cer1_T)) +
  geom_violin(linewidth = 0.2) +
  geom_point(size = 0.1, shape = 16) + 
  geom_signif(comparisons = list(c("non-polarized", "polarized")), test = "ks.test", textsize = 2.5, tip_length = 0) +
  scale_y_continuous(breaks = c(0,90,180), limits = c(0,200))+
  ylab("polarization angle") + 
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        text = element_text(color = 'black', size = 8),
        axis.title.x = element_blank(),
        axis.text = element_text(colour = "black", size = 8),
        strip.text = element_text(colour = "black", size = 8),
        axis.line = element_line(colour = "black", linewidth = 0.48),
        axis.ticks = element_line(color = "black", linewidth = 0.48),
        legend.position = "none",
        panel.background = element_rect(fill = NA)) + 
  force_panelsizes(rows = unit(2.5, "cm"),
                   cols = unit(1, "cm"))




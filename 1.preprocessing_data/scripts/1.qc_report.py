#!/usr/bin/env python
# coding: utf-8

# # Whole Image Quality Control Report

# ## Import libraries

# In[1]:


suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(platetools))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(RColorBrewer))


# ## Set paths and variables

# In[2]:


# Paths to CSV files to generate QC report
path_to_outlier_results <- file.path("./qc_results/qc_outliers.csv")
path_to_cp_qc_results <- file.path("./qc_results/localhost231120090001/Image.csv")

# Output path for bar chart
barchart_output_file <- file.path(paste0("./qc_figures/localhost231120090001_fov_per_well_chart.png"))

# Output path for platemap
platemap_output_file <- file.path(paste0("./qc_figures/localhost231120090001_platemap_fov_per_well.png"))


# ## Load in outlier and quality control results

# In[3]:


# Read in CSV files
outlier_df <- read.csv(path_to_outlier_results)
qc_df <- read.csv(path_to_cp_qc_results)

# Add 'Metadata_Well' and 'Metadata_Site' columns using filename from 'FileName_OrigActin'
qc_df <- qc_df %>%
  mutate(Metadata_Well = gsub("_", "", str_extract(FileName_OrigActin, "_([^f]+)")),
         Metadata_Site = str_extract(FileName_OrigActin, "f.{2}"))
         
head(qc_df)


# ## Generate data frame removing poor quality sites and counting the number of sites left per well

# In[4]:


# Load in platemap and add cell type to include in plate map figure
platemap_df <- read.csv("../metadata/localhost231120090001_platemap.csv")

# Find the number of FOVs that are being removed due to being poor images
unique_combinations <- outlier_df[c('Metadata_Well', 'Metadata_Site')] %>% distinct()
cat("The number of image sets being removed is", nrow(unique_combinations), "out of 960 image sets\n")

# Remove the image sets from the main data frame to generate bar chart
result_df <- anti_join(qc_df, unique_combinations, by = c('Metadata_Well', 'Metadata_Site'))

# Group by 'Metadata_Well' and count the occurrences of 'Metadata_Site'
counts_per_well <- result_df %>% group_by(Metadata_Well) %>% summarise(count = n())

# Select only relevant columns
platemap_subset <- platemap_df[, c("well_position", "cell_type")]

# Add cell type to counts per well
counts_per_well <- merge(counts_per_well, platemap_subset, 
                         by.x = "Metadata_Well", by.y = "well_position", all.x = TRUE)

dim(counts_per_well)
head(counts_per_well)


# ## Generate bar chart to show how many FOVs were processed after QC

# In[5]:


# Extract the first letter from 'Metadata_Well' to create a color palette
unique_starting_letters <- unique(substr(counts_per_well$Metadata_Well, 1, 1))

# Create a color palette based on unique starting letters
color_palette <- colorRampPalette(brewer.pal(8, "Set1"))(length(unique_starting_letters))

# Map each unique starting letter to a color in a named vector
color_dict <- setNames(color_palette, unique_starting_letters)

# Create a new column 'Color' in the data frame based on the starting letter
counts_per_well$Color <- color_dict[substr(counts_per_well$Metadata_Well, 1, 1)]

# Increase the figure size to extend the chart horizontally
options(repr.plot.width=14, repr.plot.height=6)

# Create a bar chart using ggplot2 with the 'Color' column for colors
fov_chart <- ggplot(counts_per_well, aes(x = Metadata_Well, y = count, fill = Color)) +
  geom_bar(stat = 'identity') +
  labs(x = 'Metadata_Well', y = 'Number of FOVs', title = 'Count of FOVs per Well for plate localhost231120090001') +
  scale_fill_identity() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save plot to qc_figures
ggsave(
        barchart_output_file,
        fov_chart,
        dpi = 500,
        height = 6,
        width = 10
    )

# Display the plot in the notebook
print(fov_chart)


# ## Generate plate map to show distribution of FOVs left after QC

# In[6]:


fov_platemap <- platetools::raw_map(
    data = as.character(counts_per_well$count),
    well = counts_per_well$Metadata_Well,
    plate = 96,
    size = 8
    ) +

    ggtitle(paste("Platemap of FOV Count Per Well in plate localhost231120090001")) +
    theme(plot.title = element_text(size = 10, face = "bold")) +
    geom_point(aes(shape = counts_per_well$cell_type)) +
    scale_shape_discrete(name = "Cell Type") +
    scale_fill_manual(
    name = "FOV Count",
    values = c("12" = "#ffffff", "13" = "#c5d7f0", "14" = "#688cce", "15" = "#4a5ecc", "16" = "#2424d1") # Assign colors manually
  )
    ggsave(
    platemap_output_file,
    fov_platemap,
    dpi = 500,
    height = 3.5,
    width = 6
    )

# Display the plot in the notebook
print(fov_platemap)


# In[7]:


# Filter wells with less than 16 FOVs
filtered_counts <- counts_per_well %>% filter(count < 16)

# Group by cell_type and count the occurrences
count_by_cell_type <- filtered_counts %>% group_by(cell_type) %>% summarise(well_count = n())

# View the result
count_by_cell_type


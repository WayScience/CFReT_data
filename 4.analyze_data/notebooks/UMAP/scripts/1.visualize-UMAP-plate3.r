suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))

install.packages(c("gridExtra", "grid", "ggarrange", "patchwork"))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(grid))


# Set directory and file structure
umap_dir <- file.path("results")
plate3_name <- "localhost230405150001"
umap_file <- "UMAP_localhost230405150001_sc_feature_selected.tsv.gz"
umap_path <- file.path(umap_dir, umap_file)


output_fig_dir <- file.path("figures")
umap_prefix <- "UMAP_"

# Define output figure paths

output_umap_file <- file.path(
    output_fig_dir,
    paste0(umap_prefix, plate3_name)
)

print(output_umap_file)


# Load in the umap data
df <- readr::read_tsv(
    umap_path,
    col_types = readr::cols(
        .default = "d",
        "Metadata_WellRow" = "c",
        "Metadata_WellCol" = "c",
        "Metadata_heart_number" = "c",
        "Metadata_treatment" = "c",
        "Metadata_ImageNumber" = "c",
        "Metadata_cell_type" = "c",
        "Metadata_heart_failure_type" = "c",
        "Metadata_Plate" = "c",
        "Metadata_Well" = "c",
        "Metadata_Site" = "c",
        "Metadata_Cytoplasm_Parent_Cells" = "c",
        "Metadata_Cytoplasm_Parent_Nuclei" = "c",
        "Metadata_Cells_Number_Object_Number" = "c",
        "Metadata_Nuclei_Number_Object_Number" = "c"
    )
)
cell_count_df <- df %>%
    dplyr::group_by(Metadata_Well) %>%
    dplyr::count()

umap_cp_df <- df %>%
    dplyr::left_join(cell_count_df, by = "Metadata_Well")


# Generating UMAPs 
options(repr.plot.width = 12, repr.plot.height = 10)  # Adjust width and height as desired
output_file <- paste0(output_umap_file, "_treatment.png")

treatment_gg <- ggplot(umap_cp_df, aes(x = UMAP0, y = UMAP1)) +
  geom_point(aes(color = Metadata_treatment), size = 0.5, alpha = 0.4) +  # Adjusted alpha value
  theme_bw() +
  scale_color_manual(name = "Treatment",
                     values = c("DMSO" = "#BA5A31", "drug_x" = "#69DC9E", "TGFRi" = "#DC69B4")) +
  facet_wrap(~ Metadata_cell_type) +
  guides(colour = guide_legend(override.aes = list(size = 4))) + 
  coord_fixed(1.1)

ggsave(output_file, treatment_gg, dpi = 500, height = 7, width = 10)

print(treatment_gg)


# Generating UMAP comparing failing and healthy with DMSO and Drug_X
options(repr.plot.width = 12, repr.plot.height = 10)  # Adjust width and height as desired
output_file <- paste0(output_umap_file, "_treatment_dmso_drugx.png")

treatment_gg <- ggplot(umap_cp_df, aes(x = UMAP0, y = UMAP1)) +
  geom_point(aes(color = Metadata_treatment), size = 0.4, alpha = 0.9) +
  theme_bw() +
  scale_color_manual(
    name = "Treatment",
    values = c("DMSO" = "#BA5A31", "drug_x" = "#69DC9E")
  ) +
  facet_grid(~ Metadata_cell_type) +
  guides(colour = guide_legend(override.aes = list(size = 4))) + 
  coord_fixed(1.1)


ggsave(output_file, treatment_gg, dpi = 500, height = 7, width = 10)


print(treatment_gg)


# Filter the data frames
failing_drug_x_df <- umap_cp_df %>% filter(Metadata_cell_type == "failing" & Metadata_treatment == "drug_x")
failing_DMSO_df <- umap_cp_df %>% filter(Metadata_cell_type == "failing" & Metadata_treatment == "DMSO")
healthy_DMSO_df <- umap_cp_df %>% filter(Metadata_cell_type == "healthy" & Metadata_treatment == "DMSO")
failing_TGFRi_df <- umap_cp_df %>% filter(Metadata_cell_type == "failing" & Metadata_treatment == "TGFRi")


print(paste0("failing", " + ", "drug_x"))
print(paste0("failing", "+", "DMSO"))
print(paste0("healthy", " + ", "DMSO"))
print(paste0("failing", " + ", "TGFRi"))


# UMAP Plot healthy w/ DMSO and failing w/ treatment
# plot 1 -> failing w/ Drug_X
# plot 2 -> healthy w/ DMSO
options(repr.plot.width = 12, repr.plot.height = 10)  # Adjust width and height as desired
output_file <- paste0(output_umap_file, "_treatment_health_vs_failing.png")

# plot 
plot1 <- ggplot(failing_drug_x_df, aes(x = UMAP0, y = UMAP1)) +
  geom_point(aes(color = Metadata_treatment), size = 0.6, alpha = 0.9) +
  theme_bw() +
  scale_color_manual(
    name = "Treatment",
    values = c("drug_x" = "#69DC9E")
  ) +
  guides(colour = guide_legend(override.aes = list(size = 4))) +
  labs(title = "Failing") +
  ylim(min(umap_cp_df$UMAP1), max(umap_cp_df$UMAP1))  # Set the same y-axis limits as umap_cp_df

plot2 <- ggplot(healthy_DMSO_df, aes(x = UMAP0, y = UMAP1)) +
  geom_point(aes(color = Metadata_treatment), size = 0.6, alpha = 0.9) +
  theme_bw() +
  scale_color_manual(
    name = "Treatment",
    values = c("DMSO" = "#BA5A31")
  ) +
  guides(colour = guide_legend(override.aes = list(size = 4))) +
  labs(title = "Healthy") +
  ylim(min(umap_cp_df$UMAP1), max(umap_cp_df$UMAP1))  # Set the same y-axis limits as umap_cp_df

# merging both plots into one figure
merged_plot <- grid.arrange(
  plot1, plot2, ncol = 2,
  widths = c(1, 1), heights = c(1, 1)
)

# saving image
ggsave(output_file, merged_plot, dpi = 500, height = 12, width = 10)


# Set the size of each plot
# UMAP Plot healthy w/ DMSO and failing w/ treatment
output_file <- paste0(output_umap_file, "_failing_DMSO_and_DrugX_w_healthy_DMSO.png")

# combine all dfs
combined_df <- rbind(
  failing_drug_x_df %>% mutate(Group = "Failing with Drug_X"),
  healthy_DMSO_df %>% mutate(Group = "Healthy with DMSO"),
  failing_DMSO_df %>% mutate(Group = "Failing with DMSO")
)

# Plot the combined data with facets
merged_plot <- ggplot(combined_df, aes(x = UMAP0, y = UMAP1, color = Metadata_treatment)) +
  geom_point(size = 0.6, alpha = 0.9) +
  theme_bw() +
  scale_color_manual(
    name = "Treatment",
    values = c("drug_x" = "#69DC9E", "DMSO" = "#BA5A31")
  ) +
  facet_grid(~ Group, scales = "free") + 
  theme(aspect.ratio = 1)


# Saving image
ggsave(output_file, merged_plot, dpi = 500, height = 5, width = 12)

print(merged_plot)


# UMAP merged plots (from above)
output_file <- paste0(output_umap_file, "_failing_DMSO_and_DrugX_w_healhy_DMSO_Merged.png")

# Create a combined data frame with the "Group" column
combined_df <- rbind(
  failing_drug_x_df %>% mutate(Group = "failing + drug_x"),
  healthy_DMSO_df %>% mutate(Group = "healthy + DMSO"),
  failing_DMSO_df %>% mutate(Group = "failing + DMSO")
)

# Plot the combined data without facets
merged_plot <- ggplot(combined_df, aes(x = UMAP0, y = UMAP1, color = Group)) +
  geom_point(size = 1.0 , alpha = 0.5) +  # Adjust alpha value to make points more transparent
  labs(title = "UMAP of single-cells demonstrating a shift from failing to\nhealthy cells after failing cells are treated with drug_x") + 
  theme_bw(base_size = 22) +  # Set the base font size to 12
  scale_color_manual(
    name = "Group",
    values = c("failing + drug_x" = "#69DC9E", "failing + DMSO" = "#BA5A31", "healthy + DMSO" = "#8269dc")
  ) +
  guides(colour = guide_legend(override.aes = list(size = 4))) +
  ylim(min(umap_cp_df$UMAP1), max(umap_cp_df$UMAP1)) +  # Set the same y-axis limits as umap_cp_df
  theme(
    legend.position = c(0.15, 0.88),  # Adjust the legend position (x, y)
    legend.background = element_rect(fill = "white", color = "black"),  # Add a white background to the legend
    legend.key = element_rect(color = "white"),  # Make legend key (color boxes) white

  )
# print and save plot
print(merged_plot)
ggsave(output_file, merged_plot, dpi = 500, height = 12, width = 12)


# combine all dfs with TGFRi
combined_df <- rbind(
  failing_TGFRi_df %>% mutate(Group = "Failing with TGFRi"),
  healthy_DMSO_df %>% mutate(Group = "Healthy with DMSO"),
  failing_DMSO_df %>% mutate(Group = "Failing with DMSO")
)

# UMAP merged plots (from above)
output_file <- paste0(output_umap_file, "_failing_DMSO_and_TGFRi_w_healthy_DMSO_Merged.pdf")

# Create a combined data frame with the "Group" column
combined_df <- rbind(
  failing_TGFRi_df %>% mutate(Group = "failing + TGFRi"),
  healthy_DMSO_df %>% mutate(Group = "healthy + DMSO"),
  failing_DMSO_df %>% mutate(Group = "failing + DMSO")
)

# Plot the combined data without facets
merged_TGFRi_plot <- ggplot(combined_df, aes(x = UMAP0, y = UMAP1, color = Group)) +
  geom_point(size = 1.0 , alpha = 0.5) +  # Adjust alpha value to make points more transparent
  labs(title = "UMAP of single-cells demonstrating a shift from failing to\nhealthy cells after failing cells are treated with TGFRi") + 
  theme_bw(base_size = 22) +  # Set the base font size to 12
  scale_color_manual(
    name = "Group",
    values = c("failing + TGFRi" = "#69DC9E", "failing + DMSO" = "#BA5A31", "healthy + DMSO" = "#8269dc")
  ) +
  guides(colour = guide_legend(override.aes = list(size = 4))) +
  ylim(min(umap_cp_df$UMAP1), max(umap_cp_df$UMAP1)) +  # Set the same y-axis limits as umap_cp_df
  theme(
    legend.position = c(0.15, 0.88),  # Adjust the legend position (x, y)
    legend.background = element_rect(fill = "white", color = "black"),  # Add a white background to the legend
    legend.key = element_rect(color = "white"),  # Make legend key (color boxes) white

  )
# print and save plot
print(merged_TGFRi_plot)
ggsave(output_file, merged_TGFRi_plot, dpi = 500, height = 12, width = 12)


# UMAP merged plots (from above)
output_file <- paste0(output_umap_file, "_failing_DMSO_and_DrugX_w_healhy_DMSO_Merged.png")

# Create a combined data frame with the "Group" column
combined_df <- rbind(
  failing_drug_x_df %>% mutate(Group = "failing + drug_x"),
  healthy_DMSO_df %>% mutate(Group = "healthy + DMSO"),
  failing_DMSO_df %>% mutate(Group = "failing + DMSO")
)

# Plot the combined data without facets
merged_plot <- ggplot(combined_df, aes(x = UMAP0, y = UMAP1, color = Group)) +
  geom_point(size = 1.0 , alpha = 0.5) +  # Adjust alpha value to make points more transparent
  labs(title = "UMAP of single-cells demonstrating a shift from failing to\nhealthy cells after failing cells are treated with drug_x") + 
  theme_bw(base_size = 22) +  # Set the base font size to 12
  scale_color_manual(
    name = "Group",
    values = c("failing + drug_x" = "#69DC9E", "failing + DMSO" = "#BA5A31", "healthy + DMSO" = "#8269dc")
  ) +
  guides(colour = guide_legend(override.aes = list(size = 4))) +
  ylim(min(umap_cp_df$UMAP1), max(umap_cp_df$UMAP1)) +  # Set the same y-axis limits as umap_cp_df
  theme(
    legend.position = c(0.15, 0.88),  # Adjust the legend position (x, y)
    legend.background = element_rect(fill = "white", color = "black"),  # Add a white background to the legend
    legend.key = element_rect(color = "white"),  # Make legend key (color boxes) white

  )
# print and save plot
print(merged_plot)
ggsave(output_file, merged_plot, dpi = 500, height = 12, width = 12)


# UMAP merged plots (from above)
output_file <- paste0(output_umap_file, "_failing_DMSO_and_DrugX_w_healthy_DMSO_Merged_remove_edge.png")

# Define the edge buffer size
edge_buffer <- 250

# Create a combined data frame with the "Group" column
combined_df <- rbind(
failing_drug_x_df %>% mutate(Group = "failing + drug_x"),
healthy_DMSO_df %>% mutate(Group = "healthy + DMSO"),
failing_DMSO_df %>% mutate(Group = "failing + DMSO")
)

# Filter rows that are not close to the image edges
combined_df <- combined_df %>%
filter(
Metadata_Cells_Location_Center_X > edge_buffer,
Metadata_Cells_Location_Center_X < (1104 - edge_buffer),  
Metadata_Cells_Location_Center_Y > edge_buffer,
Metadata_Cells_Location_Center_Y < (1105 - edge_buffer) 
)

# Plot the combined data without facets
merged_plot <- ggplot(combined_df, aes(x = UMAP0, y = UMAP1, color = Group)) +
geom_point(size = 1.0, alpha = 0.5) +
labs(title = "UMAP of single-cells demonstrating a shift from failing to\nhealthy cells after failing cells are treated with drug_x") + 
theme_bw(base_size = 22) +
scale_color_manual(
name = "Group",
values = c("failing + drug_x" = "#69DC9E", "failing + DMSO" = "#BA5A31", "healthy + DMSO" = "#8269dc")
) +
guides(colour = guide_legend(override.aes = list(size = 4))) +
ylim(min(umap_cp_df$UMAP1), max(umap_cp_df$UMAP1)) +
theme(
legend.position = c(0.15, 0.88),  # Adjust the legend position (x, y)
legend.background = element_rect(fill = "white", color = "black"),
legend.key = element_rect(color = "white")
)
# print and save plot
print(merged_plot)
ggsave(output_file, merged_plot, dpi = 500, height = 12, width = 12)


# Only extracting DMSO only entries
dmso_df <- subset(umap_cp_df, Metadata_treatment == "DMSO")
dmso_failing_df <- subset(dmso_df, Metadata_cell_type == "failing")
dmso_healthy_df <- subset(dmso_df, Metadata_cell_type == "healthy")


# plotting UMAP failing and non-failing DMSO only 
options(repr.plot.width = 12, repr.plot.height = 10)  # Adjust width and height as desired
output_file <- paste0(output_umap_file, "_treatment_F_H_dmso_.png")

treatment_gg <- ggplot(umap_cp_df, aes(x = UMAP0, y = UMAP1)) +
  geom_point(aes(color = Metadata_treatment), size = 0.4, alpha = 0.9) +
  theme_bw() +
  scale_color_manual(
    name = "Treatment",
    values = c("DMSO" = "#BA5A31")
  ) +
  facet_grid(~ Metadata_cell_type) +
  guides(colour = guide_legend(override.aes = list(size = 4))) + 
  coord_fixed(1.1)


ggsave(output_file, treatment_gg, dpi = 500, height = 7, width = 10)


print(treatment_gg)


# plotting UMAP failing and non-failing DMSO only 
options(repr.plot.width = 12, repr.plot.height = 10)  # Adjust width and height as desired
output_file <- paste0(output_umap_file, "_treatment_F_H_dmso_cleaned.png")

treatment_gg <- ggplot(dmso_df, aes(x = UMAP0, y = UMAP1)) +
  geom_point(aes(color = Metadata_treatment), size = 0.4, alpha = 0.9) +
  theme_bw() +
  scale_color_manual(
    name = "Treatment",
    values = c("DMSO" = "#BA5A31")
  ) +
  facet_grid(~ Metadata_cell_type) +
  guides(colour = guide_legend(override.aes = list(size = 4))) + 
  coord_fixed(1.2)

ggsave(output_file, treatment_gg, dpi = 500, height = 7, width = 10)

print(treatment_gg)


# output name 
output_file <- paste0(output_umap_file, "_treatment_F_H_dmso_cleaned_merged.png")

# Create a combined data frame with the "Group" column
combined_df <- rbind(
  dmso_failing_df %>% mutate(Group = "failing + DMSO"),
  dmso_healthy_df %>% mutate(Group = "healthy + DMSO")
)

merged_plot <- ggplot(combined_df, aes(x = UMAP0, y = UMAP1, color = Group)) +
  geom_point(size = 1.0, alpha = 0.5) +  # Adjust alpha value to make points more transparent
  theme_bw(base_size = 20) +  # Set the base font size to 20 (adjust as desired)
  scale_color_manual(
    name = "Group",
    values = c("failing + DMSO" = "#BA5A31", "healthy + DMSO" = "#8269dc")
  ) +
  guides(colour = guide_legend(override.aes = list(size = 4))) +
  ylim(min(umap_cp_df$UMAP1), max(umap_cp_df$UMAP1))

print(merged_plot)
ggsave(output_file, merged_plot, dpi = 500, height = 12, width = 12)


# Remove decimal from nuclei object number to avoid issues with merging
umap_cp_df <- umap_cp_df %>%
  mutate(Metadata_Nuclei_Number_Object_Number = gsub("\\.0", "", as.character(Metadata_Nuclei_Number_Object_Number)))

# Load in parquet file with nearest neighbors per single cell (not normalized)
neighbors_df <- arrow::read_parquet(
    "../../../3.process_cfret_features/data/single_cell_profiles/localhost230405150001_sc_annotated.parquet"
)

desired_columns <- c("Cells_Neighbors_NumberOfNeighbors_Adjacent", 
                     "Metadata_Well", 
                     "Metadata_Site", 
                     "Metadata_Nuclei_Number_Object_Number")

neighbors_df <- neighbors_df[, desired_columns]

neighbors_df <- neighbors_df %>%
  mutate(Metadata_Nuclei_Number_Object_Number = as.character(Metadata_Nuclei_Number_Object_Number))

# Assuming you want to merge based on the specified columns
merge_columns <- c("Metadata_Well", "Metadata_Site", "Metadata_Nuclei_Number_Object_Number")

# Merge number of neighbors based on cell to umap dataframe
merged_data <- merge(neighbors_df, umap_cp_df, by = merge_columns)

umap_cell_neighbors_figure <- 
  ggplot(merged_data, aes(x = UMAP0, y = UMAP1)) +
  geom_point(aes(color = Cells_Neighbors_NumberOfNeighbors_Adjacent), alpha = 0.5) +
  theme_bw() +
  scale_color_continuous(
    name = "Number of\nCell Neighbors", 
    low = "light blue",
    high = "blue"
  ) +
  labs(title = "UMAP of Plate 3 Morphology Space Comparing\nBetween Number of Neighbors per Cell", x = "UMAP0", y = "UMAP1") +
  theme(
    # Increase title size
    plot.title = element_text(size = 20),
    
    # Increase axis text size
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    
    # Increase legend size
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
  ) +
  coord_fixed(1.1)

# saving image
output_file <- file.path(paste0(output_umap_file, "_number_neighbors_cells.png"))
ggsave(output_file, umap_cell_neighbors_figure, dpi = 500, height = 10, width = 10)


print(umap_cell_neighbors_figure)


# Set directory and file structure
umap_dir <- file.path("results")
umap_file <- "UMAP_localhost230405150001_fs_filtered.tsv.gz"
umap_path <- file.path(umap_dir, umap_file)

# Load in the umap data
filtered_df <- readr::read_tsv(
    umap_path,
    col_types = readr::cols(
        .default = "c",
        "UMAP0" = "d",
        "UMAP1" = "d"
    )
)
cell_count_df <- filtered_df %>%
    dplyr::group_by(Metadata_Well) %>%
    dplyr::count() %>%
    dplyr::rename(Metadata_Cell_Count = n)

umap_filtered_df <- filtered_df %>%
    dplyr::left_join(cell_count_df, by = "Metadata_Well") %>%
    mutate_all(~ifelse(is.na(.), "None", .))

head(umap_filtered_df)


# Filter the data frames from the filtered data (only not clustered single cells)
filtered_failing_drug_x_df <- umap_filtered_df %>% filter(Metadata_cell_type == "failing" & Metadata_treatment == "drug_x")
filtered_failing_DMSO_df <- umap_filtered_df %>% filter(Metadata_cell_type == "failing" & Metadata_treatment == "DMSO")
filtered_healthy_DMSO_df <- umap_filtered_df %>% filter(Metadata_cell_type == "healthy" & Metadata_treatment == "DMSO")

# UMAP merged plots (from above)
output_file <- paste0(output_umap_file, "_filtered_failing_DMSO_and_DrugX_w_healhy_DMSO_Merged.png")

# Create a combined data frame with the "Group" column
filtered_combined_df <- rbind(
  filtered_failing_drug_x_df %>% mutate(Group = "failing + drug_x"),
  filtered_healthy_DMSO_df %>% mutate(Group = "healthy + DMSO"),
  filtered_failing_DMSO_df %>% mutate(Group = "failing + DMSO")
)

# Plot the combined data without facets
merged_filtered_plot <- ggplot(filtered_combined_df, aes(x = UMAP0, y = UMAP1, color = Group)) +
  geom_point(size = 2.0 , alpha = 0.5) +  # Adjust alpha value to make points more transparent
  labs(title = "UMAP of non-clustered single-cells comparing\ntreatment and cell type ") + 
  theme_bw(base_size = 22) +  # Set the base font size to 12
  scale_color_manual(
    name = "Group",
    values = c("failing + drug_x" = "#69DC9E", "failing + DMSO" = "#BA5A31", "healthy + DMSO" = "#8269dc")
  ) +
  guides(colour = guide_legend(override.aes = list(size = 4))) +
  ylim(min(filtered_combined_df$UMAP1), max(filtered_combined_df$UMAP1)) +  # Set the same y-axis limits as umap_cp_df
  theme(
    legend.position = c(0.15, 0.88),  # Adjust the legend position (x, y)
    legend.background = element_rect(fill = "white", color = "black"),  # Add a white background to the legend
    legend.key = element_rect(color = "white"),  # Make legend key (color boxes) white

  )
# print and save plot
print(merged_filtered_plot)
ggsave(output_file, merged_filtered_plot, dpi = 500, height = 12, width = 12)


umap_filtered_cell_neighbors_figure <- 
  ggplot(umap_filtered_df, aes(x = UMAP0, y = UMAP1)) +
  geom_point(aes(color = Metadata_Neighbors_Adjacent), alpha = 0.5, size = 2) +
  theme_bw() +
  labs(title = "UMAP of Filtered Plate 3 Morphology Space Comparing\nBetween Number of Neighbors per Cell", x = "UMAP0", y = "UMAP1") +
  scale_color_discrete(name = "Number of\nNeighbors\nAdjacent")

# saving image
output_file <- file.path(paste0(output_umap_file, "_filtered_number_neighbors_cells.png"))
ggsave(output_file, umap_filtered_cell_neighbors_figure, dpi = 500, height = 8, width = 8)


print(umap_filtered_cell_neighbors_figure)



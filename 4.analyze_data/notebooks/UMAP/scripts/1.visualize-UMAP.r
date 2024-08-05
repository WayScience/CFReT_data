suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))

# Set directory and file structure
umap_dir <- file.path("results")
umap_files <- list.files(umap_dir, full.names = TRUE)
umap_files <- grep("KK22-05-198", umap_files, value = TRUE)
print(umap_files)

output_fig_dir <- file.path("figures")
umap_prefix <- "UMAP_"
plate_suffix <- "_sc_feature_selected.tsv.gz"

# Define output figure paths
output_umap_files <- list()
for (umap_file in umap_files) {
    # Manipulate the 
    plate <- stringr::str_remove(
        stringr::str_remove(
            unlist(
                strsplit(umap_file, "/")
            )[2],
            umap_prefix
        ),
        plate_suffix
    )
    output_umap_files[plate] <- file.path(
        output_fig_dir,
        paste0(umap_prefix, plate)
    )
}
        
print(output_umap_files)

plate

# Load data
umap_cp_df <- list()
for (plate in names(output_umap_files)) {
    # Find the umap file associated with the plate
    umap_file <- umap_files[stringr::str_detect(umap_files, plate)]
    
    # Load in the umap data
    df <- readr::read_tsv(
        umap_file,
        col_types = readr::cols(
            .default = "d",
            "Metadata_WellRow" = "c",
            "Metadata_WellCol" = "c",
            "Metadata_heart_number" = "c",
            "Metadata_treatment" = "c",
            "Metadata_dose" = "d",
            "Metadata_dose_unit" = "c",
            "Metadata_ImageNumber" = "c",
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

    umap_cp_df[[plate]] <- df %>%
        dplyr::left_join(cell_count_df, by = "Metadata_Well")
}

# Order the dose column low to high
dose_unit = "uM"
for (plate in names(umap_cp_df)) {
    # Pull the associated dataframe
    df <- umap_cp_df[[plate]]
    
    # Ensure doses are numeric and add dose unit
    numeric_doses <- as.numeric(df$Metadata_dose)
    
    # Append the dose unit
    df$Metadata_dose <- paste0(numeric_doses, dose_unit)
    
    # Extract unique dose order and sort, then recompile
    dose_order <- paste0(
        sort(unique(numeric_doses), na.last = TRUE),
        dose_unit
    )
    
    # Perform the factor ordering
    df$Metadata_dose <- factor(df$Metadata_dose, levels = dose_order)
    
    # Update the original dataframe in the list
    umap_cp_df[[plate]] <- df
}


for (plate in names(umap_cp_df)) {
    # Treatment UMAP
    output_file <- output_umap_files[[plate]]
    output_file <- paste0(output_file, "_treatment.png")
    
    treatment_gg <- (
        ggplot(umap_cp_df[[plate]], aes(x = UMAP0, y = UMAP1))
        + geom_point(
            aes(color = Metadata_treatment), size = 0.4, alpha = 0.7
        )
        + theme_bw()
        + scale_color_manual(
            name = "Treatment",
            values = c("DMSO" = "#BA5A31", "drug_x" = "#69DC9E")
        )
    )
    
    ggsave(output_file, treatment_gg, dpi = 500, height = 6, width = 6)
    
    # Dose UMAP
    output_file <- output_umap_files[[plate]]
    output_file <- paste0(output_file, "_dose_facet.png")
    
    umap_dose_gg <- (
        ggplot(umap_cp_df[[plate]], aes(x = UMAP0, y = UMAP1))
        + geom_point(
            aes(color = n), size = 0.4, alpha = 0.7
        )
        + theme_bw()
        + facet_wrap("~Metadata_dose")
        + theme(
            strip.background = element_rect(colour = "black", fill = "#fdfff4")
        )
        + scale_color_continuous(name = "Number of\nsingle cells\nper well")
    )
    
    ggsave(output_file, umap_dose_gg, dpi = 500, height = 6, width = 6)
    
}

# Plot with multiple hearts
plate_with_multiple_hearts <- c("localhost220512140003_KK22-05-198")

for (plate in plate_with_multiple_hearts) {
    # Heart UMAP
    output_file <- output_umap_files[[plate]]
    output_file <- paste0(output_file, "_hearts.png")
    
    heart_gg <- (
        ggplot(umap_cp_df[[plate]], aes(x = UMAP0, y = UMAP1))
        + geom_point(
            aes(color = Metadata_heart_number), size = 0.4, alpha = 0.7
        )
        + theme_bw()
    )
    
    ggsave(output_file, heart_gg, dpi = 500, height = 6, width = 6)
    }

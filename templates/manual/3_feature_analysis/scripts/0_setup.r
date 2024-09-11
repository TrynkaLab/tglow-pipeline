

# Load the core tglow functions
tglow.core <- c(list.files(paste0("/software/teamtrynka/tglow-core/", "r-core"), full.names=T),
list.files(paste0("/software/teamtrynka/installs/tglow-core/", "r-core"), full.names=T))

for (file in tglow.core) {
  source(file)
}

# Global parameters
suffix     <- "v1"
prefix     <- "tglow_panel_"
path       <- paste0("../../2_feature_extraction/output/", prefix, suffix)


# Sheet with sample metadata in standardized format
sample.sheet <- "http://link/to/google/sheet"

# Sheet with cell and feature filters following standardized format
filter.sheet <- "http://link/to/google/sheet"
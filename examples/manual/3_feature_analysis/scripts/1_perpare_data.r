project <- ""
setwd(paste0("/lustre/scratch125/humgen/projects/cell_activation_tc/projects/",project,"/3_feature_analysis/scripts"))

# Read project variables
source("0_setup.r")

library(googlesheets4)

# Read tflow data folder
output <- tglow.read.dir(path, pat.img = "_images.tsv")

# Don't take along ID features in the QC
output$features$analyze <- T
output$features$analyze[grep("Number", output$features$category)] <- F
output$features$analyze[grep("Parent", output$features$category)] <- F
output$features$analyze[grep("Children", output$features$category)] <- F
output$features$analyze[grep("ObjectNumber", output$features$measurement)] <- F
output$features$analyze[grep("Object_Number", output$features$measurement)] <- F
output$features$analyze[output$features$type == "character"] <- F

# Samplesheet
sample.meta <- read_sheet(sample.sheet)

sample.meta$batch_well <- paste0(sample.meta$plate, "_", sample.meta$well_name)
output$meta$plate_well <- paste0(output$meta$Metadata_plate, "_", output$meta$Metadata_well)

imgs.to.keep    <- output$meta$plate_well %in% sample.meta$batch_well
output          <- tglow.filter.img.apply(output, imgs.to.keep)

# Add the sample level meta to the cell metadata
output$meta           <- merge(output$meta, sample.meta, by.x="plate_well", by.y = "batch_well", all.x=T, sort=F)
rownames(output$meta) <- output$meta$ImageNumber_Global

save(output, file=paste0("../output/raw_", prefix, suffix, ".RData"))

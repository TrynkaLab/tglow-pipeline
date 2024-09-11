project <- ""
setwd(paste0("/lustre/scratch125/humgen/projects/cell_activation_tc/projects/",project,"/3_feature_analysis/scripts"))

# Read project variables
source("0_setup.r")

library(cowplot)
library(googlesheets4)
library(matrixStats)
library(irlba)
library(magick)
library(Rfast)
library(ggrastr)

# ------------------------------------------------------------------------------
# Image level QC
# ------------------------------------------------------------------------------
load(file=paste0("../output/raw_", prefix, suffix, ".RData"))

# Select features to run PCA on
features <- c(grep("ImageQuality_", colnames(output$meta), value=T),
              grep("Correlation_", colnames(output$meta), value=T),
              grep("Intensity_", colnames(output$meta), value=T))

# Remove features with zero variance, then center and scale
tmp   <- as.matrix(output$meta[,features])
vars  <- Rfast::colVars(tmp, na.rm=T)
tmp   <- scale(tmp[, vars!=0])

# Set NA's to zero, this is an ugly hack
tmp[is.na(tmp)] <- 0
tmp <- na.omit(tmp)

#-------------------------------------------------------------------------------
# Run PCA per plate to identify outlier images
plots          <- list()
outlier.images <- rep(FALSE, nrow(tmp))
npc            <- 5

pc1            <- 1
pc2            <- 2

# Determine qc grouping
cur.meta          <- output$meta[rownames(tmp),]
cur.meta$qc_group <- paste0(output$meta$Metadata_plate, <QC GROUP HERE>)

for (plate in unique(cur.meta[,"qc_group"])) {
  cur <- cur.meta[,"qc_group"] == plate
  pca <- prcomp(tmp[cur,])

  outlier <- rep(0, sum(cur))
  for (i in 1:npc) {
    #cur.pc.scaled <- scale(pca$x[,i])
    cur.pc.scaled <- tglow.mod.zscore(pca$x[,i])
    outlier[abs(cur.pc.scaled) > 3.5] <- outlier[abs(cur.pc.scaled) > 3.5]+ 1
  }
  outlier <- outlier > 0
  outlier.images[cur][outlier] <- T

  plots[[plate]] <- theme.plain(tglow.plot.xy(tglow.mod.zscore(pca$x[,pc1]),
                                        tglow.mod.zscore(pca$x[,pc2]),
                                        xlab=paste0("PC", pc1),
                                        ylab=paste0("PC", pc2),
                                        do.lm=F,
                                        col=outlier) +
                                  scale_color_viridis_d(name="outlier") +
                                  ggtitle(plate) +
                                  geom_hline(yintercept=c(-3.5, 3.5), lty=2) +
                                  geom_vline(xintercept=c(-3.5, 3.5), lty=2),
                                legend=F)
}

# PCA plot
pdf(width=15, height=4*(length(unique(cur.meta$qc_group))/6), file=paste0("../output/plots/",prefix,"pca_outlier_image_outlier_detection_",suffix,".pdf"))
plot_grid(plotlist=plots, ncol=6)
dev.off()

# How many images are removed
sum(outlier.images)
ncells.pre <- nrow(output$cells)

# Apply filter
output <- tglow.filter.img.apply(output, !outlier.images)

# How many cells are removed
sum(ncells.pre - nrow(output$cells))
 
# Cleanup enviroment
rm(cur, pca, tmp, outlier, plots, plate, cur.pc.scaled)

# ------------------------------------------------------------------------------
# QC
# ------------------------------------------------------------------------------
# Apply QC filters
feature.filters <- data.frame(read_sheet(filter.sheet, sheet="feature_filters"))
cell.filters    <- data.frame(read_sheet(filter.sheet, sheet="cell_filters"))

# Filter features
filtered.features <- tglow.filter.features.calc(output, feature.filters)
output            <- tglow.filter.features.apply(output, filtered.features)

# Filter cells
filtered.cells  <- tglow.filter.cells.calc(output, cell.filters)
cells.removed   <- nrow(output$cells) -  colSums(filtered.cells)
cells.removed
output          <- tglow.filter.cells.apply(output, filtered.cells)
 
# Save the QC'ed output for later re-use
save(output, filtered.cells, filtered.features, cells.removed, file=paste0("../output/",prefix,"simple_qced_data_cache_",suffix,".RData"))
 
# ------------------------------------------------------------------------------
# More strict QC
# ------------------------------------------------------------------------------
load(paste0("../output/",prefix,"simple_qced_data_cache_",suffix,".RData"))

scaled <- scale(output$cells[,output$features$analyze])
scaled <- scaled[,colSums(is.na(scaled)) == 0]
scaled <- scaled[,colVars(as.matrix(scaled))!=0]

pcs    <- irlba::prcomp_irlba(scaled, n=25, center=F, scale=F)
um     <- uwot::umap(pcs$x)
save(pcs, um,file=paste0("../output/",prefix,"simple_qced_data_cache_umap_",suffix,".RData"))

# Filter for modified z-score
output$meta$qc_group <- paste0(output$meta$Metadata_plate, "_", <QC GROUP HERE>)

grouping             <- output$meta[output$cells$Image_ImageNumber_Global, "qc_group"]

# Outlier features
outlier.featureset   <- c("cell_AreaShape_Eccentricity",
                        "cell_AreaShape_Area",
                        "cell_AreaShape_MajorAxisLength",
                        "cell_AreaShape_MinorAxisLength",
                        "cell_AreaShape_Compactness",
                        "cell_AreaShape_Extent")

outlier.featureset <- c(outlier.featureset, grep("cell_AreaShape_Zernike", output$features$id, value=T))
outlier.featureset <- c(outlier.featureset, grep("nucl_AreaShape_Zernike", output$features$id, value=T))

data              <- output$cells[,outlier.featureset]
fixed.outliers    <- !filter.mod.z.perc(data, thresh=3.5, thresh2=1, grouping=grouping)

test <- tglow.pca.outliers(output,
                           grouping=grouping,
                           features=output$features$analyze,
                           pc.thresh=0.75,
                           pc.max=25,
                           method="mod.z")
filtered.cells.z <- !test$outliers & !fixed.outliers

# ------------------------------------------------------------------------------
# UMAPs

fcols <- list()
fcols[["outlier"]]         <- !filtered.cells.z
fcols[["population"]]      <- output$meta[output$cells$Image_ImageNumber_Global, "population"]
fcols[["plate"]]           <- output$meta[output$cells$Image_ImageNumber_Global, "plate"]
# <MORE VARIABLES HERE>


pdf(width=6, height=5, file=paste0("../output/plots/",prefix,"all_cells_umaps_categorical_",suffix,".pdf"))

for (fcol in names(fcols)) {
  p1 <- tglow.plot.xy(um[,1],
                um[,2],
                do.lm=F,
                xlab="umap 1",
                ylab="umap 2",
                col=fcols[[fcol]],
                raster=T,
                alpha=0.25) + scale_color_viridis_d(name=fcol)
  plot(heme.plain(p1))
}

dev.off()

fcols <- list()
fcols[["Log2 major axis"]]        <- log2(output$cells$cell_AreaShape_MajorAxisLength)
# <MORE VARIABLES HERE>

pdf(width=6, height=5, file=paste0("../output/plots/",prefix,"all_cells_umaps_numeric_",suffix,".pdf"))
for (fcol in names(fcols)) {
  p1 <- tglow.plot.xy(um[,1],
                um[,2],
                do.lm=F,
                xlab="umap 1",
                ylab="umap 2",
                col=fcols[[fcol]],
                raster=T,
                alpha=0.25) + scale_color_viridis_c(name=fcol)
  plot(theme.plain(p1))
}
dev.off()

# -------------------------------------------------------------------------------
# Save results
nrow(output$cells) - sum(filtered.cells.z)

output.f <- tglow.filter.cells.apply(output, filtered.cells.z)
save(output.f, file=paste0("../output/",prefix,"qced_data_cache_",suffix,".RData"))
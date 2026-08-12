

I would like to add an additional process.

Currently cellprofiler produces either a per well zip file that has the features, which forms the final output of the pipeline or the same data, as an unzipped folder (if cpr_no_zip=true). I would like to change this to aggregate the contents of the zip/folder into a per plate parquet file. The zip file/folder is a archive with the following folder structure:

unzip 250307_185211-V_B02.zip 
Archive:  250307_185211-V_B02.zip
  inflating: features/250307_185211-V/B/2/tglow_cat_2d_v4_cell.txt  
  inflating: features/250307_185211-V/B/2/tglow_cat_2d_v4_cyto.txt  
  inflating: features/250307_185211-V/B/2/tglow_cat_2d_v4_Experiment.txt  
  inflating: features/250307_185211-V/B/2/tglow_cat_2d_v4_Image.txt  
  inflating: features/250307_185211-V/B/2/tglow_cat_2d_v4_ki67.txt  
  inflating: features/250307_185211-V/B/2/tglow_cat_2d_v4_lipids.txt  
  inflating: features/250307_185211-V/B/2/tglow_cat_2d_v4_memb.txt  
  inflating: features/250307_185211-V/B/2/tglow_cat_2d_v4_mitoNetwork.txt  
  inflating: features/250307_185211-V/B/2/tglow_cat_2d_v4_mito.txt  
  inflating: features/250307_185211-V/B/2/tglow_cat_2d_v4_nucl.txt  
  inflating: features/250307_185211-V/B/2/tglow_cat_2d_v4_Object relationships.txt  

An example is in 250307_185211-V_B02.zip. These zip files are organized plate/row/col/zipfile. 

where tglow_cat_2d_v4 is a prefix which changes between runs. _cell.txt describes the main objects. _Experiment.txt is a log that can be ignored. _Image contains the image level data, this is important. _Object relationships.txt contains a description of how _cell relates to child objects, this is encoded in the actual files, so you can ignore this file as well. The remaining files _<suffix>.txt then describe child objects whose parent is cell. I would like to arregate these child objects by mean, median or sum (configurable parameter). If it is a categorical variable and they are all the same, it should be set to that value, if they are different, put an NA. Child objects can be matched to the parent by the column "Parent_cell", but make it a configurable parameter. When averaging the child objects, you should add the <suffix> that described the type of child object to the column name seperated by an _. I have attached an R function that reads a zip file as context for how you could make the python version. One tricky thing is that the cell id's need to be globally unqiue. By default the cellprofiler ones are not. As you are running per plate and I want to keep the Id as short of a string as possible, you can maybe assign a unique index for each plate like P1, P2, P3 etc, as you cannot copy 1:1 the R function, which runs over the whole dataset. This needs to be consistent and robust to the order nextflow executes. 

To enable this you will need to do the following:

- Figure out how to collect the per well cellprofiler output into a plate level channel, add a unique index for each plate which can be used to make a globally unique id.
- Make a new process "concat_cellprofiler" which takes a channel of plate, well level zip files, outputs a parquet file for the object level data and for the image level data.
- Make a python script that reads the contents of the zip file/folder over the folder tree and produces a parquet per plate for the _cell objects and the image objects


#-------------------------------------------------------------------------------
#' Read a cell level fileset type B
#'
#'
#' @description Reads a .zip file with cellprofiler features, each file other then
#' _Image, _Experiment and _Object Relationships are assumed to be an object
#' Will match objects on order with an appropriate matching strategy
#'
#' @param prefix Path prefix to fileset
#' @param return.feature.meta Should the dataframe with feature metadata be added
#' @param add.global.id Should extra id columns be added that are globally unique
#' @param merging.strategy How to consolidate 1:many relationships between cell: children. Accepted values: 'mean', 'none'
#' @param pat.exp The suffix pattern to identify the experiment file
#' @param pat.img The suffix pattern to identify the image level data
#' @param pat.cells The suffix pattern to identify the cell level data
#' @param pat.orl The suffix pattern to identify object relationships
#' @param pat.others The pattern to use to extract the child object names from the filename. First regex group is used as object name
#' @param na.rm Should NA's be removed when applying merging.strategy
#' @param skip.orl Should object relationships be read (not used, can be quite large)
#' @param skip.children List of child object names to skip during reading.
#' @param verbose Should I be chatty?
#' @param fileset.id The global fileset id to add if add.global.id = TRUE
#'
#' @returns list with data frames:
#' - cells (cell level features)
#' - meta (image level features)
#' - objectRelations
#' - features (optional)
#' Output is NULL if no cells are detected
#'
#' @export
read_cellprofiler_fileset_b <- function(prefix,
                                        return.feature.meta = F,
                                        add.global.id = T,
                                        merging.strategy = "mean",
                                        parent.col = "Parent_cell",
                                        pat.exp = ".*Experiment.txt",
                                        pat.img = ".*Image.txt",
                                        pat.cells = ".*cell.txt",
                                        pat.orl = ".*Object relationships.txt",
                                        pat.others = "^.*_([a-zA-Z]+\\d*).txt$",
                                        na.rm = F,
                                        skip.children = NULL,
                                        verbose = F,
                                        fileset.id = NULL) {
  if (add.global.id) {
    if (is.null(fileset.id)) {
      stop("fileset.id cannot be NULL if add.global.id=T")
    }
    # assign("FILESET_ID", FILESET_ID + 1, envir = .GlobalEnv)
    global.prefix <- paste0("FS", fileset.id)
  }

  index <- unzip(prefix, list = T)
  index$FileName <- basename(index$Name)

  tmpdir <- tempdir()

  # Clean up the tmpdir
  unlink(paste0(tmpdir, "/features"), recursive = T)

  # Unzip into tmp folder
  unzip(prefix, exdir = tmpdir)

  cells <- data.table::fread(paste0(tmpdir, "/", index[grep(pat.cells, index$FileName), "Name"]), data.table = F, showProgress = FALSE)
  img <- data.table::fread(paste0(tmpdir, "/", index[grep(pat.img, index$FileName), "Name"]), data.table = F, showProgress = FALSE)

  if (nrow(cells) == 0) {
    warning("No cells detected for ", index[grep(pat.cells, index$FileName), "Name"], " returning NULL.")
    return(NULL)
  }

  exclude <- c(
    grep(pat.cells, index$FileName),
    grep(pat.img, index$FileName),
    grep(pat.orl, index$FileName),
    grep(pat.exp, index$FileName)
  )

  index <- index[!seq_len(nrow(index)) %in% exclude, ]

  colnames(cells) <- paste0("cell_", colnames(cells))
  rownames(cells) <- paste0(cells$cell_ImageNumber, "_", cells$cell_ObjectNumber)

  children <- list()
  if (nrow(index) > 0) {
    index$object <- gsub(pat.others, "\\1", index$FileName)

    for (i in seq_len(nrow(index))) {
      obj <- index[i, "object"]
      
      # Optionally skip adding this child object
      if (obj %in% skip.children) {
        next()
      }
      
      cur <- data.table::fread(paste0(tmpdir, "/", index[i, "Name"]), data.table = T, showProgress = FALSE)

      # Remove these columns from the merging strategy
      exclude.cols <- c("Group.1", grep("ObjectNumber", colnames(cur), value = T), grep("Number_Object_Number", colnames(cur), value = T))

      # If there are no cols, return NA
      if (nrow(cur) == 0) {
        # next()
        # If the file is empty, set these columns to NA
        cur <- as.data.frame(cur)
        colnames(cur) <- paste0(obj, "_", colnames(cur))
        cells[, c(colnames(cur), paste0(obj, "_QC_Object_Count"))] <- NA
        warning(paste0(obj, " assay for ", index[i, "Name"], " is empty. Returning NA for these cols."))
        
      } else if (merging.strategy == "mean") {
        if (verbose) cat("[DEBUG] ", as.character(index[i, ]), "\n")

        if (!parent.col %in% colnames(cur)) {
          stop(paste0("Parent col: '", parent.col, "' not found for child object '", obj, "'. Either update the files, or skip reading the child object with skip.children"))
        }

        cur <- cur[as.logical(cur[[parent.col]] != 0), ]
        colnames(cur) <- paste0(obj, "_", colnames(cur))
        selector <- paste0(cur[[paste0(obj, "_ImageNumber")]], "_", cur[[paste0(obj, "_", parent.col)]])
        counts <- table(selector)
        cur$Group.1 <- selector

        if (verbose) cat("[DEBUG] NA's in selector", sum(is.na(selector)), "\n")
        if (verbose) cat("[DEBUG] Slice: ", head(cur$Group.1), "\n")

        # Calculate the mean per group
        tmp <- as.data.frame(cur[, lapply(.SD, mean, na.rm = na.rm), by = Group.1, .SDcols = colnames(cur)[!colnames(cur) %in% exclude.cols]])
        rownames(tmp) <- tmp$Group.1

        # Add the object count as a sanity check
        tmp[, paste0(obj, "_QC_Object_Count")] <- counts[tmp$Group.1]
        tmp <- tmp[, !colnames(tmp) %in% exclude.cols]

        # Assign the columns to the output matrix
        cells[selector, colnames(tmp)] <- tmp[selector, ]
      } else {
        stop("Only valid merging strategy is 'mean'")
      }
    }
  }
  # Clean up the tmpdir
  unlink(paste0(tmpdir, "/features"), recursive = T)


  # Standardize ID's across filesets into the following format
  # FS#I#O#
  # where FS = file set, I = image within fileset and O = object within image
  # This makes it easier to match across data with a unique ID
  if (add.global.id) {
    # Cell level information
    #-----------
    cells <- add_global_ids(cells, fileset.id)

    # IMG image level information
    #-----------
    img[, "ImageNumber_Global"] <- paste0(global.prefix, "_I", img[, "ImageNumber"])
  }

  if (return.feature.meta) {
    feature.meta <- tglowr::get_feature_meta_from_names(colnames(cells))
  }

  out.list <- list(cells = cells, meta = img)

  if (return.feature.meta) {
    out.list <- c(out.list, list(features = feature.meta))
  }

  if (length(children) > 0) {
    out.list <- c(out.list, list(children = children))
  }

  return(out.list)
}
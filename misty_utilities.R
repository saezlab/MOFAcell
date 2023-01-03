# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

future::plan(future::multisession)

#' Catalog of MISTy utilities
#' 
run_misty_seurat <- function(visium.slide,
                             # Seurat object with spatial transcriptomics data.
                             view.assays,
                             # Named list of assays for each view.
                             view.features = NULL,
                             # Named list of features/markers to use.
                             # Use all by default.
                             view.types,
                             # Named list of the type of view to construct
                             # from the assay.
                             view.params,
                             # Named list with parameters (NULL or value)
                             # for each view.
                             spot.ids = NULL,
                             # spot IDs to use. Use all by default.
                             out.alias = "results",
                             # folder name for output
                             bypass.intra = FALSE
                             # Shall you model intraview?
) {
  
  mistyR::clear_cache()
  
  # Extracting geometry
  geometry <- GetTissueCoordinates(visium.slide,
                                   cols = c("row", "col"), scale = NULL
  )
  
  # Extracting data
  view.data <- map(view.assays,
                   extract_seurat_data,
                   geometry = geometry,
                   visium.slide = visium.slide
  )
  
  # Constructing and running a workflow
  build_misty_pipeline(
    view.data = view.data,
    view.features = view.features,
    view.types = view.types,
    view.params = view.params,
    geometry = geometry,
    spot.ids = spot.ids,
    out.alias = out.alias,
    bypass.intra = bypass.intra
  )
}


# Extracts data from an specific assay from a Seurat object
# and aligns the IDs to the geometry
extract_seurat_data <- function(visium.slide,
                                assay,
                                geometry) {
  print(assay)
  data <- GetAssayData(visium.slide, assay = assay) %>%
    as.matrix() %>%
    t() %>%
    as_tibble(rownames = NA)
  
  return(data %>% slice(match(rownames(.), rownames(geometry))))
}

# Filters data to contain only features of interest
filter_data_features <- function(data,
                                 features) {
  if (is.null(features)) features <- colnames(data)
  
  return(data %>% rownames_to_column() %>%
           select(rowname, all_of(features)) %>% rename_with(make.names) %>%
           column_to_rownames())
}

# Builds views depending on the paramaters defined
create_default_views <- function(data,
                                 view.type,
                                 view.param,
                                 view.name,
                                 spot.ids,
                                 geometry) {
  
  mistyR::clear_cache()
  
  view.data.init <- create_initial_view(data)
  
  if (!(view.type %in% c("intra", "para", "juxta"))) {
    view.type <- "intra"
  }
  
  if (view.type == "intra") {
    data.red <- view.data.init[["intraview"]]$data %>%
      rownames_to_column() %>%
      filter(rowname %in% spot.ids) %>%
      select(-rowname)
  } else if (view.type == "para") {
    view.data.tmp <- view.data.init %>%
      add_paraview(geometry, l = view.param)
    
    data.ix <- paste0("paraview.", view.param)
    data.red <- view.data.tmp[[data.ix]]$data %>%
      mutate(rowname = rownames(data)) %>%
      filter(rowname %in% spot.ids) %>%
      select(-rowname)
  } else if (view.type == "juxta") {
    view.data.tmp <- view.data.init %>%
      add_juxtaview(
        positions = geometry,
        neighbor.thr = view.param
      )
    
    data.ix <- paste0("juxtaview.", view.param)
    data.red <- view.data.tmp[[data.ix]]$data %>%
      mutate(rowname = rownames(data)) %>%
      filter(rowname %in% spot.ids) %>%
      select(-rowname)
  }
  
  if (is.null(view.param) == TRUE) {
    misty.view <- create_view(
      paste0(view.name),
      data.red
    )
  } else {
    misty.view <- create_view(
      paste0(view.name, "_", view.param),
      data.red
    )
  }
  
  return(misty.view)
}

# Builds automatic MISTy workflow and runs it
build_misty_pipeline <- function(view.data,
                                 view.features,
                                 view.types,
                                 view.params,
                                 geometry,
                                 spot.ids = NULL,
                                 out.alias = "default",
                                 bypass.intra =  FALSE) {
  
  # Adding all spots ids in case they are not defined
  if (is.null(spot.ids)) {
    spot.ids <- rownames(view.data[[1]])
  }
  
  # First filter the features from the data
  view.data.filt <- map2(view.data, view.features, filter_data_features)
  
  # Create initial view
  views.main <- create_initial_view(view.data.filt[[1]] %>%
                                      rownames_to_column() %>%
                                      filter(rowname %in% spot.ids) %>%
                                      select(-rowname))
  
  # Create other views
  view.names <- names(view.data.filt)
  
  all.views <- pmap(list(
    view.data.filt[-1],
    view.types[-1],
    view.params[-1],
    view.names[-1]
  ),
  create_default_views,
  spot.ids = spot.ids,
  geometry = geometry
  )
  
  pline.views <- add_views(
    views.main,
    unlist(all.views, recursive = FALSE)
  )
  
  
  # Run MISTy
  run_misty(pline.views, out.alias, cached = FALSE, bypass.intra = bypass.intra)
}

# Optimizer
# rename all views so they correspond to the pseudo headers or just remove the parameters from the header and from the file names

#' Function to get optimal results
#' Here, first we identify which l parameter is the best for each target
#' Then we create an optimal directory and copy specific information
#' 
#' @param out_dir_name = out_alias used in MISTy runs
#' @param ls = l parameters used in the MISTy runs
#' @return a folder out_dir_name_optim with optimal info
get_optimal = function(out_dir_name,
                       ls){
  
  system(paste0("mkdir ", out_dir_name,"_optim"))
  
  l <- ls
  
  perf <- l %>% map_dfr(function(p) {
    performance <- read_delim(paste0(out_dir_name, "_",p, "/performance.txt"),
                              delim = " ", col_types = cols()
    ) %>% distinct()
    
    #selection criterion maximum R2 performance improvement per marker
    performance %>% arrange(target) %>%
      mutate(impr = multi.R2 - intra.R2) %>% 
      dplyr::select(target,impr) %>% 
      column_to_rownames("target") %>% t %>% as.data.frame
  })
  
  # For each target get the maximum of improvement
  # distinct ls
  optimal.l = colnames(perf) %>% enframe(name = NULL, value = "target") %>% 
    mutate(l = apply(perf, 2, which.max) %>% (l)[.])
  
  ##save optimal results
  
  write.table(optimal.l,
              file = paste0(out_dir_name,"_optim/optim_l.txt"),
              col.names = T, row.names = F, quote = F, sep = "\t")
  
  # Copy the relevant documents to the new directory, deleting the l parameter info
  walk2(optimal.l$target, optimal.l$l, function(.x, .y) {
    files <- list.files(paste0(out_dir_name,"_", .y, "/"), 
                        paste0("importances_", .x, '*'), 
                        full.names = TRUE)
    
    files %>% walk(~file.copy(., paste0(out_dir_name,"_optim/", 
                                        str_replace(last(str_split(., "/")[[1]]), 
                                                    "([\\._][0-9]+)+\\.txt", ".txt"))))
  })
  
  #very suboptimal
  walk2(optimal.l$target, optimal.l$l, function(.x, .y) {
    performance <- read_delim(paste0(out_dir_name,"_", .y, "/performance.txt"),
                              delim = " ", col_types = cols()
    ) %>% distinct()
    
    coeff <- read_delim(paste0(out_dir_name,"_", .y, "/coefficients.txt"),
                        delim = " ", col_types = cols()
    ) %>% distinct()
    
    
    if(!file.exists(paste0(out_dir_name,"_optim/performance.txt"))){
      write(colnames(performance) %>% paste0(collapse=" "), 
            paste0(out_dir_name,"_optim/performance.txt"))
    }
    
    if(!file.exists(paste0(out_dir_name,"_optim/coefficients.txt"))){
      write(str_remove_all(colnames(coeff) %>% paste0(collapse=" "), "([\\._][0-9]+)+"), 
            paste0(out_dir_name,"_optim/coefficients.txt"))
    }
    
    write(performance %>% filter(target == .x) %>% unlist %>% unname %>% paste0(collapse=" "), 
          paste0(out_dir_name,"_optim/performance.txt"), append = T)
    
    write(coeff %>% filter(target == .x) %>% unlist %>% unname %>% paste0(collapse=" "), 
          paste0(out_dir_name,"_optim/coefficients.txt"), append = T)
    
  })
  
}









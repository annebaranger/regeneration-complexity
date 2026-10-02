plot_name="GER02"
tar_load(path_plot)
folders=c("regen","hesitation")
path_trees=file.path(path_plot,
                     paste0(plot_name,"-processing"),
                     "trees",folders)
is.segmented=sum(file.exists(path_trees))


if(is.segmented>0){
  tree_files <- list.files(path_trees, pattern = "\\.txt$", full.names = FALSE)
  # to_check=as.logical(rowSums(sapply(
  #   paste0("000", list_dbh_check_ger10, ".txt"),  
  #   function(x) grepl(pattern = x, x = tree_files))))
  tree_files_check <- tree_files #[to_check]
  
  # Storage
  DBHs_checked <- data.frame(
    File = character(),
    DBH = numeric(),
    R2  = numeric(),
    H_approx=numeric(),
    slice_height = numeric(),
    slice_thickness = numeric(),
    stringsAsFactors = FALSE
  )
  
  files_hand <- character()
  
  # Parameters to try
  slice_heights   <- c(1.3, 0.9, 0.6, 1.5)
  slice_thickness <- c(0.1, 0.25, 0.5)   # tweak as needed
  
  for (tree_file in tree_files_bis) {
    
    path_tree=file.path(path_trees, tree_file)[file.exists(file.path(path_trees, tree_file))]
    
    tree_check <- tryCatch(
      read_tree_pc(path =path_tree ),
      error = function(e) {
        message("❌ Error reading: ", tree_file, " | ", e$message)
        return(NULL)
      }
    )
    if (is.null(tree_check)) {
      files_hand <- c(files_hand, tree_file)
      next
    }
    
    H_approx=quantile(tree_check$Z,probs=0.999)-min(tree_check$Z)
    
    if(H_approx>4){
      D_out = diameter_slice_pc(
        pc = tree_check,
        slice_height = 1.3,
        slice_thickness = 0.1
      )
      if(D_out$diameter>0.5|is.na(D_out$diameter)|D_out$R2>0.0001){
        adjust = TRUE
      }else{
        DBHs_checked <- rbind(
          DBHs_checked,
          data.frame(
            File = tree_file,
            DBH = D_out$diameter,
            R2  = D_out$R2,
            H_approx=H_approx,
            slice_height = 1.3,
            slice_thickness = 0.1,
            stringsAsFactors = FALSE
          )
        )
        adjust=FALSE
      }
    }else{
      DBHs_checked <- rbind(
        DBHs_checked,
        data.frame(
          File = tree_file,
          DBH = "small",
          R2  = NA,
          H_approx=H_approx,
          slice_height = NA,
          slice_thickness = NA,
          stringsAsFactors = FALSE
        )
      )
      adjust=FALSE
    }
    
    if(adjust){
      cat("\n=============================\n")
      cat("Adjusting file:", tree_file, "\n")
      accepted <- FALSE
      
      for (t in slice_thickness) {
        if (accepted) break
        
        for (h in slice_heights) {
          
          D_out <- tryCatch(
            diameter_slice_pc(
              pc = tree_check,
              slice_height = h,
              slice_thickness = t,
              plot = FALSE,
              functional = FALSE
            ),
            error = function(e) {
              message("❌ Error at h=", h, ", t=", t, ": ", e$message)
              return(NULL)
            }
          )
          if (is.null(D_out)) next
          if(D_out$diameter>0.4|is.na(D_out$diameter)|D_out$R2>0.001) next
          
          
          cat("Try: height =", h,
              "| thickness =", t,
              "| diameter =", D_out$diameter,
              "| R2 =", D_out$R2, "\n")
          
          diameter_slice_pc(
            pc = tree_check,
            slice_height = h,
            slice_thickness = t,
            plot = TRUE,
            functional = FALSE
          )
          ans <- tolower(trimws(readline("Does the fit look good? (y/n): ")))
          
          if (ans %in% c("y", "yes")) {
            DBHs_checked <- rbind(
              DBHs_checked,
              data.frame(
                File = tree_file,
                DBH = D_out$diameter,
                R2  = D_out$R2,
                H_approx=H_approx,
                slice_height = h,
                slice_thickness = t,
                stringsAsFactors = FALSE
              )
            )
            accepted <- TRUE
            break
          }
        }
      }
      
      if (!accepted) {
        H_out=tree_height_pc(sample_frac(tree_check,0.2),plot=TRUE)

        cat("Try: height =", H_out$h,
            "| file =", tree_file)
        
        ans_1 <- tolower(trimws(readline("Is it regen ? (y/n): ")))
        
        if (ans_1 %in% c("y", "yes")) {
          DBHs_checked <- rbind(
            DBHs_checked,
            data.frame(
              File = tree_file,
              DBH = "unknown",
              R2  = "unknown",
              H_approx=H_out$h,
              slice_height = NA,
              slice_thickness = NA,
              stringsAsFactors = FALSE
            )
          )
        } else{
          files_hand <- c(files_hand, tree_file)
        }
      }
    }
  }
}

write.csv(DBHs_checked,paste0("output/dbh_seg_",plot_name,".csv"))
save(files_hand,file=paste0("output/not_reg_",plot_name,".RData"))

dtm=tar_read(dtm_GER02)
bbox=tar_read(bbox_GER02)
if(is.segmented>0){
  dbh_seg=read.csv(paste0("output/dbh_seg_",plot_name,".csv")) %>% 
    select(-any_of("X")) %>%
    mutate(across(-any_of("File"), as.numeric)) %>% 
    filter(DBH<0.12|is.na(DBH)) %>% 
    mutate(file_long=case_when(grepl("hesit",File)~file.path(path_plot,
                                                             paste0(plot_name,"-processing"),
                                                             "trees","hesitation",File),
                               grepl("regen",File)~file.path(path_plot,
                                                             paste0(plot_name,"-processing"),
                                                             "trees","regen",File)))
  
  files=dbh_seg$file_long
  
  dtm_rast <- rast(dtm)
  
  read_xyz <- function(fn) {
    dt <- read_tree_pc(fn)
    setDT(dt)
    dt <- dt[
      is.finite(X) &
        is.finite(Y) &
        is.finite(Z)
    ]
    dt[]
  }
  
  tree_list <- lapply(files, read_xyz)
  tree_list <- Filter(Negate(is.null), tree_list)
  merged <- rbindlist( tree_list,
                       use.names = TRUE,
                       fill = TRUE)
  las_merged=LAS(merged[, c("X", "Y", "Z")])
  
  las_clipped <- clip_rectangle(las_merged,
                        xleft   = bbox$xmin,
                        ybottom = bbox$ymin,
                        xright  = bbox$xmax,
                        ytop    = bbox$ymax)
  
  # load dtm 
  dtm_rast <- rast(dtm)
  
  # normalize pc height
  las_norm <- normalize_height(las_clipped, dtm_rast)
  
  las_norm@data$X <- las_norm@data$X - min(las_norm@data$X)
  las_norm@data$Y <- las_norm@data$Y - min(las_norm@data$Y)
  
  las_norm <- las_update(las_norm)
  
  # las_clip <- filter_poi(las_norm, Z < 5)
  
  path_las=paste0("output/",plot_name,"_pcnorm_regen.rdata")
  save(las_norm,file=path_las)
  return(path_las)
}


file.segmented<-function(path_plot,
                       plot_name,
                       folders=c("regen","hesitation")){
  path_trees=file.path(path_plot,
                       paste0(plot_name,"-processing"),
                       "trees",folders)
  is.segmented=sum(file.exists(path_trees))
  
  if(is.segmented>0){
    file=paste0("output/dbh_seg_",plot_name,".csv")
  }else file=NULL
  return(file)
}

get_pc_norm_regen<-function(path_plot,
                            plot_name,
                            file.segmented,
                            dtm,
                            bbox){
  if(!is.null(file.segmented)){
    
    ## gather trees txt files
    #------------------------
    dbh_seg=read.csv(file.segmented) %>% 
      select(-any_of("X")) %>%
      mutate(across(-any_of("File"), as.numeric)) %>% 
      filter(DBH<0.12|is.na(DBH)) %>% 
      mutate(file_long=case_when(grepl("hesit",File)~file.path(path_plot,
                                                               paste0(plot_name,"-processing"),
                                                               "trees","hesitation",File),
                                 grepl("regen",File)~file.path(path_plot,
                                                               paste0(plot_name,"-processing"),
                                                               "trees","regen",File)))
    
    files=dbh_seg$file_long
    
    ## load dtm
    dtm_rast <- rast(dtm)
    
    
    ## load trees files and merge in LAS
    #-----------------------------------
    read_xyz <- function(fn) {
      dt <- read_tree_pc(fn)
      setDT(dt)
      dt <- dt[
        is.finite(X) &
          is.finite(Y) &
          is.finite(Z)
      ]
      dt[]
    }
    
    tree_list <- lapply(files, read_xyz)
    tree_list <- Filter(Negate(is.null), tree_list)
    merged <- rbindlist( tree_list,
                         use.names = TRUE,
                         fill = TRUE)
    las_merged=LAS(merged[, c("X", "Y", "Z")])
    
    las_clipped <- clip_rectangle(las_merged,
                                  xleft   = bbox$xmin,
                                  ybottom = bbox$ymin,
                                  xright  = bbox$xmax,
                                  ytop    = bbox$ymax)
  
    
    # normalize pc height
    las_norm <- normalize_height(las_clipped, dtm_rast)
    
    las_norm@data$X <- las_norm@data$X - min(las_norm@data$X)
    las_norm@data$Y <- las_norm@data$Y - min(las_norm@data$Y)
    
    las_norm <- las_update(las_norm)
    
    # las_clip <- filter_poi(las_norm, Z < 5)
    
    path_las=paste0("output/",plot_name,"_pcnorm_regen.rdata")
    save(las_norm,file=path_las)
  }else{path_las=NULL}
  return(path_las)
}
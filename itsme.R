library(ITSMe)

path_trees="C:/Users/amob2/OneDrive - University of Cambridge/2. FLF project/ger10-processing/GER10-regeneration-seg/"


### Batch process height/diamter ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Tree height:
# Plot the tree point clouds (check for outliers)
# Specify dtm and r in case lower part of tree is not sampled
Hs <- plot_tree_height_pcs(
  PCs_path = "C:/Users/amob2/OneDrive - University of Cambridge/2. FLF project/ger10-processing/GER10-regeneration-seg/",
  extension = ".txt",
  OUT_path = "C:/Users/amob2/OneDrive - University of Cambridge/2. FLF project/ger10-processing/GER10-regeneration-seg/"
)

# DBH:
# try out different thresholdR2 and slice_thickness values when default values fail
DBHs <- plot_dbh_fit_pcs(
  PCs_path = "C:/Users/amob2/OneDrive - University of Cambridge/2. FLF project/ger10-processing/GER10-regeneration-seg/",
  extension = ".txt",
  OUT_path = "C:/Users/amob2/OneDrive - University of Cambridge/2. FLF project/ger10-processing/GER10-regeneration-seg/",
  thresholdR2 = 0.0025, slice_thickness = 0.2
)

save(DBHs,file = file.path(path_trees,"DBHs.Rdata"))
save(Hs[[c("File","Heights")]],file=file.path(path_trees,"Hs.Rdata"))


### Index manual check ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%
list_dbh_check_ger10<-c(16,21,33,36,38,39,41,43,45,47,50,53,56,57,58,63,64,68,69,77,78,86,93,94,97,98,103,113,114,117,118,119,125,126,127,132,134,135,145,146)


### Process by hand ####
#%%%%%%%%%%%%%%%%%%%%%%%

# Files to process
tree_files <- list.files(path_trees, pattern = "\\.txt$", full.names = FALSE)
to_check=as.logical(rowSums(sapply(
  paste0("000", list_dbh_check_ger10, ".txt"),  
  function(x) grepl(pattern = x, x = tree_files))))
tree_files_check <- tree_files[to_check]

# Storage
DBHs_checked <- data.frame(
  File = character(),
  DBH = numeric(),
  R2  = numeric(),
  slice_height = numeric(),
  slice_thickness = numeric(),
  stringsAsFactors = FALSE
)

files_hand <- character()

# Parameters to try
slice_heights   <- c(1.3, 0.9, 0.6, 1.5)
slice_thickness <- c(0.1, 0.25, 0.5)   # tweak as needed

for (tree_file in tree_files_check) {
  
  cat("\n=============================\n")
  cat("File:", tree_file, "\n")
  
  tree_check <- tryCatch(
    read_tree_pc(path = file.path(path_trees, tree_file)),
    error = function(e) {
      message("❌ Error reading: ", tree_file, " | ", e$message)
      return(NULL)
    }
  )
  if (is.null(tree_check)) {
    files_hand <- c(files_hand, tree_file)
    next
  }
  
  accepted <- FALSE
  
  for (t in slice_thickness) {
    if (accepted) break
    
    for (h in slice_heights) {
      
      D_out <- tryCatch(
        diameter_slice_pc(
          pc = tree_check,
          slice_height = h,
          slice_thickness = t,
          plot = TRUE,
          functional = FALSE
        ),
        error = function(e) {
          message("❌ Error at h=", h, ", t=", t, ": ", e$message)
          return(NULL)
        }
      )
      if (is.null(D_out)) next
      
      cat("Try: height =", h,
          "| thickness =", t,
          "| diameter =", D_out$diameter,
          "| R2 =", D_out$R2, "\n")
      
      ans <- tolower(trimws(readline("Does the fit look good? (y/n): ")))
      
      if (ans %in% c("y", "yes")) {
        DBHs_checked <- rbind(
          DBHs_checked,
          data.frame(
            File = tree_file,
            DBH = D_out$diameter,
            R2  = D_out$R2,
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
    files_hand <- c(files_hand, tree_file)
  }
}

### Data check in total ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%
Height_tot=data.frame(
  File = Hs$File,
  H = Hs$Heights,
  stringsAsFactors = FALSE
) 

DBH_tot=data.frame(
  File = DBHs$File[!to_check],
  DBH = DBHs$DBHs[!to_check],
  R2  = DBHs$R2s[!to_check],
  slice_height = 1.3,
  slice_thickness = 0.6,
  stringsAsFactors = FALSE
) 

DBH_recheck=DBH_tot %>% filter(is.nan(DBH))
DBH_tot=DBH_tot%>% filter(!is.nan(DBH))
files_hand=c(files_hand,DBH_recheck$File)

### Manual estimation ####
#%%%%%%%%%%%%%%%%%%%%%%%%%

# Make R2 character so we can store "est"
DBHs_checked$R2 <- as.character(DBHs_checked$R2)

ask_dbh <- function(prompt = "Enter estimated DBH (in same units as ITSMe output): ") {
  repeat {
    x <- readline(prompt)
    x <- gsub(",", ".", x)   # in case of comma decimals
    if (!is.na(suppressWarnings(as.numeric(x)))) {
      return(as.numeric(x))
    }
    cat("Please enter a numeric value.\n")
  }
}

for (tree_file in DBH_recheck$File) {
  
  cat("\n=============================\n")
  cat("Manual handling for:", tree_file, "\n")
  
  tree_check <- tryCatch(
    read_tree_pc(path = file.path(path_trees, tree_file)),
    error = function(e) {
      message("❌ Error reading: ", tree_file, " | ", e$message)
      return(NULL)
    }
  )
  if (is.null(tree_check)) next
  
  dbh_est <- ask_dbh("Type estimated DBH and press Enter: ")
  
  DBHs_checked <- rbind(
    DBHs_checked,
    data.frame(
      File = tree_file,
      DBH = dbh_est,
      R2  = "est",
      slice_height = NA_real_,
      slice_thickness = NA_real_,
      stringsAsFactors = FALSE
    )
  )
}

DBHs_checked= DBHs_checked %>% 
  mutate(dead=if_else(DBH==0,T,F))

### Data gathering ####
#%%%%%%%%%%%%%%%%%%%%%%
DBH_all=bind_rows(DBH_tot %>% mutate(dead=FALSE,R2=as.character(R2)),
                  DBHs_checked,
                  )
ger10_inventory = Height_tot %>% 
  left_join(DBH_all) %>%
  filter(!is.nan(DBH)) %>% 
  mutate(adult=DBH>0.12,
         regen_est_12=DBH<0.12&H>1.30,
         regen_unest_12=DBH<0.12&H<1.30,
         regen_est_10=DBH<0.1&H>1.30,
         regen_unest_10=DBH<0.1&H<1.30,
         regen_est_7.5=DBH<0.075&H>1.30,
         regen_unest_7.5=DBH<0.075&H<1.30) 

# Save results
write.csv(DBHs_checked, file.path(path_trees, "DBHs_checked.csv"), row.names = FALSE)
writeLines(files_hand, file.path(path_trees, "files_to_handle_manually.txt"))
write.csv(ger10_inventory, file.path(path_trees, "GER10_inventory.csv"), row.names = FALSE)







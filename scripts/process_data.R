
# setup -------------------------------------------------------------------

library(tidyverse)

# Path to the id of the file (without the top-level google sheets path):

gs_aaron <- "1AnpepLklUc9GzxCWhiJ7USdv2BIrP-8Bb9hw-dmyqio"
gs_jared <- "1YcuPuBNz_1wSAQ4WWCyvKXqGoPodmnpnLUXIPaVq_IA"

# Read in all of the sheets and assign each:

list(
  "sticky_traps",
  "pitfall_traps",
  "treatments"
) %>% 
  set_names() %>% 
  map(
    \(.sheet) {
      lst(gs_aaron, gs_jared) %>% 
        set_names(
          names(.) %>% 
            str_replace("gs_", "")
        ) %>%
        imap(
          \(.student, .name) {
            out_frame <- 
              file.path(
                "https://docs.google.com/spreadsheets/d", 
                .student
              ) %>% 
              googlesheets4::read_sheet(sheet = .sheet)
            if(nrow(out_frame) > 1) {
              out_frame %>% 
                mutate(
                  observer = .name
                )
            } else {
              NULL
            }
          }
        ) %>% 
        list_rbind() %>% 
        distinct()
    }
  ) %>% 
  list2env(envir = .GlobalEnv)

# Remove the assignments we will no longer need:

rm(gs_aaron, gs_jared)

# process treatment data --------------------------------------------------

treatment_proc <-
  treatments %>% 
  mutate(
    date = as_date(deploy_date)
  ) %>% 
  select(
    transect_id,
    date,
    treatment
  ) %>% 
  distinct()

# Remove the assignment we will no longer need:

rm(treatments)

# process trap data -------------------------------------------------------

trap_data <-
  
  # Named list of trap data:
  
  lst(sticky_traps, pitfall_traps) %>% 
  
  # imap allows us to iterate across the data and name of the list:
  
  imap(
    \(.data, .name) {
      .data %>% 
        
        # Rename the trap_id column so that it's the same across sticky and
        # pitfall traps:
        
        rename(
          transect_id = matches("_id")
        ) %>% 
        
        # Subset to records at distance == 0:
        
        filter(
          str_detect(transect_id, "light")
        ) %>% 
        
        # Add and edit columns:
        
        mutate(
          transect_id = 
            transect_id %>% 
            
            # Remove the portion of the sticky trap ID that identifies the
            # location and side of the trap:
            
            str_remove("-l.*$") %>% 
            
            # There were some issues with habitat designation:
            
            str_replace("lands", "land") %>% 
            
            # Spell out grassland and forest:
            
            str_replace("-g-", "-grassland-") %>% 
            str_replace("-f-", "-forest-"),
          
          # Add a trap type column based on the name of the dataset (and remove
          # the "_traps" suffix):
          
          trap_type = str_remove(.name, "_traps"),
          
          # Change date-time column to a date:
          
          date = as_date(deploy_date),
        ) %>% 
        
        # Count the number of arthropods by event and trap_type:
        
        summarize(
          
          # Only count arthropods that are greater or equal to 2 mm:
          
          count = sum(count[length >= 2], na.rm = TRUE),
          .by = 
            c(
              transect_id, 
              date, 
              trap_type,
              observer
            )
        )
    }
  ) %>% 
  
  # Combine into a single data frame:
  
  list_rbind() %>% 
  
  # Join with treatment data:
  
  left_join(
    treatment_proc,
    by = join_by(transect_id, date)
  ) %>% 
  
  # Remove NA counts:
  
  drop_na(count) %>%
  
  # Re-arrange the columns:
  
  select(
    transect_id:date,
    treatment,
    trap_type:count
  )

# Remove unnecessary assignments:

rm(pitfall_traps, sticky_traps)

# write to file -----------------------------------------------------------

trap_data %>% 
  write_rds("data/trap_data.rds")

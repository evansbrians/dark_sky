
# setup -------------------------------------------------------------------

library(tidyverse)

# Path to spreadsheet (second part is the id of the file):

gs_url <-
  file.path(
    "https://docs.google.com/spreadsheets/d",
    "1AnpepLklUc9GzxCWhiJ7USdv2BIrP-8Bb9hw-dmyqio"
  )

# Read in all of the sheets and assign each:

googlesheets4::sheet_names(gs_url) %>% 
  set_names() %>% 
  map(
    ~ googlesheets4::read_sheet(gs_url, sheet = .x)
  ) %>% 
  list2env(.GlobalEnv)

# Remove the assignment we will no longer need:

rm(gs_url)

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
  )

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
          
          n = sum(count[length >= 2]),
          .by = 
            c(
              transect_id, 
              date, 
              trap_type
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
  
  # Re-arrange the columns:
  
  select(
    transect_id:date,
    treatment,
    trap_type:n
  )

# Remove unnecessary assignments:

rm(pitfall_traps, sticky_traps)

# write to file -----------------------------------------------------------

trap_data %>% 
  write_rds("data/trap_data_aaron.rds")

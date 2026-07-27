
#DATASET 3  -->  NO DAP available

# Author:        Anna Peiró
# Date:          2026-07-22 
# ------------------------------
# Description:  Input of doses (mGy) from a database designated as 'dataset 3'.
# ------------------------------
# Modification:       Missing age in data considered (but not recommended; no source PDF available for input in this script).
# Modification date:  2026-07-23
# --------------
#
#///  INPUT
#   - Data frame (df_input) containing information on age, sex and ExamCode. Details below.
#
#///  OUTPUT
#   - Optional: 200 files containing the 200 realizations.
#   - df_ds3_doses_final.rds : includes the procedures, along with the mean and median of the 200 realizations.
# --------------
# Make sure that main_table_AM and df_input have consistent column types in the variables used for matching and calculations (factor, character, numeric, etc.)
# --------------

library(dplyr)
library(tidyr)
library(purrr)

setwd("Your working directory")

# Load df_input
#Columns must include: PatientID, doe or CIncN (as procedure identification), age, sex, exam code.
 #age categories for the exam code simulation if needed (name in df_input and script: "Age_class")
    #cat 1 -> < 4 months
    #cat 2 -> 4 months to 90 months
    #cat 3 -> >= 91 months
 #age categories for the simulation of the doses (same categories as in the main table). (name in df_input and script: "age")
    #cat 1  -> < 4 months
    #cat 5  -> 4 months to 30 months
    #cat 10 -> 31 months to 90 months
    #cat 15 -> 91 months to 210 months
    #cat 30 -> >= 211 months
 #sex: 1 -> Male; 2 -> Female 
    #If sex is missing -> asssign sex 1 to under 10 years old and simulate the rest.
 #exam codes ("Final_ExamCode") in two digits except 111, 112, 131 and 184


# Load main table (main_table_AM: "sex", "age", "Final_ExamCode", "AM_dose") 
  #with AM doses from which we will sample a dose using parameters (sex, age, Final_ExamCode)

# Load PDFs in case of missings in sex and/or exam code:

PDF_sex <- data.frame(ExamCode = c(rep(c(111,112,12,13,131,14,15,16,17,18,184,21,22,23,24,31,32,33,41,43),2)),
                      Sex = c(rep("2",20), rep("1",20)), # 2 -> Female; 1 -> Male
                      Freq = c(0.5872,0.5,0.5864,0.3158,0.5061,0.36,0.3778,0.3014,0.2577,0.3372,0.4242,0.4204,0.4444,0.3077,
                               0.5228,0.4237,0.5736,0.4115,1, 0.5462,0.4128,0.5,0.4136,0.6842,0.4939,0.64,0.6222,0.6986,
                               0.7423,0.6628,0.5758,0.5796,0.5556,0.6923,0.4772,0.5763,0.4264,0.5885,0,0.4538))
PDF_exam <- data.frame(Age_class = c(rep(c(1,2,3),20)),
                       ExamCode = c(rep(111,3),rep(112,3),rep(12,3),rep(13,3),rep(131,3),rep(14,3),rep(15,3),rep(16,3),rep(17,3),
                                    rep(18,3),rep(184,3),rep(21,3),rep(22,3),rep(23,3),rep(24,3),rep(31,3),rep(32,3),rep(33,3),
                                    rep(41,3),rep(43,3)),
                       Freq = c(0.0294,0.1499,0.124,0,0,0.0013,0.1917,0.2235,0.0395,0.0406,0.0073,0.007,0.1373,0.0202,0.0128,
                                0.0363,0.0784,0.0778,0.0406,0.0187,0.0332,0.0009,0.001,0.0376,0.1313,0.004,0.0016,0.0803,0.0642,
                                0.0577,0.0009,0.0119,0.0092,0.0699,0.0657,0.0759,0.0665,0.118,0.0752,0.0026,0.0048,0.0064,0.0613,
                                0.0665,0.0561,0,0.0068,0.0319,0.0009,0.0126,0.0985,0.006,0.0344,0.1868,0,0,0.0003,0.1036,0.112,0.0673))

# Build look-up tables for simulation 

lookup_table_sex <- PDF_sex %>%
  group_by(ExamCode) %>%
  arrange(desc(Freq)) %>% 
  mutate(
    prob = Freq / sum(Freq),
    cum_prob = cumsum(prob) 
  ) %>%
  ungroup()

lookup_table_ec <- PDF_exam %>%
  group_by(Age_class) %>%
  arrange(desc(Freq)) %>% 
  mutate(
    prob = Freq / sum(Freq),
    cum_prob = cumsum(prob) 
  ) %>%
  ungroup()

# Function for AM dose sampling - Expected to be repeated the 200 realizations

AM_dose_ds3 <- function(df, df_source) {
  df_source_clean <- df_source %>%
    filter(!is.na(AM_dose)) %>%
    group_by(sex, age, Final_ExamCode) %>%
    summarise(candidatos = list(AM_dose), .groups = "drop")
  
  df %>%
    left_join(df_source_clean, by = c("sex", "age", "Final_ExamCode")) %>%
    rowwise() %>% 
    mutate(
      AM_dose = if (is.null(candidatos) || length(candidatos) == 0) {
        NA_real_
      } else {
        sample(candidatos, 1)
      }
    ) %>%
    ungroup() %>%
    select(-candidatos) 
}


df0 <- df_input

# we save every realization of the doses in a separated matrix to calculate summaries
n_iter  <- 200
n_rows  <- nrow(df0)      
AM_doses <- matrix(NA_real_, nrow = n_rows, ncol = n_iter)

for(i in 1:200){
  
  print(paste0(i, " out of ", 200))
  
  df <- df0
  
  #exam code
  
  df_missEC <- df[which(is.na(df$Final_ExamCode)),]
  df_missEC <- df_missEC[!duplicated(df_missEC[c("PatientID","doe")]),] # change "doe" to "CIncN" if needed
  
  df_missEC <- df_missEC %>%
    group_by(Age_class) %>%
    mutate(
      rand_val = runif(n()), 
      EC_sim = {
        current_age <- first(Age_class)
        ref_data <- lookup_table_ec %>% 
          filter(Age_class == current_age)
        
        if (nrow(ref_data) == 0 || is.na(current_age)) {
          rep(NA_real_, n())
        } else {
          cortes <- c(0, ref_data$cum_prob)
          indices <- findInterval(rand_val, cortes)
          indices <- pmax(1, indices) 
          ref_data$ExamCode[indices]
        }
      }
    ) %>%
    ungroup()
  
  df <-  left_join(df, df_missEC[,c("PatientID","doe","EC_sim")], # change "doe" to "CIncN" if needed
                    by = c("PatientID","doe"))
  
  df <- df %>%
    mutate(
      Final_ExamCode = coalesce(Final_ExamCode, EC_sim) # Note: column types may need verification
    ) %>%
    select(-EC_sim)
  
  #sex
  #if missing -> asssign sex 1 to under 10 years old
  
  df_missSex <- df[which(is.na(df$sex)),]
  df_missSex <- df_missSex[,c("PatientID","Final_ExamCode")]
  df_missSex <- df_missSex[which(!duplicated(df_missSex$PatientID)),]  # independently from the procedures, the Sex will be the same for each patient
  
  df_missSex <- df_missSex %>%
    group_by(Final_ExamCode) %>%
    mutate(
      rand_val = runif(n()), 
      
      Sex_sim = {
        current_code <- first(Final_ExamCode)
        ref_data <- lookup_table_sex %>% 
          filter(ExamCode == current_code)
        
        if(nrow(ref_data) == 0 || is.na(current_code)) {
          rep(NA_character_, n())
        } else {
          cortes <- c(0, ref_data$cum_prob)
          indices <- findInterval(rand_val, cortes)
          indices <- pmax(1, indices)
          ref_data$Sex[indices]
        }
      }
    ) %>%
    ungroup()
  
  df <-  left_join(df, df_missSex[,c("PatientID","Sex_sim")], by = "PatientID")
  
  df <- df %>%
    mutate(
      sex = coalesce(sex, Sex_sim) # Note: column types may need verification
    ) %>%
    select(-Sex_sim)
  
  #Dose calculation  
  
  df_dose <- AM_dose_ds3(df, main_table_AM)
  
  ## Collect the doses into a matrix
  AM_doses[, i] <- df_dose[["AM_dose"]]
  
  #Optional at this point: Save the 200 files containing the 200 realizations.
  
}

AM_doses_df <- as.data.frame(AM_doses)
names(AM_doses_df) <- paste0("AM_dose_sim_", seq_len(n_iter))

AM_mean <- data.frame(`AM_dose(mGy)_mean` = apply(AM_doses_df, 1, mean, na.rm = T))
AM_median <- data.frame(`AM_dose(mGy)_median` = apply(AM_doses_df, 1, median, na.rm = T))


df_ds3_doses_final <- cbind(df0[,c("PatientID", "doe")],  # change "doe" to "CIncN" if needed
                            AM_median, AM_mean)

saveRDS(df_ds3_doses_final,"df_ds3_doses_final.rds")




# Load necessary libraries
library("data.table")
library("tidyverse")
library("lubridate")
library("mltools")
library("ggplot2")
library("ggpubr")
library("ggforce")
library("cowplot")
library("DescTools")
library("parallelMap")
library("iml")
library("dplyr")
library("MLmetrics")
library("foreach")
library("doParallel")
library("writexl")

# Set system language to English
Sys.setenv(LANG = "en")

# Set random number generator kind
RNGkind("L'Ecuyer-CMRG")
seed_nr <- 123

# Set working directory
setwd("data_dir")

# Read data
bakteriaemi_data <- fread(file = "data.csv")

# Data preprocessing
colnames(bakteriaemi_data)[colnames(bakteriaemi_data) == "alder"] <- "Age"
colnames(bakteriaemi_data)[colnames(bakteriaemi_data) == "kon"] <- "Sex"

# Handle missing age values
missing_age <- which(bakteriaemi_data$Age <= 0 | bakteriaemi_data$Age > 120)
bakteriaemi_data$Age[missing_age] <- NA

# Convert date format
bakteriaemi_data$TestDate <- lubridate::dmy(bakteriaemi_data$afsendt)
bakteriaemi_data$afsendt <- NULL

# Read contaminants
contaminants <- read_excel("contaminants.xlsx", col_names = FALSE)
cont_idx <- which(bakteriaemi_data$Bakterie %in% contaminants$...1[which(contaminants$...2 == 0)])
bakteriaemi_data$Bakterie[cont_idx] <- ""

# Read additional age data
new_ages <- fread("age_data.csv")
colnames(new_ages)[4] <- "ID"
colnames(new_ages)[2] <- "TestDate"
colnames(new_ages)[3] <- "Age"
new_ages$Age <- as.numeric(gsub(",", ".", new_ages$Age))
new_ages <- new_ages[,-c(1)]  # Remove unnecessary columns

# Merge age data with bakteriaemi_data
bakteriaemi_data <- merge(bakteriaemi_data, new_ages, by = c("ID", "TestDate"))
bakteriaemi_data$Age <- as.numeric(gsub(",", ".", bakteriaemi_data$Age))
bakteriaemi_data <- bakteriaemi_data[!(bakteriaemi_data$Age < 18),]

# Read bacteria information
bacteria_info <- read_xlsx("pathogen_info.xlsx")

# Remove redundant columns
removable_names <- c("Blandingsflora", "Tarmpatogene bakterier ikke påvist", "Ingen vækst af bakterier og svampe")
bakteriaemi_data <- bakteriaemi_data[,(setdiff(colnames(bakteriaemi_data), removable_names))]

# Merge columns with the same bacteria info
bakteriaemi_data1 <- bakteriaemi_data
cols_to_merge_list <- c()
for (i in 5:ncol(bakteriaemi_data1)) {
  binfo_idx <- which(colnames(bakteriaemi_data1)[i] == bacteria_info$pathogen...2)
  if (sum(bacteria_info$...4 == bacteria_info$...4[binfo_idx]) > 1 &
      !(any(cols_to_merge_list == colnames(bakteriaemi_data1)[i]))) {
    cols_to_merge <- bacteria_info$pathogen...2[which(bacteria_info$...4 == bacteria_info$...4[binfo_idx])]
    bakteriaemi_data1 <- cbind(bakteriaemi_data1, FALSE)
    colnames(bakteriaemi_data1)[ncol(bakteriaemi_data1)] <- bacteria_info$...4[binfo_idx]
    bakteriaemi_data1[,(bacteria_info$...4[binfo_idx]) :=
                        apply(bakteriaemi_data1[, ..cols_to_merge], 1, function(x) { any(x) })]
    cols_to_merge_list <- c(cols_to_merge, cols_to_merge_list)
  }
}
# Remove duplicates in bakteriaemi_data_wide1
bakteriaemi_data_wide1[,(cols_to_merge_list):=NULL]
bacteria_names <- colnames(bakteriaemi_data_wide1)[5:ncol(bakteriaemi_data_wide1)]
bacteria_names <- bacteria_names[-which(bacteria_names %in% c("Age"))]
bakteriaemi_data_wide1$BSIClass_new <- apply(bakteriaemi_data_wide1[,..bacteria_names],1,any)
bakteriaemi_data_wide1$BSIClass_new <- factor(bakteriaemi_data_wide1$BSIClass_new,
                                              levels = c(FALSE,TRUE),
                                              labels = c("noBSI","BSI"))
bakteriaemi_data_wide1$BSIClass <- NULL
colnames(bakteriaemi_data_wide1)[which(colnames(bakteriaemi_data_wide1)=="BSIClass_new")] <- "BSIClass"
bakteriaemi_data_wide <- bakteriaemi_data_wide1
bakteriaemi_data_wide1 <- NULL

# Read biochemical data
biochemi_data <- fread(file = "biochemical_data.csv")

# Parse date time in biochemical data
datetime1 <- parse_date_time(biochemi_data$SVARTID, orders = c("dmy_HMS", "mdy_HMS"))
datetime2 <- parse_date_time(biochemi_data$PRVTAG_TID, orders = c("dmy_HMS", "mdy_HMS"))

# Handle problem entries
problem_entries <- biochemi_data$PRVTAG_TID[is.na(datetime2)]
print(problem_entries)

datetime2[is.na(datetime2)] <- dmy(biochemi_data$PRVTAG_TID[is.na(datetime2)])
problem_entries <- biochemi_data$SVARTID[is.na(datetime1)]
print(problem_entries)
datetime1[is.na(datetime1)] <- dmy(biochemi_data$SVARTID[is.na(datetime1)])

# Calculate time difference
time_difference <- datetime2 - datetime1
days_difference <- as.numeric(time_difference, units = "days")
days_difference_abs_largerthanone <- abs(days_difference[abs(days_difference)>1])

# Create a histogram with a log-scaled x-axis
hist(
  days_difference_abs_largerthanone,
  main = "Days Between Two Dates",
  xlab = "Days Difference",
  col = "skyblue",
  border = "black",
  breaks = 100000,
  log = "x"  # Set log argument to "x" for a log-scaled x-axis
)

# Filter out entries with a time difference larger than seven days
larger_than_seven_days <- which(abs(days_difference)>7)
biochemi_data <- biochemi_data[-larger_than_seven_days,]

# Remove unnecessary columns
biochemi_data$ANALYSEKODE_LABKA <- NULL
biochemi_data$REQUESTERCODE <- NULL
biochemi_data$REPLY <- gsub(",", ".", biochemi_data$REPLY)
biochemi_data$REPLY <- as.numeric(biochemi_data$REPLY)
biochemi_data$INVESTIGATION_NAME <- as.factor(biochemi_data$INVESTIGATION_NAME)
biochemi_data$SVARTID <- NULL

# Check for units consistency
varswithmultipleunits <- c()
k <- 1
multiunits <- list()
invnames <- as.character(unique(biochemi_data$INVESTIGATION_NAME))
for (i in 1:length(invnames)){
  idx <- which(biochemi_data$INVESTIGATION_NAME==invnames[i])
  if (length(unique(biochemi_data$REPLYUNIT[idx]))>1){
    print(invnames[i])
    varswithmultipleunits <- c(varswithmultipleunits,invnames[i])
    multiunits[[k]] <- unique(biochemi_data$REPLYUNIT[idx])
    k <- k+1
  }
}

# Remove units from REPLYUNIT column
biochemi_data$REPLYUNIT <- NULL

# Display table for INVESTIGATION_NAME
table(biochemi_data[,c("INVESTIGATION_NAME")])

# Remove the INVESTIGATION_NAME column
biochemi_data$INVESTIGATION_NAME <- NULL

# Parse the PRVTAG_TID column into a new column named parsed_date
library(lubridate)
biochemi_data$parsed_date <- parse_date_time(biochemi_data$PRVTAG_TID, orders = c("dmy_HMS", "mdy_HMS"))
# If parsing failed, try to parse using dmy format for the date part only
biochemi_data$parsed_date[is.na(biochemi_data$parsed_date)] <- 
  parse_date_time(biochemi_data$PRVTAG_TID[is.na(biochemi_data$parsed_date)], orders = "dmy")

# Extract date from parsed_date and store it in a new column named date_only
biochemi_data$date_only <- date(biochemi_data$parsed_date)

# Assuming datetime1 is a column containing date-time values (a gap of 7 days is considered)
biochemi_data$date_only <- biochemi_data$date_only + ddays(7)

# Remove the parsed_date column if it's no longer needed
biochemi_data$parsed_date <- NULL
dim(biochemi_data)
biochemi_data <- biochemi_data[-which(biochemi_data$date_only>ending_date | biochemi_data$date_only<starting_date),]
dim(biochemi_data)
biochemi_data$PRVTAG_TID <- NULL
summary(biochemi_data$date_only)
biochemi_data_wide <- dcast(biochemi_data,ID + date_only ~ ANALYSEKODE_MAIDS, value.var = "REPLY",fun.aggregate = function(x) mean(x,na.rm=T))
dim(biochemi_data_wide)
colnames(biochemi_data_wide)[which(colnames(biochemi_data_wide)=="date_only")] <- "TestDate"

# LPR data
LPR_data = fread(file = "LPR_data.csv")
LPR_data$STAMSGH <- NULL
LPR_data$STAMAFD <- as.character(LPR_data$STAMAFD) # hospital department/section the test sample was taken
LPR_data$SORID_SGH <- NULL
LPR_data$SORID_AFD <- as.character(LPR_data$SORID_AFD) # hospital department/section the test sample was taken
LPR_data$PATIENTTYPE <- as.character(LPR_data$PATIENTTYPE) # inpatient, outpatient etc. (could be relevant)
INDDATO <- lubridate::date(lubridate::dmy_hms(LPR_data$INDDATO_DATO_TID))
any(is.na(INDDATO))
all(LPR_data$INDDATO_DATO_TID[which(is.na(INDDATO))]=="")
rm(INDDATO)
LPR_data$INDDATO_DATO_TID <- lubridate::date(lubridate::dmy_hms(LPR_data$INDDATO_DATO_TID)) # extract dates of admission

LPR_data$UDDATO_DATO_TID <- lubridate::date(lubridate::dmy_hms(LPR_data$UDDATO_DATO_TID)) # extract dates of discharge
UDDATO <- LPR_data$UDDATO_DATO_TID
any(is.na(UDDATO))
all(LPR_data$UDDATO_DATO_TID[which(is.na(UDDATO))]=="")
rm(UDDATO)
LPR_data <- LPR_data[-which(is.na(LPR_data$INDDATO_DATO_TID)),]

summary(LPR_data$INDDATO_DATO_TID) 
LPR_data$INDDATO_DATO_TID[which(LPR_data$INDDATO_DATO_TID<"2008-01-01")] # some are not too old but still don't match with the period of the other datasets so excluded
if (filter_starting_date){
  LPR_data <- LPR_data[-which(LPR_data$INDDATO_DATO_TID<starting_date),] # eligible period: 2012-2020
}

summary(LPR_data$UDDATO_DATO_TID) 
LPR_data$UDDATO_DATO_TID[which(LPR_data$UDDATO_DATO_TID>ending_date)]# eligible period: 2012-2020
LPR_data <- LPR_data[-which(LPR_data$UDDATO_DATO_TID>ending_date),]

common_IDs <- intersect(bakteriaemi_data_wide$ID,biochemi_data_wide$ID) # common IDs
ttl_common_IDs <- intersect(LPR_data$ID,common_IDs) 
bakteriaemi_data_wide <- bakteriaemi_data_wide[which(bakteriaemi_data_wide$ID %in% ttl_common_IDs),]
biochemi_data_wide <- biochemi_data_wide[which(biochemi_data_wide$ID %in% ttl_common_IDs),]
LPR_data <- LPR_data[which(LPR_data$ID %in% ttl_common_IDs),]

# keep only records with common IDs in all three datasets
dim(bakteriaemi_data_wide)
dim(biochemi_data_wide)
dim(LPR_data)
LPR_data_df <- as.data.frame(table(LPR_data$ID))
records_dist_plt1 <- ggplot(bakteriaemi_data_wide, aes(x=TestDate)) + 
  geom_freqpoly(binwidth = 10)+
  ylab("counts") +
  xlab("date") +
  theme_article(base_size = 10)
ggsave(filename = "bakteriaemi_data_wide_dist_plt.pdf",plot = records_dist_plt1)
fwrite(bakteriaemi_data_wide,"bacteremia_wide.csv")
bakteriaemi_data_wide$most_common_pathogens <- bakteriaemi_data_wide$`Staphylococcus aureus` | bakteriaemi_data_wide$`E. coli` | bakteriaemi_data_wide$`Enterococcus faecium`
any(duplicated(bakteriaemi_data_wide, by=c("ID","TestDate")))
any(duplicated(biochemi_data_wide, by=c("ID","TestDate")))
dim(biochemi_data_wide)
biochemi_data_wide_uniq <- unique(biochemi_data_wide)
dim(biochemi_data_wide_uniq)
rm(biochemi_data_wide_uniq)

full_dt1 <- merge(x=biochemi_data_wide,y=bakteriaemi_data_wide, by = c("ID","TestDate"), all = TRUE)
table(full_dt1$BSIClass, useNA = "always")
full_dt1$PrevAdmissionRate <- NA

# there are duplicated rows in LPR data to remove
LPR_data_dup_rows <- which(duplicated(LPR_data))
length(LPR_data_dup_rows)
LPR_data <- LPR_data[-LPR_data_dup_rows,]
any(LPR_data$UDDATO_DATO_TID<starting_date)

# keep only records with common IDs in all three datasets
dim(bakteriaemi_data_wide)
dim(biochemi_data_wide)
dim(LPR_data)

LPR_data <- unique(LPR_data,by = c("ID", "INDDATO_DATO_TID", "UDDATO_DATO_TID"))
dim(LPR_data)

hist(as.numeric(LPR_data$UDDATO_DATO_TID-LPR_data$INDDATO_DATO_TID))

seemingly_invalid_entries <- which(as.numeric(LPR_data$UDDATO_DATO_TID-LPR_data$INDDATO_DATO_TID)>365) # patients with longer than a year hospitalization
length(seemingly_invalid_entries)

LPR_data <- LPR_data[-seemingly_invalid_entries,]
dim(LPR_data)
library(dplyr)

# Group the data by INDDATO_DATO_TID and ID, and find the row with the highest UDDATO_DATO_TID
result <- LPR_data %>%
  group_by(INDDATO_DATO_TID, ID) %>%
  filter(UDDATO_DATO_TID == max(UDDATO_DATO_TID)) %>%
  ungroup()

# Check for duplicated rows by INDDATO_DATO_TID and ID
duplicated_rows <- result[duplicated(result[, c("INDDATO_DATO_TID", "ID")]), ]

# Print the duplicated rows
if (nrow(duplicated_rows) > 0) {
  print(duplicated_rows)
} else {
  print("No duplicated rows")
}

result <- as.data.table(result)
rm(LPR_data) 
LPR_data <- result
rm(result)

dim(LPR_data)

full_dt1 <- full_dt1[-which(!full_dt1$ID %in% LPR_data$ID),]
any(!full_dt1$ID %in% LPR_data$ID)
dim(full_dt1)

library(doParallel)
library(foreach)

# Choose an appropriate number of cores (e.g., half of available cores)
num_cores <- floor(detectCores() / 2)
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Export necessary data to worker nodes
clusterExport(cl, c("LPR_data", "full_dt1", "lookback_period"))

# Initialize a list to collect results
result_list <- foreach(i = 1:nrow(full_dt1), .combine = "c") %dopar% {
  ID_idx <- which(LPR_data$ID == full_dt1$ID[i])
  
  tempMatrix <- sum(LPR_data$UDDATO_DATO_TID[ID_idx] < full_dt1$TestDate[i] & (LPR_data$UDDATO_DATO_TID[ID_idx] + lookback_period) > full_dt1$TestDate[i])
  
  # Do other things if you want
  
  return(tempMatrix)
}

# Stop cluster
stopCluster(cl)
sample_at_start <- which(full_dt1$TestDate==starting_date)
early_samples <- which(full_dt1$TestDate<(as.Date(starting_date)+lookback_period))

lookback_early <- full_dt1$TestDate[early_samples]-as.Date(starting_date)
lookback_early[lookback_early==0] <- 1
lookback_early <- as.numeric(lookback_early)
# previous admission here includes outpatients
full_dt1$PrevAdmissionRate <- as.numeric(result_list) # accumulated number of days within the last lookback_period before the test date
full_dt1$PrevAdmissionRate[-early_samples] <- full_dt1$PrevAdmissionRate[-early_samples]/lookback_period # we normalize it to the lookback_period so that is comparable if we choose different lookback_period
full_dt1$PrevAdmissionRate[early_samples] <- full_dt1$PrevAdmissionRate[early_samples]/lookback_early
if (length(sample_at_start)>0){
  full_dt1$PrevAdmissionRate[sample_at_start] <- NA
}

full_dt1$PrevInfectionRate <- NA
# Set the number of cores to use for parallel processing (half of available cores)
num_cores <- max(1, floor(detectCores() / 2))
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Initialize a vector to store calculated values
modifiedPIR <- numeric(nrow(full_dt1))

# Parallel loop using foreach
finalMatrix <- foreach(i = 1:nrow(full_dt1), .combine = "c") %dopar% {
  # Find the current ID and TestDate
  currentID <- full_dt1$ID[i]
  currentTestDate <- full_dt1$TestDate[i]

  # Find indices of rows with the same ID that are within the lookback_period
  ID_idx_previous <- which(full_dt1$ID == currentID & full_dt1$TestDate < currentTestDate & (full_dt1$TestDate + lookback_period) > currentTestDate)

  # Find indices of rows with "BSI" or "noBSI" in the BSIClass column within the ID_idx_previous subset
  ID_idx2 <- ID_idx_previous[full_dt1$BSIClass[ID_idx_previous] %in% c("BSI", "noBSI")]
  ID_idx3 <- ID_idx_previous[which(full_dt1$BSIClass[ID_idx_previous] == "BSI")]

  return(length(ID_idx3) / length(ID_idx2))
}

# Stop parallel processing
stopCluster(cl)

# Convert the calculated values to a numeric vector and store them as modifiedPIR
full_dt1$modifiedPIR <- as.numeric(finalMatrix)

# Detect the number of available CPU cores
cores = detectCores()

# Create a parallel cluster with a specified number of cores
cl <- makeCluster(cores[1] - 1) # Avoid overloading your computer

# Register the parallel cluster for parallel computations
registerDoParallel(cl)

# Initialize a vector to store calculated values
finalMatrix <- foreach(i = 1:nrow(full_dt1), .combine = "c") %dopar% {
  # Find indices of rows with matching ID in the full_dt1 data frame
  ID_idx1 <- which(full_dt1$ID == full_dt1$ID[i])

  # Find indices of rows with "BSI" in the BSIClass column within the ID_idx1 subset
  ID_idx2 <- which(full_dt1$BSIClass[ID_idx1] == "BSI")

  # Calculate the number of previous infections within a specific time range
  tempMatrix = sum(
    full_dt1$TestDate[ID_idx1[ID_idx2]] < full_dt1$TestDate[i] &
      (full_dt1$TestDate[ID_idx1[ID_idx2]] + lookback_period) > full_dt1$TestDate[i]
  )

  # Return the calculated value for the current iteration
  tempMatrix
}

# Stop the parallel cluster
stopCluster(cl)

# Convert the calculated values to a numeric vector and store them as PrevInfectionRate
full_dt1$PrevInfectionRate <- as.numeric(finalMatrix)

# Normalize PrevInfectionRate values based on conditions
full_dt1$PrevInfectionRate[-early_samples] <- full_dt1$PrevInfectionRate[-early_samples] / lookback_period
full_dt1$PrevInfectionRate[early_samples] <- full_dt1$PrevInfectionRate[early_samples] / lookback_early

# Set PrevInfectionRate values to NA for specific rows
if (length(sample_at_start) > 0) {
  full_dt1$PrevInfectionRate[sample_at_start] <- NA
}

# Get the dimensions of the modified full_dt1 data frame
dim(full_dt1)

# plot a random patient data over time
# Choose a random ID with multiple data points in a given period
# Get IDs with multiple rows
multi_id <- full_dt1 %>%
  group_by(ID) %>%
  filter(n() > 1) %>%
  distinct(ID)

# Randomly select one of the IDs with multiple rows
set.seed(seed_nr)
id <- sample(multi_id$ID, 1)

# Create a subset of the data for the selected ID
id_df <- full_dt1[full_dt1$ID == id, ]

# set.seed(seed_nr)
# id <- sample(unique(full_dt1$ID), 1)
# id_df <- full_dt1[full_dt1$ID == id, ]
id_df$Sex <- as.numeric(factor(id_df$Sex, labels = c(0,1), levels = c("Female","Male")))
id_df$BSIClass <- as.numeric(factor(id_df$BSIClass, labels = c(0,1), levels = c("noBSI","BSI")))

# Select numeric variables to normalize
numeric_vars <- colnames(id_df[,-c("ID", "TestDate")])

# Normalize the numeric variables
id_df_norm <- id_df %>%
  mutate(across(all_of(numeric_vars), scale))

# Pivot the data to long format
id_df_long <- pivot_longer(id_df_norm, cols = -c(ID, TestDate), names_to = "Variable")

# Plot the data
example_plt <- ggplot(id_df_long, aes(x = TestDate, y = value, color = Variable)) + 
  geom_line(size = 1) +
  scale_color_discrete(breaks = c(0, 0.2, 0.4, 0.8, 1))
labs(title = paste("ID:", id), x = "Date", y = "Value") +
  theme_minimal()

rm(example_plt)
LPR_data$STAMAFD <- as.factor(LPR_data$STAMAFD)
LPR_data$SORID_AFD <- NULL # not useful as discussed in the meeting 23-06-2022

LPR_data$PATIENTTYPE <- as.factor(LPR_data$PATIENTTYPE)
LPR_data$PATIENTTYPE <- factor(LPR_data$PATIENTTYPE, 
                               levels = c("0", "2", "A", "F", "I", "S"),
                               labels = c("inpatient","outpatient","others","others","others","others"))

dim(LPR_data)
# Create a new data table with unique rows from LPR_data
LPR_data_unique <- unique(LPR_data)
dim(LPR_data_unique)
rm(LPR_data_unique)
setkey(LPR_data,ID)

full_dt1[,(bacteria_names):=NULL]
dim(full_dt1)
biochemi_data <- NULL
biochemi_data_wide_uniq <- NULL
gc()
fwrite(LPR_data,paste0("LPR_data_2010_2020_cleaned.csv"))
rm(LPR_data)
gc()
library(foreach)
library(doParallel)
library(tidyr)
library(magrittr)

# Find duplicated rows
duplicated_rows <- full_dt1[duplicated(full_dt1), ]
any(is.na(full_dt1$Sex))

# Group by ID and impute missing Sex with available Sex values for the same ID
full_dt1 <- full_dt1 %>%
  group_by(ID) %>%
  mutate(Sex = ifelse(is.na(Sex), unique(na.omit(Sex)), Sex)) %>%
  ungroup()

full_dt2 <- full_dt1
library(dplyr)

# convert TestDate to date format
full_dt2$TestDate <- as.Date(full_dt2$TestDate)

# sort by ID and TestDate
full_dt2 <- full_dt2 %>% arrange(ID, TestDate)

# calculate AgeDiff for each ID
full_dt2 <- full_dt2 %>% 
  group_by(ID) %>% 
  mutate(AgeDiff = difftime(TestDate, first(TestDate), units = "days")) %>% 
  ungroup()

as.numeric(full_dt2$AgeDiff)

# calculate Age for missing values based on AgeDiff and Age at the available TestDate
full_dt2 <- full_dt2 %>% 
  group_by(ID) %>% 
  mutate(Age = ifelse(is.na(Age), Age[!is.na(Age)][1] + as.numeric(AgeDiff)/365, Age)) %>% 
  ungroup()

# remove AgeDiff column
full_dt2$AgeDiff <- NULL

dim(full_dt2)
gc()

# rows with all NA for biochemical data are removed
biochemical_cols <- which(!colnames(full_dt2) %in% c("ID","TestDate","Sex","Age","BSIClass","PrevAdmissions","PrevInfections"))

all_NA_biochemical_cols <- which(apply(full_dt2[,biochemical_cols],1,function(x) {all(is.na(x))}))
gc()

setDT(full_dt2)
setwd("mydir_results")

full_dt2$Sex <- as.numeric(full_dt2$Sex)-1
colnames(full_dt2)[which(colnames(full_dt2)=="NA")] <- "Sodium"

if (temporal_feature_engineering){
  full_dt2$year_day <- as.numeric(full_dt2$year_day)
}

naAge_idx <- which(is.na(full_dt2$Age))
if (length(naAge_idx)>1){
  full_dt2 <- full_dt2[-naAge_idx,]
}

below18 <- which(full_dt2$Age<18)
if (length(below18)>1){
  full_dt2 <- full_dt2[-below18,]
}

dim(full_dt2)

full_dt2$Sex <- factor(full_dt2$Sex, levels = c(0,1), labels = c("F","M"))

full_dt2$GLUF <- as.numeric(full_dt2$GLUF)
full_dt2$GLUF[which(full_dt2$GLUF=="NaN")] <- NA

full_dt2$ALAT_binary <- as.factor(case_when(full_dt2$ALAT>=10 & full_dt2$ALAT<=45 & full_dt2$Sex=="F" ~ "nrm",
                                            full_dt2$ALAT>=10 & full_dt2$ALAT<=70 & full_dt2$Sex=="M" ~ "nrm",
                                            is.na(full_dt2$ALAT) ~ "NA",
                                            TRUE                      ~ "abnrm"))

full_dt2$ALB_binary <- as.factor(case_when(full_dt2$ALB>=36 & full_dt2$ALB<=48 & full_dt2$Age<40 ~ "nrm",
                                           full_dt2$ALB>=36 & full_dt2$ALB<=45 & full_dt2$Age>=40 & full_dt2$Age<70 ~ "nrm",
                                           full_dt2$ALB>=34 & full_dt2$ALB<=45 & full_dt2$Age>=70 ~ "nrm",
                                           is.na(full_dt2$ALB) ~ "NA",
                                           TRUE                      ~ "abnrm"))

full_dt2$APTT_binary <- as.factor(case_when(full_dt2$APTT>=25 & full_dt2$APTT<=37 ~ "nrm",
                                            is.na(full_dt2$APTT) ~ "NA",
                                            TRUE                      ~ "abnrm"))

full_dt2$BASO_binary <- as.factor(case_when(full_dt2$BASO<=0.1 ~ "nrm",
                                            is.na(full_dt2$BASO) ~ "NA",
                                            TRUE                      ~ "abnrm"))

full_dt2$BASOPOC_binary <- as.factor(case_when(full_dt2$BASOPOC<=0.1 ~ "nrm",
                                               is.na(full_dt2$BASOPOC) ~ "NA",
                                               TRUE                      ~ "abnrm"))

full_dt2$CA_binary <- as.factor(case_when(full_dt2$CA>2.15 & full_dt2$CA<=2.51 ~ "nrm",
                                          is.na(full_dt2$CA) ~ "NA",
                                          TRUE                      ~ "abnrm"))

full_dt2$CAI_binary <- as.factor(case_when(full_dt2$CAI>1.18 & full_dt2$CAI<=1.32 ~ "nrm",
                                           is.na(full_dt2$CAI) ~ "NA",
                                           TRUE                      ~ "abnrm"))
full_dt2$CAIPOC_binary <- as.factor(case_when(full_dt2$CAIPOC>1.18 & full_dt2$CAIPOC<=1.32 ~ "nrm",
                                              is.na(full_dt2$CAIPOC) ~ "NA",
                                              TRUE                      ~ "abnrm"))

full_dt2$CHOL_binary <- as.factor(case_when(full_dt2$CHOL<=5 ~ "nrm",
                                            is.na(full_dt2$CHOL) ~ "NA",
                                            TRUE                      ~ "abnrm"))

full_dt2$CL_binary <- as.factor(case_when(full_dt2$CL>=98 & full_dt2$CL<=106 ~ "nrm",
                                          is.na(full_dt2$CL) ~ "NA",
                                          TRUE                      ~ "abnrm"))

full_dt2$CLPOC_binary <- as.factor(case_when(full_dt2$CLPOC>=98 & full_dt2$CLPOC<=106 ~ "nrm",
                                             is.na(full_dt2$CLPOC) ~ "NA",
                                             TRUE                      ~ "abnrm"))

full_dt2$CREA_binary <- as.factor(case_when(full_dt2$CREA>=50 & full_dt2$CREA<=90 & full_dt2$Sex=="F" ~ "nrm",
                                            full_dt2$CREA>=60 & full_dt2$CREA<=105 & full_dt2$Sex=="M" ~ "nrm",
                                            is.na(full_dt2$CREA) ~ "NA",
                                            TRUE                      ~ "abnrm"))

full_dt2$CRP_binary <- as.factor(case_when(full_dt2$CRP<=10 ~ "nrm",
                                           is.na(full_dt2$CRP) ~ "NA",
                                           TRUE                      ~ "abnrm"))

full_dt2$EOS_binary <- as.factor(case_when(full_dt2$EOS<=0.05 ~ "nrm",
                                           is.na(full_dt2$EOS) ~ "NA",
                                           TRUE                      ~ "abnrm"))

full_dt2$ERY_binary <- as.factor(case_when(full_dt2$ERY>=3.94 & full_dt2$ERY<=5.16 & full_dt2$Sex=="F" ~ "nrm",
                                           full_dt2$ERY>=4.25 & full_dt2$ERY<=5.71 & full_dt2$Sex=="M" ~ "nrm",
                                           is.na(full_dt2$ERY) ~ "NA",
                                           TRUE                      ~ "abnrm"))

full_dt2$FERRITIN_binary <- as.factor(case_when(full_dt2$FERRITIN>=12 & full_dt2$FERRITIN<=300 ~ "nrm",
                                                is.na(full_dt2$FERRITIN) ~ "NA",
                                                TRUE                      ~ "abnrm"))

full_dt2$GGT_binary <- as.factor(case_when(full_dt2$GGT>=10 & full_dt2$GGT<=45 & full_dt2$Sex=="F" & full_dt2$Age<=40 ~ "nrm",
                                           full_dt2$GGT>=10 & full_dt2$GGT<=75 & full_dt2$Sex=="F" & full_dt2$Age>40 ~ "nrm",
                                           full_dt2$GGT>=10 & full_dt2$GGT<=80 & full_dt2$Sex=="M" & full_dt2$Age<=40 ~ "nrm",
                                           full_dt2$GGT>=15 & full_dt2$GGT<=115 & full_dt2$Sex=="M" & full_dt2$Age>40 ~ "nrm",
                                           is.na(full_dt2$GGT) ~ "NA",
                                           TRUE                      ~ "abnrm"))

full_dt2$GLU_binary <- as.factor(case_when(full_dt2$GLU>=4.2 & full_dt2$GLU<=6.3 ~ "nrm",
                                           is.na(full_dt2$GLU) ~ "NA",
                                           TRUE                      ~ "abnrm"))

full_dt2$GLUF_binary <- as.factor(case_when(full_dt2$GLUF>=4.2 & full_dt2$GLUF<=6.1 ~ "nrm",
                                            is.na(full_dt2$GLUF) ~ "NA",
                                            TRUE                      ~ "abnrm"))

full_dt2$GLUPOC_binary <- as.factor(case_when(full_dt2$GLUPOC>=4.2 & full_dt2$GLUPOC<=6.3 ~ "nrm",
                                              is.na(full_dt2$GLUPOC) ~ "NA",
                                              TRUE                      ~ "abnrm"))

full_dt2$HAPTO_binary <- as.factor(case_when(full_dt2$HAPTO>=0.35 & full_dt2$HAPTO<=1.85 & full_dt2$Age<=50 ~ "nrm",
                                             full_dt2$HAPTO>=0.47 & full_dt2$HAPTO<=2.05 & full_dt2$Age>50 ~ "nrm",
                                             is.na(full_dt2$HAPTO) ~ "NA",
                                             TRUE                      ~ "abnrm"))

full_dt2$HB_binary <- as.factor(case_when(full_dt2$HB>=7.3 & full_dt2$HB<=9.5 & full_dt2$Sex=="F" ~ "nrm",
                                          full_dt2$HB>=8.3 & full_dt2$HB<=10.5 & full_dt2$Sex=="M" ~ "nrm",
                                          is.na(full_dt2$HB) ~ "NA",
                                          TRUE                      ~ "abnrm"))

full_dt2$HBA1C_binary <- as.factor(case_when(full_dt2$HBA1C<=48 ~ "nrm",
                                             is.na(full_dt2$HBA1C) ~ "NA",
                                             TRUE                      ~ "abnrm"))

full_dt2$HBPOC_binary <- as.factor(case_when(full_dt2$HBPOC>=7.3 & full_dt2$HBPOC<=9.5 & full_dt2$Sex=="F" ~ "nrm",
                                             full_dt2$HBPOC>=8.3 & full_dt2$HBPOC<=10.5 & full_dt2$Sex=="M" ~ "nrm",
                                             is.na(full_dt2$HBPOC) ~ "NA",
                                             TRUE                      ~ "abnrm"))

full_dt2$HDL_binary <- as.factor(case_when(full_dt2$HDL>=1 ~ "nrm",
                                           is.na(full_dt2$HDL) ~ "NA",
                                           TRUE                      ~ "abnrm"))

full_dt2$INR_binary <- as.factor(case_when(full_dt2$INR<=1.2 ~ "nrm",
                                           is.na(full_dt2$INR) ~ "NA",
                                           TRUE                      ~ "abnrm"))

full_dt2$JERN_binary <- as.factor(case_when(full_dt2$JERN>=9 & full_dt2$JERN<=34 ~ "nrm",
                                            is.na(full_dt2$JERN) ~ "NA",
                                            TRUE                      ~ "abnrm"))

full_dt2$K_binary <- as.factor(case_when(full_dt2$K>=3.5 & full_dt2$K<=44 ~ "nrm",
                                         is.na(full_dt2$K) ~ "NA",
                                         TRUE                      ~ "abnrm"))

full_dt2$KPOC_binary <- as.factor(case_when(full_dt2$KPOC>=3.5 & full_dt2$KPOC<=44 ~ "nrm",
                                            is.na(full_dt2$KPOC) ~ "NA",
                                            TRUE                      ~ "abnrm"))

full_dt2$LDH_binary <- as.factor(case_when(full_dt2$LDH>=105 & full_dt2$LDH<=205 & full_dt2$Age<=70 ~ "nrm",
                                           full_dt2$LDH>=115 & full_dt2$LDH<=255 & full_dt2$Age<=70 ~ "nrm",
                                           is.na(full_dt2$LDH) ~ "NA",
                                           TRUE                      ~ "abnrm"))

full_dt2$LDL_binary <- as.factor(case_when(full_dt2$LDL<=3 ~ "nrm",
                                           is.na(full_dt2$LDL) ~ "NA",
                                           TRUE                      ~ "abnrm"))

full_dt2$LEU_binary <- as.factor(case_when(full_dt2$LEU>=3.5 & full_dt2$LEU<=8.8 ~ "nrm",
                                           is.na(full_dt2$LEU) ~ "NA",
                                           TRUE                      ~ "abnrm"))

full_dt2$LEUPOC_binary <- as.factor(case_when(full_dt2$LEUPOC>=3.5 & full_dt2$LEUPOC<=8.8 ~ "nrm",
                                              is.na(full_dt2$LEUPOC) ~ "NA",
                                              TRUE                      ~ "abnrm"))

full_dt2$LYMFO_binary <- as.factor(case_when(full_dt2$LYMFO>=1 & full_dt2$LYMFO<=3.5 ~ "nrm",
                                             is.na(full_dt2$LYMFO) ~ "NA",
                                             TRUE                      ~ "abnrm"))

full_dt2$LYMFOPOC_binary <- as.factor(case_when(full_dt2$LYMFOPOC>=1 & full_dt2$LYMFOPOC<=3.5 ~ "nrm",
                                                is.na(full_dt2$LYMFOPOC) ~ "NA",
                                                TRUE                      ~ "abnrm"))

full_dt2$MG_binary <- as.factor(case_when(full_dt2$MG>=0.71 & full_dt2$MG<=0.94 ~ "nrm",
                                          is.na(full_dt2$MG) ~ "NA",
                                          TRUE                      ~ "abnrm"))

full_dt2$MONO_binary <- as.factor(case_when(full_dt2$MONO>=0.2 & full_dt2$MONO<=0.8 ~ "nrm",
                                            is.na(full_dt2$MONO) ~ "NA",
                                            TRUE                      ~ "abnrm"))

full_dt2$MONOPOC_binary <- as.factor(case_when(full_dt2$MONOPOC>=0.2 & full_dt2$MONOPOC<=0.8 ~ "nrm",
                                               is.na(full_dt2$MONOPOC) ~ "NA",
                                               TRUE                      ~ "abnrm"))

full_dt2$Sodium_binary <- as.factor(case_when(full_dt2$Sodium>=137 & full_dt2$Sodium<=144 ~ "nrm",
                                              is.na(full_dt2$Sodium) ~ "NA",
                                              TRUE                      ~ "abnrm"))

full_dt2$NAPOC_binary <- as.factor(case_when(full_dt2$NAPOC>=137 & full_dt2$NAPOC<=144 ~ "nrm",
                                             is.na(full_dt2$NAPOC) ~ "NA",
                                             TRUE                      ~ "abnrm"))

full_dt2$NEUTRO_binary <- as.factor(case_when(full_dt2$NEUTRO>=1.6 & full_dt2$NEUTRO<=5.9 ~ "nrm",
                                              is.na(full_dt2$NEUTRO) ~ "NA",
                                              TRUE                      ~ "abnrm"))

full_dt2$NEUTROPOC_binary <- as.factor(case_when(full_dt2$NEUTROPOC>=1.6 & full_dt2$NEUTROPOC<=5.9 ~ "nrm",
                                                 is.na(full_dt2$NEUTROPOC) ~ "NA",
                                                 TRUE                      ~ "abnrm"))

full_dt2$PROCAL_binary <- as.factor(case_when(full_dt2$PROCAL<=0.5 ~ "nrm",
                                              is.na(full_dt2$PROCAL) ~ "NA",
                                              TRUE                      ~ "abnrm"))

full_dt2$THROM_binary <- as.factor(case_when(full_dt2$THROM>=145 & full_dt2$THROM<=390 ~ "nrm",
                                             is.na(full_dt2$THROM) ~ "NA",
                                             TRUE                      ~ "abnrm"))

full_dt2$TRANS_binary <- as.factor(case_when(full_dt2$TRANS>=1.91 & full_dt2$TRANS<=3.26 ~ "nrm",
                                             is.na(full_dt2$TRANS) ~ "NA",
                                             TRUE                      ~ "abnrm"))

full_dt2$TRIG_binary <- as.factor(case_when(full_dt2$TRIG<=2 ~ "nrm",
                                            is.na(full_dt2$TRIG) ~ "NA",
                                            TRUE                      ~ "abnrm"))

full_dt2$TSH_binary <- as.factor(case_when(full_dt2$TSH>=0.4 & full_dt2$TSH<=4.8 ~ "nrm",
                                           is.na(full_dt2$TSH) ~ "NA",
                                           TRUE                      ~ "abnrm"))
# Identify binary features
binary_feats <- colnames(full_dt2)[grep("binary", names(full_dt2))]

# Define malignancy scores based on binary features
full_dt2$biochemical_abnormality_score <- apply(full_dt2[, ..binary_feats], 1, function(x) (sum(x == "abnrm", na.rm = TRUE) / length(binary_feats)) * 100)
full_dt2$biochemical_abnormality_score_NAexcluded <- apply(full_dt2[, ..binary_feats], 1, function(x) (sum(x == "abnrm", na.rm = TRUE) / sum(x %in% c("abnrm", "nrm"), na.rm = TRUE)) * 100)

# Select specific binary features for modified score
selected_binaryfeats <- paste0(c("CRP", "LEU", "PROCAL", "THROM", "ALAT", "CREA"), "_binary")
full_dt2$modified_biochemical_abnormality_score <- apply(full_dt2[, ..selected_binaryfeats], 1, function(x) (sum(x == "abnrm", na.rm = TRUE) / sum(x %in% c("abnrm", "nrm"), na.rm = TRUE)) * 100)

# Remove redundant binary features
full_dt2 <- full_dt2[, -..binary_feats]

# Extended filtering based on CRP value
if (extended_filtering){
  exclude_CRP_idx <- which(full_dt2$CRP < 20 | is.na(full_dt2$CRP))
  full_dt2 <- full_dt2[-exclude_CRP_idx,]
}

# Combine same features with different names
features_to_combine <- c("CAI", "BASO", "GLU", "HB", "LEU", "LYMFO", "MONO", "Sodium", "K", "CL")
for (feature in features_to_combine){
  full_dt2[[feature]] <- apply(cbind(full_dt2[[feature]], full_dt2[[paste0(feature, "POC")]]), 1, function(x) mean(x, na.rm = TRUE))
  full_dt2[[paste0(feature, "POC")]] <- NULL
}

# Temporal feature engineering
if (temporal_feature_engineering){
  full_dt2[, c("weekend", "week", "day", "year_day")] <- NULL
}

# Index of samples with BSI
BSI_idx <- which(full_dt2$BSIClass == "BSI")

# Identify mostly NA features
mostly_NA_feats <- colnames(full_dt2)[apply(full_dt2[BSI_idx, ], 2, function(x) sum(is.na(x)) / length(x) > 0.9)]

# Load required libraries
library("gtsummary")
library("dplyr")
library("forcats")

# Generate summary table
full_dt2[,-c("ID")] %>%      
  mutate(BSIClass = factor(BSIClass) %>% forcats::fct_explicit_na(na_level = "Unknown")) %>%
  tbl_summary(
    type = all_continuous() ~ "continuous2",
    missing = "no",
    statistic = list(
      all_continuous() ~ c("{median} ({p25}, {p75})", "{N_miss} ({p_miss}%)"),
      TestDate = "{range}"
    )
  ) %>%
  modify_table_body(
    mutate,
    label = ifelse(label == "N missing (% missing)",
                   "Unknown",
                   label)
  )

# Multiple imputation
if (do_multiple_imputation){
  library("mice")
  formulas <- make.formulas(full_dt_full)
  meth <- make.method(full_dt_full)
  imputedDataset <- parlmice(full_dt_full, method = meth, formulas = formulas, m = 5, n.core = 4, n.imp.core = 2)
  imputedDS <- mice::complete(imputedDataset, action = "long", include = FALSE)
  imputedDS <- imputedDS[, -c(1, 2)]
  fwrite(imputedDS, paste0("7DBSIdata_multimputed_", lookback_period, "days.csv"))
}

# Test the algorithm on a subset of the data on the local computer
if (test_on_subset){
  set.seed(seed_nr)
  subset_idx <- sort(sample(x = 1:nrow(full_dt2), size = 10000, replace = FALSE), decreasing = FALSE)
  full_dt2 <- full_dt2[subset_idx,]
}

# Define columns of interest for rolling means
cols <- c("Sodium", "ALAT", "ALB", "APTT", "BASO", "CA", "CAI", "CHOL", "CL", "CREA", "CRP", "EOS", "ERY", "FERRIT", "GGT", "GLU", "GLUF", "HAPTO", "HB", "HBA1C", "HDL", "INR", "JERN", "K", "LDH", "LDL", "LEU", "LYMFO", "MG", "MONO", "NEUTRO", "PROCAL", "THROM", "TRANS", "TRIG", "TSH")

# Create new columns for rolling means
for (j in cols) {
  full_dt2[, paste0(j, "_base") := NA_real_]
}

# Replace NAs with "Unknown"
full_dt2$BSIClass <- forcats::fct_explicit_na(full_dt2$BSIClass, "Unknown")

# Define the columns of interest
biochem_cols <- which(colnames(full_dt2) %in% cols)
biochem_base_cols <- paste0(cols, "_base")

# Create a unique list of IDs in full_dt2
unique_ids <- unique(full_dt2$ID)

# Loop over each unique ID
for (id in unique_ids){
  # Get the indices of rows with the current ID
  idx <- which(full_dt2$ID == id)
  
  # Loop over each row with the current ID
  for (i in idx){
    # Check if there are any rows with the same ID and BSIClass is not BSI within the lookback period
    idx1 <- which((full_dt2$BSIClass[idx] != "BSI") & (full_dt2$TestDate[idx] < full_dt2$TestDate[i]) & (full_dt2$TestDate[i] < (full_dt2$TestDate[idx] + lookback_period)))
    
    # If there are rows, calculate the mean of the biochemical variables
    if (length(idx1) > 0){
      # Subset the data to rows with idx1
      data_subset <- full_dt2[idx[idx1], ..biochem_cols]
      
      # Calculate the mean using sapply
      mean_values <- sapply(data_subset, mean, na.rm = TRUE)
      
      # Assign the mean values to the corresponding columns in the current row
      full_dt2[i, (biochem_base_cols) := as.list(as.numeric(mean_values))]
    }
  }
}

# Write the final result to a CSV file
fwrite(full_dt2, paste0("7DMerged_BSI_data_new_", lookback_period, "days.csv"))
full_dt2_binary <- full_dt2[-which(full_dt2$BSIClass == "Unknown"),]
# Summary table for baseline values of supplemental data
full_dt2_binary[,-c("ID")] %>%      
  mutate(BSIClass = factor(BSIClass) %>% forcats::fct_explicit_na(na_level = "Unknown")) %>%
  tbl_summary(
    type = all_continuous() ~ "continuous2",
    missing = "no",
    statistic = list(
      all_continuous() ~ c("{median} ({p25}, {p75})", "{N_miss} ({p_miss}%)"),
      TestDate = "{range}"
    )
  ) %>%
  modify_table_body(
    mutate,
    label = ifelse(label == "N missing (% missing)",
                   "Unknown",
                   label)
  )

# Most common pathogens analysis
baktr_plt <- bakteriaemi_data_wide[, -c("TestDate", "Sex", "Age", "most_common_pathogens")]
dim(baktr_plt)

baktr_plt <- baktr_plt[which(baktr_plt$ID %in% full_dt2$ID), ]
dim(baktr_plt)

baktr_plt <- baktr_plt[, -1]
baktr_plt <- baktr_plt[-which(baktr_plt$BSIClass == "noBSI"), ]
baktr_plt <- baktr_plt[, -("BSIClass")]

baktr_plt1 <- apply(baktr_plt, 2, function(x) as.numeric(factor(x, levels = c(FALSE, TRUE), labels = c(0, 1))) - 1)
pathogen_perc <- (sort(apply(baktr_plt1, 2, sum)) / nrow(baktr_plt1)) * 100

pathogen_perc_df <- as.data.frame(pathogen_perc)
pathogen_perc_df$pathogen <- rownames(pathogen_perc_df)
write_xlsx(pathogen_perc_df, path = paste0("pathogen_perc", seed_nr, ".xlsx"))

coinfection_perc <- apply(baktr_plt1, 1, sum)
round(prop.table(table(coinfection_perc)), 4)

any(is.na(baktr_plt1))

library("corrplot")
library("RColorBrewer")

pathogens_to_exclude <- names(which(apply(baktr_plt1, 2, var) == 0))
baktr_plt2 <- as.data.frame(baktr_plt1)
baktr_plt2 <- baktr_plt2[, -which(apply(baktr_plt1, 2, var) == 0)]

library("lares")

corr_cross(
  baktr_plt2,
  max_pvalue = 0.05,
  top = 10
)

all_NA_cols <- which(apply(full_dt2_binary, 1, function(x) {all(is.na(x))}))
total_patient_idx <- unique(full_dt2_binary$ID)

set.seed(seed_nr)
dev_idx <- sort(sample(total_patient_idx, round(0.8 * length(total_patient_idx))), decreasing = FALSE)
test_idx <- sort(total_patient_idx[which((total_patient_idx %in% dev_idx) == FALSE)], decreasing = FALSE)

test_instance_idx <- which(full_dt2_binary$ID %in% test_idx)
dev_instance_idx <- which(full_dt2_binary$ID %in% dev_idx)

any(dev_idx %in% test_idx)
any(test_idx %in% dev_idx)
any(dev_instance_idx %in% test_instance_idx)
any(test_instance_idx %in% dev_instance_idx)

if (do_data_dicotomize){
  categ_feat <- colnames(full_dt2_binary)[grep("binary", colnames(full_dt2_binary))]
  full_dt2_binary <-
    one_hot(
      dt = full_dt2_binary,
      cols = categ_feat,
      sparsifyNAs = FALSE,
      naCols = TRUE,
      dropCols = TRUE,
      dropUnusedLevels = TRUE
    )
  categ_feat_NA <- colnames(full_dt2_binary)[grep("_NA", colnames(full_dt2_binary))]

  full_dt2_binary[, categ_feat_NA] <- NULL
}

full_dt_tosave <- full_dt2_binary
any(is.na(full_dt_tosave))

devset_noFS_noIMPU <- full_dt_tosave[dev_instance_idx,]
testset_noFS_noIMPU <- full_dt_tosave[test_instance_idx,]

devset_noFS_noIMPU <- devset_noFS_noIMPU[, -c("ID", "TestDate")]
testset_noFS_noIMPU <- testset_noFS_noIMPU[, -c("ID", "TestDate")]

devset_noFS_noIMPU <- as.data.frame(devset_noFS_noIMPU)
testset_noFS_noIMPU <- as.data.frame(testset_noFS_noIMPU)

fwrite(devset_noFS_noIMPU, paste0("7DBSIdata_devset_", lookback_period, "days_noFS_noIMPU.csv"))
fwrite(testset_noFS_noIMPU, paste0("7DBSIdata_testset_", lookback_period, "days_noFS_noIMPU.csv"))

full_dt_tosave[is.na(full_dt_tosave)] <- -1
any(is.na(full_dt_tosave))

full_dt_tosave$Sex <- as.numeric(full_dt_tosave$Sex) - 1
full_dt_tosave$BSIClass <- as.numeric(full_dt_tosave$BSIClass) - 1

devset <- full_dt_tosave[dev_instance_idx,]
testset <- full_dt_tosave[test_instance_idx,]

# Save data including ID and TestDate
fwrite(devset, paste0("7DBSIdata_devset_", lookback_period, "days.csv"))
fwrite(testset, paste0("7DBSIdata_testset_", lookback_period, "days.csv"))

# Print summary of clinical (dev) data
print('Summary of clinical (dev) data:')
table(devset$BSIClass)
length(unique(devset$ID))
table(devset$BSIClass)

class_proportions <- prop.table(table(devset$BSIClass))

# Print summary of clinical (test) data
print('Summary of clinical (test) data:')
table(testset$BSIClass) 
length(unique(testset$ID))
table(testset$BSIClass)

class_proportions <- prop.table(table(testset$BSIClass))

library("report")  # To report packages that have been used
oursessionreport <- report(sessionInfo())
summary(oursessionreport)

oursessionreport_df <- as.data.frame(oursessionreport)
write_xlsx(oursessionreport_df, path = paste0("sessioninfo", seed_nr, ".xlsx"))
summary(as.data.frame(oursessionreport))

if (write_logs_to_file){
  closeAllConnections()  # Close connection to log file
}

# Find and remove constant clinical features
fwrite(devset[,-c("ID", "TestDate", "most_common_pathogens")], paste0("7DBSIdata_devset_nopatho_", lookback_period, "days.csv"))
fwrite(testset[,-c("ID", "TestDate", "most_common_pathogens")], paste0("7DBSIdata_testset_nopatho_", lookback_period, "days.csv"))

#############################################
######## mumpup data analysis example #######
#############################################


library(dplyr)
library(vegan)


  #### FUNDAMENTAL STEP BY STEP PROCESS ####

# preprocessing can be simplified by piping or by using the preprocess function


  ## PRE-PROCESS DATA ##

# 1. read data

myData<-read.csv("GCMS_data.csv", stringsAsFactors=F, check.names=F)

# 2. extract abundance values

myData_ea <- extractAbundance(myData)

# 3. split data into control samples and non control samples

myData_controls <- get_control(myData_ea)
myData_non_controls <- get_non_control(myData_ea)

# 4. filter/subset non control samples

myData_subset <- gen_subset_filter(myData_non_controls, "", "", "", "", "M", "", "", "")

# 5. remove abundances found in selected controls from non control samples

myData_subset <- remove_control(myData_controls, myData_subset, "KI_2016_C1_T6")

# 6. zero abundances for compounds which only occur in a single sample

myData_subset <- gen_zero_singles(myData_subset)

# 7. trim retention time

myData_subset <- gen_subset_select_rt(myData_subset, 10, 20)

# 8. write subset to file

export_file(myData_subset)


  ## ANALYSE DATA ##

# 1. get numerical only columns of data and perform wisconsin + sqrt transformations

myData_ws <- gen_prepared_data(myData_subset)

# 2. get numerical only columns of data and perform metaMDS transformation

myData_metaMDS <- gen_metaMDS(myData_subset)

# 3. execute adonis

  # one field 
myData_adonis_1 <- execute_adonis(myData_ws, myData_subset, "body")
  # two fields
myData_adonis_2 <- execute_adonis(myData_ws, myData_subset, "body", "age")

# 4. execute pairwise adonis


# 5. plot data
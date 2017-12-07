#############################################
##### functions for mumpup data analysis ####
#############################################


library(dplyr)
library(vegan)

##### PREPROCCESSING #####

# extract abundance values

extractAbundance <- function(data){
  data_extract <- data.frame(data, check.names = F)
  for(i in 26:length(data_extract)){
    if(class(data_extract[,i])!="integer"){
      for(j in 1:nrow(data)){
        if(data_extract[j,i]!="0")
          data_extract[j,i] <- convert(data_extract[j,i])
      }
    }
  }
  data_extract[,26:length(data_extract)]<- sapply(data_extract[,26:length(data_extract)], as.numeric)
  return(data_extract)
}

convert <- function (st) {
  E<-"E"
  value <-regexpr(E,st)[1]
  first<-substr(st,value-3,value-3)
  dec<-substr(st,value-1,value-1)
  numzero<-substr(st,value+1,value+1)
  return((as.numeric(first)*10+as.numeric(dec))*10^(as.numeric(numzero)-1))
}

# get control samples

get_control <- function(data){
  controls <- subset(data, (is.na(data[,5])) & (is.na(data[,6])) & (is.na(data[,7])))
  return(controls)
}

# get non-control samples

get_non_control <- function(data){
  non_controls <- subset(data, (!is.na(data[,5])) | (!is.na(data[,6])) | (!is.na(data[,7])))
}


# filter samples which have variable values which appear in each parameter - missing parameter includes all

gen_subset_filter <- function(data, locations, years, bodies, ages, sexes, nitrogen, mumpups, pouch){
  if(missing(locations)){
    locations <- unique(select(data, where))[,1]
  }
  if(missing(years)){
    years <- unique(select(data, year))[,1]
  }
  if(missing(bodies)){
    bodies <- unique(select(data, body))[,1]
  }
  if(missing(ages)){
    ages <- unique(select(data, age))[,1]
  }
  if(missing(sexes)){
    sexes <- unique(select(data, sex))[,1]
  }
  if(missing(nitrogen)){
    nitrogen <- unique(select(data, Nitrogen_error))[,1]
  }
  if(missing(mumpups)){
    mumpups <- unique(select(data, mumpup_pair))[,1]
  }
  if(missing(pouch)){
    pouch <- unique(select(data, pouches))[,1]
  }
  
  data_subset <- filter(data, where %in% locations,
                        year %in% years,
                        body %in% bodies,
                        age %in% ages,
                        sex %in% sexes,
                        nitrogen %in% Nitrogen_error,
                        mumpup_pair %in% mumpups,
                        pouches %in% pouch)
  return(data_subset)
}




#remove control samples from a subset (flexible with multiple inputs) 

remove_control <- function(data_control, data_non_control, w,...){
  input<<-c(w,...)
  for(i in 1:length(input)){
    for(a in 1:nrow(data_control)){
      if(match(input[i], data_control[[i]], nomatch=0)!="0"){
        index<- match(input[i],control[[1]])
        yearofindex<-data_control[index,3]
        locationofindex<-data_control[index,2]
        
        vector<-logical(length=length(data_control)-25)
        for(j in 26:length(data_control)){
          if(data_control[index,j]!=0){
            vector[j-25] <- TRUE
          } else{
            vector[j-25] <- FALSE
          }
        }
        
        for(k in 26:length(data_non_control)){
          if(isTRUE(vector[k-25])){
            for(m in 1:nrow(data_non_control)){
              if(data_non_control[m,2]==locationofindex & data_non_control[m,3]==yearofindex){
                data_non_control[m,k]<-0
              }
            }
          }
        }
      }
    }
    
  } 
  return(data_non_control)
}


# zero abundances for compounds which only occur in a single observation

gen_zero_singles <- function(data){
  data_zeros <- cbind(data[,26:length(data)])
  for(i in 1:length(data_zeros)){
    if(sum(data_zeros[,i] != 0) == 1)
      data_zeros[,i][data_zeros[,i] > 0] <- 0
  }
  data_zeros <- cbind(data[,1:25], data_zeros)
  return(data_zeros)
}


# transform data to relative abundance

gen_relative_abundance <- function(data){
  data_relative <- cbind(data[,26:length(data)])
  sums <- apply(data_relative, 1, sum)
  data_relative <- sweep(data_relative, 1, sums, "/")
  data_relative <- cbind(data[,1:25], data_relative)
  return(data_relative)
}


# select only columns with retention time in range

gen_subset_select_rt <- function(data, bound_lower, bound_upper){
  data_retentions <- cbind(data[,26:length(data)])
  column_names <- as.double(colnames(data_retentions))
  columns <- (column_names >= bound_lower) & (column_names <= bound_upper)
  data_retentions_subset <<- data_retentions[columns]
  data_subset <- cbind(data[,1:25], data_retentions_subset)
  return(data_subset)
}

# Creates csv file
filename_Safe <- function(string) {
  safeString <- gsub("[^[:alnum:]]", "_", string)
  safeString <- gsub("_+", "_", safeString)
  safeString
}

export_file <- function(data,subset){
  if(all.equal(data$where,subset$where)){
    where_string <- "ALL"
  } else{
    where_string<-filename_Safe(toString(levels(factor(subset$where))))
  }
  if(all.equal(data$year,subset$year)){
    year_string <- "ALL"
  } else{
    year_string<-filename_Safe(toString(levels(factor(subset$year))))
  }
  if(all.equal(data$body,subset$body)){
    body_string <- "ALL"
  } else{
    body_string<-filename_Safe(toString(levels(factor(subset$body))))
  }
  if(all.equal(data$age,subset$age)){
    age_string <- "ALL"
  } else{
    age_string<-filename_Safe(toString(levels(factor(subset$age))))
  }
  if(all.equal(data$sex,subset$sex)){
    sex_string <- "ALL"
  } else{
    sex_string<-filename_Safe(toString(levels(factor(subset$sex))))
  }
  if(all.equal(data$mumpup_pair,subset$mumpup_pair)){
    mumpup_pair_string <- "ALL"
  } else{
    mumpup_pair_string<-filename_Safe(toString(levels(factor(subset$mumpup_pair))))
  }
  if(all.equal(data$nitrogen_error,subset$nitrogen_error)){
    Nitrogen_error_string <- "ALL"
  } else{
    Nitrogen_error_string<-filename_Safe(toString(levels(factor(subset$Nitrogen_error))))
  }
  if(all.equal(data$pouches,subset$pouches)){
    pouches_string <- "ALL"
  } else{
    pouches_string<-filename_Safe(toString(levels(factor(subset$pouches))))
  }
  
  input_string<-filename_Safe(paste(input))
  mintime_string<-filename_Safe(names(data_retentions_subset)[1])
  maxtime_string<-filename_Safe(toString(names(data_retentions_subset)[length(data_retentions_subset)]))
  
  mystring<- paste("Where",where_string,"Year",year_string,"Body",body_string,"Age",age_string,"Sex",sex_string,"Nitrogen_error",Nitrogen_error_string,"Mumpup_pair",mumpup_pair_string,"Pouches",pouches_string, "Control", input_string, "Time", mintime_string, maxtime_string)
  write.csv(subset, paste0(filename_Safe(mystring),".csv"),row.names=F)
}



##### ANALYSIS #####

# get numerical columns and perform wisconsin + sqrt transforms

gen_prepared_data <- function(data){
  data_prepared <- cbind(data)
  data_prepared <- data_prepared %>% 
    gen_subset_numerical() %>%
    gen_wisconsin_sqrt()
  return(data_prepared)
}

# get numerical columns and perform metaMDS transform

gen_metaMDS <- function(data){
  data_transformed <- cbind(data)
  data_transformed <- data_transformed %>%
    gen_subset_numerical() %>%
    metaMDS()
  return(data_transformed)
}

# select only numerical columns for analysis

gen_subset_numerical <- function(data){
  data_subset <- data[,26:length(data)]
  return(data_subset)
}

# Wisconsin transform of sqare root of data

gen_wisconsin_sqrt <- function(data){
  data_transformed <- cbind(data)
  data_transformed <- data_transformed %>%
    sqrt() %>%
    wisconsin()
  return(data_transformed)
}

# execute adonis
 
 execute_adonis <- function(data_transformed, data, field_1){
   dt <- cbind(data_transformed)
   d <- cbind(data)
   s <- field_1
   ad <- adonis(dt~d$s)
   return(ad)
 } 

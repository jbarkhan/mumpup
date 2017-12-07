#############################################
##### functions for mumpup data analysis ####
#############################################


library(dplyr)
library(vegan)

##### PREPROCCESSING #####

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

#remove control samples from a subset (flexible with multiple inputs) 

remove_control <- function(data,w,...){
  temp<-data[,1]
  input<-c(w,...)
  for(i in 1:length(input)){
    if(input[i]!=""){
      if(match(input[i], temp, nomatch=0)!="0"){
        index<- match(input[i], temp, nomatch=0)
        data <- data[-index,]
      }
    } 
  }
  return(data)
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
  data_retentions_subset <- data_retentions[columns]
  data_subset <- cbind(data[,1:25], data_retentions_subset)
  return(data_subset)
}



##### ANALYSIS #####

# get numerical columns and perform wisconsin + sqrt transforms

gen_prepared_data <- function(data){
  data_prepared <- cbind(data)
  data_prepared <- data_prepared %>% 
    gen_wisconsin_sqrt() %>%
    gen_subset_numerical()
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
  data_subset <- select_if(data, is.numeric)
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
 
 execute_adonis <- function(data_transformed, data, field_1, field_2){
   if(missing(field_2)){
     temp <- adonis(data_transformed ~ data $ field_1)
     return(temp)
   }
   temp <- adonis(data_transformed ~ data $ field_1 * data $ field_2)
   return(temp)
 } 

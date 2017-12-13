#############################################
##### functions for mumpup data analysis ####
#############################################


library(dplyr)
library(vegan)
myData<-read.csv("GCMS_data.csv", stringsAsFactors=F, check.names=F)

#PRE-PROCESSING FUNCTION (If perform this, not required to perform all other functions below)

preprocess <- function(data, locations, year, bodies, ages, sexes, nitrogen, mumpup, pouch, lower_bound, upper_bound, remove){ #Insert "" if there's no input to the variable
  data_processed <- data %>% extractAbundance() %>% gen_subset_filter(locations, year, bodies, ages, sexes, nitrogen, mumpup, pouch)
  data_processed_control<-get_control(data_processed)
  data_processed_non_control<- get_non_control(data_processed)
  data_processed<-remove_control(data_processed_control, data_processed_non_control, remove)
  data_processed<- data_processed %>% gen_zero_singles() %>% gen_relative_abundance() %>% gen_subset_select_rt(lower_bound, upper_bound)
  return(data_processed)
}

#THE FOLLOWING FUNCTIONS BELOW MUST BE PERFORMED IN ORDER LABELLED 1 to 9 (If wish to do it seperately)

# 1. extract abundance values from original dataset
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

# helper function for extractAbundance (Don't need to run this)
convert <- function (st) {
  value <-regexpr("E",st)[1]
  first<-substr(st,value-3,value-3)
  dec<-substr(st,value-1,value-1)
  numzero<-substr(st,value+1,value+1)
  return((as.numeric(first)*10+as.numeric(dec))*10^(as.numeric(numzero)-1))
}

# 2. filter samples which have variable values which appear in each parameter - missing parameter or empty string includes all
gen_subset_filter <- function(data, locations=c(...), years=c(...), bodies=c(...), ages=c(...), sexes=c(...), nitrogen=c(...), mumpups=c(...), pouch=c(...)){
  if(locations==""){
    locations <- unique(select(data, where))[,1]
  } else if(all(match(locations,levels(factor(data$where)), nomatch=0))==FALSE){
    return("ERROR")
  }
  if(years==""){
    years <- unique(select(data, year))[,1]
  } else if(all(match(years,levels(factor(data$year)), nomatch=0))==FALSE){
    return("ERROR")
  }
  if(bodies==""){
    bodies <- unique(select(data, body))[,1]
  } else if(all(match(bodies,levels(factor(data$body)), nomatch=0))==FALSE){
    return("ERROR")
  }
  if(ages==""){
    ages <- unique(select(data, age))[,1]
  } else if(all(match(ages,levels(factor(data$age)), nomatch=0))==FALSE){
    return("ERROR")
  }
  if(sexes==""){
    sexes <- unique(select(data, sex))[,1]
  } else if(all(match(sexes,levels(factor(data$sex)), nomatch=0))==FALSE){
    return("ERROR")
  }
  if(nitrogen==""){
    nitrogen <- unique(select(data, Nitrogen_error))[,1]
  } else if(all(match(nitrogen,levels(factor(data$Nitrogen_error)), nomatch=0))==FALSE){
    return("ERROR")
  }
  if(mumpups==""){
    mumpups <- unique(select(data, mumpup_pair))[,1]
  } else if(all(match(mumpups,levels(factor(data$mumpup_pair)), nomatch=0))==FALSE){
    return("ERROR")
  }
  if(pouch==""){
    pouch <- unique(select(data, pouches))[,1]
  } else if(all(match(pouch,levels(factor(data$pouches)), nomatch=0))==FALSE){
    return("ERROR")
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

# 3. get control samples after obtain subset of data
get_control <- function(data){
  controls <- subset(data, (is.na(data[,5])) & (is.na(data[,6])))
  return(controls)
}

# 4. get non-control samples after obtain subset of data
get_non_control <- function(data){
  non_controls <- subset(data, (!is.na(data[,5])) | (!is.na(data[,6])))
}

# 5. remove control samples from a subset (flexible with multiple inputs) 
remove_control <- function(data_control, data_non_control, remove=c(...)){ #data_control (object from running get_control), data_non_control (object from running get_non_control), ... (Optional Arguments for names of control sample))
  input<<- remove
  
  if(remove==""){
    return(data_non_control)
  }
  
  for(i in 1:length(remove)){
    found <- FALSE
    for(j in 1:nrow(data_control) && found==FALSE){
      if(match(remove[i], data_control[[j]], nomatch=0)!="0"){
        index<- match(remove[i],data_control[[1]]) #index of the row related to the input name
        yearofindex<-data_control[index,3] 
        locationofindex<-data_control[index,2] 
        for(k in 26:length(data_control)){
          if(data_control[index,k]!=0){
            for(m in 1:nrow(data_non_control)){
              if(data_non_control[m,2]==locationofindex & data_non_control[m,3]==yearofindex){
                data_non_control[m,k]<-0 
              }
            }#end m loop
          } 
        }#end k loop
        found<- TRUE #set to true when the name from the input has been found in the"sample names(control)"
      }
    }#end j loop
    if(found==FALSE){
      return("ERROR")
    }
  }#end i loop
  return(data_non_control)
}

# 6. zero abundances for compounds which only occur in a single observation after remove controls
gen_zero_singles <- function(data){
  data_zeros <- cbind(data[,26:length(data)])
  zerovector <-vector()
  for(i in 1:length(data_zeros)){
    if(isTRUE(sum(data_zeros[1:NROW(data_zeros),i] != 0) == 1 || sum(data_zeros[1:NROW(data_zeros),i] != 0) == 0))
      zerovector<-c(zerovector,i)
  }
  data_zeros <- data_zeros[1:NROW(data_zeros),-zerovector]
  data_zeros <- cbind(data[,1:25], data_zeros)
  return(data_zeros)
}

# 7. transform data to relative abundance after remove zero abundances from previous step
gen_relative_abundance <- function(data){
  data_relative <- cbind(data[,26:length(data)])
  sums <- apply(data_relative, 1, sum)
  data_relative <- sweep(data_relative, 1, sums, "/")
  data_relative <- cbind(data[,1:25], data_relative)
  return(data_relative)
}

# 8. select only columns with retention time in range after transform data to relative abundance

gen_subset_select_rt <- function(data, bound_lower, bound_upper){
  data_retentions <- cbind(data[,26:length(data)])
  column_names <- as.double(colnames(data_retentions))
  columns <- (column_names >= bound_lower) & (column_names <= bound_upper)
  data_retentions_subset <<- data_retentions[columns]
  data_subset <- cbind(data[,1:25], data_retentions_subset)
  return(data_subset)
}

# 9. Creates csv file
filename_Safe <- function(string) {
  safeString <- gsub("[^[:alnum:]]", "_", string)
  safeString <- gsub("_+", "_", safeString)
  safeString
}

export_file <- function(subset){
  factorwhere<- levels(factor(c("DMM", "KI", "OI", "ZOO")))
  factoryear<- levels(factor(c("2015", "2016")))
  factorbody<- levels(factor(c("A", "B", "B/E", "Bel", "E", "Fb", "Fh", "H", "M", "Mammae", "N", "Tn")))
  factorage<- levels(factor(c("ad", "1", "2", "3", "SAM")))
  factorsex<- levels(factor(c("F", "M")))
  factormumpuppair<- levels(factor(c("Y", "N")))
  factorpouches<- levels(factor(c("Y", "N")))
  factornitrogenerror<- levels(factor(c("Y", "N")))
  
  if(all(levels(factor(subset$where))==factorwhere)){
    where_string<- "ALL"
  } else{
    where_string<-filename_Safe(toString(levels(factor(subset$where))))
  }
  
  if(all(levels(factor(subset$year))==factoryear)){
    year_string<-"ALL"
  }else{
    year_string<-filename_Safe(toString(levels(factor(subset$year))))
  }
  
  if(all(levels(factor(subset$body))==factorbody)){
    body_string<-"ALL"
  }else{
    body_string<-filename_Safe(toString(levels(factor(subset$body))))
  }
  
  if(all(levels(factor(subset$age))==factorage)){
    age_string<-"ALL"
  }else{
    age_string<-filename_Safe(toString(levels(factor(subset$age))))
  }
  
  if(all(levels(factor(subset$sex))==factorsex)){
    sex_string<-"ALL"
  }else{
    sex_string<-filename_Safe(toString(levels(factor(subset$sex))))
  }
  
  if(all(levels(factor(subset$mumpup_pair))==factormumpuppair)){
    mumpup_pair_string<-"ALL"
  }else{
    mumpup_pair_string<-filename_Safe(toString(levels(factor(subset$mumpup_pair))))
  }
  
  if(all(levels(factor(subset$Nitrogen_error))==factornitrogenerror)){
    Nitrogen_error_string<-"ALL"
  }else{
    Nitrogen_error_string<-filename_Safe(toString(levels(factor(subset$Nitrogen_error))))
  }
  
  if(all(levels(factor(subset$pouches))==factorpouches)){
    pouches_string<-"ALL"
  }else{
    pouches_string<-filename_Safe(toString(levels(factor(subset$pouches))))
  }
  
  input_string<-filename_Safe(toString(paste(input)))
  mintime_string<-filename_Safe(names(data_retentions_subset)[1])
  maxtime_string<-filename_Safe(toString(names(data_retentions_subset)[length(data_retentions_subset)]))
  
  mystring<- paste0("Where ",where_string," Year ",year_string," Body ",body_string," Age ",age_string," Sex ",sex_string," Nitrogen_error ",Nitrogen_error_string," Mumpup_pair ",mumpup_pair_string," Pouches ",pouches_string, " Control ", input_string, " Time ", mintime_string, maxtime_string)
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

execute_adonis <- function(data_transformed, data, field_1, field_2) {
  if(missing(field_2)){
    i <- grep(field_1, colnames(data))
    ad <- adonis(data_transformed ~ data[,i])
    return(ad)
  }
  i <- grep(field_1, colnames(data))
  j <- grep(field_2, colnames(data))
  ad <- adonis(data_transformed ~ data[,i] * data[,j])
  return(ad)
} 


# pairwise adonis function

pairwise.adonis <- function(x,factors, sim.function = 'vegdist', sim.method = 'bray', p.adjust.m ='holm')
{
  co = combn(unique(as.character(factors)),2)
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  
  
  for(elem in 1:ncol(co)){
    if(sim.function == 'daisy'){
      library(cluster); x1 = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],metric=sim.method)
    } else{x1 = vegdist(x[factors %in% c(co[1,elem],co[2,elem]),],method=sim.method)}
    
    ad = adonis(x1 ~ factors[factors %in% c(co[1,elem],co[2,elem])] );
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  sig = c(rep('',length(p.adjusted)))
  sig[p.adjusted <= 0.05] <-'.'
  sig[p.adjusted <= 0.01] <-'*'
  sig[p.adjusted <= 0.001] <-'**'
  sig[p.adjusted <= 0.0001] <-'***'
  
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted,sig)
  print("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1")
  return(pairw.res)
} 

#### PLOTTING ####
# Need to check if more than two fields will ever be required
gen_plot <- function(data_mds, data, field1, field2){
  
  i <- grep(field1, colnames(data))
  j <- grep(field2, colnames(data))
  
  data[,i] = factor(data[,i]) # Data needs to be refactored so 
  data[,j] = factor(data[,j]) # that values can be grouped
  
  plt1<-plot(data_mds, 
             display="sites",
             type= "n", 
             ylim=c(-2, 1.5), 
             cex.axis=1.2, 
             ylab="Dimension 1", 
             xlab="Dimension 2",
             cex.lab=1.3)
  cols1 <- c("black","red")
  points (plt1$sites, 
          pch=c(16, 17)[as.numeric(data[,j])],
          col=cols1[data[,i]], 
          cex=1.3)
  legend("bottomleft", 
         legend=c("Adult female", "Adult male", "Pup female", "Pup male"),
         pch=c(16, 17), col=c("black", "black","red", "red"), 
         cex=1.3)
  
  ordiellipse(data_mds, data[,i], col="black", show.groups='A', lwd=2.5)
  ordiellipse(data_mds, data[,i], col="red", show.groups='P', lwd=2.5)
}
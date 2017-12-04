library(dplyr)

# filter samples which have variable values which appear in each parameter - empty list includes all values

  gen_subset_filter <- function(data, locations, years, bodies, ages, sexes, nitrogen, mumpups, pouch){
    if(length(locations) == 0){
      locations <- unique(select(data, where))
    }
    if(length(years) == 0){
      years <- unique(select(data, year))
    }
    if(length(bodies) == 0){
      bodies <- unique(select(data, body))
    }
    if(length(ages) == 0){
      ages <- unique(select(data, age))
    }
    if(length(sexes) == 0){
      sexes <- unique(select(data, sex))
    }
    if(length(nitrogen) == 0){
      nitrogen <- unique(select(data, nitrogen_error))
    }
    if(length(mumpups) == 0){
      mumpups <- unique(select(data, mumpup_pair))
    }
    if(length(pouch) == 0){
      pouch <- unique(select(data, pouches))
    }

    data_subset <- filter(data, where %in% locations,
                          year %in% years,
                          body %in% bodies,
                          age %in% ages,
                          sex %in% sexes,
                          nitrogen %in% nitrogen_error,
                          mumpup_pair %in% mumpups,
                          pouches %in% pouch)
    return(data_subset)
  }


# select only columns with retention time in range

  gen_subset_select_rt <- function(data, bound_lower, bound_upper){
    data_subset <- select(data, num_range(bound_lower:bound_upper))
    return(data_subset)
}


# select only numerical columns for analysis

  gen_subset_numerical <- function(data){
    data_subset <- select_if(data, is.numeric)
    return(data_subset)
  }

cept <- function(x,except='gn_sp') # form merging control tables
  return(names(x)[!names(x) %in% except])



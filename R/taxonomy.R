
gs.check <- function(genus.species, sep='[ _\\.]'){
  gs.fmt.txt <- paste("'Genus",sep,"species' format",sep='')
  gs.pattern <- paste('^[A-Z][a-z]+',sep,'[a-z]+$',sep='')
  gs.ct <- sum(grepl(pattern=gs.pattern,x=genus.species))
  if(all(gs.ct))
   print(paste("all names appear to be of the",gs.fmt.txt))
  else
   stop(paste(gs.ct, "of of names are of the",gs.fmt.txt))
}

gs.rename <- function(df, alias.table.path, retro=FALSE, update.gn_sp=FALSE){ #uses rownames
  lookup <- read.delim(alias.table.path, stringsAsFactors=FALSE)

  .us2spc <- function(n){ gsub(x=n,'_',' ')}
  .spc2us <- function(n){ gsub(x=n,' ','_')}

  gs.check(rownames(df)) 
  
  print("replacing spaces or periods in (Genus?species names) with underscores")
  rownames(df) <- sub(x=rownames(df),'[ \\.]','_') 
  
  rn.us.ct <- sum(grepl(pattern='_',x=rownames(df))) #this check may be redundant now
  if(rn.us.ct != nrow(df)) stop(paste('number of rownames containing underscores (',rn.us.ct,') does not equal dataframe size (',nrow(df),')',sep=''))
 
  #will need to check for existing alpha rownames and offer to create them
  for(i in 1:nrow(lookup)){ 
   nw <- .us2spc(lookup$new.name[i])
   ld <- .us2spc(lookup$old.name[i])
   fnd=ld; rpl=nw; if(retro){ fnd=nw; rpl=ld;};#print(fnd);print(rpl)
   find <- match(fnd, .us2spc(rownames(df))); 
   if(length(find)>1) stop('found more than one match')#this should never happen!
   if(!is.na(find)){
    new.row.names<- .spc2us(rpl)
    if(sum(grepl(new.row.names, rownames(df)))>0){ 
       warning(paste('COULD NOT PEFORM DUPLICATE RENAME',fnd,'TO PREEXISTING',rpl));
       
     next()}else{
     print(paste('finding and renaming',fnd,'to',rpl));  
     rownames(df)[find] <- new.row.names
   }}
  }
  if(update.gn_sp){ df$gn_sp <- gsub(x=.us2spc(rownames(df)), " ","_")} #to matchup with phy in comparative}
  return(df)
}

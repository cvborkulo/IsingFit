netlogo <- function(x)
{
  if (is.character(x))
  {
    x[] <- paste0('"',x,'"')
  }
  if (is.vector(x))
  {
    return(paste("\n",paste(paste0("[",paste(x, collapse = " "),"]"),collapse="\n"),"\n\n"))
  } else if (is.matrix(x))
  {
    return(paste("[\n",paste(apply(x,1,function(s)paste0("[",paste(s, collapse = " "),"]")),collapse="\n"),"\n]")  )
  } else stop("Object not supported")
  
}
# 
# write(netlogo(x),file=paste0(name,'txt'))

# zo maken dat de weiadj en thresholds als aparte objecten eruit komen die als txt-file opgeslagen kunnen worden.
# Ook zo maken dat deze functie individueel opgeroepen kan worden. Bijv om de symptoomnamen in file te krijgen.
exportNetLogo.Risingfit <- function(object,....)
{
  if (is.character(object))
  {
    object[] <- paste0('"',object,'"')
  }
  if (is.vector(object))
  {
    return(paste("\n",paste(paste0("[",paste(object, collapse = " "),"]"),collapse="\n"),"\n\n"))
  } else if (is.matrix(object))
  {
    return(paste("[\n",paste(apply(object,1,function(s)paste0("[",paste(s, collapse = " "),"]")),collapse="\n"),"\n]")  )
  } else stop("Object not supported")
  
}
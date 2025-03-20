# Version 1.5 (October 2016)

wb.stats <- function(file, nodes=character(0), node.statistics=c("mean", "sd", "MC error", "2.5%", "median", "97.5%", "start", "sample"), return.as.list=length(node.statistics)>1)
{
  wblog <- scan(file, sep='\n', what="", quiet=T)
  tab.found <- grep('\t', wblog)
  wblog <- wblog[tab.found]
  
  # Remove node stats header line(s) -- if found twice, ignore lines below second occurence
  median.found <- grep('\\smedian', wblog, perl=T)
  if (length(median.found) > 1) wblog <- wblog[seq(median.found[2]-1)]
  wblog <- wblog[-median.found[1]]
  
  # Exhaustive list of node stats returned by WinBUGS
  wb.node.statistics <- c("mean", "sd", "MC error", "2.5%", "median", "97.5%", "start", "sample")
  
  # remove lines with comments (possible if output was modified)
  found.comment <- grep('#', wblog)
  if (length(found.comment)) wblog <- wblog[-found.comment]
  
  wblog <- gsub('\\s+', ' ', wblog, perl=T)
  wblog <- sub('\\s', '', wblog, perl=T) # remove leading space
  wblog <- wblog[nchar(wblog)>0] # drop empty lines
    # drop lines with /odds ratio/
    slashes <- grep("/", wblog)
    if (length(slashes)) wblog <- wblog[-slashes]
  wblog <- matrix(unlist(strsplit(wblog, ' ')), byrow=T, nrow=length(wblog))
  nodes.names <- wblog[,1]
  wblog <- matrix(suppressWarnings(as.numeric(wblog[,-1])), nrow=nrow(wblog))
  colnames(wblog) <- wb.node.statistics
  rownames(wblog) <- nodes.names
  
  # Keep only required node statistics
  w <- match(tolower(node.statistics), tolower(wb.node.statistics))
  wblog <- wblog[,w,drop=F]
  
  if (length(nodes))
  {
    nodes.names0 <- sub('\\[.*', '', nodes.names, perl=T)
    w <- !is.na(match(nodes.names0, nodes)) | !is.na(match(nodes.names, nodes))
    wblog <- wblog[w,,drop=F]
  }
  
  
  if (return.as.list)
  {
    stats <- sort(colnames(wblog))
    n.stats <- length(stats)

    dim <- rownames(wblog)
    dim <- sub("\\[.*", "", dim, perl=T)
    dimensions <- unique(sort(dim))
    n.dim <- length(dimensions)
  
    out <- list()
  
    for (i in seq(n.dim))
    {
      d <- dimensions[i]
      w <- which(dim==d)
    
      my.out <- list()
      my.i <- 0
      for (s in stats)
      {
        my.i <- my.i + 1
        my.out[[my.i]] <- wblog[w, s]
      } 
    
      if (n.stats == 1)
      {
        my.out <- unlist(my.out)
      }
      else
      {
        names(my.out) <- stats
      }
    
      out[[i]] <- my.out
    }
  
    names(out) <- dimensions
  }
  else
  {
    out <- wblog
  }
  
  out
} # end of wb.stats

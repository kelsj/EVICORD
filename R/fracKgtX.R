#' fracKgtX.R
#' for a vector v, returns fraction of values greater than x

fracKgtX = function(v,x){ return(sum(v>x)/length(v))}

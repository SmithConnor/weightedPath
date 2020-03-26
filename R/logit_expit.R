#' @title logit function
#' @param x numeric
#' @export

logit = function(x){
  if(any(x <= 0) | any(x >= 1)){
    error("x must be between zero and 1")
  }

  return(log(x) - log(1-x))
}


#' @title Expit function
#' @param x numeric
#' @export

expit = function(x){
  return(1/(1+exp(-x)))
}

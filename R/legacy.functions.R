## Legacy function
ffirst.lagm <- function(PKindex,...) {
  ffirst.lag(PKindex,...,Tlag=TRUE,MMe=TRUE)
}

## Legacy function
ffirst.nolagm <- function(PKindex,...) {
  ffirst.nolag(PKindex,...,Tlag=FALSE,MMe=TRUE)  ### original was 'ffirst.lag(PKindex,...,Tlag=FALSE,MMe=TRUE)'; I wondered if it's correct.
}

## Legacy function
sfirst.lag <- function(...) {
  sfirst.nolag(...,Tlag=TRUE,MMe=FALSE)
}


## Legacy function
sfirst.lagm <- function(...) {
  sfirst.nolag(...,Tlag=TRUE,MMe=TRUE)
}

## Legacy function
sfirst.nolagm <- function(...) {
  sfirst.nolag(...,Tlag=FALSE,MMe=TRUE)
}

## Legacy function
fzero.nolagm <- function(PKindex,...) {
  fzero.nolag(PKindex,...,MMe=TRUE)
}

## Legacy function
szero.nolagm <- function(...) {
  szero.nolag(...,MMe=TRUE)
}

## Legacy function
finfu.mm <- function(PKindex,...) {
  finfu1(PKindex,...,MMe=TRUE)
}

## Legacy function
sinfu.mm <- function(...) {
  sinfu1(...,MMe=TRUE)
}

## Legacy function
fbolus.mm <- function(PKindex,...) {
  fbolus1(PKindex,...,MMe=TRUE)
}

## Legacy function
sbolus.mm <- function(...) {
  sbolus1(...,MMe=TRUE)
}

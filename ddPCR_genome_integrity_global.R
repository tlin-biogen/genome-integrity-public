# Purpose: functions to be used in the shiny app globally.

process_csv <- function(csv.temp){

  if(names(csv.temp)[1] != "row.names"){
    csv <- csv.temp
  } else{
    names(csv.temp) <- c(names(csv.temp)[-1], "ToBeRemoved")
    csv.temp$ToBeRemoved <-NULL
    csv <- csv.temp
  }
  return(csv)
  
}


# target: the count of target molecule; AD = # of accepted droplet
# size = how many trials
# seq(n):  1....n out of m trials are a "hit"
polysolver <- function(m, AD, target, n) {
  p <- as.polynomial(c(-1, dbinom(seq(n), size = m, prob = 1 / AD) * AD / target))
  p
  pz <- solve(p)
  good_pz <- Re(pz[which(abs(Im(pz)) < 1e-6)])
  root1 <- good_pz[good_pz >= 0 & good_pz <= 1]
  stopifnot(length(root1) == 1)
  return(root1)
}


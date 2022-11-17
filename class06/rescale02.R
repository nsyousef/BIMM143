rescale02 <- function(x, na.rm=TRUE) {
  rng <- range(x, na.rm=na.rm)
  (x - rng[1]) / (rng[2] - rng[1])
}

# In R console, run source("rescale02.R") to use this function in the console.
# Can also run source("rescale02.R") in a file if you want to use this function
# there.


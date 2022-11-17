add <- function(x, y=1) {
  x + y
}

add(x=4)

add(x=4, y=10)

add(c(1, 2, 5, 10), y=10)
#add(c(1, 2, 5, 10), 10, 20) # gives error since too many arguments

add(c(1, 2, 5, 10), y="string") # also gives an error since arg is not right type


# Our second function: rescale()

rescale <- function(x) {
  rng <- range(x)
  (x - rng[1]) / (rng[2] - rng[1])
}

# Test on a simple example
rescale(1:10)

rescale(c(1, 4, 10, 20))

rescale(c(1, NA, 4, 10, 20)) #NA is a missing data point, like null I think
# We want to figure out how we want to handle this input
# We want it to give back NA only for the slot that it is in

# The fixed function (NA issue fixed):
rescale01 <- function(x) {
  rng <- range(x, na.rm=TRUE)
  (x - rng[1]) / (rng[2] - rng[1])
}

rescale01(c(1, NA, 4, 10, 20))

rescale01(c(1, 4, "word"))

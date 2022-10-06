# Central limit theorem functions

## Means

#' Calculate the geometric mean of a numeric vector
#'
#' @param x a numeric vector
#' @param na.rm  TRUE or FALSE
#'
#' @return the geometric mean
#' @export
#'
#' @examples
#' a <- c(1,2,3,4,5,5,6); geo_mean(a)
geo_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
}

#' Bootstrapped upper confidence limit
#'
#' Just a wrapper around Hmisc.cl.boot() to extract upper cl
#' @param x a numeric vector
#' @param ... arguments to be passed to smean.cl.boot()
#'
#' @return
#' @export
#'
#' @examples
#' df <- data.frame(y = seq(1:10), x = c(51,52,53,54,55,55,54,53,52,99), g = c(rep("a",5),rep("b",5)))
#' library(tidyverse);library(magrittr)
#' df %>% group_by(g) %>% dplyr::summarise(u = upper_cl_boot(x,B=100))
upper_cl_boot = function(x, ...){
  y <- Hmisc::smean.cl.boot(x, ...)[3]
  y
}

#' Bootstrapped lower confidence limit
#'
#' Just a wrapper around Hmisc.cl.boot() to extract lower cl
#' @param x a numeric vector
#' @param ... arguments to be passed to smean.cl.boot()
#'
#' @return
#' @export
#'
#' @examples
#' df <- data.frame(y = seq(1:10), x = c(51,52,53,54,55,55,54,53,52,99), g = c(rep("a",5),rep("b",5)))
#' library(tidyverse);library(magrittr)
#' df %>% group_by(g) %>% dplyr::summarise(u = lower_cl_boot(x,B=100))
lower_cl_boot = function(x, ...){
  y <- Hmisc::smean.cl.boot(x, ...)[2]
  y
}

## Modes

#' Calculate the mode of a numeric or character vector
#'
#' @param x a numeric or character vector
#' @param na.rm TRUE or FALSE
#'
#' @return value with the highest number of occurrences in x
#' @export
#'
#' @examples
#' a <- c(1,2,3,4,5,5,6); get_mode(a)
#' b <- c("dog","cat","dog","bird"); get_mode(b)
#' d <- c(1,2,3,4,NA,NA,NA,5,5,6); get_mode(d) # also works with NA
get_mode <- function(x, na.rm = FALSE) {
  uniqx <- unique(x)
  uniqx[which.max(tabulate(match(x, uniqx)))]
}

#' Estimate the mode of a numeric distribution
#'
#' @param x a numeric vector
#'
#' @return estimated value being the maximum value in the density function of a numeric vector
#' @export
#'
#' @examples
#' a <- c(1,2,3,4,5,5,6); estimate_mode(a)
estimate_mode <- function(x) {
  d <- density(x, na.rm = TRUE)
  d$x[which.max(d$y)]
}

## Variances

#' CV
#' calculates CV in the classic way assuming a normal distribution. sensitive to skewed data and outliers
#' @param x a vector of numeric values
#' @param na.rm default = TRUE
#'
#' @return the coefficient of variation
#' @export
#'
#' @examples
#' n <- 100
#' sd <- 0.2
#' u <- 0.5
#' set.seed(123)
#' y <- rnorm(n = n, mean = u, sd = sd)
#' cv(y)
cv <- function(x, na.rm = T){
  sd(x,na.rm = na.rm)/mean(x, na.rm = na.rm)
}

#' Relative standard deviation
#'
#' @param x a vector of numeric values
#' @param na.rm
#'
#' @return
#' @export
#'
#' @examples
#' n <- 100
#' sd <- 0.2
#' u <- 0.5
#' set.seed(123)
#' y <- rnorm(n = n, mean = u, sd = sd)
#' rsd(y)
rsd <-  function(x, na.rm = T){

  abs(sd(x, na.rm = na.rm)/mean(x, na.rm = na.rm))
}


#' Robust CV based on quantiles
#' As per Arachige et al 2020 arXiv:1907.01110v3
#' @param x a vector of numeric values
#' @param na.rm default = TRUE
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' n <- 100
#' sd <- 0.2
#' u <- 0.5
#' set.seed(123)
#' y <- rnorm(n = n, mean = u, sd = sd)
#' rcvq(y)
rcvq <- function(x, na.rm = TRUE, ...){
  0.75 * (IQR(c(x), na.rm = TRUE)/median(x, na.rm = TRUE))
}


#' Robust CV based on MAD
#' As per Arachige et al 2020 arXiv:1907.01110v3
#' @param x a vector of numeric values
#' @param na.rm TRUE or FALSE
#' @param constant a constant value you supply
#' @param ... other arguments to be passed
#'
#' @return
#' @export
#'
#' @examples
#' n <- 100
#' sd <- 0.2
#' u <- 0.5
#' set.seed(123)
#' y <- rnorm(n = n, mean = u, sd = sd)
#' rcvm(y)
rcvm <- function(x, na.rm = T, constant, ...){
  mad(x, na.rm = TRUE, constant = 1.4826)/median(x, na.rm = TRUE)
}

## Outliers

#' Remove outliers
#'
#' @param x a numeric vector
#' @param na.rm TRUE or FALSE
#' @param ... other arguments to be passed
#'
#' @return
#' @export
#'
#' @examples
#' x <- c(1,2,3,4,5,6,5,4,3,2,1,99); remove_outliers(x)
#' df <- data.frame(y = seq(1:10), x = c(51,52,53,54,55,55,54,53,52,99), g = c(rep("a",5),rep("b",5)))
#' library(tidyverse);library(magrittr)
#' df %<>% dplyr::group_by(g) %>% dplyr::summarise(x = remove_outliers(x))
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

# Functions to summarise data

#' Calculate closest value in a vector to the specified value
#'
#' @param x a numeric vector
#' @param sv a specific value for which you want to find the closest value in a vector
#'
#' @return
#' @export
#'
#' @examples
#' x <- c(1,2,3,4,5,6,5,4,3,2,1,99); sv = 100; closest(x, sv = sv)
closest <- function(x, sv, na.rm = T){
  if(na.rm)
    x <- x[!is.na(x)]
  a <- which(abs(x - sv) == min(abs(x - sv)))
  b <- x[a]
  b
}


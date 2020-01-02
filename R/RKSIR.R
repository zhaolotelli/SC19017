#' @title GKernel
#' @description Radial Gaussian Kernel
#' @param x the first vector
#' @param z the second vector
#' @param sigma bandwidth
#' @return a numeric value
#' \item{k}{value of k(x, z)}
#' @examples
#' \dontrun{
#' x = c(1, 0); z = c(0, 1)
#' GKernel(x, z)
#' }
#' @export
GKernel <- function(x, z, sigma = 2) 
{
  k <- exp(-sum((x - z)^2)/(2 * sigma^2))
  return(k)
}

#' @title AGKernel
#' @description Additive Gaussian Kernel
#' @param x the first vector
#' @param z the second vector
#' @param sigma bandwidth
#' @return a numeric value
#' \item{k}{value of k(x, z)}
#' @examples
#' \dontrun{
#' x = c(1, 0); z = c(0, 1)
#' AGKernel(x, z)
#' }
#' @export
AGKernel <- function(x, z, sigma = 2) 
{
  k <- sum(exp(-(x - z)^2/(2 * sigma^2)))
  return(k)
}

#' @title PKernel
#' @description Second Order Polynomial Kernel
#' @param x the first vector
#' @param z the second vector
#' @param sigma bandwidth
#' @return a numeric value
#' \item{k}{value of k(x, z)}
#' @examples
#' \dontrun{
#' x = c(1, 0); z = c(0, 1)
#' PKernel(x, z)
#' }
#' @export
PKernel <- function(x, z, sigma = 1) {
  k <- (1 + sum(x * z)/sigma^2)^2
  return(k)
}

#' @title KX
#' @description Compute kernel matrix
#' @param X data matrix
#' @param kernel kernelfunction, "AG" means \code{\link{AGKernel}}, "G" means \code{\link{GKernel}}, "P" means \code{\link{PKernel}}
#' @param sigma bandwidth
#' @param center logic value, centralize or not
#' @return a \code{n} times \code{n} matrix
#' \item{K}{Gram matrix of the raw data}
#' @examples
#' \dontrun{
#' X <- matrix(1:10, 2, 5)
#' KX(X)
#' }
#' @export
KX <- function(X, kernel = "AG", sigma = 2, center = TRUE) 
{
  if(kernel == "AG") {
    kernelfunc <- AGKernel
  }
  if(kernel == "P") {
    kernelfunc <- PKernel
  }
  if(kernel == "G") {
    kernelfunc <- GKernel
  }
  n <- nrow(X)
  d <- ncol(X)
  K <- matrix(0, nrow = n, ncol = n)
  for(i in 1:n) {
    for(j in 1:n) {
      K[i, j] <- kernelfunc(X[i, ], X[j, ], sigma)
    }
  }
  if(center) {
    K <- (diag(1, n, n) - 1/n * matrix(1, n, n)) %*% K %*% (diag(1, n, n) - 
                                                              1/n * matrix(1, n, n))
  }
  return(K)
}

#' @title KX2
#' @description Compute kernel matrix of X1 with X2
#' @param X1,X2 data matrixs
#' @param kernel kernelfunction, "AG" means \code{\link{AGKernel}}, "G" means \code{\link{GKernel}}, "P" means \code{\link{PKernel}}
#' @param sigma bandwidth
#' @param center logic value, centralize or not
#' @return a kernel matrix
#' \item{K}{Gram matrix of the raw data}
#' @examples
#' \dontrun{
#' X1 <- matrix(1:20, 4, 5)
#' X2 <- matrix(1:10, 2, 5)
#' KX2(X1, X2)
#' }
#' @export
KX2 <- function(X1, X2, kernel = "AG", sigma = 2, center = TRUE) 
{
  if(kernel == "AG") {
    kernelfunc <- AGKernel
  }
  if(kernel == "P") {
    kernelfunc <- PKernel
  }
  if(kernel == "G") {
    kernelfunc <- GKernel
  }
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  K <- matrix(0, nrow = n1, ncol = n2)
  for(i in 1:n1) {
    for(j in 1:n2) {
      K[i, j] <- kernelfunc(X1[i, ], X2[j, ], sigma)
    }
  }
  if(center) {
    K <- (diag(1, n1, n1) - 1/n1 * matrix(1, n1, n1)) %*% K %*% (diag(1, n2, n2) - 
                                                                   1/n2 * matrix(1, n2, n2))
  }
  return(K)
}

#' @title Jy
#' @description Compute the slice matrix J
#' @param y response vector
#' @param numeric logic value, numeric response or categorial response
#' @param h if response is numeric, give the number of slices for \code{y}
#' @return a slice matrix
#' \item{J}{slice matrix of the response}
#' @examples
#' \dontrun{
#' y <- sample(1:15, 15)
#' Jy(y, h = 5)
#' }
#' @export
Jy <- function(y, numeric = TRUE, h = NULL) 
{
  n <- length(y)
  J <- matrix(0, nrow = n, ncol = n)
  if(numeric) {
    # when y is numerical
    qs <- seq(1/h, 1, 1/h)
    ps <- quantile(y, qs)
    ysorted <- y
    for(i in 2:h) {
      ysorted[which(y <= ps[i] & y > ps[i - 1])] <- i
    }
    ysorted[which(y <= ps[1])] <- 1
    t <- as.vector(table(ysorted))
    for(i in 1:n) {
      for(j in 1:n) {
        if(ysorted[i] == ysorted[j]) {
          J[i, j] <- 1/t[ysorted[i]]
        }
      }
    }
  }
  else {
    # when y is categorical
    yf <- as.factor(y)
    ns <- summary(yf)
    for(i in 1:n) {
      for(j in 1:n) {
        if(y[i] == y[j]) {
          J[i, j] <- 1/ns[y[i]]
        }
      }
    }
  }
  return(J)
}

#' @title Ei
#' @description Compute the eigensystem of RKSIR
#' @param Kx kernel matrix, better be the result of \code{\link{KX}}
#' @param J slice matrix, better be the result of \code{\link{Jy}}
#' @param s regularization parameter
#' @return a eigen analysis result
#' \item{ei}{the eigensystem of RKSIR}
#' @examples
#' \dontrun{
#' X <- matrix(rnorm(30), 15, 2)
#' Kx <- KX(X)
#' y <- sample(1:15, 15)
#' J <- Jy(y, h = 5)
#' Ei(Kx, J, s = 0.1)
#' }
#' @export
Ei <- function(Kx, J, s) 
{
  n <- ncol(Kx)
  A <- t(Kx) %*% J %*% Kx
  B <- t(Kx) %*% Kx + n^2 * s * diag(1, n, n)
  Sigma <- solve(B) %*% A
  ei <- eigen((Sigma + t(Sigma))/2)
  return(ei)
}

#' @title Vtr
#' @description Compute the e.d.r. direction on the training data
#' @param Xtrain training data
#' @param Xsple if data is too large, can choose some sample of training data
#' @param p dimensions of e.d.r. directions
#' @param ei eigen analysis result of \code{\link{Ei}}
#' @param kernel kernel functions, "AG" means \code{\link{AGKernel}}, "G" means \code{\link{GKernel}}, "P" means \code{\link{PKernel}}
#' @param sigma bandwidth
#' @return a matrix with column vector being e.d.r. direction on the training set
#' \item{Vtr}{column vectors are the e.d.r. directions of training data}
#' @examples
#' \dontrun{
#' X <- matrix(rnorm(30), 15, 2)
#' Kx <- KX(X)
#' y <- sample(1:15, 15)
#' J <- Jy(y, h = 5)
#' ei <- Ei(Kx, J, s = 0.1)
#' Vtr(X, X, p = 2, ei)
#' }
#' @export
Vtr <- function(Xtrain, Xsple = Xtrain, p, ei, kernel = "AG", sigma = 2) 
{
  Ktr <- KX2(Xtrain, Xsple, kernel, sigma, center = TRUE)
  alpha <- as.matrix(ei$vectors[, 1:p])
  Vtr <- apply(Ktr, 1, function(o) {
    t(alpha) %*% o
  })
  if(p == 1) {
    return(t(t(Vtr)))
  }
  return(t(Vtr))
}

#' @title Vte
#' @description Compute the e.d.r. direction on the testing data
#' @param Xtest testing data
#' @param Xtrain training data
#' @param Xsple if data is too large, can choose some sample of training data
#' @param p dimensions of e.d.r. directions
#' @param ei eigen analysis result of \code{\link{Ei}}
#' @param kernel kernel functions, "AG" means \code{\link{AGKernel}}, "G" means \code{\link{GKernel}}, "P" means \code{\link{PKernel}}
#' @param sigma bandwidth
#' @return a matrix with column vector being e.d.r. direction on the testing set
#' \item{Vtr}{column vectors are the e.d.r. directions of testing data}
#' @examples
#' \dontrun{
#' Xtrain <- matrix(rnorm(30), 15, 2)
#' Xtest <- matrix(rnorm(20), 10, 2)
#' Kx <- KX(Xtrain)
#' y <- sample(1:15, 15)
#' J <- Jy(y, h = 5)
#' ei <- Ei(Kx, J, s = 0.1)
#' Vte(Xtest, Xtrain, Xtrain, p = 2, ei)
#' }
#' @export
Vte <- function(Xtest, Xtrain, Xsple, p, ei, kernel = "AG", sigma = 2) 
{
  if(kernel == "AG") {
    kernelfunc <- AGKernel
  }
  if(kernel == "P") {
    kernelfunc <- PKernel
  }
  if(kernel == "G") {
    kernelfunc <- GKernel
  }
  n1 <- nrow(Xtrain)
  n2 <- nrow(Xsple)
  m <- nrow(Xtest)
  Ktr <- KX2(Xtrain, Xsple, kernel, sigma, center = FALSE)
  Kte.raw <- t(apply(Xtest, 1, function(o) {
    apply(Xsple, 1, function(oo) {
      kernelfunc(o, oo, sigma)
    })
  }))
  Kte <- (Kte.raw - 1/n1 * matrix(1, nrow = m, ncol = n1) %*% Ktr) %*% 
    (diag(1, n2, n2) - 1/n2 * matrix(1, n2, n2))
  alpha <- as.matrix(ei$vectors[, 1:p])
  Vte <- apply(Kte, 1, function(o) {
    t(alpha) %*% o
  })
  if(p == 1) {
    return(t(t(Vte)))
  }
  return(t(Vte))
}

#' @title split
#' @description Split data randomly into training set and testing set
#' @param y response vector, must be 1,2,...
#' @param ns a vector which contains the length of each class in the training data
#' @return a vector
#' \item{train}{a vector which contains the number of data in the training set}
#' @examples
#' \dontrun{
#' y <- c(rep(1, 7), rep(2, 10), rep(3, 12))
#' ns <- c(3, 4, 5)
#' split(y, ns)
#' }
#' @export
split <- function(y, ns = NULL) {
  n <- length(y)
  train <- NULL
  for(i in 1:length(ns)) {
    train <- c(train, sample((1:n)[which(y == i)], ns[i]))
  }
  return(train)
}
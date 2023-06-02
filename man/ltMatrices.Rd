\name{ltMatrices}
\alias{ltMatrices}
\alias{Tcrossprod}
\alias{Crossprod}
\alias{diagonals}
\alias{Mult}
\alias{solve.ltMatrices}
\alias{as.array.ltMatrices}
\alias{as.array.syMatrices}
\alias{diagonals.ltMatrices}
\alias{diagonals.syMatrices}
\alias{diagonals.matrix}
\alias{diagonals<-}
\alias{diagonals<-.ltMatrices}
\alias{diagonals<-.syMatrices}
\alias{diagonals.integer}
\alias{Lower_tri}
\alias{chol.syMatrices}
\alias{aperm.ltMatrices}
\alias{adddiag}
\alias{chol2cov}
\alias{invchol2chol}
\alias{chol2invchol}
\alias{invchol2cov}
\alias{invchol2pre}
\alias{chol2pre}
\alias{Dchol}
\alias{invcholD}
\alias{chol2cor}
\alias{invchol2cor}
\alias{invchol2pc}
\alias{chol2pc}
\alias{vectrick}
\alias{standardize}
\alias{destandardize}
\title{
  Multiple Lower Triangular Matrices
}
\description{
  A class representing multiple lower triangular matrices and some methods.
}
\usage{
ltMatrices(object, diag = FALSE, byrow = FALSE, names = TRUE)
\S3method{as.array}{ltMatrices}(x, symmetric = FALSE, \dots)
\S3method{as.array}{syMatrices}(x, \dots)
\S3method{diagonals}{ltMatrices}(x, \dots)
\S3method{diagonals}{syMatrices}(x, \dots)
\S3method{diagonals}{matrix}(x, \dots)
\S3method{diagonals}{integer}(x, \dots)
diagonals(x) <- value
\S3method{diagonals}{ltMatrices}(x) <- value
\S3method{diagonals}{syMatrices}(x) <- value
\S3method{solve}{ltMatrices}(a, b, transpose = FALSE, \dots)
\S3method{chol}{syMatrices}(x, \dots)
\S3method{aperm}{ltMatrices}(a, perm, is_chol = FALSE, \dots)
Lower_tri(x, diag = FALSE, byrow = attr(x, "byrow"))
Mult(x, y, transpose = FALSE)
Tcrossprod(x, diag_only = FALSE)
Crossprod(x, diag_only = FALSE)
chol2cov(x)
invchol2chol(x)
chol2invchol(x)
invchol2cov(x)
invchol2pre(x)
chol2pre(x)
Dchol(x, D = 1 / sqrt(Tcrossprod(x, diag_only = TRUE)))
invcholD(x, D = sqrt(Tcrossprod(solve(x), diag_only = TRUE)))
chol2cor(x)
invchol2cor(x)
chol2pc(x)
invchol2pc(x)
vectrick(C, S, A, transpose = c(TRUE, TRUE))
standardize(chol, invchol)
destandardize(chol = solve(invchol), invchol, score_schol)
}
\arguments{
  \item{object}{a \code{matrix} representing the lower triagular elements of
                \eqn{N} lower triangular matrix, each of dimension \eqn{J \times
                J}. Dimensions of \code{object} depend on \code{diag}: With
                diagonal elements, \code{object} is a \eqn{J(J+1)/2 \times N} matrix,
                otherwise, the number of rows is \eqn{J(J - 1) / 2}.
}
  \item{diag}{logical, \code{object} contains diagonal elements if
              \code{TRUE}, otherwise unit diagonal elements are assumed.
}
  \item{byrow}{logical, \code{object} represents matrices in row-major order
               if \code{TRUE} or, otherwise, in column-major order.
}
  \item{names}{logical or character vector of length \eqn{J}.
}
  \item{symmetric}{logical, object is interpreted as a symmetric matrix if
\code{TRUE}.
}
  \item{diag_only}{logical, compute diagonal elements of crossproduct only
if \code{TRUE}.
}
  \item{x,chol,invchol}{object of class \code{ltMatrices} or \code{syMatrices} (for \code{chol}).
}
  \item{value}{a matrix of diagonal elements to be assigned (of dimension \eqn{J \times N}).
}
  \item{a}{object of class \code{ltMatrices}.
}
  \item{perm}{a permutation of the covariance matrix corresponding to \code{a}.
}
  \item{is_chol}{a logical indicating if \code{a} is the Cholesky of the
covariance (\code{chol = TRUE}) of the precision matrix.
}
  \item{D}{a matrix (of dimension \eqn{J \times N}) of diagonal elements to be multiplied with.
}
  \item{y}{matrix with \eqn{J} rows.
}
  \item{b}{matrix with \eqn{J} rows.
}
  \item{C}{an object of class \code{ltMatrices}.}
  \item{S}{an object of class \code{ltMatrices} or a matrix with \eqn{J^2} rows 
           representing multiple \eqn{J x J} matrices (columns of vec operators).}
  \item{A}{an object of class \code{ltMatrices}.}
  \item{transpose}{a logical of length two indicating if \code{A} or \code{B} shall be transposed in \code{vectrick}.
                   For \code{solve}, this argument being true computes
                   \code{solve(t(a), b)} (in absence of a \code{t()} method for \code{ltMatrices}
                   objects).}
  \item{score_schol}{score matrix for a standardized \code{chol} object.}
  \item{\dots}{additional arguments, currently ignored.
}
}
\details{
\code{ltMatrices} interprets a matrix as lower triangular elements of
multiple lower triangular matrices. The corresponding class can be used to
store such matrices efficiently. Matrix multiplications, solutions to linear
systems, explicite inverses, and crossproducts can be computed based on such
objects. Details can be found in the \code{lmvnorm_src} package vignette.
}
\value{
The constructor \code{ltMatrices} returns objects of class \code{ltMatrices}
with corresponding methods.
}
\seealso{\code{vignette("lmvnorm_src", package = "mvtnorm")}}
\examples{

  J <- 4L
  N <- 2L
  dm <- paste0("d", 1:J)
  xm <- paste0("x", 1:N)
  (C <- ltMatrices(matrix(runif(N * J * (J + 1) / 2), 
                          ncol = N, dimnames = list(NULL, xm)), 
                   diag = TRUE, names = dm))

  ## dimensions and names
  dim(C)
  dimnames(C)
  names(C)

  ## subset
  C[,2:3]

  ## multiplication
  y <- matrix(runif(N * J), nrow = J)
  Mult(C, y)

  ## solve
  solve(C)
  solve(C, y)

  ## tcrossprod
  Tcrossprod(C)

  ## convert to matrix
  as.array(solve(C[1,]))[,,1]

}
\keyword{matrix}

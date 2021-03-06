\name{SPC}
\alias{SPC}
\title{Perform sparse principal component analysis}
\description{
Performs sparse principal components analysis by applying PMD to a data
matrix with lasso ($L_1$) penalty on the columns and no penalty on the rows.
}
\usage{
SPC(x, sumabsv=4, niter=20, K=1, orth=FALSE, trace=TRUE, v=NULL,
center=TRUE, cnames=NULL, vpos=FALSE, vneg=FALSE, compute.pve=TRUE)
}
\arguments{
  \item{x}{Data matrix of dimension $n x p$, which can contain NA
    for missing values. We are interested in finding sparse principal
    components of dimension $p$.}
  \item{sumabsv}{How sparse do you want v to be? This is the sum of
    absolute values of elements of v. It must be between 1 and square
    root of number of columns of data. The smaller it is, the sparser v
    will be.}
  \item{niter}{How many iterations should be performed. It is best
    to run at least 20 of so. Default is 20.}
  \item{K}{The number of factors in the PMD to be returned; default is
    1.}
  \item{v}{The first right singular vector(s) of the data. (If missing
    data is present, then the missing values are imputed before the
    singular vectors are calculated.) v is used as the initial value for
    the iterative PMD($L_1$, $L_1$) algorithm. If x is large, then this
    step can be time-consuming; therefore, if PMD is to be run multiple
    times, then v should be computed once and saved.}
  \item{trace}{Print out progress as iterations are performed? Default
    is TRUE.}
  \item{orth}{If TRUE, then
    use method of Section 3.2 of Witten, Tibshirani and Hastie (2008) to
    obtain multiple sparse principal components. Default is FALSE.}
  \item{center}{Subtract out mean of x? Default is TRUE}
  \item{cnames}{An optional vector containing a name for each column.}
  \item{vpos}{Constrain the elements of v to be positive? TRUE or
    FALSE.}
  \item{vneg}{Constrain the elements of v to be negative? TRUE or
    FALSE.}
  \item{compute.pve}{Compute percent variance explained? Default
    TRUE. If not needed, then choose FALSE to save time.}
}
\details{
PMD(x,sumabsu=sqrt(nrow(x)), sumabsv=3, K=1) and SPC(x,sumabsv=3, K=1) give the
same result, since the SPC method is simply PMD with an L1 penalty on
the columns and no penalty on the rows.


In Witten, Tibshirani, and Hastie (2008), two methods are presented for
obtaining multiple factors for SPC. The methods are as follows:

(1) If one has already obtained factors $k-1$ factors then oen can
compute residuals by subtracting out these factors. Then $u_k$ and $v_k$ can be obtained
by applying the SPC/PMD algorithm to the residuals.

(2) One can require that $u_k$ be orthogonal to $u_i$'s with $i<k$; the
method is slightly more complicated, and is explained in WT&H(2008).

Method 1 is performed by running SPC with option orth=FALSE (the
default) and Method 2 is performed using option orth=TRUE. Note that
Methods 1 and 2 always give identical results for the first component,
and often given quite similar results for later components.
}
\value{
  \item{u}{u is output. If you asked for multiple factors then each
    column of u is a factor. u has dimension nxK if you asked for K factors.}
  \item{v}{v is output. These are the sparse principal components. If you asked for multiple factors then each
    column of v is a factor. v has dimension pxK if you asked for K
    factors.}
  \item{d}{d is output; it is the diagonal of the matrix $D$ in the
    penalized matrix decomposition. In the case of the rank-1
    decomposition, it is given in the formulation $||X-duv'||_F^2$ subject to $||u||_1 <= sumabsu$,
    $||v||_1 <= sumabsv$. Computationally, $d=u'Xv$ where $u$ and $v$
    are the sparse factors output by the PMD function and $X$ is the
    data matrix input to the PMD function.}
  \item{prop.var.explained}{A vector containing the proportion of
    variance explained by the first 1, 2, ..., K sparse principal
    components obtaineds. Formula for proportion of variance explained
    is on page 20 of Shen & Huang (2008), Journal of Multivariate
    Analysis 99: 1015-1034.}
  \item{v.init}{The first right singular vector(s) of the data; these
    are returned to save on computation time if PMD will be run again.}
  \item{meanx}{Mean of x that was subtracted out before SPC was performed.}
}
\references{Witten, DM and Tibshirani, R and T Hastie (2008) A penalized
  matrix decomposition, with applications to
  sparse principal components and canonical correlation
  analysis. Submitted. <http://www-stat.stanford.edu/~dwitten>}
\author{Daniela M. Witten and Robert Tibshirani}
\examples{
        
# A simple simulated example
set.seed(1)
u <- matrix(c(rnorm(50), rep(0,150)),ncol=1)
v <- matrix(c(rnorm(75),rep(0,225)), ncol=1)
x <- u\%*\%t(v)+matrix(rnorm(200*300),ncol=300)
# Perform Sparse PCA - that is, decompose a matrix w/o penalty on rows
# and w/ L1 penalty on columns
# First, we perform sparse PCA and get 4 components, but we do not
# require subsequent components to be orthogonal to previous components
out <- SPC(x,sumabsv=3, K=4)
print(out,verbose=TRUE)
# We could have selected sumabsv by cross-validation, using function SPC.cv
# Now, we do sparse PCA using method in Section 3.2 of WT&H(2008) for getting
# multiple components - that is, we require components to be orthogonal
out.orth <- SPC(x,sumabsv=3, K=4, orth=TRUE)
print(out.orth,verbose=TRUE)
par(mfrow=c(1,1))
plot(out$u[,1], out.orth$u[,1], xlab="", ylab="")
# Note that the first components w/ and w/o orth option are identical,
# since the orth option only affects the way that subsequent components
# are found
print(round(t(out$u)\%*\%out$u,4)) # not orthogonal
print(round(t(out.orth$u)\%*\%out.orth$u,4)) # orthogonal

# Use SPC.cv to choose tuning parameters:
cv.out <- SPC.cv(x)
print(cv.out)
plot(cv.out)
out <- SPC(x, sumabsv=cv.out$bestsumabsv)
print(out)
# or we could do
out <- SPC(x, sumabsv=cv.out$bestsumabsv1se)
print(out)
                   
}


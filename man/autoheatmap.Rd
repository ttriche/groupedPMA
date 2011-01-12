\name{autoheatmap}
\alias{autoheatmap}
\title{Automated heatmap production for plotting cross-correlated predictors}
\description{
  This function automates the production of heatmaps from supervised or
  unsupervised canonical correlation and principal component analyses.
}
\usage{autoheatmap(dat, design, x.covs, z.covs, y.cov=NULL, pen.x, pen.z, ...)}
\details{
  autoheatmap is a convenience function for quickly dissecting rich datasets 
  based on orthogonal decomposition of the cross-covariance matrices between 
  various assays or covariates.  Further details are provided in the vignette.
}
\examples{

data(AML)

# methylation changes at homeobox gene promoters X mutation
HOXloci <- featureNames(AML)[which(fData(AML)$HOX==TRUE)]
CCA.HOX <- autoheatmap(AML, pen.x=0.5, pen.z=0.3, 
                       covs=colnames(AML.colorcoded), 
                       colorbar=AML.colorcoded,
                       loci=HOXloci)
AML.mutations[CCA.HOX$u<0] # DNMT3A, FLT3, NPM1

}

#### calculate over-dispersed genes for a gene x cell matrix (used for the C1 dataset)
# counts.nodups: normalized gene (row) x cell (column) matrix
sel.by.cv <- function(counts.nodups) {
  library(statmod); library(fastICA)
  ed <- counts.nodups
  winsorize <- function (x, fraction=0.05) {
    if(length(fraction) != 1 || fraction < 0 ||
       fraction > 0.5) {
      stop("bad value for 'fraction'")
    }
    lim <- quantile(x, probs=c(fraction, 1-fraction))
    x[ x < lim[1] ] <- lim[1]
    x[ x > lim[2] ] <- lim[2]
    x
  }
  wed <- t(apply(ed, 1, winsorize, fraction=2/ncol(ed))) 
  dim(wed)
  means = rowMeans(wed); vars = apply(wed,1,var); cv2 <- vars/means^2
  useForFit <- means >= unname( quantile( means[ which( cv2 > .3 ) ], .95 ) ) 
  fit <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/means[useForFit] ),cv2[useForFit] )
  xg <- exp(seq( min(log(means[means>0])), max(log(means), na.rm=T), length.out=1000 ))
  afit <- fit$coef["a1tilde"]/means+fit$coef["a0"]
  vfit <- fit$coef["a1tilde"]/xg+fit$coef["a0"]
  varFitRatio <- vars/(afit*means^2)
  varorder <- order(varFitRatio,decreasing=T)
  return(varorder)
}
#' @rdname mvnpdf
#' @export
mvnpdfoptim <- function(x, mean =  rep(0, nrow(x)),
                        varcovM = diag(nrow(x)), Log=TRUE){

  n <- ncol(x)
  p <- nrow(x)
  x0 <- x-mean
  Rinv = backsolve(chol(varcovM),x=diag(p)) # équivalent de solve mais un peu plus rapide (choleski)
  xRinv <- matrix(apply(X=x0, MARGIN=2, FUN=crossprod, y=Rinv)) # crossprod :
  logSqrtDetvarcovM <- sum(log(diag(Rinv)))

  quadform <- apply(X=xRinv, MARGIN=2, FUN=crossprod)
  y <- (-0.5*quadform + logSqrtDetvarcovM -p*log(2*pi)/2)

  if(!Log){
    y <- exp(y)
  }

  return(y)
}

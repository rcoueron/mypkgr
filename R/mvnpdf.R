#' mvnpdf
#'
#' densite d’une loi normale multivariee sur R^p en n points.
#'
#' @param x x : une matrice, a n colonnes (les observations) et p lignes
#' @param mean un vecteur de moyennes
#' @param varcovM une matrice de variance-covariance
#' @param Log un paramètre logique valant TRUE par défaut
#'
#' @return renvoie une liste contenant la matrice x ainsi qu’un vecteur des images
#' des points de x par la fonction de densite de la variable aleatoire de loi normale multivariée consideree.
#'
#' @export
#'
#' @examples
#' mvnpdf(x=matrix(1.96),log=FALSE)
#' dnorm(1.96)
#'
mvnpdf <- function(x, mean =  rep(0, nrow(x)),
                   varcovM = diag(nrow(x)), Log = TRUE) {
  n <- ncol(x)
  p <- nrow(x)
  x0 <- x - mean
  # si inverse dans la boucle prend du temps
  Rinv <- solve(varcovM)
  LogDetvarcovM <- log(det(varcovM))
  # si log = TRUE on applique le logarithme au résultat
  # meilleur precision avec le log car avec exp on trouvera 0
  y <- NULL
  for (j in 1:n) {
    yj <- - p/2 * log(2*pi) - 0.5 * LogDetvarcovM -
      0.5 * t(x0[, j]) %*% Rinv %*% x0[, j]
    y <- c(y, yj)
  }
  if (!Log) {
    y <- exp(y)
  }
  res <- return(list(x=x,y=y))
  class(res) <- "mvnpdf"
  return(res)
}

#' Plot of the mvnpdf function
#'
#' @param x an object of class \code{mvnpdf} resulting from "????"
#' \code{mvnpdf()} function.
#' @param ... graphical parameters passed to \code{plot()} "????"
#'
#' @return Nothing is returned, only a plot is given
#' @export
#' @importFrom graphics plot
#' @examples
#' pdfvalues <- mvnpdf(x=matrix(seq(-3, 3, by=0.1), nrow "????")
#' plot(pdfvalues)
plot.mvnpdf <- function(x,...){
  graphics::plot(x$x, x$y, type="l", ...)
}


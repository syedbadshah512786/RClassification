#' Leave-one-out (Lachenbruch "Holdout") Procedure for equal and unequal variance
#'
#' @param x is a 1st Group
#' @param y is a 2nd group
#' @param i is such a variable which you went to test whether it from 1st or 2nd Group
#' @param var is use to select a method whether it has for equal or unequal variance (Var="Equal" for equal and Var="NotEqual" for unequal variance)
#' @param g is use to define a group i.e the test value is belong from 1 group or from 2nd group (g=1 for 1st Group and g=2 for 2nd group)
#'
#' @author Syed Hammad Raza
#' @export

s.holdout <-   function(x,y,i,var,g)
{
  if(var=="Equal")
  {
    if(g==1)
    {
      cat("Remember: You test a",i,"Observation of 1st Group when variance is Equal \n")
      xo <- x[i,]
      x1 <- x[-i,]
      x2 <- y

      n1<-nrow(x1)
      n2<-nrow(x2)

      x1_bar <- colMeans(x1)
      cat("Mean of 1st Population is: \n")
      print(x1_bar)

      x2_bar <- colMeans(x2)
      cat("Mean of 2nd Population is: \n")
      print(x2_bar)

      sigma_x1 <- cov(x1)
      cat("Variance of 1st Population is: \n")
      print(sigma_x1)

      sigma_x2 <- cov(x2)
      cat("Variance of 2nd Population is: \n")
      print(sigma_x2)

      Sp <- (((n1-1) / (n1-1+n2-1))*sigma_x1) + (((n2-1) / (n1-1+n2-1))*sigma_x2)
      cat("Pooled Variance is given as: \n")
      print(Sp)

      Sp.inv <- solve(Sp)
      cat("Inverse of Pooled Variance is given as: \n")
      print(Sp.inv)

      a_hat <- t(x1_bar - x2_bar) %*% Sp.inv
      cat("a_hat is given as: \n")
      print(a_hat)

      y1_bar <- a_hat %*% x1_bar
      cat("Y1_bar is given as: \n")
      print(y1_bar)

      y2_bar <- a_hat %*% x2_bar
      cat("Y2_bar is given as: \n")
      print(y2_bar)

      m_hat <- 1/2*(y1_bar + y2_bar)
      cat("m_hat is given as: \n")
      print(m_hat)

      y_hat <- a_hat %*% xo
      cat("Y_hat is given as: \n")
      print(y_hat)

      if (y_hat >= m_hat) {cat( "The",i,"Observation of 1st Group is Correctly classified")} else {cat(" The",i,"Observation of 1st Group is Misclassified") }
    }
    if(g==2)
    {
      cat("Remember: You test a",i,"Observation of 2nd Group when variance is Equal \n")
      xo <- y[i,]
      x1 <- x
      x2 <- y[-i,]

      n1<-nrow(x1)
      n2<-nrow(x2)

      x1_bar <- colMeans(x1)
      cat("Mean of 1st Population is: \n")
      print(x1_bar)

      x2_bar <- colMeans(x2)
      cat("Mean of 2nd Population is: \n")
      print(x2_bar)

      sigma_x1 <- cov(x1)
      cat("Variance of 1st Population is: \n")
      print(sigma_x1)

      sigma_x2 <- cov(x2)
      cat("Variance of 2nd Population is: \n")
      print(sigma_x2)

      Sp <- (((n1-1) / (n1-1+n2-1))*sigma_x1) + (((n2-1) / (n1-1+n2-1))*sigma_x2)
      cat("Pooled Variance is given as: \n")
      print(Sp)

      Sp.inv <- solve(Sp)
      cat("Inverse of Pooled Variance is given as: \n")
      print(Sp.inv)

      a_hat <- t(x1_bar - x2_bar) %*% Sp.inv
      cat("a_hat is given as: \n")
      print(a_hat)

      y1_bar <- a_hat %*% x1_bar
      cat("Y1_bar is given as: \n")
      print(y1_bar)

      y2_bar <- a_hat %*% x2_bar
      cat("Y2_bar is given as: \n")
      print(y2_bar)

      m_hat <- 1/2*(y1_bar + y2_bar)
      cat("m_hat is given as: \n")
      print(m_hat)

      y_hat <- a_hat %*% xo
      cat("Y_hat is given as: \n")
      print(y_hat)

      if (y_hat >= m_hat) {cat( "The",i,"Observation of 2nd Group is Misclassified")} else {cat(" The",i,"Observation of 2nd Group is Correctly classified") }
    }
  }
  if(var=="NotEqual")
  {
    if(g==1)
    {
      cat("Remember: You test a",i,"Observation of 1st Group when variance is Not Equal \n")
      x1 <- x[-i,]
      x2 <- y
      xo  <- x[i,]

      n1 <- nrow(x1)
      n2 <- nrow(x2)

      x1_bar <- colMeans(x1)
      cat("Mean of 1st Population is: \n")
      print(x1_bar)

      x2_bar <- colMeans(x2)
      cat("Mean of 2nd Population is: \n")
      print(x2_bar)

      sigma_x1 <- cov(x1)
      cat("Variance of 1st Population is: \n")
      print(sigma_x1)

      sigma_x2 <- cov(x2)
      cat("Variance of 2nd Population is: \n")
      print(sigma_x2)

      sigma_x1.inv <- solve(sigma_x1)
      cat("Inverse Variance of 1st Population is: \n")
      print(sigma_x1.inv)

      sigma_x2.inv <- solve(sigma_x2)
      cat("Inverse Variance of 2nd Population is: \n")
      print(sigma_x2.inv)

      det.sigma_x1 <- det(sigma_x1)
      cat("Determinent of variance of 1st Population is: \n")
      print(det.sigma_x1)

      det.sigma_x2 <- det(sigma_x2)
      cat("Determinent of variance of 2nd Population is: \n")
      print(det.sigma_x2)

      RHS <- 1/2 * log(det.sigma_x1 / det.sigma_x2) + 1/2 * (t(x1_bar) %*% sigma_x1.inv %*% x1_bar - t(x2_bar) %*% sigma_x2.inv %*% x2_bar)
      cat("The RHS of the Formula for Unequal Variance is or it is the value of K: \n")
      print(RHS)

      LHS <- -1/2 * xo %*% (sigma_x1.inv - sigma_x2.inv) %*% xo + (t(x1_bar) %*% sigma_x1.inv - t(x2_bar) %*% sigma_x2.inv) %*% xo
      cat("The LHS of the Formula for Unequal Variance is: \n")
      print(LHS)

      if(LHS >= RHS){cat( "The",i,"Observation of 1st Group is Correctly classified")} else {cat(" The",i,"Observation of 1st Group is Misclassified") }
    }
    if(g==2)
    {
      cat("Remember: You test a",i,"Observation of 2nd Group when variance is Not Equal \n")
      x1  <- x
      x2  <- y[-i,]
      xo  <- y[i,]

      n1 <- nrow(x1)
      n2 <- nrow(x2)

      x1_bar <- colMeans(x1)
      cat("Mean of 1st Population is: \n")
      print(x1_bar)

      x2_bar <- colMeans(x2)
      cat("Mean of 2nd Population is: \n")
      print(x2_bar)

      sigma_x1 <- cov(x1)
      cat("Variance of 1st Population is: \n")
      print(sigma_x1)

      sigma_x2 <- cov(x2)
      cat("Variance of 2nd Population is: \n")
      print(sigma_x2)

      sigma_x1.inv <- solve(sigma_x1)
      cat("Inverse Variance of 1st Population is: \n")
      print(sigma_x1.inv)

      sigma_x2.inv <- solve(sigma_x2)
      cat("Inverse Variance of 2nd Population is: \n")
      print(sigma_x2.inv)

      det.sigma_x1 <- det(sigma_x1)
      cat("Determinent of variance of 1st Population is: \n")
      print(det.sigma_x1)

      det.sigma_x2 <- det(sigma_x2)
      cat("Determinent of variance of 2nd Population is: \n")
      print(det.sigma_x2)

      RHS <- 1/2 * log(det.sigma_x1 / det.sigma_x2) + 1/2 * (t(x1_bar) %*% sigma_x1.inv %*% x1_bar - t(x2_bar) %*% sigma_x2.inv %*% x2_bar)
      cat("The RHS of the Formula for Unequal Variance is or it is the value of K: \n")
      print(RHS)

      LHS <- -1/2 * xo %*% (sigma_x1.inv - sigma_x2.inv) %*% xo + (t(x1_bar) %*% sigma_x1.inv - t(x2_bar) %*% sigma_x2.inv) %*% xo
      cat("The LHS of the Formula for Unequal Variance is: \n")
      print(LHS)

      if(LHS <= RHS){cat( "The",i,"Observation of 2nd Group is Correctly classified")} else {cat(" The",i,"Observation of 2nd Group is Misclassified") }
    }
  }
}

#' Classification of Multivariate Data Procedure for equal and unequal variance
#'
#' @param x is a 1st Group
#' @param y is a 2nd group
#' @param i is such a variable which you went to test whether it from 1st or 2nd Group
#' @param var is use to select a method whether it has for equal or unequal variance (Var="Equal" for equal and Var="NotEqual" for unequal variance)
#' @param g is use to define a group i.e the test value is belong from 1 group or from 2nd group (g=1 for 1st Group and g=2 for 2nd group)
#'
#' @author Syed Hammad Raza
#' @export

CF.method <-   function(x,y,i,var,g)
{
  if(var=="Equal")
  {
    if(g==1)
    {
      cat("Remember: You test a",i,"Observation of 1st Group when variance is Equal \n")
      xo <- x[i,]
      x1 <- x
      x2 <- y

      n1<-nrow(x1)
      n2<-nrow(x2)

      x1_bar <- colMeans(x1)
      cat("Mean of 1st Population is: \n")
      print(x1_bar)

      x2_bar <- colMeans(x2)
      cat("Mean of 2nd Population is: \n")
      print(x2_bar)

      sigma_x1 <- cov(x1)
      cat("Variance of 1st Population is: \n")
      print(sigma_x1)

      sigma_x2 <- cov(x2)
      cat("Variance of 2nd Population is: \n")
      print(sigma_x2)

      Sp <- (((n1-1) / (n1-1+n2-1))*sigma_x1) + (((n2-1) / (n1-1+n2-1))*sigma_x2)
      cat("Pooled Variance is given as: \n")
      print(Sp)

      Sp.inv <- solve(Sp)
      cat("Inverse of Pooled Variance is given as: \n")
      print(Sp.inv)

      a_hat <- t(x1_bar - x2_bar) %*% Sp.inv
      cat("a_hat is given as: \n")
      print(a_hat)

      y1_bar <- a_hat %*% x1_bar
      cat("Y1_bar is given as: \n")
      print(y1_bar)

      y2_bar <- a_hat %*% x2_bar
      cat("Y2_bar is given as: \n")
      print(y2_bar)

      m_hat <- 1/2*(y1_bar + y2_bar)
      cat("m_hat is given as: \n")
      print(m_hat)

      y_hat <- a_hat %*% xo
      cat("Y_hat is given as: \n")
      print(y_hat)

      if (y_hat >= m_hat) {cat( "The",i,"Observation of 1st Group is Correctly classified")} else {cat(" The",i,"Observation of 1st Group is Misclassified") }
    }
    if(g==2)
    {
      cat("Remember: You test a",i,"Observation of 2nd Group when variance is Equal \n")
      xo <- y[i,]
      x1 <- x
      x2 <- y

      n1<-nrow(x1)
      n2<-nrow(x2)

      x1_bar <- colMeans(x1)
      cat("Mean of 1st Population is: \n")
      print(x1_bar)

      x2_bar <- colMeans(x2)
      cat("Mean of 2nd Population is: \n")
      print(x2_bar)

      sigma_x1 <- cov(x1)
      cat("Variance of 1st Population is: \n")
      print(sigma_x1)

      sigma_x2 <- cov(x2)
      cat("Variance of 2nd Population is: \n")
      print(sigma_x2)

      Sp <- (((n1-1) / (n1-1+n2-1))*sigma_x1) + (((n2-1) / (n1-1+n2-1))*sigma_x2)
      cat("Pooled Variance is given as: \n")
      print(Sp)

      Sp.inv <- solve(Sp)
      cat("Inverse of Pooled Variance is given as: \n")
      print(Sp.inv)

      a_hat <- t(x1_bar - x2_bar) %*% Sp.inv
      cat("a_hat is given as: \n")
      print(a_hat)

      y1_bar <- a_hat %*% x1_bar
      cat("Y1_bar is given as: \n")
      print(y1_bar)

      y2_bar <- a_hat %*% x2_bar
      cat("Y2_bar is given as: \n")
      print(y2_bar)

      m_hat <- 1/2*(y1_bar + y2_bar)
      cat("m_hat is given as: \n")
      print(m_hat)

      y_hat <- a_hat %*% xo
      cat("Y_hat is given as: \n")
      print(y_hat)

      if (y_hat >= m_hat) {cat( "The",i,"Observation of 2nd Group is Misclassified")} else {cat(" The",i,"Observation of 2nd Group is Correctly classified") }
    }
  }
  if(var=="NotEqual")
  {
    if(g==1)
    {
      cat("Remember: You test a",i,"Observation of 1st Group when variance is Not Equal \n")
      x1 <- x
      x2 <- y
      xo  <- x[i,]

      n1 <- nrow(x1)
      n2 <- nrow(x2)

      x1_bar <- colMeans(x1)
      cat("Mean of 1st Population is: \n")
      print(x1_bar)

      x2_bar <- colMeans(x2)
      cat("Mean of 2nd Population is: \n")
      print(x2_bar)

      sigma_x1 <- cov(x1)
      cat("Variance of 1st Population is: \n")
      print(sigma_x1)

      sigma_x2 <- cov(x2)
      cat("Variance of 2nd Population is: \n")
      print(sigma_x2)

      sigma_x1.inv <- solve(sigma_x1)
      cat("Inverse Variance of 1st Population is: \n")
      print(sigma_x1.inv)

      sigma_x2.inv <- solve(sigma_x2)
      cat("Inverse Variance of 2nd Population is: \n")
      print(sigma_x2.inv)

      det.sigma_x1 <- det(sigma_x1)
      cat("Determinent of variance of 1st Population is: \n")
      print(det.sigma_x1)

      det.sigma_x2 <- det(sigma_x2)
      cat("Determinent of variance of 2nd Population is: \n")
      print(det.sigma_x2)

      RHS <- 1/2 * log(det.sigma_x1 / det.sigma_x2) + 1/2 * (t(x1_bar) %*% sigma_x1.inv %*% x1_bar - t(x2_bar) %*% sigma_x2.inv %*% x2_bar)
      cat("The RHS of the Formula for Unequal Variance is or it is the value of K: \n")
      print(RHS)

      LHS <- -1/2 * xo %*% (sigma_x1.inv - sigma_x2.inv) %*% xo + (t(x1_bar) %*% sigma_x1.inv - t(x2_bar) %*% sigma_x2.inv) %*% xo
      cat("The LHS of the Formula for Unequal Variance is: \n")
      print(LHS)

      if(LHS >= RHS){cat( "The",i,"Observation of 1st Group is Correctly classified")} else {cat(" The",i,"Observation of 1st Group is Misclassified") }
    }
    if(g==2)
    {
      cat("Remember: You test a",i,"Observation of 2nd Group when variance is Not Equal \n")
      x1  <- x
      x2  <- y
      xo  <- y[i,]

      n1 <- nrow(x1)
      n2 <- nrow(x2)

      x1_bar <- colMeans(x1)
      cat("Mean of 1st Population is: \n")
      print(x1_bar)

      x2_bar <- colMeans(x2)
      cat("Mean of 2nd Population is: \n")
      print(x2_bar)

      sigma_x1 <- cov(x1)
      cat("Variance of 1st Population is: \n")
      print(sigma_x1)

      sigma_x2 <- cov(x2)
      cat("Variance of 2nd Population is: \n")
      print(sigma_x2)

      sigma_x1.inv <- solve(sigma_x1)
      cat("Inverse Variance of 1st Population is: \n")
      print(sigma_x1.inv)

      sigma_x2.inv <- solve(sigma_x2)
      cat("Inverse Variance of 2nd Population is: \n")
      print(sigma_x2.inv)

      det.sigma_x1 <- det(sigma_x1)
      cat("Determinent of variance of 1st Population is: \n")
      print(det.sigma_x1)

      det.sigma_x2 <- det(sigma_x2)
      cat("Determinent of variance of 2nd Population is: \n")
      print(det.sigma_x2)

      RHS <- 1/2 * log(det.sigma_x1 / det.sigma_x2) + 1/2 * (t(x1_bar) %*% sigma_x1.inv %*% x1_bar - t(x2_bar) %*% sigma_x2.inv %*% x2_bar)
      cat("The RHS of the Formula for Unequal Variance is or it is the value of K: \n")
      print(RHS)

      LHS <- -1/2 * xo %*% (sigma_x1.inv - sigma_x2.inv) %*% xo + (t(x1_bar) %*% sigma_x1.inv - t(x2_bar) %*% sigma_x2.inv) %*% xo
      cat("The LHS of the Formula for Unequal Variance is: \n")
      print(LHS)

      if(LHS <= RHS){cat( "The",i,"Observation of 2nd Group is Correctly classified")} else {cat(" The",i,"Observation of 2nd Group is Misclassified") }
    }
  }
}



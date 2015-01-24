\name{wage.rates}
\docType{data}
\alias{wage.rates}
\title{ Wage Rates of 753 Women}
\description{ Wage rates of 753 married white women with left censoring 
}
\usage{data(wage.rates)}
\format{ 
A data frame with 753 observed wage rates of married white women in 1975. Some wage rates are set equal to zero, this means that those wives did not work in 1975, therefore, these observations are considered left censored at zero 
}

\encoding{latin1}

\source{ 
Mroz, T.A. 1987. "The sesitivity of an empirical model of married women´s hours of work to economic and statistical assumptions". Full data set is available at 
}

\examples{
##Load the data
data(wage.rates)

##Setting y, x and cc
y <- wage.rates$wage
x <- cbind(1,wage.rates$age,wage.rates$educ,wage.rates$kidslt6,wage.rates$kidsge6)
cc<- (wage.rates$inlf==0)+0
p <- ncol(x)
n <- nrow(x)

## Fits a left censored Normal model and calculates Generalized Cook Distance
Normal <- em.cens(cc,-x,-y,dist="Normal",diagnostic=TRUE,typediag=1)

#Graphics the diagnostic measures
GD <- Normal$measure$GD
plot(GD,main="Generalized Cook Distance")
abline(h=2*(p+1)/n,col="red",lty=2)

GDbeta <- Normal$measure$GDbeta
plot(GDbeta, main="Generalized Cook Distance - Influence in Beta")
abline(h=2*p/n,col="red",lty=2)

GDsigma2 <- Normal$measure$GDsigma2
plot(GDsigma2,main="Generalized Cook Distance - Influence in Sigma2")
abline(h=2/n,col="red",lty=2)


## Fits a left censored Student-t model, does not calculate any diagnostic measures
T <- em.cens(cc,-x,-y,nu=5,dist="T")

}
\keyword{datasets}




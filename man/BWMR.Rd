\name{BWMR}
\alias{BWMR}
\title{
  Inference the causality based on BWMR method
}
\description{
Estimate the causal effect between the exposure and the outcome.
}
\usage{
BWMR(gammahat, Gammahat, sigmaX, sigmaY)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gammahat}{
     SNP-exposure effects.
}
  \item{Gammahat}{
     SNP-outcome effects.
}
  \item{sigmaX}{
     Standard errors of SNP-exposure effects.
}
  \item{sigmaY}{
     Standard errors of SNP-outcome effects.
}
}
\details{
  \code{BWMR} obtain the causal effect based on summary statistics.
}
\value{
  \item{mu_beta}{
     Estimate of parameter \code{beta}.
}
  \item{se_beta}{
     Estimate of the standard error of parameter \code{beta}.
}
  \item{P_value}{
     P_value.
}
  \item{plot1}{
     Plot of Data with Standard Error Bar.
}
  \item{plot2}{
     Trace Plot of Logarithm of Approximate Data Likelihood.
}
  \item{plot3}{
     Estimate of Weight of Each Data Point.
}
  \item{plot4}{
     Plot of Weighted Data and Its Regression Result.
}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Jia Zhao
}
\seealso{
}
\examples{
  library(BWMR)
  library(ggplot2)

  data(ExampleData)

  MRres <- BWMR(ExampleData$beta.exposure, ExampleData$beta.outcome, ExampleData$se.exposure, ExampleData$se.outcome)

  beta.hat <- MRres$beta
  beta.se <- MRres$se_beta
  P_value <- MRres$P_value

  # Plot1: Plot of Data with Standard Error Bar
  plot1 <- MRres$plot1
  plot1

  # Plot2: Trace Plot of Logarithm of Approximate Data Likelihood
  plot2 <- MRres$plot2
  plot2

  # Plot3: Estimate of Weight of Each Data Point
  plot3 <- MRres$plot3
  plot3

  # Plot4: Plot of Weighted Data and Its Regression Result
  plot4 <- MRres$plot4
  plot4
}

\name{mccimm}
\alias{mccimm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Constructing Confidence Intervals for Moderated Mediation with 5,000,000 Monte Carlo Simulated Samples
}
\description{
Percentile and bias-corrected confidence intervals for moderating mediation created by 5,000,000 Monte Carlo Simulated samples. Can also create confidence intervals for a function.
}
\usage{
mccimm(estP="Ctemp1.txt", varcovP="Ctemp2.txt", NoSP, LstdZ="NA", LstdW="NA",
                   La1="NA", Lz1="NA", Lw1="NA", Lzw1="NA",
                   La2="NA", Lz2="NA", Lw2="NA", Lzw2="NA",
                   La3="NA", Lz3="NA", Lw3="NA", Lzw3="NA",
                   La4="NA", Lz4="NA", Lw4="NA", Lzw4="NA",
                   Sfunction="NULL")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{estP}{
File with estimated parameters
}
  \item{varcovP}{
File with estimated variance-covariance matrix of estimated parameters (TECH3 in Mplus).
}
  \item{NoSP}{
Number of defined additional parameters.
}
  \item{LstdZ}{
Location of the standard deviation of moderator Z in the additional parameters.
}
  \item{LstdW}{
Location of the standard deviation of moderator w in the additional parameters.
}
  \item{La1}{
Location of direct path a1 in the additional parameters.
}
  \item{Lz1}{
Location of moderating effect z1 in the additional parameters.
}
  \item{Lw1}{
Location of moderating effect w1 in the additional parameters.
}
  \item{Lzw1}{
Location of 3-way interaction path zw1 in the additional parameters.
}
  \item{La2}{
Location of direct path a2 in the additional parameters.
}
  \item{Lz2}{
Location of moderating effect z2 in the additional parameters.
}
  \item{Lw2}{
Location of moderating effect w2 in the additional parameters.
}
  \item{Lzw2}{
Location of 3-way interaction path zw2 in the additional parameters.
}
  \item{La3}{
Location of direct path a3 in the additional parameters.
}
  \item{Lz3}{
Location of moderating effect z3 in the additional parameters.
}
  \item{Lw3}{
Location of moderating effect w3 in the additional parameters.
}
  \item{Lzw3}{
Location of 3-way interaction path zw3 in the additional parameters.
}
  \item{La4}{
Location of direct path a4 in the additional parameters.
}
  \item{Lz4}{
Location of moderating effect z4 in the additional parameters.
}
  \item{Lw4}{
Location of moderating effect w4 in the additional parameters.
}
  \item{Lzw4}{
Location of 3-way interaction path zw4 in the additional parameters.
}
  \item{Sfunction}{
User defined function in quotations using a, z, w, zw.
}
}
\details{
	By specifying the location in the additional parameters for direct paths (a1, a2, a3, and a4); the moderating effects from moderator Z to the direct paths (z1, z2, z3, and z4); the moderating effects from moderator W to the direct paths (w1, w2, w3, and w4), the interaction effects between W and Z on the direct paths (three-way interactions zw1, zw2, zw3, and zw4), as well as the standard deviation  of Z (stdz) and/or standard deviation of W (stdw), confidence intervals and p-values of the mediating effects or moderated mediating effects can be estimated by Monte Carlo simulation.
	The above general model has covered many models in Hayes' templates. Specify the location of the parameters in the New/Additional Parameters defined in the NEW function in Mplus. Only need to specify the location of paths and standard deviation of moderator in the model (by adding an "L" before the parameter). Moderated direct paths can also specify as a separate model using the above general model.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{

}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

One Moderator Z
Model 1: La1, Lz1, LstdZ
Model 4: La1, La2
Model 5: La1, La2 (mediating effect) or La1, Lz1, LstdZ (moderated direct effect)
Model 6: La1, La2, La3
Model 7: La1, Lz1, La2, LstdZ
Model 8: La1, Lz1, La2, LstdZ (moderated mediating effect) or La1, Lz1, LstdZ (moderated direct effect)
Model 14: La1, La2, Lz2, LstdZ
Model 15: La1, La2, Lz2, LstdZ (moderated mediating effect) or La1, Lz1, LstdZ (moderated direct effect)
Model 58: La1, Lz1, La2, Lz2, LstdZ
Model 59: La1, Lz1, La2, Lz2, LstdZ (moderated mediating effect) or La1, Lz1, LstdZ (moderated direct effect)

Two Moderators Z & W
Model 2: La1, Lz1, Lw1, LstdZ, LstdW
Model 3: La1, Lz1, Lw1, Lzw1, LstdZ, LstdW
Model 9: La1, Lz1, Lw1, La2, LstdZ, LstdW
Model 10: La1, Lz1, Lw1, La2, LstdZ, LstdW (moderated mediating effect) or La1, Lz1, Lw1, LstdZ, LstdW (moderated direct effect)
Model 11: La1, Lz1, Lw1, Lzw1, La2, LstdZ, LstdW
Model 12: La1, Lz1, Lw1, Lzw1, La2, LstdZ, LstdW (moderated mediating effect) or La1, Lz1, Lw1, Lzw1, LstdZ, LstdW (moderated direct effect)
Model 13: La1, Lz1, Lw1, Lzw1, La2, LstdZ, LstdW (moderated mediating effect) or La1, Lz1, LstdZ (moderated direct effect)
Model 16: La1, La2, Lz2, Lw2, LstdZ, LstdW
Model 17: La1, La2, Lz2, Lw2, LstdZ, LstdW (moderated mediating effect) or La1, Lz1, Lw1, LstdZ, LstdW (moderated direct effect)
Model 18: La1, La2, Lz2, Lw2, Lzw2, LstdZ, LstdW
Model 19: La1, La2, Lz2, Lw2, Lzw2, LstdZ, LstdW (moderated mediating effect) or La1, Lz1, Lw1, Lzw1, LstdZ, LstdW (moderated direct effect)
Model 20: La1, La2, Lz2, Lw2, Lzw2, LstdZ, LstdW (moderated mediating effect) or La1, Lz1, LstdZ (moderated direct effect)
Model 21: La1, Lz1, La2, Lw2, LstdZ, LstdW
Model 22: La1, Lz1, La2, Lw2, LstdZ, LstdW (moderated mediating effect) or La1, Lz1, LstdZ (moderated direct effect)
Model 28: La1, Lz1, La2, Lw2, LstdZ, LstdW (moderated mediating effect) or La1, Lz1, LstdZ (moderated direct effect)
Model 29: La1, Lz1, La2; Lw2, LstdZ, LstdW (moderated mediating effect) or La1, Lz1, Lw1, LstdZ, LstdW (moderated direct effect)
Model 60: La1, Lz1, Lw1, La2, Lw2, LstdZ, LstdW
Model 61: La1, Lz1, Lw1, La2, Lw2, LstdZ, LstdW (moderated mediating effect) or La1, Lw1, LstdW (moderated direct effect)
Model 62: La1, Lz1, Lw1, La2, Lw2, LstdZ, LstdW (moderated mediating effect) or La1, Lz1, LstdZ (moderated direct effect)
Model 63: La1, Lz1, Lw1, La2, Lw2, LstdZ, LstdW (moderated mediating effect) or La1, Lw1, LstdW (moderated direct effect)
Model 64: La1, Lz1, La2, Lz2, Lw2, LstdZ, LstdW
Model 65: La1, Lz1, La2, Lz2, Lw2, LstdZ, LstdW (moderated mediating effect) or La1, Lz1, LstdZ (moderated direct effect)
Model 66: La1, Lz1, La2, Lz2, Lw2, LstdZ, LstdW (moderated mediating effect) or La1, Lw1, LstdW (moderated direct effect)
Model 67: La1, Lz1, La2, Lz2, Lw2, LstdZ, LstdW (moderated mediating effect) or La1, Lz1, Lw1, LstdZ, LstdW (moderated direct effect)
Model 68: La1, Lz1, Lw1, Lzw1, La2, Lz2, LstdZ, LstdW
Model 69: La1, Lz1, Lw1, Lzw1, La2, Lz2, LstdZ, LstdW (moderated mediating effect) or La1, Lz1, Lw1, Lzw1, LstdZ, LstdW (moderated direct effect)
Model 70: La1, Lz1, La2, Lz2, Lw2, Lzw2, LstdZ, LstdW
Model 71: La1, Lz1, La2, Lz2, Lw2, Lzw2, LstdZ, LstdW (moderated mediating effect) or La1, Lz1, Lw1, Lzw1, LstdZ, LstdW (moderated direct effect)
Model 72: La1, Lz1, Lw1, Lzw1, La2, Lz2, Lw2, Lzw2, LstdZ, LstdW
Model 73: La1, Lz1, Lw1, Lzw1, La2, Lz2, Lw2, Lzw2, LstdZ, LstdW (moderated mediating effect) or La1, Lz1, Lw1, Lzw1, LstdZ, LstdW (moderated direct effect)
Model 75: La1, Lz1, Lw1, La2, Lz2, Lw2, LstdZ, LstdW
Model 76: La1, Lz1, Lw1, La2, Lz2, Lw2, LstdZ, LstdW (moderated mediating effect) or La1, Lz1, Lw1, LstdZ, LstdW (moderated direct effect)

For single-level moderated mediation models, a1, a2, a3, and a4 are the direct effects whereas z1, z2, z3, z4 and w1, w2, w3, w4 are two-way interaction effects between the moderator and the predictor on a criterion.  Three-way interaction effects between the moderators on the effects from the predictor on a criterion are represented by zw1, zw2, zw3, and zw4.
When there are three-way interactions, the difference in simple slope tests are provided based on Dawson (2014).
For cross-level moderated mediation models, a1, a2, a3, and a4 are the mean of random slopes where z1, z2, z3, z4 and w1, w2, w3, w4 are the effects of the Level 2 moderators on the random slopes.
The general model actually goes beyond the above models because more than two mediating effects (maximum four) can be specified.  Please note that when there are two or more moderating effects, Index of Moderated Mediation is not defined (same for PROCESS).

EXAMPLE A: (Model 7 with 7 parameters in New/Additional Parameters)
mccimm(estP="Ctemp1.txt", varcovP="Ctemp2.txt", NoSP=7, La1=2, Lz1=3, La2=4, LstdZ=1)

EXAMPLE B: Generate Confidence Intervals with a function
mccimm(estP="Ctemp1.txt", varcovP="Ctemp2.txt", NoSP=4, La1=1, Lz1=2, Lw1=3, Lzw1=4, Sfunction="(a1+z1)*zw1")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.

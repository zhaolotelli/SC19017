#include <Rcpp.h>
using namespace Rcpp;

//' @title MetropolisC
//' @description A Metropolis sampler using Rcpp
//' @param sigma the standard deviation of normal distribution
//' @param x0 the first element of the sampler
//' @param N the number of samples
//' @return a random sample of size \code{n}
//' \item{x}{random walk samples}
//' \item{k}{rejections}
//' @examples
//' \dontrun{
//' rwC <- rw_MetropolisC(1,25,2000)$x
//' plot(rwC, type="l", xlab=bquote(sigma == 1), ylab="X", ylim=range(rwC))
//' }
//' @export
// [[Rcpp::export]]
List rw_MetropolisC(double sigma, double x0, int N) 
{
  //Metropolis Randomwalk using C
  NumericVector x(N);
  x[0] = x0;
  double u, y;
  int k = 0;
  for (int i = 1; i < N; i++) 
  {
    y = rnorm(1, x[i-1], sigma)[0];
    u = runif(1)[0];
    if (u <= exp(-(abs(y) - abs(x[i-1])))) 
    {
      x[i] = y; 
    }
    else 
    {
      x[i] = x[i-1];
      k++;
    }
  }
  return List::create(Named("x", x), Named("k", k));
}

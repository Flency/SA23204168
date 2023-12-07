#include <Rcpp.h>
using namespace Rcpp;

//' @title Compute the maximum l2 norm for each column of the matrix
//' @description Compute the maximum l2 norm for each column of the matrix
//' @param x the object matrix
//' @return the maximum l2 norm for each column of the matrix
//' @examples
//' \dontrun{
//' x <- matrix(c(1, 2, 4, 9), 2, 2)
//' print(lmax(x))
//' }
//' @export
// [[Rcpp::export]]
double lmax(NumericMatrix x){
  int nr = x.nrow();
  int nc = x.ncol();
  double y, y1;
  y = 0.0;
  NumericVector ycol(nr);
     for(int j=0; j<nc; j++){
        ycol = x.column(j);
        y1 = sum(ycol * ycol);
        if(y1 > y) y = y1;
    }
  y = sqrt(y);
  return y;
}

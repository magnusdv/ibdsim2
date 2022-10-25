#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix recombine(NumericMatrix strand1, NumericMatrix strand2, NumericVector cross) {
	int nc = cross.size();
  if(nc == 0) 
    return strand1;
  
  if(is_true(any(cross < 0))) stop("Inadmissible crossover point");

  int nr1 = strand1.nrow(), nr2 = strand2.nrow();
  NumericMatrix gamete(nr1 + nr2 + nc, 2);

  int ic, i1 = 0, i2 = 1, i3 = -1;
  int k = 0;
  double cx1, cx2;

  // Rprintf("\nStarting recombine()\n"); Rf_PrintValue(strand1); Rf_PrintValue(strand2);
  
  for(ic = 0; ic < nc - 1; ic += 2) {
    cx1 = cross[ic];
    cx2 = cross[ic+1];
    
    // Copy strand1 until before break point
    for(k = i1;  k < nr1 && strand1(k, 0) < cx1; k++)
      gamete(++i3, _) = strand1(k, _);

    i1 = k;
    
    // Skip down strand2 until cx1
    for(k = i2-1; k < nr2 && strand2(k,0) <= cx1; k++)  {}

    i2 = k;
    
    // Insert break point S1 --> S2
    gamete(++i3, 0) = cx1;
    gamete(i3, 1) = strand2(i2-1,1);

    // Copy strand2 until before break point
    for(k = i2;  k < nr2 && strand2(k,0) < cx2; k++)
      gamete(++i3, _) = strand2(k, _);

    i2 = k;
    
    // Skip down strand1 until cx2
    for(k = i1-1; k < nr1 && strand1(k,0) <= cx2; k++) {}

    i1 = k;
    
    // Insert break point S2 --> S1
    gamete(++i3, 0) = cx2;
    gamete(i3, 1) = strand1(i1-1,1);
    
    // Rprintf("\nCrossed at %f and %f. Result:\n", cx1, cx2); Rf_PrintValue(gamete);
  }
  
  if(ic == nc - 1) {
    cx1 = cross[ic];
    
    // Copy strand1 until before break point
    for(k = i1;  k < nr1 && strand1(k, 0) < cx1; k++)
      gamete(++i3, _) = strand1(k, _);

    i1 = k;
    
    // Skip down strand2 until cx1
    for(k = i2-1; k < nr2 && strand2(k,0) <= cx1; k++)  {}

    i2 = k;
    
    // Insert break point S1 --> S2
    gamete(++i3, 0) = cx1;
    gamete(i3, 1) = strand2(i2-1,1);

    // Copy strand2 until end
    for(k = i2;  k < nr2; k++)
      gamete(++i3, _) = strand2(k, _);
    
    // Rprintf("Final cross S1 --> S2 at %f. Result:\n", cx1); Rf_PrintValue(gamete);
  }
  else {
    // Copy strand1 until end
    for(k = i1;  k < nr1; k++)
      gamete(++i3, _) = strand1(k, _);
  }
  
  return gamete(Range(0, i3), _);
}


// [[Rcpp::export]]
NumericVector sort_dbl_C(NumericVector x) {
  NumericVector y(x);
  std::sort(y.begin(), y.end());
  return y;
} 


// [[Rcpp::export]]
double sample_12_C() {
  double r = R::rnorm(0,1);
  if(r < 0) {
    return(1);
  }
  else {
    return(2);
  }
}


// [[Rcpp::export]]
NumericVector sample_int_C(double n, double size) {
  return ceiling(Rcpp::runif(size, 0, n));
}


// [[Rcpp::export]]
NumericVector convert_pos_C(NumericVector pos, 
                            NumericVector mapFrom, 
                            NumericVector mapTo,
                            double extValue) {
  // Output vector
  NumericVector res(pos.length());
    
  double src, rate;
  int idx = 0;
  int mapN = mapFrom.length();
  
  for(int i=0; i < res.length(); ++i) {
    
    src = pos[i];
    
    // Values outside of map
    if(src < mapFrom[0]) {
      res[i] = 0;
    }
    else if(src > mapFrom[mapN-1]) {
      res[i] = extValue;
    }
    else {
      // findInterval(src, mapFrom). (Extract 1 because cpp is 0-indexed.)
      idx = std::distance(mapFrom.begin(), std::upper_bound(mapFrom.begin(), mapFrom.end(), src)) - 1;
      
      // Rprintf("\nidx = %d. mapto[idx+1] = %f. mapto[idx] = %f", idx, mapTo[idx + 1], mapTo[idx]);
      
      // Rate of change
      rate = (mapTo[idx + 1] - mapTo[idx]) / (mapFrom[idx + 1] - mapFrom[idx]);
      
      // Extrapolate
      res[i] = mapTo[idx] + (src - mapFrom[idx]) * rate;
    }
  }
  return res;
}


// [[Rcpp::export]]
NumericMatrix build_allelemat_C(NumericVector pos, List haplolist) { 
  // Each entry of haplolist: matrix with 2 columns (breaks - allele)
  
  int idx = 0;
  
  NumericMatrix am(pos.length(), haplolist.length());
  //if(is.null(haplo)) return(rep(0, length(posvec)))
  
  for(int j=0; j < haplolist.length(); ++j) {
    NumericMatrix h = haplolist[j]; 
    NumericVector breaks = h(_, 0);
    NumericVector als = h(_, 1);
    for(int i=0; i < pos.length(); ++i) {
      if(breaks.length() == 1) {
        idx = 0;
      }
      else {
        idx = std::distance(breaks.begin(), std::upper_bound(breaks.begin(), breaks.end(), pos[i])) - 1;
      }
      am(i, j) = als[idx];
    }
  }
  return am;
}
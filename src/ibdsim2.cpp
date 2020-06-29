#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix recombine(NumericMatrix strand1, NumericMatrix strand2, NumericVector cross) {
	int nc = cross.size();
  if(nc == 0) return strand1;
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



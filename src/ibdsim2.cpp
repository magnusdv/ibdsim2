#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix recombine(NumericMatrix strand1, NumericMatrix strand2, NumericVector cross) {
	int nc = cross.size();
  if(nc == 0) return strand1;
  if(is_true(any(cross < 0))) stop("Inadmissible crossover point");

  int nr1 = strand1.nrow(), nr2 = strand2.nrow();
  NumericMatrix gamete(nr1 + nr2 + nc, 2);

  int ic, i1=0, i2 = 1, i3 = -1;
  double cx1, cx2;

  for(ic = 0; ic < nc - 1; ic += 2) {
    cx1 = cross[ic];
    cx2 = cross[ic+1];
    //Rprintf("\nStarting loop involving crosses %f and %f.\nPresent result:\n", cx1, cx2); Rf_PrintValue(gamete);

    // Copy strand1 until before break point
    for(i1 = i1;  i1 < nr1 && strand1(i1, 0) < cx1; i1++)
      gamete(++i3, _) = strand1(i1, _);

    // Skip down strand2 until cx1
    for(i2 = i2-1; i2 < nr2 && strand2(i2,0) <= cx1; i2++)  {}

    // Insert break point S1 --> S2
    gamete(++i3, 0) = cx1;
    gamete(i3, 1) = strand2(i2-1,1);

    // Copy strand2 until before break point
    for(i2 = i2;  i2 < nr2 && strand2(i2,0) < cx2; i2++)
      gamete(++i3, _) = strand2(i2, _);

    // Skip down strand1 until cx2
    for(i1 = i1-1; i1 < nr1 && strand1(i1,0) <= cx2; i1++) {}

    // Insert break point S2 --> S1
    gamete(++i3, 0) = cx2;
    gamete(i3, 1) = strand1(i1-1,1);
  }
  if(ic == nc - 1) {
    cx1 = cross[ic];
    // Rprintf("Final cross S1 --> S2 at %f.\nPresent result:\n", cx1); Rf_PrintValue(gamete);

    // Copy strand1 until before break point
    for(i1 = i1;  i1 < nr1 && strand1(i1, 0) < cx1; i1++)
      gamete(++i3, _) = strand1(i1, _);

    // Skip down strand2 until cx1
    for(i2 = i2-1; i2 < nr2 && strand2(i2,0) <= cx1; i2++)  {}

    // Insert break point S1 --> S2
    gamete(++i3, 0) = cx1;
    gamete(i3, 1) = strand2(i2-1,1);

    // Copy strand2 until end
    for(i2 = i2;  i2 < nr2; i2++)
      gamete(++i3, _) = strand2(i2, _);
  }
  else {
    // Copy strand1 until end
    for(i1 = i1;  i1 < nr1; i1++)
      gamete(++i3, _) = strand1(i1, _);
  }
  return gamete(Range(0, i3), _);
}



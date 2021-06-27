//
//  Auxiliaries.cpp
//
//
//  Initially created by Arne Pommerening on 17.01.2013. Common functions outsourced by Henrike Haebel on 26.04.2018.
//  Copyright 2018 Philodendron International. All rights reserved.
//

#include <Rcpp.h>
#include <cmath>
#include "Auxiliaries.h"
using namespace Rcpp;
using namespace std;


/*
 * Calculates distance between trees with periodic boundary conditions (Illian et al., 2008, p. 184).
 */
// [[Rcpp::export]]
double getEuclideanDistance(double xmax, double ymax, double x1, double y1, double x2, double y2) {
  double dx = fabs(x1 - x2);
  double dy = fabs(y1 - y2);
  dx = min(dx, xmax - dx);
  dy = min(dy, ymax - dy);
  double dz = sqrt(dx * dx + dy * dy);
  return dz;
}


/**
 * Impulse function as in Illian et al. (2008), p. 435.
 */
// double impulse(NumericVector abdn, double distance, double y) {
//   return pow(y, abdn[0]) * exp(-abdn[2] * pow(distance, 2) / pow(y, abdn[1]));
  // return pow(y, abdn[0]) * exp(-pow(distance, 2) /  pow(y, abdn[1]));
  // return y * pow(abdn[0], 2) / PI * exp(-pow(abdn[0], 2) * pow(distance, 2));
//}
 
 double impulse(NumericVector abdn, double distance, double y) {
   double g = pow(y, abdn[0]);
   if(distance > (y / 2))
     g *= exp(-abdn[2] * pow(distance, 2) / pow(y, abdn[1]));
   return g;
 }

 
 /**
 * Critical cumulative 5-years RGR function.
 */
// [[Rcpp::export]]
NumericVector pcrit(NumericVector m, NumericVector par) {
  return par[0] * pow(m, par[1]);
//  return par[0] * exp(par[1] * m);
}



//
//  GrowthInteraction.cpp 
//
//
//  This is code related to interaction and growth integrating both models, GI and SN.
//  Synthesised by Henrike Haebel on 24/04/2018 based on 
//  Arne Pommerening on 11/03/2017. Last modified by him on 07/12/20.
//  Copyright 2018 Philodendron International. All rights reserved.
//

#include <Rcpp.h>
#include <cmath>
#include "Auxiliaries.h"
#include <string>
using namespace Rcpp;
using namespace std;



/********************************************************************************************/
/* Interaction                                                                              */
/********************************************************************************************/


/*
 * Calculates interaction for SN model simulation.
 */
Rcpp::NumericVector calcTreeInteractionSN(NumericVector abdn, DataFrame xData, double xmax, double ymax, double super, bool superO) {
  vector<double> x = as< vector<double> > (xData["x"]);
  vector<double> y = as< vector<double> > (xData["y"]);
  vector<double> mark = as< vector<double> > (xData["mark"]);
  int n = x.size();
  NumericVector ci(n);
  for (int i = 0; i < n; i++) {
    ci[i] = 0; 
    for (int j = n - 1; j > -1; j--) {
      if (i != j) {
        double distance = getEuclideanDistance(xmax, ymax, x[i], y[i], x[j], y[j]);
        if(superO == true) {
          double superRule = (mark[i] / 2 + mark[j] / 2) / distance;
          if(superRule < super) 
            ci[i] += impulse(abdn, distance, mark[j]) / pow(mark[i], abdn[0]); 
        }
        else
          ci[i] += impulse(abdn, distance, mark[j]) / pow(mark[i], abdn[0]); 
        // else Rcout << "superRule: " << superRule << std::endl;
      }
    }
    ci[i] = exp(-abdn[3] / ci[i]); // /= pow(dbh[i], abdn[0]) + ci[i];
  }
  return 1 - ci;
}

/*
 * Calculates SN field values for a grid.
 */
// [[Rcpp::export]]
  NumericVector calcShotNoiseField(NumericVector abdn, NumericVector x, NumericVector y, DataFrame xData, double xmax, double ymax) {
  vector<double> xx = as< vector<double> > (xData["x"]);
  vector<double> yy = as< vector<double> > (xData["y"]);
  vector<double> mark = as< vector<double> > (xData["mark"]);
  int n = x.size();
  int m = xx.size();
  NumericVector ci(n);
  for (int i = 0; i < n; i++) {
    ci[i] = 0;
    for (int j = 0; j < m; j++) {
      double distance = getEuclideanDistance(xmax, ymax, x[i], y[i], xx[j], yy[j]);
      ci[i] += impulse(abdn, distance, mark[j]);
    }
  }
  return ci;
}

/*
 * Calculates SN field values for regeneration without superorganisms.
 */
// [[Rcpp::export]]
NumericVector calcShotNoiseFieldForRegen(NumericVector abdn, DataFrame xData, double xmax, double ymax) {
  vector<double> x = as< vector<double> > (xData["x"]);
  vector<double> y = as< vector<double> > (xData["y"]);
  vector<double> mark = as< vector<double> > (xData["mark"]);
  int n = x.size();
  NumericVector ci(n);
  for (int i = 0; i < n; i++) {
    ci[i] = 0;
    for (int j = n - 1; j > -1; j--) {
      if (i != j) {
        double distance = getEuclideanDistance(xmax, ymax, x[i], y[i], x[j], y[j]);
        ci[i] += impulse(abdn, distance, mark[j]);
      }
    }
  }
  return ci;
}

/*
 * Calculates interaction for SN regression (parameter estimation).
 */
// [[Rcpp::export]]
NumericVector estTreeInteraction(NumericVector abdn, DataFrame xData, double super, bool superO) {
  vector<double> x = as< vector<double> > (xData["x"]);
  vector<double> y = as< vector<double> > (xData["y"]);
  vector<double> mark = as< vector<double> > (xData["cw"]);
  vector<double> year = as< vector<double> > (xData["year"]);
  vector<int> plotno = as< vector<int> > (xData["plotno"]);
  vector<double> xmax = as< vector<double> > (xData["xmax"]);
  vector<double> ymax = as< vector<double> > (xData["ymax"]);
  int n = mark.size();
  NumericVector ci(n);
  for (int i = 0; i < n; i++) {
    ci[i] = 0; 
    for (int j = n - 1; j > -1; j--) {
      if ((i != j) && (year[i] == year[j]) && (plotno[i] == plotno[j])) {
        double distance = getEuclideanDistance(xmax[i], ymax[i], x[i], y[i], x[j], y[j]);
        if(superO == true) {
          double superRule = (mark[i] / 2 + mark[j] / 2) / distance;
          if(superRule < super) 
            ci[i] += impulse(abdn, distance, mark[j]) / pow(mark[i], abdn[0]); // +=
        }
        else
          ci[i] += impulse(abdn, distance, mark[j]) / pow(mark[i], abdn[0]); // +=
      }
    }
    ci[i] = exp(-abdn[3] / ci[i]); 
    // ci[i] = exp(-abdn[3] * ci[i]); 
    // ci[i] = pow(exp(ci[i]), -abdn[3]); 
    // ci[i] = ci[i] * abdn[3]; 
    // ci[i] = exp(1 / (-abdn[3] * ci[i])); 
  }
  return 1 - ci;
}




/*********************************************************************************************/
/* Relative growth rate RGR                                                                 */
/********************************************************************************************/

/**
 * Annual RGR model simulation.
 */
// [[Rcpp::export]]
  NumericVector simulateRelativeGrowthRate(DataFrame xData, NumericVector abdn, double xmax, double ymax, NumericVector potGrowth, double super, bool superO) {
  NumericVector y = as< NumericVector > (xData["mark"]);
  int n = y.size();
  NumericVector g(n);
    // g = estimateGrowthPotential(y, param) * abdn[3] * calcTreeInteraction(abdn, xData, xmax, ymax, minteraction);
    g = potGrowth * calcTreeInteractionSN(abdn, xData, xmax, ymax, super, superO);
    // Rcout << "Selected routine: multiplicative" <<  " Flag: " << mgrowth.substr(0, 5) << std::endl;
  return g;
}

/**
 * Annual RGR regression.
 */
// [[Rcpp::export]]
NumericVector estimateRelativeGrowthRate(DataFrame xData, NumericVector param, NumericVector abdn, double super, bool superO) {
  NumericVector y = as< NumericVector > (xData["cw"]);
  int n = y.size();
  NumericVector g(n);
    // g = estimateGrowthPotential(y, param) * abdn[3] * estTreeInteraction(abdn, xData, minteraction);
    g = pcrit(y, param) *  estTreeInteraction(abdn, xData, super, superO);
    // g = estimateGrowthPotential(y, param) * abdn[2] * estTreeInteraction(abdn, xData, minteraction);
    //Rcout << "Selected routine: multiplicative" <<  " Flag: " << mgrowth.substr(0, 5) << std::endl;
  return g;
}  



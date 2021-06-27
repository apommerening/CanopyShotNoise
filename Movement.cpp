//
//  Movement.cpp
//  
//
//  Created by Arne Pommerening on 26/10/2020. Modified on 07/12/2020.
//  Copyright 2020 Philodendron International. All rights reserved.
//
//

#include <Rcpp.h>
#include <cmath>
#include "Auxiliaries.h"
#include <vector>
#include <random>
// [[Rcpp::plugins("cpp11")]]

using namespace Rcpp;
using namespace std;

int seed = 1; //std::chrono::system_clock::now().time_since_epoch().count();
static std::default_random_engine gen(seed);
static std::uniform_real_distribution<double> u(0, 1);


double calcEuclideanDistance(double xmax, double ymax, double x1, double y1, double x2, double y2) {
  double dx = fabs(x1 - x2);
  double dy = fabs(y1 - y2);
  // dx = min(dx, xmax - dx);
  // dy = min(dy, ymax - dy);
  double dz = sqrt(dx * dx + dy * dy);
  return dz;
}

int findIndexOfSmallestElement(vector<double> vec) {
    int size = vec.size();
        int index = 0; // vec[0];
        double myValue = 999;
        for(int i = 0; i < size; i++) {
            if((vec[i] < myValue) & (vec[i] > 0)) {
                index = i;
                myValue = vec[i];
            }
        }
    return index;
}

/*
 * Calculates interaction for SN model simulation for one tree.
 */
double calcTreeInteractionOne(NumericVector abdn, NumericVector x, NumericVector y, NumericVector mark, double xmax, double ymax, int index, double super, bool superO) {
  int n = x.size();
  double ci;
  ci = 0;
  for(int j = n - 1; j > -1; j--) {
    if (index != j) {
      double distance = getEuclideanDistance(xmax, ymax, x[index], y[index], x[j], y[j]);
      if(superO == true) {
        double superRule = (mark[index] / 2.0 + mark[j] / 2.0) / distance;
        if(superRule < super)
          ci += impulse(abdn, distance, mark[j]); // / pow(mark[index], abdn[0]);
      }
      else
        ci += impulse(abdn, distance, mark[j]); // +=
    }
  }
  return ci;
}


/**
 * Determines the movement of tree crowns in one year.
 */
// [[Rcpp::export]]
DataFrame moveCrowns(int year, double maxMov, int maxAttempts, NumericVector x, NumericVector y, NumericVector mark, double xmax, double ymax, 
                     NumericVector param, double super, bool superO){
                     // NumericVector param, double super, bool superO, double crossArea, double mort, bool silent){
  int n = x.size();
  // NumericVector dummyX = clone(x); // Copy of original coordinates
  // NumericVector dummyY = clone(y);
  // NumericVector xx = clone(x); // New coordinates
  // NumericVector yy = clone(y);
  for(int j = 0; j < n; j++) {
    // if(!silent)
    //   Rcout << "Year: " << year << " Tree #: " << j + 1 << " n: " << n << " crossArea: " << crossArea << " superRule: " << super << " mortThresh: " << mort << std::endl;
    vector<double> ciList(maxAttempts + 1, 999), xList(maxAttempts + 1, 999), yList(maxAttempts + 1, 999);
    xList[0] = x[j];
    yList[0] = y[j];
    ciList[0] = calcTreeInteractionOne(param, x, y, mark, xmax, ymax, j, super, superO);
    NumericVector dummyX = clone(x);
    NumericVector dummyY = clone(y);
    for(int k = 1; k < (maxAttempts + 1); k++) {
      double dist = 999;
      double ux = 0, uy = 0;
      while(dist > maxMov) {
        // ux = (u(gen) * (x[j] + maxMov - x[j])) + x[j];
        // uy = (u(gen) * (y[j] + maxMov - y[j])) + y[j];
        ux = u(gen) * maxMov + x[j];
        uy = u(gen) * maxMov + y[j];
        dist = calcEuclideanDistance(xmax, ymax, x[j], y[j], ux, uy);
        // double a = u(gen);
        // Rcout << "(int)(u(gen): " << a << " ux: " << ux << " uy: " << uy << " x[j]: " << x[j] << " y[j]: " << y[j] << " dist: " << dist << std::endl;
      }
      if(ux > xmax)
        ux = ux - xmax;
      if(uy > ymax)
        uy = uy - ymax;
      xList[k] = ux;
      yList[k] = uy;
      x[j] = ux;
      y[j] = uy;
      ciList[k] = calcTreeInteractionOne(param, x, y, mark, xmax, ymax, j, super, superO);
      // Rcout << "Year: " << year << " Tree #: " << j << " ciList[k]: " << ciList[k] << " dist: " << dist << " param[1]: " << param[1]  << std::endl;
      x = clone(dummyX);
      y = clone(dummyY);
    }
    int z = findIndexOfSmallestElement(ciList); 
    // Rcout << "Year: " << year << " Tree #: " << j << " z: " << z << std::endl;
    x[j] = xList[z];
    y[j] = yList[z];
    // xx[j] = xList[z];
    // yy[j] = yList[z];
  }
  return DataFrame::create( Named("x") = x, Named ("y") = y);
  // return DataFrame::create( Named("x") = xx, Named ("y") = yy);
}

  

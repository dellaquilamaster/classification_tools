#ifndef CLASSIFICATION_TOOLS_HIDDEN_H
#define CLASSIFICATION_TOOLS_HIDDEN_H

//
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <map>
#include <random>
//
#include <TCanvas.h>
#include <TPad.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TRandom3.h>
//

//
typedef struct classification_tools_pattern_ {
  std::vector <double> xy_values;
  double weight;
} classification_tools_pattern;
//

//
typedef struct classification_tools_data_ {
  unsigned int num_patterns;
  unsigned int num_x;
  unsigned int num_y;
  std::vector <std::string> x_names;
  std::vector <std::string> y_names;
  std::vector <double> x_costs;
  std::vector <double> y_weights;
  std::map<std::string, classification_tools_pattern> patterns;
  void (* bp_function)(double *, double *);
  //
  bool is_pattern_weights;
  bool is_y_weights;
  bool is_x_costs;
} classification_tools_data;
//

#endif

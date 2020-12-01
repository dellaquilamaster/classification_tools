#include <string>
#include "classification_tools.h"

void calcola(double * features, double * output) {
  output[0]=(0.000103439)*features[0]*features[1];
}

int main (int argc, char ** argv)
{
  std::string data_file("/media/dellaquila/Seagate_SAS/Enrico/NHANES/202003_dysglycemia/dataset06.tst");
  std::string xnames_file("/media/dellaquila/Seagate_SAS/Enrico/NHANES/202003_dysglycemia/dataset05.xnames");
  std::string ynames_file("/media/dellaquila/Seagate_SAS/Enrico/NHANES/202003_dysglycemia/dataset05.ynames");
  std::string weight_file("/media/dellaquila/Seagate_SAS/Enrico/NHANES/202003_dysglycemia/dataset06.wtst");
  
  //
  void * data = create_classification_data(data_file.c_str(), xnames_file.c_str(), ynames_file.c_str(), 2, 1);
  add_pattern_weights(data, weight_file.c_str());
  set_function(data, calcola);
  //
  
  //
  double threshold=0.5;
  //
  
  //
  draw_roc(data,"try.png");
  draw_confusion_matrix(data,&threshold,"try1.png");
  draw_fuzzy_values(data,"32","try2.png");
  //
}

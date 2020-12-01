#include "classification_tools.h"
#include "classification_tools_hidden.h"

//____________________________________________________
void * create_classification_data (const char * file_name, const char * file_x_names, const char * file_y_names, unsigned int num_x, unsigned int num_y, const char * unique_id_file)
{
  //
  std::ifstream file_data(file_name);
  if(!file_data.is_open()) {
    return 0;
  }
  //
  //
  std::ifstream file_data_xnames(file_x_names);
  if(!file_data_xnames.is_open()) {
    return 0;
  }
  //
  //
  std::ifstream file_data_ynames(file_y_names);
  if(!file_data_ynames.is_open()) {
    return 0;
  }
  //
  
  //
  std::vector <std::string> unique_id;
  //
  
  //reading unique id file
  if(unique_id_file)
  {
    //
    std::ifstream file_data_unique_id(unique_id_file);
    if(!file_data_unique_id.is_open()) {
      return 0;
    }
    //
    
    //
    while (!file_data_unique_id.eof())
    {     
      std::string LineRead;
      std::getline(file_data_unique_id, LineRead);
            
      if(LineRead.empty()) continue;
      if(LineRead.find_first_not_of(' ') == std::string::npos) continue;
      
      std::istringstream LineStream(LineRead);
      
      std::string unique_id_data;
      
      LineStream>>unique_id_data;
      
      unique_id.push_back(unique_id_data);
    }
    //
  }
  //end of reading unique id file
  //
  
  //
  //creation of data
  classification_tools_data * pdata = new classification_tools_data;
  //
  pdata->is_pattern_weights=false;
  pdata->is_y_weights=false;
  pdata->is_x_costs=false;
  //
  pdata->num_patterns=0;
  pdata->num_x=num_x;
  pdata->num_y=num_y;
  //
  
  //
  unsigned int file_line=0;
  //
  
  //
  while (!file_data.eof())
  {
    //
    file_line++;
    //
    
    std::string LineRead;
    std::getline(file_data, LineRead);

    if(LineRead.empty()) continue;
    if(LineRead.find_first_not_of(' ') == std::string::npos) continue;

    std::istringstream LineStream(LineRead);
    
    std::string this_unique_id(unique_id.size()==0 ? Form("%u", file_line) : unique_id[file_line]);
    std::vector<double> * pxy_values = &(pdata->patterns[this_unique_id].xy_values);
    pdata->patterns[this_unique_id].weight=1;
    
    for(unsigned int i_feature = 0; i_feature<pdata->num_x; i_feature++)
    {
      double data_read;
      LineStream>>data_read;
      //
      pxy_values->push_back(data_read);
      //
    } 
    for(unsigned int i_output = 0; i_output<pdata->num_y; i_output++)
    {
      double data_read;
      LineStream>>data_read;
      //
      pxy_values->push_back(data_read);
      //      
    }
    pdata->num_patterns++;
  }
  //
  
  //reading x_names
  //
  while (!file_data_xnames.eof())
  {
    //
    file_line++;
    //
    
    std::string LineRead;
    std::getline(file_data_xnames, LineRead);

    if(LineRead.empty()) continue;
    if(LineRead.find_first_not_of(' ') == std::string::npos) continue;

    std::istringstream LineStream(LineRead);
    
    std::string name_data;
    
    LineStream>>name_data;
    
    pdata->x_names.push_back(name_data);
    pdata->x_costs.push_back(1);
  }
  //end of reading x_names
  
  //reading y_names
  //
  while (!file_data_ynames.eof())
  {
    //
    file_line++;
    //
    
    std::string LineRead;
    std::getline(file_data_ynames, LineRead);

    if(LineRead.empty()) continue;
    if(LineRead.find_first_not_of(' ') == std::string::npos) continue;

    std::istringstream LineStream(LineRead);
    
    std::string name_data;
    
    LineStream>>name_data;
    
    pdata->y_names.push_back(name_data);
    pdata->y_weights.push_back(1);
  }
  //end of reading x_names
  
  //
  return (void *)pdata;
  //
}

//____________________________________________________
int add_pattern_weights(void * pdata, const char * weights_file)
{
  std::ifstream file_data_weights(weights_file);
  if(!file_data_weights.is_open()) {
    return -1;
  }
  //
  
  std::vector <double> pattern_weights;
  
  int n_read=0;
  
  //
  while (!file_data_weights.eof())
  {     
    std::string LineRead;
    std::getline(file_data_weights, LineRead);
    
    if(LineRead.empty()) continue;
    if(LineRead.find_first_not_of(' ') == std::string::npos) continue;
    
    std::istringstream LineStream(LineRead);
    
    double weight_data;
    
    LineStream>>weight_data;
    
    pattern_weights.push_back(weight_data);
    
    n_read++;
  }
  //
  
  //
  int i_weight=0;
  for(auto && a_pattern : ((classification_tools_data *)pdata)->patterns)
  {
    a_pattern.second.weight=pattern_weights[i_weight];
    i_weight++;
  }
  //
  
  //
  ((classification_tools_data *)pdata)->is_pattern_weights=true;
  //

  return n_read;
}  

//____________________________________________________
int add_feature_costs(void * pdata, const char * costs_file)
{
  std::ifstream file_data_costs(costs_file);
  if(!file_data_costs.is_open()) {
    return -1;
  }
  //
  
  std::vector <double> x_costs;
  
  int n_read=0;

  //
  while (!file_data_costs.eof())
  {     
    std::string LineRead;
    std::getline(file_data_costs, LineRead);
    
    if(LineRead.empty()) continue;
    if(LineRead.find_first_not_of(' ') == std::string::npos) continue;
    
    std::istringstream LineStream(LineRead);
    
    double cost_data;
    
    LineStream>>cost_data;
    
    x_costs.push_back(cost_data);
    
    n_read++;
  }
  //
  
  //
  int i_cost=0;
  for(auto && a_feature_cost : ((classification_tools_data *)pdata)->x_costs)
  {
    a_feature_cost=x_costs[i_cost];
    i_cost++;
  }
  //
  
  //
  ((classification_tools_data *)pdata)->is_x_costs=true;
  //

  return n_read;
}  

//____________________________________________________
int add_output_weights(void * pdata, const char * output_weights_file)
{
  //
  std::ifstream file_data_output_weights(output_weights_file);
  if(!file_data_output_weights.is_open()) {
    return -1;
  }
  //
  
  std::vector <double> y_weights;
  
  int n_read=0;

  //
  while (!file_data_output_weights.eof())
  {     
    std::string LineRead;
    std::getline(file_data_output_weights, LineRead);
    
    if(LineRead.empty()) continue;
    if(LineRead.find_first_not_of(' ') == std::string::npos) continue;
    
    std::istringstream LineStream(LineRead);
    
    double output_weight_data;
    
    LineStream>>output_weight_data;
    
    y_weights.push_back(output_weight_data);
    
    n_read++;
  }
  //
  
  //
  int i_output_weight=0;
  for(auto && a_y_weight : ((classification_tools_data *)pdata)->y_weights)
  {
    a_y_weight=y_weights[i_output_weight];
    i_output_weight++;
  }
  //
  
  //
  ((classification_tools_data *)pdata)->is_y_weights=true;
  //
  
  return n_read;
}

//____________________________________________________
void set_function(void * pdata, void (*function)(double *, double *))
{
  ((classification_tools_data *)pdata)->bp_function=function;
  
  return;
}

//____________________________________________________
double * get_pattern(void * pdata, const char * unique_id)
{
  return ((classification_tools_data *)pdata)->patterns[unique_id].xy_values.data();
}

//____________________________________________________
void draw_roc(void * pdata, const char * file_output_name, double min_threshold, double max_threshold, double threshold_step)
{
  //
  std::string file_image_prefix(file_output_name);
  file_image_prefix.assign(file_image_prefix.substr(0,file_image_prefix.find_last_of('.')));
  std::string file_image_extension(file_output_name);
  file_image_extension.assign(file_image_extension.substr(file_image_extension.find_last_of('.')));
  //
  
  //
  const unsigned int num_points=(max_threshold-min_threshold)/threshold_step;
  //
  const unsigned int dimensionality_output=((classification_tools_data *)pdata)->num_y;
  //
  double * calculated_output = new double[dimensionality_output];
  //
  TCanvas the_canvas("the_canvas","",800,600);
  TGraph roc_graph;
  //
  //loop on classifier output
  for(unsigned int i_output=0; i_output<dimensionality_output; i_output++)
  {
    std::string file_image_name(Form("%s_%02d%s", file_image_prefix.c_str(), i_output, file_image_extension.c_str()));
    //loop on threshold values
    for(unsigned int i_point=0; i_point<num_points; i_point++) 
    {
      double num_positives=0;
      double num_negatives=0;
      double tp=0;
      double tn=0;
      
      const double threshold=min_threshold+i_point*threshold_step;
      
      //
      //loop on patterns
      for(auto && a_pattern : ((classification_tools_data *)pdata)->patterns) 
      {
        bool expected_output=a_pattern.second.xy_values[((classification_tools_data *)pdata)->num_x+i_output];
        double the_weight=a_pattern.second.weight;
        //
        num_positives+=expected_output*the_weight; 
        num_negatives+=(!expected_output)*the_weight;
        //
        //evaluating model output
        ((classification_tools_data *)pdata)->bp_function(a_pattern.second.xy_values.data(),calculated_output);
        //
        
        //
        //making a prediction
        bool predicted_output=false;
        if(calculated_output[i_output]>=threshold) predicted_output=true;
        //
        
        if(predicted_output==expected_output && expected_output) tp+=the_weight;
        if(predicted_output==expected_output && !expected_output) tn+=the_weight;
      }
      //end of loop on patterns
      //
      //adding point to graph
      const double sensitivity=tp/num_positives;
      const double specificity=tn/num_negatives;      
      roc_graph.SetPoint(roc_graph.GetN(),1-specificity,sensitivity);
      //
    }
    // end of loop on threshold values
    //
    roc_graph.Draw("AL");
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.14);
    gPad->SetBottomMargin(0.12);
    gPad->SetFrameLineWidth(2);
    roc_graph.SetLineWidth(3);
    roc_graph.SetLineColor(kBlack);
    roc_graph.GetXaxis()->SetLabelSize(0.05);
    roc_graph.GetXaxis()->SetTitleSize(0.05);
    roc_graph.GetYaxis()->SetLabelSize(0.05);
    roc_graph.GetYaxis()->SetTitleSize(0.05);
    roc_graph.GetXaxis()->CenterTitle();
    roc_graph.GetYaxis()->CenterTitle();
    roc_graph.GetXaxis()->SetRangeUser(0,1);
    roc_graph.GetYaxis()->SetRangeUser(0,1);
    roc_graph.SetTitle(((classification_tools_data *)pdata)->y_names[i_output].c_str());
    roc_graph.GetYaxis()->SetNdivisions(9,4,0);
    roc_graph.GetXaxis()->SetNdivisions(9,4,0);
    roc_graph.GetYaxis()->SetTitle("sensitivity");
    roc_graph.GetXaxis()->SetTitle("1 - specificity");
    the_canvas.Print(file_image_name.c_str());
    //
  }
  //end of loop on classifier output
  
  //
  delete [] calculated_output;
  //
  
  return;
}

//____________________________________________________
void draw_confusion_matrix (void * pdata, double * thr, const char * file_output_name)
{
  //
  const unsigned int dimensionality_output=((classification_tools_data *)pdata)->num_y;
  const unsigned int dimensionality_input=((classification_tools_data *)pdata)->num_x;
  double * calculated_output = new double[dimensionality_output];
  //
  std::vector<const char *> labels;
  for(unsigned int i_class=0; i_class<dimensionality_output; i_class++)
  {
    labels.push_back(((classification_tools_data *)pdata)->y_names[i_class].c_str());
  }
  //
  
  //
  TCanvas the_canvas("the_canvas","",800,600);
  TH2D histogram ("histogram","Confusion matrix",dimensionality_output,-0.5,dimensionality_output-0.5,dimensionality_output,-0.5,dimensionality_output-0.5);
  histogram.SetCanExtend(TH1::kAllAxes);
  histogram.SetStats(0);
  //
  
  //
  //loop on classifier output
  for(unsigned int i_output=0; i_output<dimensionality_output; i_output++)
  {
    //
    double patterns_in_class=0;
    double * patterns_in_predicted_classes = new double[dimensionality_output]();
    //
    //loop on patterns
    for(auto && a_pattern : ((classification_tools_data *)pdata)->patterns) 
    {
      if(a_pattern.second.xy_values[i_output+dimensionality_input]) {
        //
        patterns_in_class+=a_pattern.second.weight;
        //
        //evaluating model output
        ((classification_tools_data *)pdata)->bp_function(a_pattern.second.xy_values.data(),calculated_output);
        //
        //making a prediction
        for(unsigned int i_predicted_output=0; i_predicted_output<dimensionality_output; i_predicted_output++)
        {
          if(calculated_output[i_predicted_output]>=thr[i_predicted_output]) 
          {
            patterns_in_predicted_classes[i_predicted_output]+=a_pattern.second.weight;
          }
        }
        //
      }
    }
    //end of loop on patterns
    //
    for(unsigned int i_predicted_output=0; i_predicted_output<dimensionality_output; i_predicted_output++)
    {
      histogram.Fill(labels[i_predicted_output],labels[i_output],patterns_in_predicted_classes[i_predicted_output]/patterns_in_class);
    }
    delete [] patterns_in_predicted_classes;
    //
  }
  //end of loop on classifier output
  
  //
  histogram.LabelsDeflate("X");
  histogram.LabelsDeflate("Y");
  histogram.LabelsOption("v");
  histogram.Draw("text");
  the_canvas.SetGrid();
  the_canvas.SetLeftMargin(0.15);
  the_canvas.SetBottomMargin(0.15);
  the_canvas.Print(file_output_name);
  //  
  
  //
  delete [] calculated_output;
  //
  
  return;
}

//____________________________________________________
void draw_fuzzy_values (void * pdata, const char * unique_id, const char * file_output_name)
{
  //
  const unsigned int dimensionality_output=((classification_tools_data *)pdata)->num_y;
  double * calculated_output = new double[dimensionality_output];
  //
  std::vector<const char *> labels;
  for(unsigned int i_class=0; i_class<dimensionality_output; i_class++)
  {
    labels.push_back(((classification_tools_data *)pdata)->y_names[i_class].c_str());
  }
  //
  
  //
  TCanvas the_canvas("the_canvas","",800,600);
  TH1D histogram ("histogram","Fuzzy values",dimensionality_output,-0.5,dimensionality_output-0.5);
  histogram.SetCanExtend(TH1::kAllAxes);
  histogram.SetStats(0);
  //
  
  //
  //evaluating model output
  ((classification_tools_data *)pdata)->bp_function(((classification_tools_data *)pdata)->patterns[unique_id].xy_values.data(),calculated_output);
  //
  //making a prediction
  for(unsigned int i_predicted_output=0; i_predicted_output<dimensionality_output; i_predicted_output++)
  {
    histogram.Fill(labels[i_predicted_output],calculated_output[i_predicted_output]);
  }
  //
  
  //
  histogram.LabelsDeflate("X");
  histogram.LabelsOption("v");
  histogram.GetYaxis()->SetTitle("output");
  histogram.GetYaxis()->CenterTitle();
  histogram.Draw("text hist");
  histogram.SetLineColor(kBlack);
  histogram.SetLineWidth(2);
  the_canvas.SetGrid();
  the_canvas.SetLeftMargin(0.15);
  the_canvas.SetBottomMargin(0.15);
  the_canvas.Print(file_output_name);
  //
  
  //
  delete [] calculated_output;
  //
  
  return; 
}

//____________________________________________________
void draw_fuzzy_values (void * pdata, const char * unique_id, double * thr, const char * file_output_name)
{
  //
  const unsigned int dimensionality_output=((classification_tools_data *)pdata)->num_y;
  double * calculated_output = new double[dimensionality_output];
  //
  std::vector<const char *> labels;
  for(unsigned int i_class=0; i_class<dimensionality_output; i_class++)
  {
    labels.push_back(((classification_tools_data *)pdata)->y_names[i_class].c_str());
  }
  //
  
  //
  TCanvas the_canvas("the_canvas","",800,600);
  TH1D histogram ("histogram","Fuzzy values",dimensionality_output,-0.5,dimensionality_output-0.5);
  histogram.SetCanExtend(TH1::kAllAxes);
  histogram.SetStats(0);
  //
  
  //
  //evaluating model output
  ((classification_tools_data *)pdata)->bp_function(((classification_tools_data *)pdata)->patterns[unique_id].xy_values.data(),calculated_output);
  //
  //making a prediction
  for(unsigned int i_predicted_output=0; i_predicted_output<dimensionality_output; i_predicted_output++)
  {
    if(calculated_output[i_predicted_output]<thr[i_predicted_output]) continue;
    histogram.Fill(labels[i_predicted_output],calculated_output[i_predicted_output]);
  }
  //
  
  //
  histogram.LabelsDeflate("X");
  histogram.LabelsOption("v");
  histogram.GetYaxis()->SetTitle("output");
  histogram.GetYaxis()->CenterTitle();
  histogram.Draw("text hist bar2");
  histogram.SetFillColor(kGreen+1);
  histogram.SetLineWidth(2);
  the_canvas.SetGrid();
  the_canvas.SetLeftMargin(0.15);
  the_canvas.SetBottomMargin(0.15);
  the_canvas.Print(file_output_name);
  //
  
  //
  delete [] calculated_output;
  //
  
  return; 
}

//____________________________________________________
void calculate_auc (void * pdata, double * result, double* ci_min, double * ci_max)
{
  //
  const double zcrit=1.959964;
  const double min_threshold=0;
  const double max_threshold=1;
  const double threshold_step=0.001;
  unsigned int num_positives_nominal=0;
  unsigned int num_negatives_nominal=0;  
  //
  const unsigned int num_points=(max_threshold-min_threshold)/threshold_step;
  //
  const unsigned int dimensionality_output=((classification_tools_data *)pdata)->num_y;
  //
  double * calculated_output = new double[dimensionality_output];
  //
  TGraph roc_graph;
  //
  //loop on classifier output
  for(unsigned int i_output=0; i_output<dimensionality_output; i_output++)
  {
    //loop on threshold values
    for(unsigned int i_point=0; i_point<num_points; i_point++) 
    {
      double num_positives=0;
      double num_negatives=0;
      double tp=0;
      double tn=0;
      
      const double threshold=min_threshold+i_point*threshold_step;
      
      //
      //loop on patterns
      for(auto && a_pattern : ((classification_tools_data *)pdata)->patterns) 
      {
        bool expected_output=a_pattern.second.xy_values[((classification_tools_data *)pdata)->num_x+i_output];
        double the_weight=a_pattern.second.weight;
        //
        num_positives+=expected_output*the_weight; 
        num_negatives+=(!expected_output)*the_weight;
        num_positives_nominal+=expected_output; 
        num_negatives_nominal+=(!expected_output);
        //
        //evaluating model output
        ((classification_tools_data *)pdata)->bp_function(a_pattern.second.xy_values.data(),calculated_output);
        //
        
        //
        //making a prediction
        bool predicted_output=false;
        if(calculated_output[i_output]>=threshold) predicted_output=true;
        //
        
        if(predicted_output==expected_output && expected_output) tp+=the_weight;
        if(predicted_output==expected_output && !expected_output) tn+=the_weight;
      }
      //end of loop on patterns
      //
      //adding point to graph
      const double sensitivity=tp/num_positives;
      const double specificity=tn/num_negatives;      
      roc_graph.SetPoint(roc_graph.GetN(),1-specificity,sensitivity);
      //
    }
    // end of loop on threshold values
    
    //
    //Calculation of values
    roc_graph.SetPoint(roc_graph.GetN(),1,0);
    const double AUC_integral=roc_graph.Integral();
    const double q0=AUC_integral*(1-AUC_integral);
    const double q1=AUC_integral/(2-AUC_integral)-AUC_integral*AUC_integral;
    const double q2=2*AUC_integral*AUC_integral/(1+AUC_integral)-AUC_integral*AUC_integral;
    const double AUC_se=sqrt((q0+(num_positives_nominal-1)*q1+(num_negatives_nominal-1)*q2)/(num_positives_nominal*num_negatives_nominal));
    //
    result[i_output]=AUC_integral;
    ci_min[i_output]=AUC_integral-AUC_se*zcrit;
    ci_max[i_output]=AUC_integral+AUC_se*zcrit;
    //
  }
  //end of loop on classifier output
  
  //
  delete [] calculated_output;
  //
  
  return;
}

//____________________________________________________
void calculate_sensitivity (void * pdata, double* thr, double* result, double* ci_min, double * ci_max)
{
  //
  const double zcrit=1.959964;
  unsigned int num_positives_nominal=0;
  unsigned int num_negatives_nominal=0;
  //
  const unsigned int dimensionality_output=((classification_tools_data *)pdata)->num_y;
  //
  double * calculated_output = new double[dimensionality_output];
  //
  //loop on classifier output
  for(unsigned int i_output=0; i_output<dimensionality_output; i_output++)
  {
    double num_positives=0;
    double num_negatives=0;
    double tp=0;
    double tn=0;
        
    //
    //loop on patterns
    for(auto && a_pattern : ((classification_tools_data *)pdata)->patterns) 
    {
      bool expected_output=a_pattern.second.xy_values[((classification_tools_data *)pdata)->num_x+i_output];
      double the_weight=a_pattern.second.weight;
      //
      num_positives+=expected_output*the_weight; 
      num_negatives+=(!expected_output)*the_weight;
      num_positives_nominal+=expected_output; 
      num_negatives_nominal+=(!expected_output);
      //
      //evaluating model output
      ((classification_tools_data *)pdata)->bp_function(a_pattern.second.xy_values.data(),calculated_output);
      //
      
      //
      //making a prediction
      bool predicted_output=false;
      if(calculated_output[i_output]>=thr[i_output]) predicted_output=true;
      //
      
      if(predicted_output==expected_output && expected_output) tp+=the_weight;
      if(predicted_output==expected_output && !expected_output) tn+=the_weight;
    }
    //end of loop on patterns
    //
    //Calculation of sensitivity
    const double sensitivity=tp/num_positives;
    const double sensitivity_se=sqrt(sensitivity*(1-sensitivity)/num_positives_nominal);
    //
    result[i_output]=sensitivity;
    ci_min[i_output]=sensitivity-sensitivity_se*zcrit;
    ci_max[i_output]=sensitivity+sensitivity_se*zcrit;
    //
  }
  //end of loop on classifier output
  
  //
  delete [] calculated_output;
  //
  
  return;
}

//____________________________________________________
void calculate_specificity (void * pdata, double* thr, double* result, double* ci_min, double * ci_max)
{
  //
  const double zcrit=1.959964;
  unsigned int num_positives_nominal=0;
  unsigned int num_negatives_nominal=0;
  //
  const unsigned int dimensionality_output=((classification_tools_data *)pdata)->num_y;
  //
  double * calculated_output = new double[dimensionality_output];
  //
  //loop on classifier output
  for(unsigned int i_output=0; i_output<dimensionality_output; i_output++)
  {
    double num_positives=0;
    double num_negatives=0;
    double tp=0;
    double tn=0;
        
    //
    //loop on patterns
    for(auto && a_pattern : ((classification_tools_data *)pdata)->patterns) 
    {
      bool expected_output=a_pattern.second.xy_values[((classification_tools_data *)pdata)->num_x+i_output];
      double the_weight=a_pattern.second.weight;
      //
      num_positives+=expected_output*the_weight; 
      num_negatives+=(!expected_output)*the_weight;
      num_positives_nominal+=expected_output; 
      num_negatives_nominal+=(!expected_output);
      //
      //evaluating model output
      ((classification_tools_data *)pdata)->bp_function(a_pattern.second.xy_values.data(),calculated_output);
      //
      
      //
      //making a prediction
      bool predicted_output=false;
      if(calculated_output[i_output]>=thr[i_output]) predicted_output=true;
      //
      
      if(predicted_output==expected_output && expected_output) tp+=the_weight;
      if(predicted_output==expected_output && !expected_output) tn+=the_weight;
    }
    //end of loop on patterns
    //
    //Calculation of specificity
    const double specificity=tn/num_negatives;
    const double specificity_se=sqrt(specificity*(1-specificity)/num_negatives_nominal);
    //
    result[i_output]=specificity;
    ci_min[i_output]=specificity-specificity_se*zcrit;
    ci_max[i_output]=specificity+specificity_se*zcrit;
    //
  }
  //end of loop on classifier output
  
  //
  delete [] calculated_output;
  //
  
  return;
}

//____________________________________________________
void calculate_optimal_threshold (void * pdata, double * result)
{
  //
  const unsigned int dimensionality_output=((classification_tools_data *)pdata)->num_y;
  //
  //loop on classifier output
  for(unsigned int i_output=0; i_output<dimensionality_output; i_output++)
  {
    calculate_optimal_threshold(pdata, ((classification_tools_data *)pdata)->y_names[i_output].c_str(), result+i_output);
  }
  //end of loop on classifier output
  
  return;
}

//____________________________________________________
void calculate_optimal_threshold (void * pdata, const char * y_name, double * result)
{  
  //
  const double min_threshold=0;
  const double max_threshold=1;
  const double threshold_step=0.001;
  //
  const unsigned int num_points=(max_threshold-min_threshold)/threshold_step;
  //
  unsigned int i_output;
  for(i_output=0; i_output<((classification_tools_data *)pdata)->num_y; i_output++)
  {
    if(strcmp(((classification_tools_data *)pdata)->y_names[i_output].c_str(), y_name)==0) {
      break; 
    }
  }
  //
  std::map<double, double> performance_to_threshold;
  //
  const unsigned int dimensionality_output=((classification_tools_data *)pdata)->num_y;
  //
  double * calculated_output = new double[dimensionality_output];
  //
  //loop on threshold values
  for(unsigned int i_point=0; i_point<num_points; i_point++) 
  {
    double num_positives=0;
    double num_negatives=0;
    double tp=0;
    double tn=0;
    
    const double threshold=min_threshold+i_point*threshold_step;
    
    //
    //loop on patterns
    for(auto && a_pattern : ((classification_tools_data *)pdata)->patterns)
    {
      bool expected_output=a_pattern.second.xy_values[((classification_tools_data *)pdata)->num_x+i_output];
      double the_weight=a_pattern.second.weight;
      //
      num_positives+=expected_output*the_weight; 
      num_negatives+=(!expected_output)*the_weight;
      //
      //evaluating model output
      ((classification_tools_data *)pdata)->bp_function(a_pattern.second.xy_values.data(),calculated_output);
      //
      
      //
      //making a prediction
      bool predicted_output=false;
      if(calculated_output[i_output]>=threshold) predicted_output=true;
      //
      
      if(predicted_output==expected_output && expected_output) tp+=the_weight;
      if(predicted_output==expected_output && !expected_output) tn+=the_weight;
    }
    //end of loop on patterns
    //
    //Calculation of specificity and sensitivity
    const double sensitivity=tp/num_positives;
    const double specificity=tn/num_negatives;      
    //
    performance_to_threshold[sensitivity+specificity]=threshold;
    //
  }
  // end of loop on threshold values
  
  //
  auto best_result = std::max_element(performance_to_threshold.begin(), performance_to_threshold.end());
  *result=best_result->second;
  //

  //
  delete [] calculated_output;
  //
  
  return;
}

//____________________________________________________
void split_data (void * pdata,
                double perc_lrn, const char * file_lrn_name,
                double perc_tst, const char * file_tst_name,
                double perc_crs, const char * file_crs_name)
{
  //
  std::vector<std::string> file_names;
  std::vector<double> perc_values;
  //
  file_names.push_back(Form("%s.lrn", file_lrn_name));
  file_names.push_back(Form("%s.tst", file_tst_name));
  file_names.push_back(Form("%s.crs", file_tst_name));
  //
  perc_values.push_back(perc_lrn);
  perc_values.push_back(perc_tst);
  perc_values.push_back(perc_tst);
  //
  
  //
  split_data(pdata, 3, perc_values.data(), (const char **) file_names.data());
  //
  
  return;
}

//____________________________________________________
void split_data (void * pdata, unsigned int num_files, double * perc_values, const char ** file_names)
{
  //
  //opening files
  std::ofstream * file_patterns = new std::ofstream[num_files];
  std::ofstream * file_unique_ids = new std::ofstream[num_files];
  std::ofstream * file_weights = new std::ofstream[num_files];
  for(unsigned int i_file=0; i_file<num_files; i_file++)
  {
    //
    std::string pattern_file_name(file_names[i_file]);
    //
    std::string file_name_prefix(pattern_file_name.substr(0,pattern_file_name.find_last_of('.')));
    //
    std::string pattern_file_extension(pattern_file_name.substr(pattern_file_name.find_last_of('.')+1));
    std::string unique_id_file_extension(Form("id%s", pattern_file_extension.c_str()));
    std::string weight_file_extension(Form("w%s", pattern_file_extension.c_str()));
    //
    file_patterns[i_file].open(pattern_file_name.c_str());
    file_unique_ids[i_file].open(Form("%s.%s", file_name_prefix.c_str(), unique_id_file_extension.c_str()));
    if(((classification_tools_data *)pdata)->is_pattern_weights) {
      file_weights[i_file].open(Form("%s.%s", file_name_prefix.c_str(), weight_file_extension.c_str()));
    }
    //
  }
  //
  
  //
  std::random_device rd_dev;
  std::mt19937_64 generator(rd_dev());
  std::uniform_real_distribution<> distr(0,1);
  //
  
  //
  //creation of separation values used to split patterns into files
  const unsigned int num_separation_values=num_files-1;
  double * separation_values = new double[num_separation_values];
  for(unsigned int i_separation=0; i_separation<num_separation_values; i_separation++)
  {
    separation_values[i_separation]=(i_separation>0 ? separation_values[i_separation-1] : 0) + perc_values[i_separation];
  }
  //
  
  //
  //loop on patterns
  for(auto && a_pattern : ((classification_tools_data *)pdata)->patterns)
  {
    const double random_number=distr(generator);
    //
    //determination of which file to use
    unsigned int i_file=0;
    //
    for(; i_file<num_separation_values; i_file++)
    {
      if(random_number<separation_values[i_file]) 
      {
        break; 
      }
    }
    //
    
    //
    //writing pattern to file
    for(auto && pattern_component : a_pattern.second.xy_values)
    {
      file_patterns[i_file] << std::setw(20) << pattern_component;
    }
    file_patterns[i_file] << std::endl;
    //
    //writing pattern id to file
    file_unique_ids[i_file] << a_pattern.first << std::endl;
    //
    if(((classification_tools_data *)pdata)->is_pattern_weights) {
      //
      file_weights[i_file] << a_pattern.second.weight << std::endl;
      //
    }
  }
  //
  
  //
  //closing files
  for(unsigned int i_file=0; i_file<num_files; i_file++)
  {
    file_patterns[i_file].close();
    file_unique_ids[i_file].close();
    if(((classification_tools_data *)pdata)->is_pattern_weights) {
      file_weights[i_file].close();
    }
  }
  //
  delete [] file_patterns;
  delete [] file_unique_ids;
  delete [] file_weights;
  //
  
  return;  
}

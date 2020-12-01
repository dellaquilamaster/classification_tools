#ifndef CLASSIFICATION_TOOLS_H
#define CLASSIFICATION_TOOLS_H

/*
 * classification_tools
 * Contributors: Daniele Dell'Aquila (ddellaquila@uniss.it), Marco Russo (marco.russo@ct.infn.it)
 * 
 * 27/11/2020
 * v1.0 
 * 
 * 
 */

//
//Constructors
//file_name: pathname of the file containing the data in the standard BP format
//file_x_names: pathname of the file containing the names associated to the features
//file_y_names: pathname of the file containing the names associated to the outputs
//num_x: number of features
//num_y: number of outputs
//unique_id_file: file containing the id associated to each pattern in the same order of the initial data file (if not given, the line number is used as unique id)
void *   create_classification_data (const char * file_name, const char * file_x_names, const char * file_y_names, unsigned int num_x, unsigned int num_y, const char * unique_id_file=0);
//
//weights_file: file containing the weights associated to each pattern in the same order of the initial data file (optional)
//costs_file: file containing the cost associated to each feature (optional)
//output_weights_file: file containing the weight associated to each output (optional)
int     add_pattern_weights(void * pdata, const char * weights_file);
int     add_feature_costs(void * pdata, const char * costs_file);
int     add_output_weights(void * pdata, const char * output_weights_file);
//

//
//Set BP classification function
void    set_function(void * pdata, void (*function)(double *, double *));
//

//
//Get the pattern associated to a given id
double  * get_pattern(void * pdata, const char * unique_id);
//

//
//Draw functions
//file_output_name: name of the output picture file. .png, .pdf, .ps, .eps.
void    draw_roc (void * pdata, const char * file_output_name, double min_threshold=0, double max_threshold=1, double threshold_step=0.001);
void    draw_confusion_matrix (void * pdata, double * thr, const char * file_output_name); //single-label confusion matrix. It is assumed that only one classifier returns a value greater than the thr value.
void    draw_fuzzy_values (void * pdata, const char * unique_id, const char * file_output_name);
void    draw_fuzzy_values (void * pdata, const char * unique_id, double * thr, const char * file_output_name); //draws only fuzzy values of significative output (i.e. > thr)
//altri grafici vedi Campobello
//

//
//Calculate functions
//result: final results are stored in this vector
//ci_min: minimum value of 95% confidence interval
//ci_min: maximum value of 95% confidence interval
void    calculate_auc (void * pdata, double * result, double* ci_min, double * ci_max);
void    calculate_sensitivity (void * pdata, double* thr, double* result, double* ci_min, double * ci_max);
void    calculate_specificity (void * pdata, double* thr, double* result, double* ci_min, double * ci_max);
//altri indicatori vedi Campobello
void    calculate_optimal_threshold (void * pdata, double * result); //calculates the optimal threshold to maximize the sum of specificity and sensitivity in each class
void    calculate_optimal_threshold (void * pdata, const char * y_name, double * result); // calculate the optimal threshold of a given classifier associated to the name y_name
//

//
//split data for BP
//perc_lrn: fraction of patterns used for the learning file
// same for tst and crs
//file_lrn_name: name used for the learning file (file_lrn_name.lrn, file_lrn_name.wlrn and file_lrn_name.idlrn are created)
// same for tst and crs
void    split_data (void * pdata,
                    double perc_lrn, const char * file_lrn_name,
                    double perc_tst, const char * file_tst_name,
                    double perc_crs, const char * file_crs_name);
//num_files: number of files produced by the splitting procedure
//perc_values: array containing the fraction of patterns assigned to each file
//file_names: array of strings containing the file names
void    split_data (void * pdata, unsigned int num_files, double * perc_values, const char ** file_names);
//

//
//methods to identify misclassified patterns
void *   create_list_misclassified(void * pdata, double * thr_min, double * thr_max); //NOTE: ragionare sulla misclassificazione
int     print_list_misclassified(void * misclassified_data);
//erase "misclassified_data" from "data"
int     erase_misclassified(void * pdata, void * misclassified_data);
//erase data misclassified according to thr_min and thr_max
int     erase_misclassified(void * pdata, double * thr_min, double * thr_max);
//

#endif

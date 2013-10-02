//header file for all functions in the package
#include <R.h>
#include <Rinternals.h>

/*function to run recursive loop to calculate all models with fully
independent structures that lie within threshold of maximum model -
to be called from "calc_PPAs_approx_fn"*/
void calc_approx_indep_recur(int nsamp, int nmod, int nmodcol, int incmods, int *totmods, int currsamp, int *r, int *currmods, int *currlocs, int *maxmods, int *maxlocs, double *max_lPPA, double curr_lPPA, double thresh, int **models_num, int **hyp, double *lPDM_int_mat, int *lPDM_int_ind, double **lPPA_mat);

//function to calculate correct column of lPDM_int_mat based on a single set of changes
int calc_corr_col(int q, int nrow, int *structure, int nsamp, int *ncomb_sub, int samp,  int sum, int *temp_index);

/*function to discriminate null and alternative models
based on less-stringent hypothesis*/
void calc_lessstring_fn(int *nsamp1, int *totmods1, int *models_num, int *hyp);

/*function to discriminate null and alternative models
based on stringent hypothesis*/
void calc_string_fn(int *nsamp1, int *totmods1, int *models_num, int *hyp);

/*function to discriminate null and alternative models based on 
less-stringent hypothesis for a single model structure*/
int calc_lessstring_single(int nsamp, int *models_num);

/*function to discriminate null and alternative models based on 
stringent hypothesis for a single model structure*/
int calc_string_single(int nsamp, int *models_num);

//function for calculating log[P'(D|M)]s for all unique models for a given set of base distributions
void calc_lPDM_fn(int *nsamp1, int *totmods1, int *models_num, double *lPDM_int_mat, double *lPDM_mat);

//function to calculate intermediate values for log[P'(D|M)]
void calc_lPDM_int_fn(int *nsamp1, int *z, double *pstar1, double *lPDM_int_mat);

//function to simplify a single unique model and return the log[P'(D|M)]
double calc_lPDMs_single(int r, int nsamp, int totmods, int nmod, int *models_num, double *lPDM_int_mat, int *ncomb_sub, int *sum, int *temp_index, int *mods, int *inds);

//function to simplify a single unique model and return the log[P'(D|M)]
//(EXACTLY THE SAME AS THE ABOVE FUNCTION ONLY FOR WHEN MODELS_NUM IS A SINGLE VECTOR)
double calc_lPDMs_single_vec(int nsamp, int nmod, int *models_num, double *lPDM_int_mat, int *ncomb_sub, int *sum, int *temp_index, int *mods, int *inds);

//function to calculate "top" set of models for use in approximation
SEXP calc_PPAs_approx_fn(SEXP nsamp1, SEXP ntotcol1, SEXP logc1, SEXP lpriornull1, SEXP lprioralt1, SEXP lPDM_int_mat1, SEXP string1, SEXP structure1, SEXP nstruct1, SEXP return_final1, SEXP uni1, SEXP uni_ind1);

//internal function to check for duplicated combinations
void check_duplicates(int start, int ncomb, int ncol, int *comb, int *ind);

//log-factorial function
double lfactorial(int n);

//factorial function
int factorial(int n);

//choose function
int choose_c(int n, int k);

/*function to run recursive loop to calculate final lPPAs for Occam's window approach without recording all models*/
void final_lPPA_recur(int nsamp, int nmod, int nmodcol, int incmods, int *totmods, int currsamp, int *r1, int m1, int *currmods, double curr_lPPA, int curr_hyp, int **models_num, int **hyp, double **lPPA_mat, int **multfact, int *uni, int *uni_ind, int string, double lpriornull, double lprioralt);

//function to generate all possible models based on an intermediate set of structures
void genfullmodels(int q, int nrow, int ncol, int *structure, int *r, int *totmods, int incmods, int **models_num);

//function to generate top subset of fully and patially dependent models based on an intermediate set of structures
void genfullmodels_approx(int q, int nrow, int ncol, int *structure, int *r, int *totmods, int incmods, int nmodcol, double *lPDM_int_mat, int *lPDM_int_ind,  int ntotcol, int **models_num, int **hyp, double **lPPA_mat, double logthresh, int *ncomb_sub);

//function to generate top subset of fully and patially dependent models based on an intermediate set of structures
void genfullmodels_approx_recur(int nhier, int currhier1, int *temp, int *tempind, int *currmods, int q, int nrow, int ncol, int *structure, int *r, int *totmods, int incmods, int nmodcol, double *lPDM_int_mat, int *lPDM_int_ind, int ntotcol, int **models_num, int **hyp, double **lPPA_mat, double logthresh, int *ncomb_sub, double curr_lPPA, double *max_lPPA);

/*function to generate all possible models based on an intermediate set of structures
(returns just normalising constants for efficient EXACT routine)*/
void genfullmodels_exact(int q, int nrow, int ncol, int *structure, int *r, int *models_num, double *norm_ls, double *norm_s, double *palt_ls, double *palt_s, double lpriornull_ls, double lprioralt_ls, double lpriornull_s, double lprioralt_s, double *lPDM_int_mat, int *ncomb_sub, int *sum, int *temp_index, int *mods, int *inds, int *prec_ind);

/*function to generate independent model structure recursively
(for use in efficient EXACT routine)*/
void genfullmodels_indep_exact(int currcol, int nsamp, int nmod, int *models_num, double *norm_ls, double *norm_s, double *palt_ls, double *palt_s, double lpriornull_ls, double lprioralt_ls, double lpriornull_s, double lprioralt_s, double *lPDM_int_mat, int *ncomb_sub, int *sum, int *temp_index, int *mods, int *inds, int *prec_ind);

/*function to generate independent model structure recursively
(for use in approximation routine)*/
void genfullmodels_indep_priors(int currcol, int nsamp, int nmod, int *models_num, int *nalt_ls, int *nalt_s);

/*function to generate all possible models based on an intermediate set of structures
(returns just prior information)*/
void genfullmodels_priors(int q, int nrow, int ncol, int *structure, int *r, int *models_num, int *nalt_ls, int *nalt_s);

//recursive function to be called from 'genfullmodels'
void genfullmodels_recur(int nhier, int currhier1, int *temp, int *tempind, int q, int nrow, int ncol, int *structure, int *r, int *totmods, int incmods, int **models_num);

//recursive function
void genfullmodels_recur_exact(int nhier, int currhier1, int *temp, int *tempind, int q, int nrow, int ncol, int *structure, int *r, int *models_num, double *norm_ls, double *norm_s, double *palt_ls, double *palt_s, double lpriornull_ls, double lprioralt_ls, double lpriornull_s, double lprioralt_s, double *lPDM_int_mat, int *ncomb_sub, int *sum, int *temp_index, int *mods, int *inds, int *prec_ind);

//recursive function to be called from 'genfullmodels_priors'
void genfullmodels_recur_priors(int nhier, int currhier1, int *temp, int *tempind, int q, int nrow, int ncol, int *structure, int *r, int *models_num, int *nalt_ls, int *nalt_s);

/*function to generate all model structures in an efficient manner
(returns an R vector of a given size, though it is not necessary to
specify the size in advance)*/
SEXP genmodels(SEXP nsamp1);

/*function to generate all model structures in an efficient manner
and calculate normalising for use in EXACT efficient search routine*/
SEXP genmodels_exact(SEXP nsamp1, SEXP ntotcol1, SEXP lpriors, SEXP lPDM_int_mat1);

/*function to generate all model structures in an efficient manner
and calculate prior information for use in efficient search routine*/
SEXP genmodels_priors(SEXP nsamp1);

//recursive function to be called from 'genmodels'
void genmodels_recur(int nsamp, int nstart, int *indexl, int *indexl_ind, int repeatind, int *comb_next, int ncomb_next, int ntotsamp, int *currstruct, int *q, int **structure, int *maxmodels, int incmodels);

//indexing function for matrices
inline int index2(int i, int j, int nrow)
{
	return (nrow*j)+i;
}

//indexing function for matrices
inline int index2_col(int i, int j, int ncol)
{
	return (ncol*i)+j;
}

//indexing function for matrices
inline int index3(int i, int j, int k, int nrow, int ncol)
{
	return (k*nrow*ncol)+(nrow*j)+i;
}

//recursive function to produce lexicographic orderings
void lexico_recur(int nsamp, int col, int totcol, int *temp_samp, int *temp_order1, int norder, int start, int *overall_order, int *ncomb_sub);

//function to calculate l[P'(D|M)] for a given distribution of bases
double lPDM_mod_fn(int *z, int ind, double pstar);

/*This code to calculate combinations was adapted from a post by scvalex at
http://compprog.wordpress.com/2007/10/17/generating-combinations-1/.
It's very elegant and saved me a headache---many thanks!*/
int next_comb(int *comb, int k, int n);

//function for calculating integer powers
int pow_int(int base, int exp);

//function to perform memory reallocation of intermediate structures if required
void realloc_maxmod(int nsamp, int *maxmodels, int incmodels, int **structure);

//function to perform memory reallocation of full model set if required
void realloc_maxmod_full(int nsamp, int *maxmodels, int incmodels, int **models_num);

//function to perform memory reallocation of full model set if required
void realloc_approx(int nmodcol, int *totmods, int incmods, int **models_num, int **hyp, double **lPPA_mat);

//function to perform memory reallocation of full model set if required
void realloc_approx_new(int nmodcol, int *totmods, int incmods, int **models_num, int **hyp, double **lPPA_mat, int **multfact);

/*bubble sort algorithm (sorts into decreasing order
and provides indicators in order to sort multiple columns)*/
void bubble_sort_dec(double *values, int *ind, int n);

/*bubble sort algorithm for INTEGERS (sorts into decreasing order
and provides indicators in order to sort multiple columns)*/
void bubble_sort_dec_int(int *values, int *ind, int n);

/*bubble sort algorithm for INTEGERS (sorts into increasing order
and provides indicators in order to sort multiple columns)*/
void bubble_sort_inc_int(int *values, int *ind, int n);

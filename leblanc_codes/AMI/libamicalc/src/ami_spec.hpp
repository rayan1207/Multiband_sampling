
#pragma once
#include "ami_base.hpp"
#include "ami_calc.hpp"
#include "linterp.h"

#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <complex>
#include <vector>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <random>
#include <experimental/filesystem>


/**
 * @class AmiSpec
 *
 *
 *
 * @brief  A spectral application of the AmiBase class.  Experimental.
 *
 * @note N/A
 *
 * @author JPF LeBlanc
 *
 * @version Revision: 0.4
 *
 * @date Date: 2020/11/03
 *
 *
 * Contact: jleblanc@mun.ca
 *
 *
 *
 *
 */

class AmiSpec
{

public:

AmiBase amibase;
NewAmiCalc ami;

double xi_cutoff=10;
double global_gamma=0.5;//0.5;
bool regulate=true;
bool use_sigma_file=false;

bool overflow_detected=false;

std::vector<double> AMI_spec_se_Im_vector;
std::vector<double> AMI_spec_se_Im_err_vector;
std::vector<double> AMI_spec_se_Re_vector;
std::vector<double> AMI_spec_se_Re_err_vector;
std::vector<double> AMI_spec_freq_vector;
std::vector<double> AMI_spec_ky_vector;
std::vector<double> AMI_spec_kx_vector;
std::vector<double> AMI_spec_freq_vector_simplified;
std::vector<double> AMI_spec_kx_vector_simplified;
std::vector<double> AMI_spec_ky_vector_simplified;
// Define structures that are not in the ami_base class
// Same as above but for W
std::vector<double> AMI_W_se_Im_vector;
std::vector<double> AMI_W_se_Im_err_vector;
std::vector<double> AMI_W_se_Re_vector;
std::vector<double> AMI_W_se_Re_err_vector;
std::vector<double> AMI_W_freq_vector;
std::vector<double> AMI_W_ky_vector;
std::vector<double> AMI_W_kx_vector;
std::vector<double> AMI_W_freq_vector_simplified;
std::vector<double> AMI_W_kx_vector_simplified;
std::vector<double> AMI_W_ky_vector_simplified;



typedef AmiBase::epsilon_t X_t ;

// linterp testing functions 
NewAmiCalc::k_vector_t remap_k(NewAmiCalc::k_vector_t &k_in);
void initialize_linterp();
std::vector< std::vector<double>::iterator > grid_iter_list;
array<int,3> grid_sizes;
int num_elements;

void initialize_W_linterp();
std::vector< std::vector<double>::iterator > W_grid_iter_list;
array<int,3> W_grid_sizes;
int W_num_elements;

// InterpMultilinear<3, double> interp_ML_REAL();
// InterpMultilinear<3, double> interp_ML_IMAG();

std::complex<double> linterp_sigma(NewAmiCalc::k_vector_t &k, std::complex<double> &X);
std::complex<double> linterp_W(NewAmiCalc::k_vector_t &k, std::complex<double> &X);


// the A_struct stores the symbolic information about A.  you will need your own function to take the A_struct and some self energy and return a value.
struct A_struct{
A_struct(AmiBase::epsilon_t eps, X_t x){
eps_=eps;
x_=x;
species_=0;
}

A_struct(AmiBase::epsilon_t eps, X_t x,AmiBase::species_t species){
eps_=eps;
x_=x;
species_=species;
}


/// uninitialized variant
A_struct(){
}

// epsilon_i either by index or by value
int eps_index=-1;
std::complex<double> eps_val;

AmiBase::epsilon_t eps_;
X_t x_; // not to be confused with xi_t.  X_t is symbolic.  xi_t is numbers. 
AmiBase::alpha_t x_alpha_;
AmiBase::alpha_t alpha_;
AmiBase::species_t species_;
AmiBase::stat_type eff_stat_;


};




typedef std::vector<A_struct> A_prod_t;


typedef AmiBase::g_struct delta_t;
typedef std::vector< delta_t > delta_prod_t;
typedef AmiBase::energy_t xi_t;

xi_t this_xi;



struct ami_spec_vars{




ami_spec_vars(AmiBase::energy_t eps, AmiBase::frequency_t freq, NewAmiCalc::k_vect_list_t k, xi_t x,  double Bta, std::complex<double> mu, double pf){
energy_= eps;
frequency_= freq;
prefactor=pf;
BETA_=Bta;
k_list_=k;
xi_list_=x;
MU_=mu;
}

ami_spec_vars(){prefactor=1.0;}


AmiBase::energy_t energy_;
AmiBase::frequency_t frequency_;
double prefactor;
double BETA_;
std::complex<double> MU_;

NewAmiCalc::k_vect_list_t k_list_;
xi_t xi_list_;

///Experimental parameter for spectral representation 
double gamma_=0;
std::vector< AmiBase::alpha_t > alpha_list;

};

typedef std::vector< ami_spec_vars > ami_spec_vars_list;



struct ami_sp_term{

ami_sp_term(AmiBase::term this_term, A_prod_t aprod, delta_prod_t dprod){
aprod_=aprod;
ami_term_=this_term;

}

//uninitialized
ami_sp_term(){}

A_prod_t aprod_;
AmiBase::term ami_term_;// prod of G's and prod of F's and a sign 
delta_prod_t dprod_;
delta_prod_t regprod_;

bool root=true;
int delta_count=0;

};

typedef std::vector<ami_sp_term> sp_terms;


struct spec_solution_set{
	
spec_solution_set(){}
AmiBase::terms ami_terms_;
sp_terms sp_terms_;
AmiBase::g_prod_t R0_;

double prefactor_;
int loops_;
int ct_count_=0;
int sigma_ct_count_=0;
int order_shift=0;

std::vector<AmiBase::alpha_t> bose_alphas_;
AmiBase::ami_parms ami_parms_;


};

typedef std::vector<spec_solution_set> sp_solution_set_vec_t;
typedef std::vector< std::vector<spec_solution_set> > sp_solution_set_matrix_t;
typedef std::vector< sp_solution_set_matrix_t > sp_gg_solution_set_matrix_t;




// initial attempt
std::complex<double> evaluate_sp_term(AmiBase::ami_parms &parms, AmiSpec::ami_sp_term &sp_term, NewAmiCalc::ext_vars &ev,   AmiBase::ami_vars &external, NewAmiCalc::k_vect_list_t &klist,   xi_t &xi_list);

// simplified usage 
std::complex<double> evaluate_sp_term(AmiBase::ami_parms &parms, AmiSpec::ami_sp_term &sp_term, AmiSpec::ami_spec_vars &vars);
ami_spec_vars construct_ami_spec_vars(AmiBase::g_prod_t &R0, double prefactor, NewAmiCalc::internal_state &state, NewAmiCalc::ext_vars &external, xi_t &xi);
void construct_ami_spec_vars_list(AmiBase::g_prod_t &R0, double prefactor, NewAmiCalc::internal_state &state, NewAmiCalc::ext_vars &external, xi_t &xi,ami_spec_vars_list &vars_list);

std::complex<double> evaluate_sp_terms(AmiBase::ami_parms &parms, AmiSpec::sp_terms &sp_terms, AmiSpec::ami_spec_vars &vars);

void evaluate_spectral_solutions(std::vector<double> &Re_results, std::vector<double> &Im_results, spec_solution_set &AMI, ami_spec_vars_list &ami_spec_eval_vars,   xi_t &xi_list);
void  construct_ami_spec_vars_list(AmiBase::g_prod_t &R0, double prefactor, NewAmiCalc::internal_state &state, NewAmiCalc::external_variable_list &external, xi_t &xi,ami_spec_vars_list &vars_list);



std::complex<double> A_eval(std::complex<double> &sigma, std::complex<double> &X, std::complex<double> &E);

std::complex<double> get_A(A_struct &A, double this_x, NewAmiCalc::k_vector_t k);
// std::complex<double> eval_Aprod(A_prod_t &Ap, xi_t &xi, NewAmiCalc::k_vect_list_t &klist, std::complex<double> &mu);
std::complex<double> eval_Aprod(A_prod_t &Ap, xi_t &xi, AmiBase::frequency_t &freq, NewAmiCalc::k_vect_list_t &klist, std::complex<double> &mu);
std::complex<double> eval_tb(double t, double tp, NewAmiCalc::k_vector_t &k, std::complex<double> &mu);


// std::complex<double> get_X(X_t &Xsym, xi_t &xi);
std::complex<double> get_X(X_t &Xsym, xi_t &xi, AmiBase::alpha_t &x_alpha_, AmiBase::frequency_t &freq);
std::complex<double> get_E(AmiBase::energy_t &ei, AmiBase::epsilon_t &eps);
void randomize_xi(xi_t &xi, int length);

// Sigma
std::complex<double> get_sigma(NewAmiCalc::k_vector_t &k, std::complex<double> &X);
void find_closest_points_in_vector(double &closest_lt,double &closest_gt,double point, std::vector<double> vec);
std::vector<double> return_simple_grid_vector(std::vector<double> &in_vector);
void read_self_energy(std::string file_name);
void read_W(std::string file_name);
//

//
std::complex<double> get_regdelta(double &sign, delta_prod_t &dprod,  AmiBase::ami_parms &parms ,AmiBase::ami_vars &external);
//


std::complex<double> construct_energy(AmiBase::alpha_t &alpha, NewAmiCalc::k_vect_list_t &klist, std::complex<double> &mu);

void generate_sp_terms(AmiBase::term &start_term, sp_terms &new_sp_terms, AmiBase::g_prod_t &R0);
void R0_to_Aprod(AmiBase::g_prod_t &R0, A_prod_t &Ap);


void generate_simple_sp_terms(AmiBase::term &start_term, sp_terms &new_sp_terms, AmiBase::g_prod_t &R0);
void generate_simple_sp_terms(AmiBase::terms &start_terms, sp_terms &full_sp_terms, AmiBase::g_prod_t &R0);


void generate_sp_terms(AmiBase::terms &start_terms, sp_terms &full_sp_terms, AmiBase::g_prod_t &R0);

// void reduce_deltas(ami_sp_term &term);
void resolve_deltas( sp_terms &sp_terms);
void resolve_deltas(ami_sp_term &sp_term);

void spawn_regulator_terms( sp_terms &sp_terms);
void spawn_regulator_terms( ami_sp_term &term_in, sp_terms &terms_out);


void replace_xi(int i, AmiBase::pole_array_t &pv, ami_sp_term &sp_term);
void update_spec_pole(AmiBase::pole_struct &source_pole, AmiBase::alpha_t &target_alpha, AmiBase::epsilon_t &target_eps);



// OPT stuff

void evaluate_old_spectral_solutions(std::vector<double> &Re_results, std::vector<double> &Im_results, NewAmiCalc::solution_set &AMI, ami_spec_vars_list &ami_eval_vars,    std::vector< std::complex<double>>  &xi_list, double &xi_cut);
std::complex<double> evaluate_old_spectral(AmiBase::ami_parms &parms, AmiBase::g_prod_t &R0, AmiBase::R_t &R_array, AmiBase::P_t &P_array, AmiBase::S_t &S_array, ami_spec_vars &external, AmiBase::g_prod_t &unique_g,  AmiBase::R_ref_t &Rref,AmiBase::ref_eval_t &Eval_list,  std::vector< std::complex<double>>  &xi_list);

std::complex<double> eval_spectral_product(std::vector<std::complex<double>> &Ei,  std::vector< std::complex<double>>  &xi, std::vector< std::complex<double>> &sigma_list);
void get_sigma_list(std::vector< AmiBase::alpha_t> &alpha_list,  NewAmiCalc::k_vect_list_t &k_list,  std::vector< std::complex<double>>  &xi,std::vector< std::complex<double>> &sigma_list);
///




AmiSpec(double xc);
AmiSpec();

// io functions for debugging
void print_delta_prod_t(delta_prod_t &delta_prod);
void print_reg_prod_t(delta_prod_t &regprod);
void print_delta_t(delta_t &delta);
void print_a_prod_t(A_prod_t &Aprod);
void print_a_struct(A_struct &A);
void print_sp_term(ami_sp_term &term);
void print_sp_terms(sp_terms &sp_terms);
void print_int_vec(std::vector<int> vec);

void assign_corners_indices(double freq_lt,double freq_gt,double kx_lt,double kx_gt,double ky_lt,double ky_gt,std::vector<double>& AMI_spec_freq_vector, std::vector<double>& AMI_spec_kx_vector,std::vector<double>& AMI_spec_ky_vector,int& lll_corner,int& llg_corner,int& lgl_corner,int& lgg_corner,int& gll_corner,int& glg_corner,int& ggl_corner,int& ggg_corner);
void assign_corners_freq_interp(double& se_Im_ll_corner,double& se_Im_lg_corner,double& se_Im_gl_corner,double& se_Im_gg_corner,double& se_Re_ll_corner,double& se_Re_lg_corner,double& se_Re_gl_corner,double& se_Re_gg_corner,std::vector<double>& AMI_spec_se_Im_vector,std::vector<double>& AMI_spec_se_Re_vector,int lll_corner,int llg_corner,int lgl_corner,int lgg_corner,int gll_corner,int glg_corner,int ggl_corner,int ggg_corner,std::complex<double> &X);




private:


};


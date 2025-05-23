
#pragma once
#include "ami_base.hpp"
 
#include <string> 
#include <sstream> 
#include <fstream>
#include <iostream>
#include <complex>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <random>
#include <experimental/filesystem>



/**
 * @class NewAmiCalc
 *
 *  
 *
 * @brief  A somewhat more targetted application of the AmiBase class.  Intended primarily as an example, but for some problems my have direct applicability. 
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

class NewAmiCalc 
{

public:

AmiBase amibase;

// random stuff
		std::random_device rd;
    std::mt19937 nac_rand_gen;
		int random_int(int min, int max);
		double random_real(double min, double max);
		void randomize_xi(std::vector<double> &xi, int length, double max);
		
		bool overflow_detected=false;
		
		


// momenta - since these are used to generate epsilon - typically not a part of AmiBase
typedef std::vector< double> k_vector_t;
typedef std::vector< k_vector_t> k_vect_list_t;

// Hopping variables for energy evaluation - lists and species - Typically not a part of this 
typedef double hopping_t;
typedef std::vector<hopping_t> hopping_list_t;
typedef int species_t;

// Similarly, the 'state' of the system is external to AmiBase

struct internal_state{

internal_state(int k_length, int dim){
internal_k_list_.assign(k_length, std::vector< double>(dim,0.0));	
dim_=dim;
order_=k_length;
if(2*order_-1>0){
t_list_.resize(2*order_-1,1);
tp_list_.resize(2*order_-1,0);
}
disp_=static_cast<AmiBase::disp_type>(0); // by default tight binding unless necessary to change 

mink_=0;
maxk_=2*M_PI;


internal_freq_size_=k_length;

}

internal_state(){}

void initialize(int k_length, int dim){
internal_k_list_.assign(k_length, std::vector< double>(dim,0.0));	
dim_=dim;
order_=k_length;
if(2*order_-1>0){
t_list_.resize(2*order_-1,1);
tp_list_.resize(2*order_-1,0);
}

disp_=static_cast<AmiBase::disp_type>(0); // by default tight binding unless necessary to change 
mink_=0;
maxk_=2*M_PI;
	
internal_freq_size_=k_length;	
	
}

k_vect_list_t internal_k_list_;
int order_;
int group_;
int num_;
int ext_;
int term_ind_;

std::complex<double> current_value, proposed_value, last_value, measure_value;
int sample_count_=0;
int norm_count_=0;
double norm_=0;
bool physical=true;

int dim_;

int internal_freq_size_;

hopping_list_t t_list_;
hopping_list_t tp_list_;
AmiBase::disp_type disp_;

double mink_,maxk_;
double Jk_=1;

};


// External variables 

struct ext_vars{

ext_vars(int dim){
	
KDIM_=dim;
// external_k_vector_.assign(dim,0.0);	
dummy_k_.assign(dim,0.0);
external_k_list_.push_back(dummy_k_);
external_freq_.resize(1);


	
}
ext_vars(int dim, double beta, std::complex<double> mu){
	
KDIM_=dim;
// external_k_vector_.assign(dim,0.0);
dummy_k_.assign(dim,0.0);
external_k_list_.push_back(dummy_k_);	
external_freq_.resize(1);

BETA_=beta;
MU_=mu;

H_=0;
gamma_=0.01;
	
}

ext_vars(int dim, double beta, std::complex<double> mu, double H){
	
KDIM_=dim;
// TODO: make this a vector of external momentum vectors 
// external_k_vector_.assign(dim,0.0);	
dummy_k_.assign(dim,0.0);
external_k_list_.push_back(dummy_k_);
external_freq_.resize(1);

BETA_=beta;
MU_=mu;

H_=H;
gamma_=0.01;
	
}

ext_vars(int dim, double beta, std::complex<double> mu, double H, double gamma){
	
KDIM_=dim;
// TODO: make this a vector of external momentum vectors 
// external_k_vector_.assign(dim,0.0);	
dummy_k_.assign(dim,0.0);
external_k_list_.push_back(dummy_k_);
external_freq_.resize(1);

BETA_=beta;
MU_=mu;

H_=H;
gamma_=gamma;
	
}

ext_vars(){
	external_freq_.resize(1);
	
}


// k_vector_t external_k_vector_;
k_vector_t dummy_k_;
k_vect_list_t external_k_list_;
	
int KDIM_;
AmiBase::frequency_t external_freq_;

double BETA_;
double H_;
std::complex<double> MU_;
double gamma_=0.01;
};

typedef std::vector< ext_vars> external_variable_list;
//
typedef std::vector<AmiBase::ami_vars> ami_vars_list;

// The solution set structure is the complete package that is passed for evaluation 

struct solution_set{
	
solution_set(AmiBase::g_prod_t test_R0,AmiBase::S_t  test_S, AmiBase::P_t test_P, AmiBase::R_t test_R,AmiBase::ami_parms test_amiparms, double prefactor){

R0_=test_R0;
S_=test_S;
P_=test_P;
R_=test_R;
ami_parms_=test_amiparms;
prefactor_=prefactor;


}

solution_set(AmiBase::g_prod_t test_R0,AmiBase::S_t  test_S, AmiBase::P_t test_P, AmiBase::R_t test_R,AmiBase::ami_parms test_amiparms, double prefactor, std::vector<AmiBase::alpha_t> bose){

R0_=test_R0;
S_=test_S;
P_=test_P;
R_=test_R;
ami_parms_=test_amiparms;
prefactor_=prefactor;

bose_alphas_=bose;

}

solution_set(){}
	
AmiBase::g_prod_t R0_;
AmiBase::S_t S_;
AmiBase::P_t P_;
AmiBase::R_t R_;
AmiBase::ami_parms ami_parms_;
double prefactor_;
int loops_;
int ct_count_=0;
int sigma_ct_count_=0;
int order_shift=0;

std::vector<AmiBase::alpha_t> bose_alphas_; // TODO: I think this is depricated 

//optimization

AmiBase::g_prod_t Unique;
AmiBase::R_ref_t Rref;
AmiBase::ref_eval_t Eval_list;
AmiBase::Pi_t Unique_poles;

// Term construction 

AmiBase::terms ami_terms_;

//

// Experimental Spectral Representation 
// std::vector<std::vector<int>> INT_XI_VEC; // which should have size of the epsilon of R0_ (R0_.size()). Entries of 1 mean it is an independent xi. 0 means it is dependent. 
// std::vector< std::vector< AmiBase::pole_struct >> XI_MAP_VEC;

// std::complex<double> SPECTRAL_PREFACTOR; // this will contain (-i pi)^(Ndeltas) As well as any integral bounds normalization...

	
};



/// These don't explicitly need to exist in ami library - but not harmful 
typedef std::vector<solution_set> solution_set_vec_t;
typedef std::vector< std::vector<solution_set> > solution_set_matrix_t;
typedef std::vector< solution_set_matrix_t > gg_solution_set_matrix_t;

// Functions 
// IO function for external variables

//TODO: The external variable list is hardcoded for what I need.  A user might want something very different. So their external_variable_list might contain completely different things.  Some thought into how to restructure this with epsilon is warranted

// _io
void read_external(std::string filename, external_variable_list &extern_list);

// At one point we thought we would write S,P,R to files along with prefactors etc. Largely this was abbandoned but technically these work.

// _io.cpp
// Test Priority: 3
void read_text_S_solutions(std::string filename, AmiBase::S_t &s_array);
void read_text_P_solutions(std::string eps_filename,std::string alpha_filename, AmiBase::P_t &p_array);
void read_text_R_solutions(std::string eps_filename,std::string alpha_filename, AmiBase::R_t &r_array, int size);
void read_text_R0(std::string alpha_filename, std::string eps_filename, AmiBase::g_prod_t &R0);
double load_prefactor(std::string filename, std::string mul_file, int order);
double load_mul(std::string filename);


// _io.cpp 
void write_S_readable(AmiBase::S_t &s_array);
void write_P_readable(AmiBase::P_t &p_array);
void write_R_readable(AmiBase::R_t &r_array);

void load_solutions(std::string top_directory, solution_set_matrix_t &AMI_MATRIX, int MAX_ORDER, double EREG);

// screen io 

void print_S(int dim, AmiBase::S_t &s_array);
void print_P( int dim, AmiBase::P_t &p_array);
void print_R( int dim, AmiBase::R_t &r_array);
void print_final( int dim, AmiBase::R_t &r_array, AmiBase::P_t &p_array, AmiBase::S_t &s_array);


void print_g_prod_array(AmiBase::g_prod_array_t g_array);
void print_pole_array(AmiBase::pole_array_t g);

void print_g_prod_info(AmiBase::g_prod_t g);
void print_g_struct_info(AmiBase::g_struct g);
void print_epsilon_info(AmiBase::epsilon_t eps);
void print_alpha_info(AmiBase::alpha_t alpha);

void print_pole_struct_info(AmiBase::pole_struct g);

void print_sign_array(AmiBase::sign_array_t signs);
void print_signs(AmiBase::sign_t signs);
void print_Pi( int dim, AmiBase::Pi_t &Pi_array);



// Complete up to here I think 

/// Evaluation
// evaluate functions 
// TODO: Can't recall why there are two. Found it cleaner for some reason to separate the real and imaginary measurements. 
/// Test Priority: 1 - however, this is a very complex function So rather than testing it directly, go into the function and find the functions it depends on and test those. I dont' think we can test the whole thing. 
void evaluate_solutions(std::vector<std::complex<double>> &results, solution_set &AMI, ami_vars_list &ami_eval_vars_list);
void evaluate_solutions(std::vector<double> &Re_results, std::vector<double> &Im_results, solution_set &AMI, ami_vars_list &ami_eval_vars);

// here 

void construct_ami_vars_list(AmiBase::g_prod_t &R0, double prefactor, internal_state &state, external_variable_list &external,ami_vars_list &vars_list);
AmiBase::ami_vars construct_ami_vars(AmiBase::g_prod_t &R0, double prefactor, internal_state &state, ext_vars &external);


// Energy and epsilon functions remain problematic
// testing priority: 4
AmiBase::energy_t construct_energy(AmiBase::g_prod_t &R0, internal_state &state, ext_vars &external);

// could force user to create an energy_class could provide a template for the class require that the class have a function called construct_energy that takes the state and external 


// THis can probably stay
k_vector_t construct_k(AmiBase::alpha_t &alpha, k_vect_list_t &k);

// These should live outside of ami - construct energy should be an external function
std::complex<double> eval_epsilon(hopping_t t, hopping_t tp, NewAmiCalc::k_vector_t k, species_t spin, std::complex<double> mu, double H, AmiBase::disp_type disp);

std::complex<double> eval_epsilon(hopping_t t, k_vector_t k, std::complex<double> mu  , AmiBase::disp_type disp );
std::complex<double> eval_epsilon(hopping_t t, k_vector_t k, species_t spin, std::complex<double> mu, double H, AmiBase::disp_type disp);


// Epsilon stuff

bool not_molecule=1;
std::vector<std::complex<double>> global_hii;

bool density_warning=true;


void read_hf(std::string pgrid_filename, std::string ptoi_filename, std::string sigma_filename);
std::vector<int> ptoi;
std::vector<double> pgrid, sigma_hf;
double hf_mu=0;
double hf_kstep;
double rs;
double max_k_hf=16;

double get_hf_sigma(double kk);
double hf_energy(double kk);

void read_hii(std::string filename, int maxval);


struct evaluation_set{
evaluation_set(){}

evaluation_set(solution_set s,	ext_vars ext){
ext_vars_=ext;
sol_=s;	
}
	
	
ext_vars ext_vars_;
solution_set sol_;
	
int ord, num,gnum;	
	
};


// Spectral Representation
bool find_spectral_pole(int xb, AmiBase::g_struct &g, AmiBase::pole_struct &pole);
void collect_spectral_poles(AmiBase::g_prod_t &gprod, AmiBase::Pi_t &pa);

// std::complex<double> evaluate_spectral(AmiBase::ami_parms &parms, AmiBase::R_t &R_array, AmiBase::P_t &P_array, AmiBase::S_t &S_array, AmiBase::ami_vars &external,AmiBase::g_prod_t &unique_g, AmiBase::Pi_t &Unique_poles, AmiBase::R_ref_t &Rref,AmiBase::ref_eval_t &Eval_list, internal_state &state, ext_vars &ext_var);

std::complex<double> evaluate_spectral(AmiBase::ami_parms &parms, AmiBase::R_t &R_array, AmiBase::P_t &P_array, AmiBase::S_t &S_array, AmiBase::ami_vars &external,AmiBase::g_prod_t &unique_g, AmiBase::Pi_t &Unique_poles, AmiBase::R_ref_t &Rref,AmiBase::ref_eval_t &Eval_list,  internal_state &state, ext_vars &ext_var, double &xi_cut);

std::complex<double> evaluate_spectral_v2(AmiBase::ami_parms &parms, AmiBase::R_t &R_array, AmiBase::P_t &P_array, AmiBase::S_t &S_array, AmiBase::ami_vars &external,AmiBase::g_prod_t &unique_g,  AmiBase::R_ref_t &Rref,AmiBase::ref_eval_t &Eval_list, std::vector<double> &xi_list, AmiBase::Pi_t &Unique_poles, double &xi_cut);

void evaluate_spectral_solutions(std::vector<double> &Re_results, std::vector<double> &Im_results, solution_set &AMI, ami_vars_list &ami_eval_vars, internal_state &state, external_variable_list &ext_var_list, std::vector<double> &xi_list, double &xi_cut);


void evaluate_simple_spectral_solutions(std::vector<double> &Re_results, std::vector<double> &Im_results, solution_set &AMI, ami_vars_list &ami_eval_vars,   std::vector<double> &xi_list, double &xi_cut);
std::complex<double> evaluate_simple_spectral(AmiBase::ami_parms &parms, AmiBase::R_t &R_array, AmiBase::P_t &P_array, AmiBase::S_t &S_array, AmiBase::ami_vars &external,AmiBase::g_prod_t &unique_g,  AmiBase::R_ref_t &Rref,AmiBase::ref_eval_t &Eval_list, std::vector<double> &xi_list);
std::complex<double> eval_spectral_product(std::vector<std::complex<double>> &Ei, std::vector<double> &xi, double &gamma);

void get_pp_comb(int length, std::vector< std::vector<int> > &ppv);
std::vector<int> toBinary(int n,int length);
void remap_external(AmiBase::ami_vars &external,AmiBase::ami_vars &this_external,AmiBase::Pi_t &Unique_poles, std::vector<int> &pp, std::vector<int> &used);

std::complex<double> get_xb_from_pole( AmiBase::pole_struct pole, AmiBase::ami_vars external);

std::complex<double> optimized_spectral_star(AmiBase::ami_parms &parms, AmiBase::SorF_t K, AmiBase::g_prod_t &unique_g, AmiBase::ref_v_t &Rref,AmiBase::ref_v_t &Eval_list, AmiBase::ami_vars external, std::vector<int> &pp);


std::complex<double> eval_spectral_product(AmiBase::g_prod_t &R0,  internal_state &state, ext_vars &external, AmiBase::ami_vars &external_xi);

std::complex<double> soften_divergence(AmiBase::ami_parms &parms, AmiBase::R_t &R_array, AmiBase::P_t &P_array, AmiBase::S_t &S_array, AmiBase::ami_vars &external,AmiBase::g_prod_t &unique_g,  AmiBase::R_ref_t &Rref,AmiBase::ref_eval_t &Eval_list);

bool get_mirror_state(AmiBase::g_prod_t &R0, AmiBase::g_struct &unique_g, internal_state &state, ext_vars &external, internal_state &mirror_state);

void construct_ami_mirror_vars_list(solution_set &AMI, internal_state &state, ext_vars &external, NewAmiCalc::ami_vars_list &vars_list);
void evaluate_mirror_solutions(std::vector<double> &Re_results, std::vector<double> &Im_results, solution_set &AMI, internal_state &state, external_variable_list &external);


// terms spectral functions
std::complex<double> evaluate_terms_simple_spectral(AmiBase::ami_parms &parms, AmiBase::terms ami_terms, AmiBase::ami_vars &external, AmiBase::g_prod_t &unique_g, AmiBase::R_ref_t &Rref,AmiBase::ref_eval_t &Eval_list, std::vector<double> &xi_list);

std::complex<double> evaluate_terms_spectral(AmiBase::ami_parms &parms, AmiBase::terms ami_terms, AmiBase::ami_vars &external, AmiBase::g_prod_t &unique_g, AmiBase::R_ref_t &Rref,AmiBase::ref_eval_t &Eval_list, std::vector<double> &xi_list,  AmiBase::Pi_t &Unique_poles, double &xi_cut);

std::complex<double> evaluate_OPT_spectral_term(AmiBase::ami_parms &parms, AmiBase::terms &ami_terms, AmiBase::ami_vars &external,  AmiBase::g_prod_t &unique_g, AmiBase::R_ref_t &Rref,AmiBase::ref_eval_t &Eval_list, std::vector<int> &pp);

std::complex<double> evaluate_sp_term(AmiBase::ami_parms &parms, AmiBase::term &ami_term, AmiBase::ami_vars &external, std::vector<double> &xi_list, double &xi_cut);


double delta_scale( std::vector<double> xi_list, double xi_cut, std::vector<int> &used);


private:

public:

NewAmiCalc(){
            nac_rand_gen.seed(rd());
            
        }

};
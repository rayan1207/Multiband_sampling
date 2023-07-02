#include "ami_spec.hpp"



void AmiSpec::evaluate_old_spectral_solutions(std::vector<double> &Re_results, std::vector<double> &Im_results, NewAmiCalc::solution_set &AMI, ami_spec_vars_list &ami_eval_vars,   std::vector< std::complex<double>> &xi_list, double &xi_cut){

overflow_detected=false;

Re_results.clear(); 
Im_results.clear();
Re_results.resize(ami_eval_vars.size(),0);
Im_results.resize(ami_eval_vars.size(),0);


for(int i=0; i<ami_eval_vars.size(); i++){
// std::cout<<i<<std::endl;
std::complex<double> calc_result(0,0);

calc_result=evaluate_old_spectral(AMI.ami_parms_, AMI.R0_, AMI.R_, AMI.P_, AMI.S_,  ami_eval_vars[i], AMI.Unique, AMI.Rref, AMI.Eval_list,  xi_list);

// std::cout<<"Overflow??? = "<< overflow_detected <<std::endl;
// std::cout<<calc_result<<std::endl;

Re_results[i]=calc_result.real();
Im_results[i]=calc_result.imag();	
}	

}	


std::complex<double> AmiSpec::evaluate_old_spectral(AmiBase::ami_parms &parms, AmiBase::g_prod_t &R0, AmiBase::R_t &R_array, AmiBase::P_t &P_array, AmiBase::S_t &S_array, ami_spec_vars &external, AmiBase::g_prod_t &unique_g,  AmiBase::R_ref_t &Rref,AmiBase::ref_eval_t &Eval_list, std::vector< std::complex<double>> &xi_list){

// amibase.overflow_detected=false;

AmiBase::ami_vars this_external(external.energy_,external.frequency_, external.BETA_, external.prefactor);	
	

std::vector< std::complex<double>> sigma_list;
get_sigma_list(external.alpha_list,external.k_list_, xi_list, sigma_list);	
	
std::complex<double> A_prod=eval_spectral_product(external.energy_, xi_list, sigma_list);


// std::cout<<"Evaluating with energies"<<std::endl;
for(int i=0; i< this_external.energy_.size(); i++){
 
this_external.energy_[i]=-xi_list[i];	
	// std::cout<<this_external.energy_[i]<<" ";
	
	
}
// std::cout<<"Precision is "<<amibase.precision_cutoff<<std::endl;
std::complex<double> normal=amibase.evaluate(	parms, R_array, P_array, S_array, this_external,unique_g, Rref, Eval_list);

std::complex<double> first=A_prod*normal;

// std::cout<<"Input" <<" "<<normal<<" "<<A_prod*normal<<" "<< amibase.overflow_detected<< std::endl;

// std::complex<double> x0= xi_list[0];
// std::complex<double> x3= xi_list[3];
// xi_list[0]=-x0;
// xi_list[3]=-x3;

// A_prod=eval_spectral_product(external.energy_, xi_list, sigma_list);

// std::cout<<"Evaluating with energies"<<std::endl;
// for(int i=0; i< this_external.energy_.size(); i++){
 
// this_external.energy_[i]=-xi_list[i];	
	// std::cout<<this_external.energy_[i]<<" ";
	
	
// }
// std::cout<<"Precision is "<<amibase.precision_cutoff<<std::endl;
// normal=amibase.evaluate(	parms, R_array, P_array, S_array, this_external,unique_g, Rref, Eval_list);

// std::cout<<"Symm" <<" "<<normal<<" "<<A_prod*normal<<" "<< amibase.overflow_detected<< std::endl;


// std::complex<double> second=A_prod*normal;







if(amibase.overflow_detected){ overflow_detected=true;}

// if(std::abs(A_prod)>1){
  // std::cout<<"A prod is huge "<< A_prod<<" "<< normal<<" "<< A_prod*normal<<std::endl;
// }

// if(std::abs(normal)>1.4){
  // std::cout<<"Normal is huge "<< A_prod<<" "<< normal<<" "<< A_prod*normal<<std::endl;
// }

// if(std::abs(A_prod*normal)>2){
  // std::cout<<"prod is huge "<< A_prod<<" "<< normal<<" "<< A_prod*normal<<std::endl;
// }
	
	
return first;//*second;//A_prod*normal;	
	
	
}

void AmiSpec::get_sigma_list(std::vector< AmiBase::alpha_t> &alpha_list, NewAmiCalc::k_vect_list_t &k_list, std::vector< std::complex<double>> &xi,std::vector< std::complex<double>> &sigma_list){
sigma_list.clear();
if( alpha_list.size()!= xi.size()){
  // std::cout<<"Alpha size "<< alpha_list.size()<<" xi size "<< xi.size()<<std::endl;
	throw std::runtime_error("Alpha  list is wrong size in get_sigma_list - exiting");
}
// std::cout<<"Alpha size "<< alpha_list.size()<<" xi size "<< xi.size()<<std::endl;
// std::cout<<"Use sigma is "<<use_sigma_file<<std::endl;

for(int i=0; i< alpha_list.size(); i++){
  // std::cout<<"Alpha list item "<<i<<std::endl;
NewAmiCalc::k_vector_t this_k=ami.construct_k(alpha_list[i], k_list);
		std::complex<double> this_sigma;
		std::complex<double> this_X=xi[i];
		std::complex<double> defsigma(0,-global_gamma);
		if(!use_sigma_file){
			// std::cout<<"Use sigma file boolean is "<< use_sigma_file<<std::endl;
			this_sigma=defsigma;// would replace this with a working sigma
		}else{
			// std::cout<<"Use sigma file boolean is "<< use_sigma_file<<std::endl;
			// this_sigma=get_sigma(this_k,this_X);
			
			this_sigma=linterp_sigma(this_k,this_X);
		}

// std::cout<<this_sigma<<std::endl;
sigma_list.push_back(this_sigma);
	// std::cout<<"Sigma list has size "<<sigma_list.size()<<std::endl;
	
}
	
  // std::cout<<"At end Sigma list has size "<<sigma_list.size()<<std::endl;
  
	return ;
}



std::complex<double> AmiSpec::eval_spectral_product(std::vector<std::complex<double>> &Ei,  std::vector< std::complex<double>>  &xi, std::vector< std::complex<double>> &sigma_list){
	
std::complex<double> output(1,0);

// the energies are the negative of the epsilon values 
for(int i=0; i< Ei.size(); i++){

std::complex<double> this_sigma=sigma_list[i];


std::complex<double> this_A=-this_sigma.imag()/(std::pow(xi[i] +Ei[i]-this_sigma.real(),2) + std::pow(this_sigma.imag(),2))/M_PI;	

output=output*this_A;
}
	

return output;	
	
}



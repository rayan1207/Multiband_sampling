#include "ami_spec.hpp"



void AmiSpec::construct_ami_spec_vars_list(AmiBase::g_prod_t &R0, double prefactor, NewAmiCalc::internal_state &state, NewAmiCalc::external_variable_list &external, xi_t &xi,ami_spec_vars_list &vars_list){
vars_list.clear();
vars_list.reserve(external.size());
for(int i=0; i<external.size(); i++){
// std::cout<<"Making ami vars on external "<< i<<std::endl;
vars_list.push_back(construct_ami_spec_vars(R0, prefactor, state, external[i], xi));
}	
	
}

void AmiSpec::evaluate_spectral_solutions(std::vector<double> &Re_results, std::vector<double> &Im_results, spec_solution_set &AMI, ami_spec_vars_list &ami_spec_eval_vars,   xi_t &xi_list){

overflow_detected=false;

// std::cout<<"Precision cutoff is "<<amibase.precision_cutoff<<std::endl;

Re_results.clear(); 
Im_results.clear();
Re_results.resize(ami_spec_eval_vars.size(),0);
Im_results.resize(ami_spec_eval_vars.size(),0);



for(int i=0; i<ami_spec_eval_vars.size(); i++){

std::complex<double> calc_result(0,0);

// calc_result=evaluate_simple_spectral(AMI.ami_parms_, AMI.R_, AMI.P_, AMI.S_,  ami_eval_vars[i], AMI.Unique, AMI.Rref, AMI.Eval_list,  xi_list);
calc_result=evaluate_sp_terms(AMI.ami_parms_, AMI.sp_terms_,ami_spec_eval_vars[i]);

// if( calc_result.real()> 1){
// std::cout<<"Returned evaluation "<<calc_result<<std::endl;
// }
// double norm=std::pow(2.0*xi_cut,xi_list.size()

double norm=1.0;//std::pow(2.0*xi_cutoff, ami_spec_eval_vars[i].xi_list_.size());

Re_results[i]=calc_result.real();//*norm;
Im_results[i]=calc_result.imag();//*norm;	
}	
	
// std::cout<<"Eval complete"<<std::endl;	
}	


std::complex<double> AmiSpec::evaluate_sp_terms(AmiBase::ami_parms &parms, AmiSpec::sp_terms &sp_terms, AmiSpec::ami_spec_vars &vars){

std::complex<double> sum(0,0);
/* 
AmiSpec::sp_terms hardcoded;
hardcoded.push_back(sp_terms[0]);
hardcoded.push_back(sp_terms[2]);

hardcoded[0].aprod_[0].alpha_={-1,1};
hardcoded[0].aprod_[0].x_={-1,0};
hardcoded[0].aprod_[0].x_alpha_={0,1};

hardcoded[0].aprod_[1].alpha_={1,0};
hardcoded[0].aprod_[1].x_={1,0};
hardcoded[0].aprod_[1].x_alpha_={0,0};

hardcoded[1].aprod_[0].alpha_={-1,1};
hardcoded[1].aprod_[0].x_={-1,0};
hardcoded[1].aprod_[0].x_alpha_={0,1};

hardcoded[1].aprod_[1].alpha_={1,0};
hardcoded[1].aprod_[1].x_={1,0};
hardcoded[1].aprod_[1].x_alpha_={0,0};


hardcoded[0].ami_term_.p_list[0].eps_={1,0};
hardcoded[0].ami_term_.p_list[0].alpha_={0,1};
hardcoded[0].ami_term_.p_list[0].x_alpha_={0,0};
hardcoded[0].ami_term_.sign=-1;

hardcoded[1].ami_term_.p_list[0].eps_={1,0};
hardcoded[1].ami_term_.p_list[0].alpha_={0,0};
hardcoded[1].ami_term_.p_list[0].x_alpha_={0,-1};
hardcoded[1].ami_term_.sign=1;

// print_sp_terms(hardcoded);


// exit(0);
sp_terms=hardcoded; */
// hardcoded[0].aprod_[0].alpha_={
  
   
 /*  
sp_terms[0].aprod_[0].alpha_={-1,1};
sp_terms[0].aprod_[0].x_={-1,0};
sp_terms[0].aprod_[0].x_alpha_={0,1};

sp_terms[0].aprod_[1].alpha_={1,0};
sp_terms[0].aprod_[1].x_={1,0};
sp_terms[0].aprod_[1].x_alpha_={0,0};

sp_terms[2].aprod_[0].alpha_={-1,1};
sp_terms[2].aprod_[0].x_={-1,0};
sp_terms[2].aprod_[0].x_alpha_={0,1};

sp_terms[2].aprod_[1].alpha_={1,0};
sp_terms[2].aprod_[1].x_={1,0};
sp_terms[2].aprod_[1].x_alpha_={0,0};


sp_terms[0].ami_term_.p_list[0].eps_={1,0};
sp_terms[0].ami_term_.p_list[0].alpha_={0,1};
sp_terms[0].ami_term_.p_list[0].x_alpha_={0,0};
sp_terms[0].ami_term_.sign=-1;

sp_terms[2].ami_term_.p_list[0].eps_={1,0};
sp_terms[2].ami_term_.p_list[0].alpha_={0,0};
sp_terms[2].ami_term_.p_list[0].x_alpha_={0,-1};
sp_terms[2].ami_term_.sign=1;  
  */

 
  
// sp_terms[0].aprod_[0].alpha_={1,0};
// sp_terms[0].aprod_[0].x_={0,-1};
// sp_terms[0].aprod_[0].x_alpha_={0,-1};

// sp_terms[0].aprod_[1].alpha_={-1,1};
// sp_terms[0].aprod_[1].x_={0,1};
// sp_terms[0].aprod_[1].x_alpha_={0,0};

// sp_terms[2].aprod_[0].alpha_={1,0};
// sp_terms[2].aprod_[0].x_={0,-1};
// sp_terms[2].aprod_[0].x_alpha_={0,-1};

// sp_terms[2].aprod_[1].alpha_={-1,1};
// sp_terms[2].aprod_[1].x_={0,1};
// sp_terms[2].aprod_[1].x_alpha_={0,0};


// sp_terms[0].ami_term_.p_list[0].eps_={0,1};
// sp_terms[0].ami_term_.p_list[0].alpha_={0,0};
// sp_terms[0].ami_term_.p_list[0].x_alpha_={0,1};
// sp_terms[0].ami_term_.sign=1;

// sp_terms[2].ami_term_.p_list[0].eps_={0,1};
// sp_terms[2].ami_term_.p_list[0].alpha_={0,1};
// sp_terms[2].ami_term_.p_list[0].x_alpha_={0,0};
// sp_terms[2].ami_term_.sign=-1;  
 
  

// std::cout<<"At evaluation the terms are "<<std::endl;
// print_sp_terms(sp_terms); 
 // exit(0);

for(int i=0; i< sp_terms.size(); i++){

// std::cout<<"On term "<<i<<std::endl;
std::complex<double> sp_result(0,0);
// if(i==2){
sp_result=evaluate_sp_term(parms, sp_terms[i], vars);
// }
// std::cout<<"gave "<<sp_result<<std::endl;

sum+=sp_result;
}
// if(std::abs(sum.real())>.1){
// std::cout<<"Terms sum to "<<sum<<std::endl;
// std::cout<<std::endl;
// }

// std::cout<<"Returning sum "<<sum<<std::endl;

return sum;


}

// This is the functioning one?
std::complex<double> AmiSpec::evaluate_sp_term(AmiBase::ami_parms &parms, AmiSpec::ami_sp_term &sp_term, AmiSpec::ami_spec_vars &vars){
	
	// print_sp_term(sp_term);
	// std::cout<<"Vars prefactor is "<<vars.prefactor<<std::endl;

// vars.frequency_.back()=vars.frequency_.back().real();

std::complex<double> output(0,0);
// if(sp_term.delta_count!=1){return output;}
// std::cout<<"Entering eval"<<std::endl;

// print_sp_term(sp_term);
// exit(0);


// std::cout<<"Entering eval"<<std::endl;
//TODO: Need some sort of warning or catch because this could be really dangerous....
// I HAVE DECIDED TO REMOVE THIS HERE - In case it is better to just regulate with Gamma anyways
// vars.frequency_.back()=vars.frequency_.back().real();

AmiSpec::ami_spec_vars eval_vars=vars;

eval_vars.frequency_.back()=-eval_vars.frequency_.back();


// std::cout<<"Evaluating for frequency "<< eval_vars.frequency_.back()<<std::endl;

AmiBase::ami_vars gprod_external(eval_vars.xi_list_, eval_vars.frequency_, eval_vars.BETA_, eval_vars.prefactor);

// AmiBase::ami_vars fprod_external(vars.xi_list_, vars.frequency_, vars.BETA_, vars.prefactor);

// for(int i=0; i< gprod_external.energy_.size(); i++){

// fprod_external.energy_[i]=-fprod_external.energy_[i];
// fprod_external.frequency_[i]=-fprod_external.frequency_[i];

// }

// AmiBase::ami_vars fprod_external(vars.xi_list_, vars.frequency_, vars.BETA_, vars.prefactor);

// TODO: is this sign swap necessary?
/// For the spectral representation, the xi are converted directly to...
// for(int i=0; i< gprod_external.energy_.size(); i++){

// gprod_external.energy_[i]=-gprod_external.energy_[i];

// }
 // std::cout<<"Entering eval A"<<std::endl;

std::complex<double> A_prod=eval_Aprod(sp_term.aprod_, eval_vars.xi_list_, eval_vars.frequency_, eval_vars.k_list_, eval_vars.MU_);
// std::cout<< "A_prod:  "<<A_prod<<std::endl;

// std::cout<<"Entering eval G"<<std::endl;
std::complex<double> gprod;
gprod=amibase.eval_gprod(parms, sp_term.ami_term_.g_list, gprod_external);

// if ((std::abs(std::real(gprod)) > amibase.precision_cutoff) ||
        // (std::abs(std::imag(gprod)) > amibase.precision_cutoff)) {
      // overflow_detected = true;
    // }


// std::cout<< "gprod:  "<<gprod<<std::endl;
	//std::cout<<"Entering eval F"<<std::endl;
  
  // for(int i=0; i< gprod_external.energy_.size(); i++){

// gprod_external.energy_[i]=-gprod_external.energy_[i];

// }

std::complex<double> fprod;
fprod=amibase.eval_fprod(parms, sp_term.ami_term_.p_list, gprod_external);
// std::cout<< "fprod:  "<<fprod<<std::endl;
std::complex<double> term_val(0,0);
std::complex<double> norm(1,0);

term_val=sp_term.ami_term_.sign*gprod*fprod;

std::complex<double> orig_value=gprod;

std::complex<double> imag(0.,1.0);

// August 26th 2022
// Ok, this is complicated.  Due to the somewhat sloppy implementation, during the evaluation stage what is written symbolically as  1/(omega +epsilon) is replaced with 1/(-omega +epsilon) which is the difference between the integrand for the spectral problem and the non-spectral problem.  This means that the prefactor of the symbolic omega, that was accounted for in the pole determination is the wrong sign.
// so instead of 1/(omega+i0+) -> 1/omega -ipi delta (omega)
// we are actually doing  1/(-omega-i0+) - > 1/(-omega) + i pi delta(omega)
// I'm not 100% sure of this. This may not be the only prefactor issue. Need more test cases. 
// if(!regulate){
//norm=std::pow(-imag*M_PI/(2.0*xi_cutoff), sp_term.delta_count);
norm=std::pow(imag*M_PI/(2.0*xi_cutoff), sp_term.delta_count);
// }
// else{

// SEE changes in here 
norm=norm*get_regdelta(sp_term.ami_term_.sign, sp_term.regprod_, parms, gprod_external);	
	
// }


// std::cout<< "norm:  "<<norm<<std::endl;

// std::cout<<term_val<<" "<<A_prod<<" "<< norm<<std::endl;
// std::cout<<"Eval xi_list size is "<<eval_vars.xi_list_.size()<<std::endl;
// exit(0);
// THE overall prefactor appears to be wrong.  Could be that this power of -1 is not needed because the spectral function absorbs the -1 already. hard to say. 
// Question is, is the sign also wrong for chi_FF AND chi_jj? both have even numbers of terms. 
// REMOVED extra -1 
// what about prefactors?
term_val=term_val*A_prod*norm*(std::pow(-1,eval_vars.xi_list_.size())); // the plus one is a fudge I think 

output+= term_val;

// if(std::abs(output.real())>0.005){
// std::cout<<sp_term.regprod_.size()<<std::endl;  
// std::cout<<output<<" "<<A_prod<<" "<<norm<<" "<<term_val<<" "<< orig_value<<std::endl;  
  
// }

// std::cout<<"Exiting eval "<<output<<std::endl;

// if ((std::abs(std::real(output)) > amibase.precision_cutoff) ||
        // (std::abs(std::imag(output)) > amibase.precision_cutoff)) {
      // overflow_detected = true;
    // }

// exit(0);
return output;



}
//
// TODO: Aug 26th 2022.  The regdelta function DOES work via the lines:


// This fixes the prefactor difference between the delta terms and the regulator - that would otherwise be the same expression. 
// double fix_factor=1.0;
// for(int m=0; m< dprod.size(); m++){
	// fix_factor=fix_factor*dprod[m].alpha_.back();
// }

// amibase.print_g_struct_info(dprod[i]);

// std::complex<double> this_reg_val;
// this_reg_val=amibase.eval_gprod(parms, this_g, external);//*sign;
//std::complex<double> this_val=(imag*M_PI/(2.0*xi_cutoff) +(1.0-exp(-std::pow(this_reg_val,2)))*this_reg_val*fix_factor) ; 
// **************
// HOWEVER.  this only works for second order and bubble diagrams where there is only a single denominator. It seems that the mu
std::complex<double> AmiSpec::get_regdelta(double &sign, delta_prod_t &dprod,  AmiBase::ami_parms &parms ,AmiBase::ami_vars &external){

std::complex<double> norm(1.0,0.0);
std::complex<double> imag(0.,1.0);

// if(dprod.size()==0){return 0;}// this will disable the non-regulator terms 
if(dprod.size()==0){return 1.0;}
for(int i=0; i< dprod.size(); i++){

AmiBase::g_prod_t this_g;
this_g.push_back(dprod[i]);


// This fixes the prefactor difference between the delta terms and the regulator - that would otherwise be the same expression. 
// double fix_factor=1.0;
// for(int m=0; m< dprod.size(); m++){
	// fix_factor=fix_factor*dprod[m].alpha_.back();
// }

// amibase.print_g_struct_info(dprod[i]);

std::complex<double> this_reg_val;
this_reg_val=amibase.eval_gprod(parms, this_g, external);//*sign;

if(this_reg_val.imag()!=0){throw std::runtime_error("The regulator has to be real and is returning non-zero imaginary part. Exiting.");}

double reg=this_reg_val.real();


// if ((std::abs(std::real(this_reg_val)) > amibase.precision_cutoff) ||
        // (std::abs(std::imag(this_reg_val)) > amibase.precision_cutoff)) {
      // overflow_detected = true;
    // }

// std::complex<double> this_val=(-imag*M_PI/(2.0*xi_cutoff) +(1.0-exp(-std::pow(this_reg_val,2)))*this_reg_val*fix_factor) ; 


//*****
// August 26th 2022, swapped the sign due to evaluation at -omega - see above. 

//****
// std::complex<double> this_val=(-imag*M_PI/(2.0*xi_cutoff) +(1.0-exp(-std::pow(this_reg_val,2)))*this_reg_val*fix_factor) ; 

// Old method working version

// std::complex<double> this_val=(imag*M_PI/(2.0*xi_cutoff)+(1.0-exp(-std::pow(this_reg_val,2)))*this_reg_val*fix_factor) ; 

std::complex<double> this_val=expm1(-std::pow(reg,2));//-(1.0-exp(-std::pow(this_reg_val,2))) ;
// std::cout<<this_reg_val<<std::endl;
// std::complex<double> this_val=-(1.0-exp(-std::pow(this_reg_val,2))) ; 



//(-imag*M_PI/(2.0*xi_cutoff) +(1.0-exp(-std::pow(this_reg_val,2)))*this_reg_val*fix_factor) ;//(1.0-exp(-std::pow(this_reg_val,2)))*this_reg_val;// (imag*M_PI/(2.0*xi_cutoff) +(1.0-exp(-std::pow(this_reg_val,2)))*this_reg_val) ;
// (imag*M_PI/(2.0*xi_cutoff) + (1.0-exp(-std::pow(this_reg_val,2)))*this_reg_val);

norm=norm*this_val;

// std::cout<<"In AMI Used this_reg_val of "<< this_reg_val<< " when sign was "<<sign<<" and fix factor was "<< fix_factor<<std::endl;


// if(norm.real()>1){
// std::cout<<"Inside Regdelta "<<i<<" "<<norm<<" "<< this_val <<" "<<this_reg_val<<std::endl;
// }	

}

	
// std::cout<<"Returning regdelta of "<< -norm<<std::endl;
// 2022 03: Removing hard coded -1 in attempt to debug 
return norm; // put overall negative at the end the overall negative 	
	
}


//TODO: Remove this...
std::complex<double> AmiSpec::evaluate_sp_term(AmiBase::ami_parms &parms, AmiSpec::ami_sp_term &sp_term, NewAmiCalc::ext_vars &ev,   AmiBase::ami_vars &external, NewAmiCalc::k_vect_list_t &klist,   xi_t &xi_list){

std::complex<double> output(0,0);


AmiBase::ami_vars gprod_external=external;

// TODO: is this sign swap necessary?
for(int i=0; i< gprod_external.energy_.size(); i++){

gprod_external.energy_[i]=-xi_list[i];

}

// std::cout<<"Entering eval A"<<std::endl;

std::complex<double> A_prod=eval_Aprod(sp_term.aprod_, xi_list, external.frequency_, klist, ev.MU_);

// std::cout<<"Entering eval G"<<std::endl;
std::complex<double> gprod;
gprod=amibase.eval_gprod(parms, sp_term.ami_term_.g_list, gprod_external);

	// std::cout<<"Entering eval F"<<std::endl;
std::complex<double> fprod;
fprod=amibase.eval_fprod(parms, sp_term.ami_term_.p_list, gprod_external);

std::complex<double> term_val(0,0);
std::complex<double> norm(0,0);

term_val=sp_term.ami_term_.sign*gprod*fprod;

std::complex<double> imag(0.,1.0);
norm=std::pow(-imag*M_PI/(2.0*xi_cutoff), sp_term.delta_count);

// std::cout<<term_val<<" "<<A_prod<<" "<< norm<<std::endl;

term_val=term_val*A_prod*norm;

output+= term_val;

// std::cout<<"Exiting eval"<<std::endl;

return output;


}




// TODO: t'=0 is hardcoded 
std::complex<double> AmiSpec::eval_Aprod(A_prod_t &Ap, xi_t &xi, AmiBase::frequency_t &freq, NewAmiCalc::k_vect_list_t &klist, std::complex<double> &mu){

	std::complex<double> output(1,0);

	// A(Sigma, X, E)  : Sigma: self-energy, X: is frequency, E: energy from k_vector
	for(int i=0; i< Ap.size(); i++){


		std::complex<double> this_X=get_X( Ap[i].x_, xi, Ap[i].x_alpha_, freq);
	
		NewAmiCalc::k_vector_t this_k=ami.construct_k(Ap[i].alpha_, klist);
    // std::cout<<"K list is "<<std::endl;
    // std::cout<<"("<<klist[0][0]<<","<<klist[0][1]<<") and ("<<klist[1][0]<<","<<klist[1][1]<<std::endl;
    // std::cout<<"This_k is "<< this_k[0]<<" "<<this_k[1]<<std::endl;
		
  if(Ap[i].eff_stat_==1){
    std::complex<double> this_E=eval_tb(1.,0., this_k, mu);
	
		std::complex<double> this_sigma;
		std::complex<double> defsigma(0,-global_gamma);
		if(!use_sigma_file){
			// std::cout<<"Use sigma file boolean is "<< use_sigma_file<<std::endl;
			this_sigma=defsigma;// would replace this with a working sigma
		}else{
			// std::cout<<"Use sigma file boolean is "<< use_sigma_file<<std::endl;
			// this_sigma=get_sigma(this_k,this_X);
			this_sigma=linterp_sigma(this_k,this_X);
		}

		output=output*A_eval(this_sigma, this_X, this_E);
  }
  
  if(Ap[i].eff_stat_==0){
    
    std::complex<double> this_W=linterp_W(this_k,this_X);
    
    output=output*this_W.imag()/M_PI;
    
    
  }



	}

	// std::cout<<"successsful return of eval_Aprod "<<output<<std::endl;
	return output;

}

// todo: probably don't need to pass A if this function takes in X and E already
std::complex<double> AmiSpec::A_eval( std::complex<double> &sigma, std::complex<double> &X, std::complex<double> &E){

std::complex<double> output(0,0);

// std::cout<<"Evaluating A with "<<sigma.real()<<" "<<sigma.imag()<<" "<< X<<std::endl;

output=-sigma.imag()/M_PI/( std::pow(X-E-sigma.real(),2) +std::pow(sigma.imag(),2) );
// std::cout<<"Evaluated A="<<output<<std::endl;
return output;

}



std::complex<double> AmiSpec::eval_tb(double t, double tp, NewAmiCalc::k_vector_t &k, std::complex<double> &mu){

	std::complex<double> output(0,0);

	for(int i=0; i<k.size();i++){

	output+=-2.0*t*cos(k[i]);

	}

	// uncomment this block if want t' in dispersion
	// double term=-4.0*tp;
	// for(int i=0; i<k.size(); i++){

		// term=term*cos(k[i]);
	// }
	// output+=term;


	output -= mu;

	return output;

}


std::complex<double> AmiSpec::construct_energy(AmiBase::alpha_t &alpha, NewAmiCalc::k_vect_list_t &klist, std::complex<double> &mu){

std::complex<double> result=0;

NewAmiCalc::k_vector_t this_k=ami.construct_k(alpha, klist);
result=eval_tb(1.,0., this_k, mu);


return result;

}

#include "ami_spec.hpp"


AmiSpec::AmiSpec(double xc){

xi_cutoff=xc;

}


AmiSpec::AmiSpec(){

}

std::complex<double> AmiSpec::get_E(AmiBase::energy_t &ei, AmiBase::epsilon_t &eps){
	if(ei.size()!= eps.size()){
	throw std::runtime_error("In A_t the epsilon_t and energy_t do not match in size - exiting");
	}
	std::complex<double> output(0,0);
	for(int i=0; i< ei.size(); i++){

		output+=ei[i]*(double)eps[i];

	}

	return output;
}

void AmiSpec::randomize_xi(xi_t &xi, int length){

// std::cout<<"xi is size "<< xi.size()<<std::endl;
// std::cout<<"sizing xi to new length "<<length<<std::endl;
xi.resize(length,0);
// std::cout<<"Randomizing xi"<<std::endl;
for(int i=0; i<xi.size(); i++){
	// std::cout<<i;
	xi[i]=ami.random_real(-xi_cutoff,xi_cutoff);
	// std::cout<<" "<<xi[i]<<std::endl;
}

// std::cout<<"Exiting"<<std::endl;

}

std::complex<double> AmiSpec::get_X(X_t &Xsym, xi_t &xi, AmiBase::alpha_t &x_alpha_, AmiBase::frequency_t &freq){
	if(Xsym.size()!= xi.size()){
		std::cout<<"("<<Xsym.size()<<" and "<< xi.size()<<")"<<std::endl;
	throw std::runtime_error("In A_t the X_t and xi_t do not match in size  - exiting");
	}
	std::complex<double> output(0,0);
	for(int i=0; i< xi.size(); i++){
    // std::cout<<"Multiplying Xsym_ "<<i<<" "<< Xsym[i]<<" "<<xi[i]<<std::endl;
		output+=(double)Xsym[i]*xi[i];

	}

	for(int i=0; i< freq.size(); i++){
    // std::cout<<"Multiplying x_alpha_ "<<i<<" "<< x_alpha_[i]<<" "<<freq[i]<<std::endl;
		output+=(double)x_alpha_[i]*freq[i];
	}


	return output;
}


//TODO: connor please do this part :)
/*std::complex<double> AmiSpec::get_sigma(NewAmiCalc::k_vector_t &k, std::complex<double> &X){


}*/


void AmiSpec::generate_simple_sp_terms(AmiBase::terms &start_terms, sp_terms &full_sp_terms, AmiBase::g_prod_t &R0){

std::cout<<"Generating terms"<<std::endl;
for(int i=0; i<start_terms.size(); i++){
	std::cout<<"On term "<< i<<std::endl;
sp_terms new_terms;

generate_simple_sp_terms(start_terms[i], new_terms, R0);

full_sp_terms.insert(full_sp_terms.end(), new_terms.begin(), new_terms.end());
	
	
}

return;

}



void AmiSpec::generate_simple_sp_terms(AmiBase::term &start_term, sp_terms &new_sp_terms, AmiBase::g_prod_t &R0){

// Generate the A product from R0

A_prod_t this_Ap;
R0_to_Aprod(R0, this_Ap);

new_sp_terms.resize(1);

new_sp_terms[0].aprod_=this_Ap;
new_sp_terms[0].ami_term_=start_term;


}




void AmiSpec::generate_sp_terms(AmiBase::terms &start_terms, sp_terms &full_sp_terms, AmiBase::g_prod_t &R0){

std::cout<<"Generating terms based on Nterms="<< start_terms.size()<<std::endl;
for(int i=0; i<start_terms.size(); i++){
	std::cout<<"On term "<< i<<std::endl;
sp_terms new_terms;

generate_sp_terms(start_terms[i], new_terms, R0);

full_sp_terms.insert(full_sp_terms.end(), new_terms.begin(), new_terms.end());
	
	
}

return;

}



void AmiSpec::generate_sp_terms(AmiBase::term &start_term, sp_terms &new_sp_terms, AmiBase::g_prod_t &R0){

// first lets identify which G's have external frequencies -
AmiBase::g_prod_t innert;
AmiBase::g_prod_t nert;

// double this_sign=start_term.sign;
ami.amibase.print_term(start_term);

// std::vector<int> pole_indices;
// ami_term.g_list
for(int i=0; i<start_term.g_list.size(); i++){

	if(start_term.g_list[i].alpha_.back()!=0){
		nert.push_back(start_term.g_list[i]);
	}else{
		innert.push_back(start_term.g_list[i]);
	}

}

// Next we want to create the product of ( a +d)*(b+d)*(c+d)... So we want to specify how many deltas are in each one
int max_deltas=nert.size();

std::vector<std::vector<int>> pp_v;
ami.get_pp_comb(max_deltas,pp_v);

sp_terms new_terms;
new_terms.resize(pp_v.size()); // one new term for each combination in pp_v

for(int i=0; i< pp_v.size(); i++){
	std::cout<<"Term:";
	for(int m=0; m< pp_v[i].size(); m++){
	std::cout<<pp_v[i][m]<<" ";

	}
	std::cout<<std::endl;
}

for(int i=0; i< pp_v.size(); i++){
	for(int m=0; m< pp_v[i].size(); m++){

	if(pp_v[i][m]==0){

		new_terms[i].dprod_.push_back(nert[m]);
		new_terms[i].delta_count++;
	}

	if(pp_v[i][m]==1){

	new_terms[i].ami_term_.g_list.push_back(nert[m]);
	}

	}

}

// Generate the A product from R0

A_prod_t this_Ap;
R0_to_Aprod(R0, this_Ap);

// std::cout<<"This ap has size "<<this_Ap<<std::endl;

// now attach the innert G's to each new term

for(int i=0; i< new_terms.size(); i++){

new_terms[i].ami_term_.g_list.insert(new_terms[i].ami_term_.g_list.end(), innert.begin(), innert.end());
new_terms[i].ami_term_.p_list= start_term.p_list;

// get delta prefactor in analytic continuation
// This is hard coded to handle only last variable as external 
double factor=1.0;
for(int m=0; m< new_terms[i].dprod_.size(); m++){
	factor=factor*new_terms[i].dprod_[m].alpha_.back();
}

new_terms[i].ami_term_.sign= start_term.sign*factor;
new_terms[i].aprod_=this_Ap;

}

// no update to the sign needed - the delta_count will tell about prefactors when evaluating (FALSE: Cases with no deltas need the prefactor in case it gets hidden in the Gprod )
// each ami_term.p_list is unchanged from the original
// What is missing is how to handle the A_prod_t - which at least at the start is defined by R0. and the same for every term.



new_sp_terms=new_terms;



}

//TODO: the deltas generated in regulator don't necessitate different xi selections
void AmiSpec::resolve_deltas( sp_terms &sp_terms){

for(int i=0; i< sp_terms.size(); i++){
	std::cout<<"Resolving deltas for term "<< i<<std::endl;
resolve_deltas(sp_terms[i]);
}

}


// TODO: if the prefactor is not +1 or -1, then need to modify the sign for the term to account for the delta rule delta(2*x)=delta(x)/|2|

// TODO: if there is a delta function that is squared, then it is the same as the single delta function but still gives the volume element due to normalization in the larger space.
void AmiSpec::resolve_deltas(ami_sp_term &sp_term){

// std::cout<<"Working on term "<<std::endl;

// print_sp_term(sp_term);

// std::cout<<"Number of deltas is "<<sp_term.dprod_.size()<<std::endl;

if(	sp_term.dprod_.size()==0){ return;}





// look at each delta. create a pole for each, and assign an index_ to it that represents which xi will be replaced.  do this for each and make sure each delta gets a unique xi to replace
// for each pole with respect to its xi values generate the actual pole
// for each pole replace the xi in every Aterm, g_prod, fermi_pole AND other remaining deltas. once done mark the delta for removal
// for removal define an empty integer vec of size delta.size(). initialize to zero. set to 1 for each resolved delta.  check at the end that all the entries are 1.  and then remove all of the deltas. otherwise throw an error.  OR remove only the entries that are resolved.

AmiBase::pole_array_t pv;
pv.resize(sp_term.dprod_.size());

int xi_size=sp_term.dprod_[0].eps_.size();
int delta_size=sp_term.dprod_.size();

std::vector<int> used, assigned;
used.resize(xi_size,0);
assigned.resize(delta_size,0);



// find any xi that appear only once. these are guaranteed to be safe.
std::vector<int> trivial_xi;

for(int i=0; i< xi_size; i++){

	int count=0;
	for(int m=0; m< sp_term.dprod_.size(); m++){
		if(sp_term.dprod_[m].eps_[i]!=0 ){

			count++;

		}
		if(count>1){break;}

	}

	if(count==1){
		trivial_xi.push_back(i);
	}



}

// std::cout<<"Found trivial xi size is "<< trivial_xi.size()<<std::endl;

// exit(0);

if(trivial_xi.size()< delta_size){
	throw std::runtime_error("Not enough trivial xi's to proceed - exiting");
	exit(0);
}


for(int i=0; i< trivial_xi.size(); i++){

	for(int m=0; m< sp_term.dprod_.size(); m++){

		if(assigned[m]==0){
			if(sp_term.dprod_[m].eps_[trivial_xi[i]]!=0 ){
			pv[m].index_=trivial_xi[i];
			used[trivial_xi[i]]=1;
			assigned[m]=1;


			}

		}

	}


}




// for(int i=0; i< xi_size; i++){

// //search through deltas to find a non-zero and then break

	// for(int m=0; m< sp_term.dprod_.size(); m++){
		// if( assigned[m]==0){
		// if(sp_term.dprod_[m].eps_[i]!=0 ){
			// pv[m].index_=i;
			// used[i]=1;
			// assigned[m]=1;
			// break;

		// }
		// }

	// }

// if( count(used.begin(), used.end(),1)==delta_size){
	// break;
// }


// }

if(count(used.begin(), used.end(),1)!= delta_size){
	throw std::runtime_error("Didn't resolve deltas correctly - exiting");
	exit(0);
}

// at this point we should have a list of xi we are going to replace. so generate the poles

 for(int i=0; i< sp_term.dprod_.size(); i++){

	pv[i].eps_.resize( sp_term.dprod_[i].eps_.size());

	pv[i].alpha_.resize( sp_term.dprod_[i].alpha_.size());

	int this_index=pv[i].index_;
	int this_prefactor=sp_term.dprod_[i].eps_[this_index];


	for( int m=0; m< pv[i].alpha_.size(); m++){
		pv[i].alpha_[m]=sp_term.dprod_[i].alpha_[m]*(-this_prefactor);
	}

	for(int m=0; m< pv[i].eps_.size(); m++){
		if(m==this_index){pv[i].eps_[m]=0;}
		else{

			pv[i].eps_[m]= sp_term.dprod_[i].eps_[m]*(-this_prefactor);
			// pv[i].eps_[m]= sp_term.dprod_[i].eps_[m]*(this_prefactor);
		}

	}


}

// now at this point we should have all of the poles accounted for
// std::cout<<std::count(used.begin(), used.end(),1)<<std::endl;

// std::cout<<"BEFORE: Printing poles "<<std::endl;

// for(int i=0; i< pv.size(); i++){

// std::cout<<i<<" x"<<pv[i].index_<<" alpha-";
// print_int_vec(pv[i].alpha_);
// std::cout<<" | eps-";
// print_int_vec(pv[i].eps_);
// std::cout<<std::endl;


// }


// so now replace poles
for(int i=0; i< pv.size(); i++){

replace_xi(i,pv, sp_term);


}

// std::cout<<"AFTER: Printing poles "<<std::endl;

// for(int i=0; i< pv.size(); i++){

// std::cout<<i<<" x"<<pv[i].index_<<" alpha-";
// print_int_vec(pv[i].alpha_);
// std::cout<<" | eps-";
// print_int_vec(pv[i].eps_);
// std::cout<<std::endl;


// }


// we should have a catch in case something goes wrong.
// maybe we should move this info instead of deleting it. 
// sp_term.regprod_=sp_term.dprod_;
sp_term.dprod_.clear();



// for(int i=0; i< sp_term.dprod_.size(); i++){
	// pv[i].eps_=sp_term.dprod_[i].eps_;
	// pv[i].alpha_=sp_term.dprod_[i].alpha_;
// }



// std::cout<<"--------------"<<std::endl<<"Resulting term "<<std::endl;

// print_sp_term(sp_term);


}





void AmiSpec::replace_xi(int i, AmiBase::pole_array_t &pv, ami_sp_term &sp_term){


AmiBase::pole_struct this_pole=pv[i];
// double prefactor=

// Sept 1st: shouldn't need to do this - currently using only trivial poles, so there should be no replacement because of this
// update other poles j>i
for(int j=0; j< pv.size(); j++){
	if(j==i){ continue;}
	update_spec_pole(pv[i], pv[j].alpha_, pv[j].eps_);
	
}

// update the fermi poles in the sp_term
std::cout<<"Updating fermi"<<std::endl;
for(int j=0; j< sp_term.ami_term_.p_list.size(); j++){

int size=sp_term.ami_term_.p_list[j].alpha_.size();
sp_term.ami_term_.p_list[j].x_alpha_.resize(size,0);

	update_spec_pole(pv[i],sp_term.ami_term_.p_list[j].x_alpha_,sp_term.ami_term_.p_list[j].eps_);

}

// update the g_prod
std::cout<<"Updating gprod"<<std::endl;
for(int j=0; j< sp_term.ami_term_.g_list.size(); j++){

	update_spec_pole(pv[i], sp_term.ami_term_.g_list[j].alpha_, sp_term.ami_term_.g_list[j].eps_);

}

// update the Aprod
// note that Aprod is defined by the X_t x_ values not eps_
std::cout<<"Updating Aprod"<<std::endl;
for(int j=0; j< sp_term.aprod_.size(); j++){

	// AmiBase::pole_struct this_pole=pv[i];
	// for(int m=0; m< this_pole.eps_.size(); m++){
		// this_pole.eps_[m]=-this_pole.eps_[m];
	// }
	
	
	// int index=pv[i].index_;
    // int prefactor=sp_term.ami_term_.g_list[j].eps_[index];
	
	// I think there is a prefactor missing 
	// TODO: This MIGHT be necessary due to the sign mismatch between spec and G and F etc. but I think this should be taken care of at evaluation not during construction. 
	// for(int m=0; m< this_pole.alpha_.size(); m++){
		// this_pole.alpha_[m]=-this_pole.alpha_[m];
	// }
	

	update_spec_pole(this_pole, sp_term.aprod_[j].x_alpha_, sp_term.aprod_[j].x_);

}





}



void AmiSpec::spawn_regulator_terms( sp_terms &in_terms){
sp_terms final_terms;
for(int i=0; i< in_terms.size(); i++){
	std::cout<<"Creating regulator terms for term "<< i<<std::endl;
  sp_terms terms_out;
spawn_regulator_terms(in_terms[i],terms_out);

std::cout<<"Regulated terms contains N="<<terms_out.size()<<std::endl;

final_terms.insert(final_terms.end(), terms_out.begin(), terms_out.end());

}

resolve_deltas(final_terms);

in_terms=final_terms;

}
void AmiSpec::spawn_regulator_terms( ami_sp_term &term_in, sp_terms &terms_out){

terms_out.clear();
// terms_out.push_back(term_in); 
 
std::cout<<"Starting to spawn regulator terms"<<std::endl; 
 
 
print_sp_term(term_in); 
 
// First we just put the term here 


ami_sp_term start_term=term_in;


// Just like with sp term part we identify innert G's that here contain no epsilons 

AmiBase::g_prod_t g_orig=start_term.ami_term_.g_list;

AmiBase::g_prod_t innert;
AmiBase::g_prod_t nert;

for(int i=0; i<start_term.ami_term_.g_list.size(); i++){

  bool zeros = std::all_of(start_term.ami_term_.g_list[i].eps_.begin(), start_term.ami_term_.g_list[i].eps_.end(), [](int j) { return j==0; });
	if(!zeros){
		nert.push_back(start_term.ami_term_.g_list[i]);
	}else{
		innert.push_back(start_term.ami_term_.g_list[i]);
	}

}

// Next we want to create the product of ( a +d)*(b+d)*(c+d)... So we want to specify how many deltas are in each one
int max_deltas=nert.size();
std::cout<<nert.size()<<std::endl;

std::vector<std::vector<int>> pp_v;
ami.get_pp_comb(max_deltas,pp_v);

sp_terms new_terms;
new_terms.resize(pp_v.size(),start_term);

// for(int i=0; i< pp_v.size(); i++){
	// std::cout<<"Term:";
	// for(int m=0; m< pp_v[i].size(); m++){
	// std::cout<<pp_v[i][m]<<" ";

	// }
	// std::cout<<std::endl;
// }

for(int i=0; i< pp_v.size(); i++){
	for(int m=0; m< pp_v[i].size(); m++){

	if(pp_v[i][m]==0){

		new_terms[i].dprod_.push_back(nert[m]);
    new_terms[i].regprod_.push_back(nert[m]);
    // new_terms[i].ami_term_.g_list.push_back(nert[m]);
		// The deltas are going to masquerade so not going to count the number of terms 
    // new_terms[i].delta_count++;
	}

  // Since we start by initializing the new term based from start term, there is no change to the number of green's functions. Just have to add the dprod and store the regulator. 
	// if(pp_v[i][m]==1){

	// new_terms[i].ami_term_.g_list.push_back(nert[m]);
	// }

	}

}

// std::cout<<std::endl;
// std::cout<<"Resulting regulated terms has "<<new_terms.size()<<std::endl;
// for(int i=0 ; i< new_terms.size(); i++){
// std::cout<<"On reg term i="<<i<<std::endl;
// print_sp_term(new_terms[i]);

// }

resolve_deltas(new_terms);
// restore original g-product from the starting term 
// TODO: Is this the best thing to do? 
for(int i=0 ; i< new_terms.size(); i++){
  new_terms[i].ami_term_.g_list=g_orig;
}

// std::cout<<std::endl;
// std::cout<<"After resolving deltas regulated terms has "<<new_terms.size()<<std::endl;

// for(int i=0 ; i< new_terms.size(); i++){
// std::cout<<"On reg term i="<<i<<std::endl;
// print_sp_term(new_terms[i]);

// }  
  
  
  
// if (nert.size()!=0){  
// exit(0);  
// }

terms_out=new_terms;
  
}


void AmiSpec::update_spec_pole(AmiBase::pole_struct &source_pole, AmiBase::alpha_t &target_alpha, AmiBase::epsilon_t &target_eps){

std::cout<<"Updating : starting with: alpha - ";
print_int_vec(target_alpha);
std::cout<<" | and eps=";
print_int_vec(target_eps);
std::cout<<std::endl;

std::cout<<"USING pole: alpha - ";
print_int_vec(source_pole.alpha_);
std::cout<<" | and eps=";
print_int_vec(source_pole.eps_);
std::cout<<std::endl;


int index=source_pole.index_;
std::cout<<"Replacing for index "<<index<<std::endl;
if(target_eps[index]==0){ 
std::cout<<"Nothing to do "<<std::endl;
return;} // nothing to do
int prefactor=target_eps[index];

for(int m=0; m< target_alpha.size(); m++){
	target_alpha[m]+=source_pole.alpha_[m]*prefactor;
}

for(int m=0; m< target_eps.size(); m++){
	if(m==index){ target_eps[m]=0;}
	else{
	target_eps[m]+=source_pole.eps_[m]*prefactor;
	}
}




std::cout<<"Update resulted in: alpha - ";
print_int_vec(target_alpha);
std::cout<<" | and eps=";
print_int_vec(target_eps);
std::cout<<std::endl;
std::cout<<"__________"<<std::endl;


	return;


}



// epsilon_i either by index or by value
// int eps_index=-1;
// std::complex<double> eps_val;

// AmiBase::epsilon_t eps_;
// X_t x_;
// AmiBase::alpha_t alpha_;
// AmiBase::species_t species_;


// };
void AmiSpec::R0_to_Aprod(AmiBase::g_prod_t &R0, A_prod_t &Ap){

// R0 is g_prod_t
// each g_prod_t has an epsilon and an alpha

Ap.clear();
Ap.resize(R0.size());

for(int i=0; i< R0.size(); i++){

X_t this_X;
this_X.resize(R0.size(),0);
this_X[i]=1;

Ap[i].alpha_=R0[i].alpha_;
Ap[i].x_alpha_.resize(R0[i].alpha_.size(),0);
Ap[i].species_=R0[i].species_;
Ap[i].eff_stat_=R0[i].eff_stat_;
Ap[i].x_=this_X;
Ap[i].eps_index=i;
// Ap[i].eps_=R0[i].eps_; // not sure this is needed
Ap[i].x_=R0[i].eps_; // DEBUGGING: This forces you to keep R0
}


}


AmiSpec::ami_spec_vars AmiSpec::construct_ami_spec_vars(AmiBase::g_prod_t &R0, double prefactor, NewAmiCalc::internal_state &state, NewAmiCalc::ext_vars &external, xi_t &xi){




//energy_t energy={-4,1,-1};
// std::cout<<"Beta value is "<<external.BETA_<<std::endl;
// std::cout<<"Frequency value is "<< external.external_freq_[0]<<std::endl;
//std::cout<<"construct spec vars 1"<<std::endl;
AmiBase::energy_t energy=ami.construct_energy(R0, state, external);
//std::cout<<"construct spec vars 2"<<std::endl;
NewAmiCalc::k_vect_list_t k_list;

k_list=state.internal_k_list_;
for(int i=0; i<external.external_k_list_.size(); i++){
k_list.push_back(external.external_k_list_[i]);
}


// std::cout<<"Got here"<<std::endl;

// the state 'order_' is actually just the internal k-length - or number of independent variables
AmiBase::frequency_t frequency;
int size=state.internal_k_list_.size();
frequency.reserve(size+1);
// std::cout<<"Frequency size is apparently "<<size<<std::endl;

for(int i=0;i<size;i++){ frequency.push_back(std::complex<double>(0,0));}

// TODO : this doesn't work with multiple external frequencies
// if(external.external_freq_.size()!=0){
// frequency.push_back(external.external_freq_[0]); // some number of external frequencies
// }

// This should address above todo: should allow multiple external frequencies.
// TODO: need a check somewhere that the frequency length matches the alpha length
for(int i=0; i< external.external_freq_.size(); i++){

frequency.push_back(external.external_freq_[i]);

}


if(frequency.size()!= R0[0].alpha_.size()){
	throw std::runtime_error(std::string("Frequency size does not match alpha:")+ std::to_string(frequency.size())+" "+std::to_string(R0[0].alpha_.size()));
}


AmiSpec::ami_spec_vars final_out(energy, frequency, k_list, xi, external.BETA_, external.MU_, prefactor);


final_out.alpha_list.resize(R0.size());
for(int i=0; i< R0.size(); i++){
	for(int j=0; j<R0[i].eps_.size();j++){
		if(R0[i].eps_[j]!=0){
			final_out.alpha_list[j]=R0[i].alpha_;
			
		}
	}
}



// final_out.BETA_=external.BETA_;
// final_out.prefactor=prefactor; //state.prefactor_;
return final_out;


}




// void generate_sp_terms(ami_term &start_term, sp_terms &new_sp_terms);
// void reduce_deltas(ami_sp_term &term);





AmiCalc::ami_vars AmiCalc::construct_ami_vars(AmiCalc::g_prod_t &R0, double prefactor, AmiCalc::internal_state &state, AmiCalc::ext_vars &external){
	
//energy_t energy={-4,1,-1};
// std::cout<<"Beta value is "<<external.BETA_<<std::endl;
// std::cout<<"Frequency value is "<< external.external_freq_[0]<<std::endl;

AmiCalc::energy_t energy=construct_energy(R0, state, external);

// the state 'order_' is actually just the internal k-length - or number of independent variables 
AmiCalc::frequency_t frequency;
frequency.reserve(state.order_+1);

for(int i=0;i<state.order_;i++){ frequency.push_back(std::complex<double>(0,0));}

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
	throw std::runtime_error("Frequency size does not match alpha");
}


AmiCalc::ami_vars final_out(energy, frequency);
final_out.BETA_=external.BETA_;
final_out.prefactor=prefactor; //state.prefactor_;
return final_out;

	
}


std::complex<double> AmiCalc::eval_epsilon(hopping_t t, AmiCalc::k_vector_t k, species_t spin, std::complex<double> mu, double H, disp_type disp){
	
	std::complex<double> output(0,0);
	// print_kvector(k);
	
if(disp==AmiCalc::tb){	
	for(int i=0; i<k.size();i++){
		// std::cout<<"i is "<<i<<std::endl;
	output+=-2.0*t*cos(k[i]);	
		// std::cout<<"Evaluated tb "<< -2.0*t*cos(k[i]) <<" with momentum "<< k[i]<<" and hopping "<< t <<std::endl;
	}
}

// Units are of rydbergs:  We used the atomic Rydberg units. Please see the attachment for the details. In this unit, the length scale is the Bohr radius a_0, and the energy scale is the Rydberg energy e^2/2a_0. Equivalently,  you may set \hbar=1, m=1/2 and e^2/2=1. This is why the dispersion becomes \epsilon_k=k^2/2, and the Coulomb replusion =8*pi/q^2. 
if(disp==AmiCalc::fp){

for(int i=0; i<k.size();i++){
		// std::cout<<"i is "<<i<<" ki is "<<k[i]<<std::endl;
	output+=std::pow(k[i],2);	
		
	}
}

// assuming that spin=0 is up, and spin=1 is down. then spin-1/2, gives -1/2 for up and 1/2 for down.

if(disp==AmiCalc::hf){
	double q=0.0;
	for(int i=0; i<k.size();i++){
		// std::cout<<"i is "<<i<<" ki is "<<k[i]<<std::endl;
	q+=std::pow(k[i],2);	
		
	}
	q=std::sqrt(q);

output=hf_energy(q);
}else{	
	
output+=H*(spin-0.5);	
	
output -= mu;

}

//TODO we don't subtract mu if the hf dispersion is given. it contains its own mu value. 
	// if(std::abs(output)<0.1){
// std::cout<<"Returning epsilon of "<< output<<std::endl;
	// }	
return -output;	
	
	
	
}


std::complex<double> AmiCalc::eval_epsilon(hopping_t t, AmiCalc::k_vector_t k, std::complex<double> mu , disp_type disp){
	
	std::complex<double> output(0,0);
	// print_kvector(k);
	
if(disp==AmiCalc::tb){	
	for(int i=0; i<k.size();i++){
		// std::cout<<"i is "<<i<<std::endl;
	output+=-2.0*t*cos(k[i]);	
		
	}
}
if(disp==AmiCalc::fp){

for(int i=0; i<k.size();i++){
		// std::cout<<"i is "<<i<<std::endl;
	output+=std::pow(k[i],2);	
		
	}
}

if(disp==AmiCalc::hf){
	double q=0.0;
	for(int i=0; i<k.size();i++){
		// std::cout<<"i is "<<i<<" ki is "<<k[i]<<std::endl;
	q+=std::pow(k[i],2);	
		
	}
	q=std::sqrt(q);

output=hf_energy(q);
}else{	
	
output -= mu;

}
	
	return -output;
}

AmiCalc::k_vector_t AmiCalc::construct_k(AmiCalc::alpha_t alpha, AmiCalc::k_vect_list_t &k){
	
AmiCalc::k_vector_t kout(k[0].size(),0);	

for(int j=0; j<kout.size(); j++){
for(int i=0;i<k.size(); i++){
	
kout[j]+= alpha[i]*k[i][j];	
	
}
}

return kout;	
	
}



AmiCalc::energy_t AmiCalc::construct_energy(AmiCalc::g_prod_t &R0, AmiCalc::internal_state &state, AmiCalc::ext_vars &external){

AmiCalc::energy_t result;	
AmiCalc::k_vect_list_t k_list;

k_list=state.internal_k_list_;
// if(external.external_k_vector_.size()!=0){
	
// print_g_prod_info(R0);	


for(int i=0; i<external.external_k_list_.size(); i++){
k_list.push_back(external.external_k_list_[i]);
}

// }

// std::cout<<"Momentum list is "<<std::endl;
// print_array(k_list);

int count=0;

result.resize(R0[0].eps_.size(),0);
for(int i=0; i< R0.size(); i++){
	for(int j=0; j<R0[i].eps_.size();j++){
		if(R0[i].eps_[j]==1){
			// std::cout<<"On energy item "<<j<<std::endl;
			// std::cout<<"t list entry is "<<state.t_list_[j]<<std::endl;
result[j]=eval_epsilon(state.t_list_[j], construct_k(R0[i].alpha_ , k_list) , R0[i].species_, external.MU_, external.H_, state.disp_);
// std::cout<<"energy "<<count<<" "<< result[j].real()<<std::endl;
count++; 
		}
	}
}	

// std::cout<<count<<" "<< R0[0].eps_.size();
if(count != R0[0].eps_.size()){
	std::cout<<count<<" "<< R0[0].eps_.size();
	throw std::runtime_error("Something wrong with epsilon");}
	
return result;
	
}

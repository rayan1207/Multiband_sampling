#include "amigraph.hpp"

// TODO: the effective order is determined by these bosonic lines. 
void AmiGraph::fold_in_bose_alphas(graph_t &g){


boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
// Step through all edges and collect the alpha and epsilon labels of the 'g' graph
for (boost::tie(ei,ei_end) = edges(g); ei != ei_end; ++ei){

if(g[*ei].g_struct_.stat_==AmiBase::Bose){

// if(g[*ei].g_struct_.stat_==AmiBase::Fermi){


if ( !(degree( source(*ei,g),g)==1 || degree( target(*ei,g),g)==1 )){
// std::cout<<"Pushing back"<<std::endl;
g[*ei].g_struct_.stat_=AmiBase::Fermi;
g[*ei].g_struct_.species_=1;
}
}

}

// std::cout<<"Finished with size "<<R0.size()<<std::endl;


}


// this is untested.  Need a labelled graph to test this.
AmiBase::g_prod_t AmiGraph::graph_to_R0(graph_t &g){

AmiBase::g_prod_t R0;




boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
// Step through all edges and collect the alpha and epsilon labels of the 'g' graph
for (boost::tie(ei,ei_end) = edges(g); ei != ei_end; ++ei){

//g[*ei].g_struct_.stat_==AmiBase::Bose){

if(g[*ei].g_struct_.stat_==AmiBase::Fermi){


if ( !(degree( source(*ei,g),g)==1 || degree( target(*ei,g),g)==1 )){
// std::cout<<"Pushing back"<<std::endl;
R0.push_back(g[*ei].g_struct_);
}
}

}

// std::cout<<"Finished with size "<<R0.size()<<std::endl;

return R0;

}

// note the duplication of the graph so that we can delete lines 

// TODO: this won't work for bare phonon(bosonic) propagators 
void AmiGraph::extract_bose_alphas(graph_t g, std::vector<AmiBase::alpha_t> &bose){

graph_t gc=g;	
	
// std::cout<<"Here with num_edges "<< num_edges(gc)<<std::endl;
// std::cout<<std::flush;	
if(num_edges(gc)>1){	
	
bose.clear();

vertex_vector_t extern_vert_list;
edge_vector_t extern_edge_list;

// std::cout<<"Finding externals"<<std::endl;
find_external_vertices(gc, extern_vert_list, extern_edge_list);

// std::cout<<"Delete legs"<<std::endl;
delete_legs(g,extern_vert_list, extern_edge_list);

boost::graph_traits<graph_t>::edge_iterator ei, ei_end;

		for( boost::tie(ei,ei_end)=edges(gc); ei!=ei_end; ++ei){

			if( g[*ei].g_struct_.stat_== AmiBase::Bose){
			bose.push_back(g[*ei].g_struct_.alpha_);
				
			}
		
		}
}

}

void AmiGraph::graph_to_R0(graph_t &g, AmiBase::g_prod_t &R0){

R0.clear();

// std::cout<<"Converting graph to R0_"<<std::endl;

// print_all_edge_info(g);


boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
// Step through all edges and collect the alpha and epsilon labels of the 'g' graph
for (boost::tie(ei,ei_end) = edges(g); ei != ei_end; ++ei){

//g[*ei].g_struct_.stat_==AmiBase::Bose){

if(g[*ei].g_struct_.stat_==AmiBase::Fermi){

// special case for bare lines - since boost graph 'degree' function returns 1 if you close a line on itself.	
/* 	if(num_edges(g)==1 && graph_type==AmiBase::density){
		std::cout<<"Triggered on "<<std::endl;
	g[*ei].g_struct_.species_=(AmiBase::species_t)g[*ei].spin;	
	R0.push_back(g[*ei].g_struct_);	
	} */

// if(graph_type==AmiBase::ENERGY){

// if ( !(degree( source(*ei,g),g)==1)){
	
// g[*ei].g_struct_.species_=(AmiBase::species_t)g[*ei].spin;
	// std::cout<<"G had spin "<<g[*ei].spin<<std::endl;
	// R0.push_back(g[*ei].g_struct_);
// }

// } else	
	
if( graph_type== AmiBase::ENERGY){
	
	if (num_edges(g)==1){
	g[*ei].g_struct_.species_=(AmiBase::species_t)g[*ei].spin;
		R0.push_back(g[*ei].g_struct_);	
		
	}else 
		if ( !(degree( source(*ei,g),g)==1)){
		g[*ei].g_struct_.species_=(AmiBase::species_t)g[*ei].spin;
		R0.push_back(g[*ei].g_struct_);
		
	}
	
}


if( graph_type== AmiBase::density || graph_type== AmiBase::Greens || graph_type== AmiBase::DOS  ){

g[*ei].g_struct_.species_=(AmiBase::species_t)g[*ei].spin;
	// std::cout<<"G had spin "<<g[*ei].spin<<std::endl;
	R0.push_back(g[*ei].g_struct_);
	
	
}else{


	if ( !(degree( source(*ei,g),g)==1 || degree( target(*ei,g),g)==1 )){
	// std::cout<<"Pushing back"<<std::endl;
	g[*ei].g_struct_.species_=(AmiBase::species_t)g[*ei].spin;
	// std::cout<<"G had spin "<<g[*ei].spin<<std::endl;
	R0.push_back(g[*ei].g_struct_);
	}
	}

}

}

if(bose_alphas_in_R0){
  
  // put the bosonic lines in here 
 for (boost::tie(ei,ei_end) = edges(g); ei != ei_end; ++ei){

if(g[*ei].g_struct_.stat_==AmiBase::Bose){

if ( !(degree( source(*ei,g),g)==1 || degree( target(*ei,g),g)==1 )){
	
	g[*ei].g_struct_.species_=(AmiBase::species_t)g[*ei].spin;
	g[*ei].g_struct_.eff_stat_=AmiBase::Bose;
  // std::transform(g[*ei].g_struct_.alpha_.begin(), g[*ei].g_struct_.alpha_.end(), g[*ei].g_struct_.alpha_.begin(),
               // std::bind(std::multiplies<int>(), std::placeholders::_1, -1));
  // g[*ei].g_struct_.alpha_=-g[*ei].g_struct_.alpha_;
	R0.push_back(g[*ei].g_struct_);
	}
	

}  
  
  
}

//now reset the epsilons to include the extra green's functions. 

reset_epsilons(R0);

}



// std::cout<<"Finished with size "<<R0.size()<<std::endl;

// return R0;

}

// internal_state(int k_length, int dim){
// internal_k_list_.assign(k_length, std::vector< double>(dim,0.0));	
// dim_=dim;
// order_=k_length;
// }



void AmiGraph::randomize_state(NewAmiCalc::internal_state &state, int k_length, int gsize){

// std::cout<<"State dimension is "<< state.dim_<<std::endl;
// if(state.dim_<2){ throw std::runtime_error("State dimension may not have been set"); }
	
state.internal_k_list_.resize(k_length, std::vector< double>(state.dim_,0.0));
state.internal_freq_size_=state.internal_k_list_.size(); //internal_k_list_.resize(k_length, std::vector< double>(state.dim_,0.0));
// std::cout<<"Set k_length to size "<<k_length<<std::endl;
// std::cout<<"Using tp="<<state.tp_list_[0]<<std::endl;
// reset all hoppings to 1
state.t_list_.clear();
state.t_list_.resize(gsize, 1);
double tpval=state.tp_list_[0];	
state.tp_list_.clear();
state.tp_list_.resize(gsize,tpval);
//
state.order_=k_length;
// state.prefactor_=prefactor;
	
for(int i=0; i< state.internal_k_list_.size(); i++){

for(int j=0; j< state.dim_ ; j++){
	// std::cout<<"On dim and internal k "<< j<< " "<<i<<std::endl;
	
	state.internal_k_list_[i][j]=random_real(state.mink_, state.maxk_);	
	// std::cout<<"Set state to "<< state.internal_k_list_[i][j]<< std::endl;
}


}	
	
		
}


void AmiGraph::randomize_state_spherical(NewAmiCalc::internal_state &state, int k_length, int gsize){
// std::cout<<"Random spherical state"<<std::endl;
// std::cout<<"State dimension is "<< state.dim_<<std::endl;
if(state.dim_!=3){ throw std::runtime_error("Spherical RNG only for 3D"); }
	
state.internal_k_list_.resize(k_length, std::vector< double>(state.dim_,0.0));

// std::cout<<"Set k_length to size "<<k_length<<std::endl;

// reset all hoppings to 1
state.t_list_.clear();
state.t_list_.resize(gsize, 1);
double tpval=state.tp_list_[0];	
state.tp_list_.clear();
state.tp_list_.resize(gsize,tpval);
//
state.order_=k_length;
// state.prefactor_=prefactor;

double Jk=1;
	
for(int i=0; i< state.internal_k_list_.size(); i++){
	
	// std::cout<<"generating kvector for i="<<i<<std::endl;
	
double mag=random_real(state.mink_,state.maxk_);
double theta=random_real(0.0, M_PI);
double phi=random_real(0.0,2.0*M_PI);

// std::cout<<"Generated r theta phi "<< mag<<" "<<theta<<" "<<phi<<std::endl;


state.internal_k_list_[i][0]=mag*std::sin(theta)*std::cos(phi);
state.internal_k_list_[i][1]=mag*std::sin(theta)*std::sin(phi);
state.internal_k_list_[i][2]=mag*std::cos(theta);
	
// kx=mag*std::sin(theta)*std::cos(phi);
// ky=mag*std::sin(theta)*std::sin(phi);
// kz=mag*std::cos(theta);	

Jk=Jk*std::pow(mag,2)*std::sin(theta);//*state.maxk_*2.0*std::pow(M_PI,2);
	
	
/* 
for(int j=0; j< state.dim_ ; j++){
	// std::cout<<"On dim and internal k "<< j<< " "<<i<<std::endl;
	
	state.internal_k_list_[i][j]=random_real(state.mink_, state.maxk_);	
	// std::cout<<"Set state to "<< state.internal_k_list_[i][j]<< std::endl;
} */


}	

state.Jk_=Jk;

// std::cout<<"Resulting Jk is "<<Jk<<std::endl;

	
		
}




void AmiGraph::randomize_state(NewAmiCalc::internal_state &state, int k_length){

// std::cout<<"State dimension is "<< state.dim_<<std::endl;
if(state.dim_<2){ throw std::runtime_error("State dimension may not have been set"); }
	
state.internal_k_list_.resize(k_length, std::vector< double>(state.dim_,0.0));

// reset all hoppings to 1
state.t_list_.clear();
state.t_list_.resize(2*k_length-1, 1);
double tpval=state.tp_list_[0];	
state.tp_list_.clear();
state.tp_list_.resize(2*k_length-1,tpval);	
//
state.order_=k_length;
// state.prefactor_=prefactor;
	
for(int i=0; i< state.internal_k_list_.size(); i++){

for(int j=0; j< state.dim_ ; j++){
	
	state.internal_k_list_[i][j]=random_real(state.mink_, state.maxk_);	
	
}


}	
	
		
}

void AmiGraph::shift_state_pi(NewAmiCalc::internal_state &state, int k_length){

//std::cout<<"State dimension is "<< state.dim_<<std::endl;
if(state.dim_!=2){ throw std::runtime_error("State dimension may not have been set"); }
	
state.internal_k_list_.resize(k_length, std::vector< double>(state.dim_,0.0));	
state.order_=k_length;
// state.prefactor_=prefactor;
	
for(int i=0; i< state.internal_k_list_.size(); i++){

for(int j=0; j< state.dim_ ; j++){
	
	state.internal_k_list_[i][j]+=M_PI/2.;	
	
}


}	
	
		
}


void AmiGraph::zero_state(NewAmiCalc::internal_state &state, int k_length){

//std::cout<<"State dimension is "<< state.dim_<<std::endl;
if(state.dim_!=2){ throw std::runtime_error("State dimension may not have been set"); }
	
state.internal_k_list_.resize(k_length, std::vector< double>(state.dim_,0.0));	
state.order_=k_length;
// state.prefactor_=prefactor;
	
for(int i=0; i< state.internal_k_list_.size(); i++){

for(int j=0; j< state.dim_ ; j++){
	
	state.internal_k_list_[i][j]=0.0;	
	
}


}	
	
		
}


void AmiGraph::shift_state_list(NewAmiCalc::internal_state &state,  NewAmiCalc::k_vect_list_t &shift){

//std::cout<<"State dimension is "<< state.dim_<<std::endl;
 if(state.internal_k_list_.size()>shift.size()){ throw std::runtime_error("Can't shift, not enough entries."); }

for(int m=0; m< state.internal_k_list_.size(); m++){
for(int j=0; j< state.dim_ ; j++){
	
	
	// std::cout<<"Shifting state value "<<state.internal_k_list_[m][j]<<" by "<<shift[m][j]<<std::endl;
	state.internal_k_list_[m][j]+=shift[m][j];	
	
}
}
	
	
		
}

// void AmiGraph::shift_state_withsign(NewAmiCalc::internal_state &state, int k_index, NewAmiCalc::k_vector_t &shift){

// // std::cout<<"State dimension is "<< state.dim_<<std::endl;
// if(state.dim_!=shift.size()){ throw std::runtime_error("State dimension may not have been set"); }


// for(int j=0; j< state.dim_ ; j++){
	
	// state.internal_k_list_[k_index][j]+=shift[j];	
	
// }

	
	
		
// }


void AmiGraph::shift_state(NewAmiCalc::internal_state &state, int k_index, NewAmiCalc::k_vector_t &shift){

//std::cout<<"State dimension is "<< state.dim_<<std::endl;
if(state.dim_!=shift.size()){ throw std::runtime_error("State dimension may not have been set"); }


for(int j=0; j< state.dim_ ; j++){
	
	state.internal_k_list_[k_index][j]+=shift[j];	
	
}

	
	
		
}



void AmiGraph::sobol_randomize_state(NewAmiCalc::internal_state &state, int k_length){

//std::cout<<"State dimension is "<< state.dim_<<std::endl;
if(state.dim_!=2){ throw std::runtime_error("State dimension may not have been set"); }
	
state.internal_k_list_.resize(k_length, std::vector< double>(state.dim_,0.0));	
state.order_=k_length;
// state.prefactor_=prefactor;

std::vector<double> sobol=sobol_random_real();
	
for(int j=0; j< state.dim_ ; j++){

		
	
for(int i=0; i< state.internal_k_list_.size(); i++){


	state.internal_k_list_[i][j]=sobol[i+j*state.internal_k_list_.size()]*2*M_PI;	
	
}


}	
	
		
}



void AmiGraph::construct_state(NewAmiCalc::internal_state &state){
	
//for(int i=0; i< state.internal_k_list_.size(); i++){



state.internal_k_list_[0]={0,0};
state.internal_k_list_[1]={M_PI,M_PI/3.0};

// for(int j=0; j< state.dim_ ; j++){
	
	// state.internal_k_list_[i][j]=random_real(2*M_PI);	
	
// }


//}	
	
		
}

void AmiGraph::print(AmiBase::ami_vars vars){

std::cout<<"Beta is "<< vars.BETA_<<std::endl;

std::cout<<"Printing energy"<<std::endl;	
for (int i=0; i< vars.energy_.size(); i++){

std::cout<<"("<<vars.energy_[i].real()<<","<<vars.energy_[i].imag()<<") ";

}	
std::cout<<std::endl;

std::cout<<"Printing frequencies"<<std::endl;	
for (int i=0; i< vars.frequency_.size(); i++){

std::cout<<"("<<vars.frequency_[i].real()<<","<<vars.frequency_[i].imag()<<") ";

}	
std::cout<<std::endl;	
	
	
}




void AmiGraph::print_state(NewAmiCalc::internal_state &state){
	
for(int i=0; i< state.internal_k_list_.size(); i++){
std::cout<< "[ ";
for(int j=0; j< state.dim_ ; j++){
	
	std::cout<< state.internal_k_list_[i][j]<<" ";	
	
}
std::cout<< "]"<<std::endl;

}	
	
}


void AmiGraph::print_kvector(NewAmiCalc::k_vector_t &k){
	
	std::cout<<"K_vector has size "<<k.size()<<std::endl;
for(int i=0; i< k.size(); i++){
std::cout<< "[ ";
	
	std::cout<< k[i]<<" ";	
	

std::cout<< "]"<<std::endl;

}	
	
}




void AmiGraph::set_amiparms(AmiBase::ami_parms &amiparms, int &graph_order, double &ereg){
	
amiparms.N_INT_=graph_order;
amiparms.E_REG_=ereg;	
		
	
}








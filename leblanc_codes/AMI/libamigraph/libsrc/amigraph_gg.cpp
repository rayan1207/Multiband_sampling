#include "amigraph.hpp"


void AmiGraph::ggm_construct_ami_sol(gg_matrix_t &ggm, double ereg, int mpi_rank){
	
if(mpi_rank==0){	
std::cout<<"-------------"<<std::endl;
std::cout<<"Constructing AMI solutions for graph_groups"<<std::endl;
std::cout<<"-------------"<<std::endl;
}
for(int ord=0; ord< ggm.size(); ord++){
	if(mpi_rank==0){	
	std::cout<<"Working on order "<< ord<<" with size "<<ggm[ord].size()<<std::endl;
	}
for(int i=0; i< ggm[ord].size(); i++){	
if(mpi_rank==0){	
std::cout<<"On graph ord num "<< ord<<" "<<i<<std::endl;
}
// print_all_edge_info(ggm[ord][i].graph_vec[0]);
	
gg_construct_ami_sol( ggm[ord][i], ereg);	
	

}
}
		// if(mpi_rank==0){	
		// std::cout<<"-------------"<<std::endl;
		// std::cout<<"Complete"<<std::endl;
		// std::cout<<"-------------"<<std::endl;	}
}

void AmiGraph::ggm_reset_epsilons(gg_matrix_t &ggm){

	if(mpi_rank==0){	
std::cout<<"-------------"<<std::endl;
std::cout<<"Fixing epsilons in ggm for spectral representation"<<std::endl;
std::cout<<"-------------"<<std::endl;
}
for(int ord=0; ord< ggm.size(); ord++){
	if(mpi_rank==0){	
	std::cout<<"Working on order "<< ord<<" with size "<<ggm[ord].size()<<std::endl;
	}
for(int i=0; i< ggm[ord].size(); i++){	



for(int m=0; m< ggm[ord][i].graph_vec.size(); m++){
	std::cout<<"Before the order is "<< graph_order(ggm[ord][i].graph_vec[m])<<std::endl;
reset_epsilons( ggm[ord][i].graph_vec[m]);	
	std::cout<<"After the order is "<< graph_order(ggm[ord][i].graph_vec[m])<<std::endl;
}

}
}
		
	
}


void AmiGraph::ggm_fold_in_bose_alphas(gg_matrix_t &ggm){
// std::cout<<"In fold function "<<std::endl;
// std::cout<<"ggm size is "<<ggm.size()<<std::endl;
// std::cout<<"mpi rank is "<< mpi_rank<<std::endl;
	if(mpi_rank==0){	
std::cout<<"-------------"<<std::endl;
std::cout<<"Jamming the bose alphas into the calculation"<<std::endl;
std::cout<<"-------------"<<std::endl;
}
for(int ord=0; ord< ggm.size(); ord++){
	if(mpi_rank==0){	
	std::cout<<"Working on order "<< ord<<" with size "<<ggm[ord].size()<<std::endl;
	}
for(int i=0; i< ggm[ord].size(); i++){	



for(int m=0; m< ggm[ord][i].graph_vec.size(); m++){
std::cout<<"Before the order is "<< graph_order(ggm[ord][i].graph_vec[m])<<std::endl;
fold_in_bose_alphas(ggm[ord][i].graph_vec[m]);
	std::cout<<"After the order is "<< graph_order(ggm[ord][i].graph_vec[m])<<std::endl;
// reset_epsilons( ggm[ord][i].graph_vec[m]);	
	
}

}
}
		
	
  
 // exit(0); 
}


void AmiGraph::ggm_syk_epsilons(gg_matrix_t &ggm){

	if(mpi_rank==0){	
std::cout<<"-------------"<<std::endl;
std::cout<<"Fixing epsilons for SYK in ggm for spectral representation"<<std::endl;
std::cout<<"-------------"<<std::endl;
}
for(int ord=0; ord< ggm.size(); ord++){
	if(mpi_rank==0){	
	std::cout<<"Working on order "<< ord<<" with size "<<ggm[ord].size()<<std::endl;
	}
for(int i=0; i< ggm[ord].size(); i++){	



for(int m=0; m< ggm[ord][i].graph_vec.size(); m++){
	
syk_epsilons( ggm[ord][i].graph_vec[m]);	
	
}

}
}
		
	
}



void AmiGraph::ggm_generate_sp_terms(gg_matrix_t &ggm, int mpi_rank){

std::cout<<"On rank "<<mpi_rank <<" with ggm size "<<ggm.size()<<std::endl;

	if(mpi_rank==0){	
std::cout<<"-------------"<<std::endl;
std::cout<<"Generating SP terms from existing terms "<<std::endl;
std::cout<<"-------------"<<std::endl;
}
for(int ord=0; ord< ggm.size(); ord++){
	if(mpi_rank==0){	
	std::cout<<"Working on order "<< ord<<" with size "<<ggm[ord].size()<<std::endl;
	}
for(int i=0; i< ggm[ord].size(); i++){	


ggm[ord][i].sss_vec.resize(ggm[ord][i].graph_vec.size());
for(int m=0; m< ggm[ord][i].graph_vec.size(); m++){
	
// reset_epsilons( ggm[ord][i].graph_vec[m]);	
std::cout<<"In graph GG - generating terms "<<std::endl;
ggm[ord][i].sss_vec[m].ami_terms_=ggm[ord][i].ss_vec[m].ami_terms_;
ggm[ord][i].sss_vec[m].R0_=ggm[ord][i].ss_vec[m].R0_;
// ami.print_g_prod_info(ggm[ord][i].sss_vec[m].R0_);
ggm[ord][i].sss_vec[m].ami_parms_=ggm[ord][i].ss_vec[m].ami_parms_;
ggm[ord][i].sss_vec[m].prefactor_=ggm[ord][i].ss_vec[m].prefactor_;
sp.generate_sp_terms(ggm[ord][i].sss_vec[m].ami_terms_, ggm[ord][i].sss_vec[m].sp_terms_, ggm[ord][i].sss_vec[m].R0_);
std::cout<<"In graph GG - done "<<std::endl;
std::cout<<"There are Nterms= "<<ggm[ord][i].sss_vec[m].sp_terms_.size()<<std::endl;
// sp.print_sp_terms(ggm[ord][i].sss_vec[m].sp_terms_);
sp.resolve_deltas(ggm[ord][i].sss_vec[m].sp_terms_);
// sp.print_sp_terms(ggm[ord][i].sss_vec[m].sp_terms_);
	std::cout<<"In graph GG - resolve deltas done  "<<std::endl;
  sp.spawn_regulator_terms(ggm[ord][i].sss_vec[m].sp_terms_);
  std::cout<<"In graph GG - spawned regulator terms complete"<<std::endl;
  // sp.print_sp_terms(ggm[ord][i].sss_vec[m].sp_terms_);
  
  // std::cout<<"BEGIN DEBUGGING"<<std::endl;
  // for(int t=0; t< ggm[ord][i].sss_vec[m].sp_terms_.size(); t++){
  
// if(ggm[ord][i].sss_vec[m].sp_terms_[t].delta_count==1){
// sp.print_sp_term(ggm[ord][i].sss_vec[m].sp_terms_[t]);
// }  

    
    
  // }
  // exit(0);
  
}

}
}
		
	
}


void AmiGraph::ggm_generate_simple_sp_terms(gg_matrix_t &ggm, int mpi_rank){

std::cout<<"On rank "<<mpi_rank <<" with ggm size "<<ggm.size()<<std::endl;

	if(mpi_rank==0){	
std::cout<<"-------------"<<std::endl;
std::cout<<"Generating SP terms from existing terms "<<std::endl;
std::cout<<"-------------"<<std::endl;
}
for(int ord=0; ord< ggm.size(); ord++){
	if(mpi_rank==0){	
	std::cout<<"Working on order "<< ord<<" with size "<<ggm[ord].size()<<std::endl;
	}
for(int i=0; i< ggm[ord].size(); i++){	


ggm[ord][i].sss_vec.resize(ggm[ord][i].graph_vec.size());
for(int m=0; m< ggm[ord][i].graph_vec.size(); m++){
	
// reset_epsilons( ggm[ord][i].graph_vec[m]);	
std::cout<<"In graph GG - generating terms "<<std::endl;
ggm[ord][i].sss_vec[m].ami_terms_=ggm[ord][i].ss_vec[m].ami_terms_;
ggm[ord][i].sss_vec[m].R0_=ggm[ord][i].ss_vec[m].R0_;
ggm[ord][i].sss_vec[m].ami_parms_=ggm[ord][i].ss_vec[m].ami_parms_;
ggm[ord][i].sss_vec[m].prefactor_=ggm[ord][i].ss_vec[m].prefactor_;
sp.generate_simple_sp_terms(ggm[ord][i].sss_vec[m].ami_terms_, ggm[ord][i].sss_vec[m].sp_terms_, ggm[ord][i].sss_vec[m].R0_);
std::cout<<"In graph GG - done "<<std::endl;
// NO DELTAS to resolve in simple case 
// sp.resolve_deltas(ggm[ord][i].sss_vec[m].sp_terms_);
	// std::cout<<"In graph GG - resolve deltas done  "<<std::endl;
}

}
}
		
	
}







void AmiGraph::ggm_construct_ami_term_sol(gg_matrix_t &ggm, double ereg, int mpi_rank){
	
if(mpi_rank==0){	
std::cout<<"-------------"<<std::endl;
std::cout<<"Constructing AMI solutions term-by-term for graph_groups"<<std::endl;
std::cout<<"-------------"<<std::endl;
}
for(int ord=0; ord< ggm.size(); ord++){
	if(mpi_rank==0){	
	std::cout<<"Working on order "<< ord<<" with size "<<ggm[ord].size()<<std::endl;
	}
for(int i=0; i< ggm[ord].size(); i++){	
if(mpi_rank==0){	
std::cout<<"On graph ord num "<< ord<<" "<<i<<std::endl;
}
// print_all_edge_info(ggm[ord][i].graph_vec[0]);
	
gg_construct_ami_term_sol( ggm[ord][i], ereg);	
	

}
}
		// if(mpi_rank==0){	
		// std::cout<<"-------------"<<std::endl;
		// std::cout<<"Complete"<<std::endl;
		// std::cout<<"-------------"<<std::endl;	}
}





void AmiGraph::ggm_opt_der(gg_matrix_t &ggm, int mpi_rank){
	
if(mpi_rank==0){	
std::cout<<"-------------"<<std::endl;
std::cout<<"Taking derivatives of OPT solutions for graph_groups"<<std::endl;
std::cout<<"-------------"<<std::endl;
}
for(int ord=0; ord< ggm.size(); ord++){
	if(mpi_rank==0){	
	std::cout<<"Working on order "<< ord<<" with size "<<ggm[ord].size()<<std::endl;
	}
for(int num=0; num< ggm[ord].size(); num++){	
if(mpi_rank==0){	
std::cout<<"On graph ord num "<< ord<<" "<<num<<std::endl;
}
// print_all_edge_info(ggm[ord][i].graph_vec[0]);
	for(int m=0; m< ggm[ord][num].ss_vec.size(); m++){
		
		ami.amibase.derivative_opt(ggm[ord][num].ss_vec[m].Unique, ggm[ord][num].ss_vec[m].Rref, ggm[ord][num].ss_vec[m].Eval_list);
		
		
	}
	

}
}




}

void AmiGraph::ggm_optimize_ami(gg_matrix_t &ggm, int mpi_rank){
	
if(mpi_rank==0){	
std::cout<<"-------------"<<std::endl;
std::cout<<"Factorizing AMI solutions for graph_groups"<<std::endl;
std::cout<<"-------------"<<std::endl;
}
for(int ord=0; ord< ggm.size(); ord++){
	if(mpi_rank==0){	
	std::cout<<"Working on order "<< ord<<" with size "<<ggm[ord].size()<<std::endl;
	}
for(int num=0; num< ggm[ord].size(); num++){	
if(mpi_rank==0){	
std::cout<<"On graph ord num "<< ord<<" "<<num<<std::endl;
}
// print_all_edge_info(ggm[ord][i].graph_vec[0]);
	for(int m=0; m< ggm[ord][num].ss_vec.size(); m++){
		
		if(ggm[ord][num].ss_vec[m].R_.size()!=0){
		ami.amibase.factorize_Rn(ggm[ord][num].ss_vec[m].R_.back(), ggm[ord][num].ss_vec[m].Unique, ggm[ord][num].ss_vec[m].Rref, ggm[ord][num].ss_vec[m].Eval_list);
		}else{
			ami.amibase.factorize_terms(ggm[ord][num].ss_vec[m].ami_terms_,  ggm[ord][num].ss_vec[m].Unique, ggm[ord][num].ss_vec[m].Rref, ggm[ord][num].ss_vec[m].Eval_list);
		}
		
		
		// TODO: This causes some issues - which is bad given that it is optional?
		// std::cout<<"Collecting poles"<<std::endl;
		ami.collect_spectral_poles(ggm[ord][num].ss_vec[m].Unique, ggm[ord][num].ss_vec[m].Unique_poles);
		// std::cout<<"Complete"<<std::endl;
		
	}
	

}
}
		// if(mpi_rank==0){	
		// std::cout<<"-------------"<<std::endl;
		// std::cout<<"Complete"<<std::endl;
		// std::cout<<"-------------"<<std::endl;	}
}



void AmiGraph::ggm_reduce_ami_terms(gg_matrix_t &ggm,double ereg, int mpi_rank, int max_try){

if(mpi_rank==0){	
std::cout<<"-------------"<<std::endl;
std::cout<<"Minimizing AMI solutions for graph_groups"<<std::endl;
std::cout<<"-------------"<<std::endl;
}

std::cout<<ggm.size()<<std::endl;
for(int ord=1; ord< ggm.size(); ord++){
	if(mpi_rank==0){	
	std::cout<<"Working on order "<< ord<<" with size "<<ggm[ord].size()<<std::endl;
	}
for(int i=0; i< ggm[ord].size(); i++){	
if(mpi_rank==0){	
std::cout<<"On graph ord num "<< ord<<" "<<i<<std::endl;
}
// print_all_edge_info(ggm[ord][i].graph_vec[0]);


// first try swapping alphas 
	for(int m=0; m< ggm[ord][i].graph_vec.size(); m++){
bool use_terms;
if(ggm[ord][i].ss_vec[m].R_.size()==0){	
	use_terms=true;
}else{use_terms=false;}
	
int terms;
if(!use_terms){		
construct_ami_sol( ggm[ord][i].graph_vec[m],ggm[ord][i].ss_vec[m], ereg);		
terms=ggm[ord][i].ss_vec[m].R_.back().size();}
if(use_terms){
construct_ami_terms_sol( ggm[ord][i].graph_vec[m],ggm[ord][i].ss_vec[m], ereg);	
terms=ggm[ord][i].ss_vec[m].ami_terms_.size();

}


int nint=alpha_size(ggm[ord][i].graph_vec[m])-1;
// graph_t temp_g, store_g;
// NewAmiCalc::solution_set ss, store_ss;

// bool replace=false;
// int reset=0;
for(int count=0; count<max_try; count++){
// reset++;
// if(reset>10){break;}
graph_t temp_g;
NewAmiCalc::solution_set ss;

temp_g=ggm[ord][i].graph_vec[m];

int first=random_int(0,nint-1);
int second=random_int(0,nint-1);
swap_alphas(temp_g, first,second);

if(!use_terms){
construct_ami_sol( temp_g, ss, ereg);	
}
if(use_terms){
	construct_ami_terms_sol( temp_g, ss, ereg);
}

if(!use_terms){	
if(ss.R_.back().size()<terms){
// reset=0;
terms=ss.R_.back().size();
ggm[ord][i].graph_vec[m]=temp_g;
ggm[ord][i].ss_vec[m]=ss;

}	
}else{

if(ss.ami_terms_.size()<terms){
// reset=0;
terms=ss.ami_terms_.size();
ggm[ord][i].graph_vec[m]=temp_g;
ggm[ord][i].ss_vec[m]=ss;

}	



}	


	
}

		
	
	} // end alpha swap 
/* 	
	// since we cannot label counter-term diagrams check if these are counterterms and if they are we will end 
bool no_counterterms=true;
for(int m=0; m< ggm[ord][i].graph_vec.size(); m++){
	if(ggm[ord][i].graph_vec[m][boost::graph_bundle].ct_count!=0 ||ggm[ord][i].graph_vec[m][boost::graph_bundle].sigma_ct_count!=0){no_counterterms=false;   }
}

if(no_counterterms){	
	bool allmin=true;
	int min=-1;
	for(int m=0; m< ggm[ord][i].graph_vec.size(); m++){
		
		if (m==0){ min=ggm[ord][i].ss_vec[m].R_.back().size();
		continue;
		}
		
		if( ggm[ord][i].ss_vec[m].R_.back().size()< min){
			allmin=false;
			min=ggm[ord][i].ss_vec[m].R_.back().size();
		}
				
	}
	
	// if they are all the same size then great. if not then for each graph with solution not the size of min, try to reduce the number of terms further by relabelling etc 
	if(!allmin){
		
		for(int m=0; m< ggm[ord][i].graph_vec.size(); m++){
			
			if( ggm[ord][i].ss_vec[m].R_.back().size() > min){
				
				bool cont=true;
				int count=0;
				do{
				int start=ggm[ord][i].ss_vec[m].R_.back().size();
				
				graph_t temp_g;
				NewAmiCalc::solution_set ss;

				temp_g=ggm[ord][i].graph_vec[m];
				
				bool success=false;
				repeated_labelling(temp_g, success);
				if(success!=true){throw std::runtime_error("Labelling failed in AMI term reduction");}
				
				construct_ami_sol( temp_g, ss, ereg);	
				
				int end=ss.R_.back().size();
				
				if (end==min){cont=false;}
				
				if(end<start){
					ggm[ord][i].ss_vec[m]=ss;
					ggm[ord][i].graph_vec[m]=temp_g;
				}
				
				count++;
				
				if(count>3){cont=false;}
				
				}while(cont);
				
				
			}
			
			
		}
		
	} // end !allmin 

} */	
	
	
	

}
}

// std::cout<<"Got here"<<std::endl;


}	


void AmiGraph::ggm_maximize_ami_terms(gg_matrix_t &ggm,double ereg, int mpi_rank, int max_try){

if(mpi_rank==0){	
std::cout<<"-------------"<<std::endl;
std::cout<<"Maximizing AMI solutions for graph_groups"<<std::endl;
std::cout<<"-------------"<<std::endl;
}

std::cout<<ggm.size()<<std::endl;
for(int ord=1; ord< ggm.size(); ord++){
	if(mpi_rank==0){	
	std::cout<<"Working on order "<< ord<<" with size "<<ggm[ord].size()<<std::endl;
	}
for(int i=0; i< ggm[ord].size(); i++){	
if(mpi_rank==0){	
std::cout<<"On graph ord num "<< ord<<" "<<i<<std::endl;
}
// print_all_edge_info(ggm[ord][i].graph_vec[0]);


// first try swapping alphas 
	for(int m=0; m< ggm[ord][i].graph_vec.size(); m++){
bool use_terms;
if(ggm[ord][i].ss_vec[m].R_.size()==0){	
	use_terms=true;
}else{use_terms=false;}
	
int terms;
if(!use_terms){		
construct_ami_sol( ggm[ord][i].graph_vec[m],ggm[ord][i].ss_vec[m], ereg);		
terms=ggm[ord][i].ss_vec[m].R_.back().size();}
if(use_terms){
construct_ami_terms_sol( ggm[ord][i].graph_vec[m],ggm[ord][i].ss_vec[m], ereg);	
terms=ggm[ord][i].ss_vec[m].ami_terms_.size();

}


int nint=alpha_size(ggm[ord][i].graph_vec[m])-1;
// graph_t temp_g, store_g;
// NewAmiCalc::solution_set ss, store_ss;

// bool replace=false;
// int reset=0;
for(int count=0; count<max_try; count++){
// reset++;
// if(reset>10){break;}
graph_t temp_g;
NewAmiCalc::solution_set ss;

temp_g=ggm[ord][i].graph_vec[m];

int first=random_int(0,nint-1);
int second=random_int(0,nint-1);
swap_alphas(temp_g, first,second);

if(!use_terms){
construct_ami_sol( temp_g, ss, ereg);	
}
if(use_terms){
	construct_ami_terms_sol( temp_g, ss, ereg);
}

if(!use_terms){	
if(ss.R_.back().size()>terms){
// reset=0;
terms=ss.R_.back().size();
ggm[ord][i].graph_vec[m]=temp_g;
ggm[ord][i].ss_vec[m]=ss;

}	
}else{

if(ss.ami_terms_.size()>terms){
// reset=0;
terms=ss.ami_terms_.size();
ggm[ord][i].graph_vec[m]=temp_g;
ggm[ord][i].ss_vec[m]=ss;

}	



}	


	
}

		
	
	} // end alpha swap 
	
	
	
	

}
}

std::cout<<"Got here"<<std::endl;


}	


/* void AmiGraph::ggm_reduce_ami_terms(gg_matrix_t &ggm,double ereg, int mpi_rank, int max_try){

if(mpi_rank==0){	
std::cout<<"-------------"<<std::endl;
std::cout<<"Minimizing AMI solutions for graph_groups"<<std::endl;
std::cout<<"-------------"<<std::endl;
}
for(int ord=0; ord< ggm.size(); ord++){
	if(mpi_rank==0){	
	std::cout<<"Working on order "<< ord<<" with size "<<ggm[ord].size()<<std::endl;
	}
for(int i=0; i< ggm[ord].size(); i++){	
if(mpi_rank==0){	
std::cout<<"On graph ord num "<< ord<<" "<<i<<std::endl;
}
// print_all_edge_info(ggm[ord][i].graph_vec[0]);
	for(int m=0; m< ggm[ord][i].graph_vec.size(); m++){
		
		
construct_ami_sol( ggm[ord][i].graph_vec[m],ggm[ord][i].ss_vec[m], ereg);	
int terms=ggm[ord][i].ss_vec[m].R_.back().size();
int nint=alpha_size(ggm[ord][i].graph_vec[m])-1;
// graph_t temp_g, store_g;
// NewAmiCalc::solution_set ss, store_ss;

// bool replace=false;
// int reset=0;
for(int count=0; count<max_try; count++){
// reset++;
// if(reset>10){break;}
graph_t temp_g;
NewAmiCalc::solution_set ss;

temp_g=ggm[ord][i].graph_vec[m];

int first=random_int(0,nint-1);
int second=random_int(0,nint-1);
swap_alphas(temp_g, first,second);

construct_ami_sol( temp_g, ss, ereg);	
	
if(ss.R_.back().size()<terms){
// reset=0;
terms=ss.R_.back().size();
ggm[ord][i].graph_vec[m]=temp_g;
ggm[ord][i].ss_vec[m]=ss;

}	
	
	
}

	
	
	
	
	}

}
}



}	 */



void AmiGraph::ggm_eo_split(gg_matrix_t &ggm, gg_matrix_t &ggm_even, gg_matrix_t &ggm_odd){

ggm_even.clear();
ggm_odd.clear();

ggm_even.resize(ggm.size());
ggm_odd.resize(ggm.size());

gg_vec_t even_temp, odd_temp;

graph_group ggeven_temp, ggodd_temp;

// irr -> even, red-> odd

for(int ord=0; ord< ggm.size(); ord++){	

if(ord==0){
ggm_odd[ord]=ggm[ord];	
}else{
even_temp.clear();
odd_temp.clear();
	
for(int group=0; group< ggm[ord].size(); group++){

for(int graph=0; graph< ggm[ord][group].graph_vec.size(); graph++){
std::cout<<"On graph "<<ord<<" "<< group<<" "<<graph<<std::endl;

// print_all_edge_info(ggm[ord][group].graph_vec[graph]);
if(is_hubW_graph(ggm[ord][group].graph_vec[graph])){
	std::cout<<"Is odd"<<std::endl;
	ggodd_temp.graph_vec.push_back(ggm[ord][group].graph_vec[graph]);
	
}else{
	std::cout<<"Is even"<<std::endl;
	ggeven_temp.graph_vec.push_back(ggm[ord][group].graph_vec[graph]);
}



}

if(ggeven_temp.graph_vec.size()!=0){
even_temp.push_back(ggeven_temp);
ggeven_temp.graph_vec.clear();	
}

if(ggodd_temp.graph_vec.size()!=0){
odd_temp.push_back(ggodd_temp);
ggodd_temp.graph_vec.clear();	
}



}

ggm_odd[ord]=odd_temp;	
odd_temp.clear();

ggm_even[ord]=even_temp;
even_temp.clear();

}


}// end ord



}




void AmiGraph::ggm_1P_split(gg_matrix_t &ggm, gg_matrix_t &ggm_irr, gg_matrix_t &ggm_red){

ggm_irr.clear();
ggm_red.clear();

ggm_irr.resize(ggm.size());
ggm_red.resize(ggm.size());

gg_vec_t irr_temp, red_temp;

graph_group ggirr_temp, ggred_temp;


for(int ord=0; ord< ggm.size(); ord++){	

if(ord==0){
ggm_irr[ord]=ggm[ord];	
}else{
irr_temp.clear();
red_temp.clear();
	
for(int group=0; group< ggm[ord].size(); group++){

for(int graph=0; graph< ggm[ord][group].graph_vec.size(); graph++){
std::cout<<"On graph "<<ord<<" "<< group<<" "<<graph<<std::endl;

// print_all_edge_info(ggm[ord][group].graph_vec[graph]);
if(is_1P_reducible_graph(ggm[ord][group].graph_vec[graph])){
	std::cout<<"Is reducible"<<std::endl;
	ggred_temp.graph_vec.push_back(ggm[ord][group].graph_vec[graph]);
	
}else{
	std::cout<<"Is irreducible"<<std::endl;
	ggirr_temp.graph_vec.push_back(ggm[ord][group].graph_vec[graph]);
}



}

if(ggirr_temp.graph_vec.size()!=0){
irr_temp.push_back(ggirr_temp);
ggirr_temp.graph_vec.clear();	
}

if(ggred_temp.graph_vec.size()!=0){
red_temp.push_back(ggred_temp);
ggred_temp.graph_vec.clear();	
}



}

ggm_red[ord]=red_temp;	
red_temp.clear();

ggm_irr[ord]=irr_temp;
irr_temp.clear();

}


}// end ord



}


void AmiGraph::ggm_split(gg_matrix_t &ggm, gg_matrix_t &ggm_uu, gg_matrix_t &ggm_ud){

ggm_uu.clear();
ggm_ud.clear();

ggm_uu.resize(ggm.size());
ggm_ud.resize(ggm.size());

gg_vec_t ud_temp, uu_temp;

graph_group gguu_temp, ggud_temp;


for(int ord=0; ord< ggm.size(); ord++){	

ud_temp.clear();
uu_temp.clear();
	
for(int group=0; group< ggm[ord].size(); group++){

for(int graph=0; graph< ggm[ord][group].graph_vec.size(); graph++){

int spins=count_spins(ggm[ord][group].graph_vec[graph]);

// std::cout<<"Spin count was "<< spins<<std::endl;

if(spins==2){
gguu_temp.graph_vec.push_back(ggm[ord][group].graph_vec[graph]);
}
if(spins==0){
ggud_temp.graph_vec.push_back(ggm[ord][group].graph_vec[graph]);	
}


}

if(gguu_temp.graph_vec.size()!=0){
uu_temp.push_back(gguu_temp);
gguu_temp.graph_vec.clear();	
}

if(ggud_temp.graph_vec.size()!=0){
ud_temp.push_back(ggud_temp);
ggud_temp.graph_vec.clear();	
}



}

ggm_uu[ord]=uu_temp;	
uu_temp.clear();

ggm_ud[ord]=ud_temp;
ud_temp.clear();

}



}

int AmiGraph::count_spins(graph_t &g){
	
graph_t gc;
gc=g;
int up=0;
int dn=0;
int spin_count;

// first delete external legs of gc
vertex_vector_t extern_vect_list;
edge_vector_t extern_edge_list;

// this incorrectly picks up closed graphs I think
find_external_vertices(gc, extern_vect_list, extern_edge_list);


// std::cout<<"Deleting external legs of size "<< extern_edge_list.size()<<std::endl;


delete_legs(gc,extern_vect_list, extern_edge_list);	

	
boost::graph_traits<graph_t>::edge_iterator ei, ei_end;	

for(boost::tie(ei,ei_end)=edges(gc); ei!=ei_end; ++ei){	

if( g[*ei].g_struct_.stat_==AmiBase::Fermi ){
if(g[*ei].spin==AmiGraph::up){ up++;}
if(g[*ei].spin==AmiGraph::dn){ dn++;}
	
}
	
}
std::cout<<" up and dn are "<< up<<" "<<dn<<std::endl;

spin_count=up-dn;


return spin_count;	
}



void AmiGraph::ggm_to_amim(gg_matrix_t &ggm, NewAmiCalc::solution_set_matrix_t &AMI_MATRIX){

AMI_MATRIX.clear();

AMI_MATRIX.resize(ggm.size());

for (int ord=0; ord< ggm.size(); ord++){

for(int i=0; i< ggm[ord].size(); i++){

for(int m=0; m< ggm[ord][i].ss_vec.size(); m++){

AMI_MATRIX[ord].push_back(ggm[ord][i].ss_vec[m]);

}
}

}	
	
	
}


// TODO: This doesn't work for AmiBase::density graphs 
void AmiGraph::ggm_label(gg_matrix_t &ggm, int min){
	

for(int ord=min; ord< ggm.size(); ord++){
	std::cout<<"Producing ggm labels at order "<<ord<< " with size "<< ggm.size()<<std::endl;
	
	for(int group=0; group<ggm[ord].size(); group++){
		// if(group%5==0){std::cout<<"On group "<< group<<std::endl;}
	for(int graph=0; graph< ggm[ord][group].graph_vec.size(); graph++){
	
	// systematic_vertex_label(ggm[ord][group].graph_vec[graph]);
print_all_edge_info(	ggm[ord][group].graph_vec[graph]);
	 std::cout<<"Labelling graph ord"<<ord<<" group "<<group<<" num "<<graph<<std::endl;
bool success=true;

if(graph_type==AmiBase::density || graph_type==AmiBase::Greens|| graph_type==AmiBase::DOS|| graph_type==AmiBase::ENERGY || graph_type==AmiBase::FORCE || graph_type==AmiBase::Pi_ppud){

// NOT SURE why we wouldn't just try the systematic case anyways. 
	
sys_label(ggm[ord][group].graph_vec[graph], success);

if(!success){ repeated_labelling(ggm[ord][group].graph_vec[graph], success);}

}else{
repeated_labelling(ggm[ord][group].graph_vec[graph], success);
}



// sys_label(ggm[ord][group].graph_vec[graph], success);
if(!success){throw std::runtime_error("Failed to label a graph");}	
		
	}
	}
}	
	
	
}

void AmiGraph::ggm_group_eff_order_CT(gg_matrix_t &ggm, int max_bub_ct, int max_sigma_ct){
	
gg_matrix_t ggmc;	
	
// gg_vec_t ggv;

ggmc.resize(ggm.size());

for(int ord=0; ord<ggm.size(); ord++){
	
gg_vec_t ggv;

for(int num_b=0; num_b<=max_bub_ct; num_b++){
	for(int num_s=0; num_s<=max_sigma_ct; num_s++){
	graph_group gg;
			for(int group=0; group<ggm[ord].size(); group++){
		// std::cout<<"On group "<<group<<std::endl;

			for(int num=0; num<ggm[ord][group].graph_vec.size(); num++){

			int order=ord;//graph_order(ggm[ord][group].graph_vec[num]);
			int bubble_ct=ggm[ord][group].graph_vec[num][boost::graph_bundle].ct_count;
			int sigma_ct=ggm[ord][group].graph_vec[num][boost::graph_bundle].sigma_ct_count;
			
			if(bubble_ct==num_b && sigma_ct==num_s){
			
			gg.graph_vec.push_back(ggm[ord][group].graph_vec[num]);
				
			}


			
			}}
	
		
if(gg.graph_vec.size()!=0){
		ggv.push_back(gg);	
		}
	
	}
	

	
}


	
			if(ggv.size()!=0){
			ggmc[ord]=ggv;
			}
	
	
}
	
	
ggm=ggmc; 
	
}

void AmiGraph::ggm_filter_max_eff_order(gg_matrix_t &ggm, int max, int max_bub_ct, int max_sigma_ct, int ct_size){
	
gg_matrix_t ggmc;	

// gg_vec_t ggv;
// graph_group gg; 

ggmc.resize(ggm.size());

// filter 

for(int ord=0; ord < ggm.size(); ord++){
	
	// std::cout<<"On ord "<<ord<<std::endl;

gg_vec_t ggv;
	for(int group=0; group<ggm[ord].size(); group++){
// std::cout<<"On group "<<group<<std::endl;

	graph_group gg; 

	for(int num=0; num<ggm[ord][group].graph_vec.size(); num++){
		
		std::cout<<  ord <<" "<<ggm[ord][group].graph_vec[num][boost::graph_bundle].ct_count<<" "<<ggm[ord][group].graph_vec[num][boost::graph_bundle].sigma_ct_count<<std::endl;

	int order=ord;//graph_order(ggm[ord][group].graph_vec[num]);
	int bubble_ct=ggm[ord][group].graph_vec[num][boost::graph_bundle].ct_count;
	int sigma_ct=ggm[ord][group].graph_vec[num][boost::graph_bundle].sigma_ct_count;
	
	int sum=order+bubble_ct+ct_size*sigma_ct+1;

	// std::cout<<"Sum is "<<sum<<std::endl;
	// std::cout<< bool(sum<=max)<<" "<< bool(bubble_ct<=max_bub_ct)<<" "<< bool(sigma_ct<= max_sigma_ct)<<std::endl;

	if(sum<=max && bubble_ct<=max_bub_ct && sigma_ct<= max_sigma_ct){

	// std::cout<<"Pushing"<<std::endl;
	gg.graph_vec.push_back(ggm[ord][group].graph_vec[num]);
		
	}

		
	}

		if(gg.graph_vec.size()!=0){
		ggv.push_back(gg);	
		}

	}

			if(ggv.size()!=0){
			ggmc[ord]=ggv;
			}

		}

ggm=ggmc; 	
	
}



void AmiGraph::ggm_filter_hfspinsusc(gg_matrix_t &ggm){
	
gg_matrix_t ggmc;	

gg_vec_t ggv;
// graph_group gg; 

ggmc.resize(ggm.size());

// filter 

for(int ord=0; ord < ggm.size(); ord++){

gg_vec_t ggv;
	for(int group=0; group<ggm[ord].size(); group++){


	graph_group gg; 

	for(int num=0; num<ggm[ord][group].graph_vec.size(); num++){

	if(!disconnected_loop(ggm[ord][group].graph_vec[num])){

	gg.graph_vec.push_back(ggm[ord][group].graph_vec[num]);
		
	}

		
	}

		if(gg.graph_vec.size()!=0){
		ggv.push_back(gg);	
		}

	}

			if(ggv.size()!=0){
			ggmc[ord]=ggv;
			}

		}

ggm=ggmc; 	
	
}

void AmiGraph::ggm_filter_tp(gg_matrix_t &ggm){

/* 
for(int ord=0; ord< ggm.size(); ord++){

 for(int group=0; group<ggm[ord].size(); group++){

	std::vector<graph_t> new_graph_vec;
	for(int num=0; num<ggm[ord][group].graph_vec.size();num++){

	if(count_tadpoles(ggm[ord][group].graph_vec[num])==0){
		
		new_graph_vec.push_back(ggm[ord][group].graph_vec[num]);
		
	}
	

	}
std::cout<<"New graph vector contains "<<new_graph_vec.size()<<" graphs while original was "<<ggm[ord][group].graph_vec.size()<<std::endl;
 }
}	 */
	
	  
gg_matrix_t ggmc;	

// gg_vec_t ggv;
// graph_group gg; 

ggmc.resize(ggm.size());

// filter 

for(int ord=0; ord < ggm.size(); ord++){

gg_vec_t ggv;
	for(int group=0; group<ggm[ord].size(); group++){

 // std::cout<<ord<<" "<<group<<std::endl;
	graph_group gg; 

	for(int num=0; num<ggm[ord][group].graph_vec.size(); num++){
// std::cout<<ord<<" "<<group<<" "<<num<<std::endl;
// std::cout<<"Tadpole count is "<<count_tadpoles(ggm[ord][group].graph_vec[num])<<std::endl;

	if(count_tadpoles(ggm[ord][group].graph_vec[num])==0){

	gg.graph_vec.push_back(ggm[ord][group].graph_vec[num]);
		
	}
// std::cout<<"Here on numb "<<num<<std::endl;
		
	}
	// std::cout<<"End numb"<<std::endl;

		if(gg.graph_vec.size()!=0){
			
		ggv.push_back(gg);	
		}

	}
	// std::cout<<"End group"<<std::endl;

			if(ggv.size()!=0){
				// std::cout<<"Assigning vector with size "<< ggv.size()<<" on ord="<<ord<<std::endl;
				// ggmc[ord].resize(ggv.size());
			ggmc[ord]=ggv;
			}

		}

ggm=ggmc; 
// return ggmc;
	// std::cout<<"Exiting"<<std::endl;
	// exit(0);  
}



void AmiGraph::ggm_filter_fock(gg_matrix_t &ggm){
	
gg_matrix_t ggmc;	

// gg_vec_t ggv;
// graph_group gg; 

ggmc.resize(ggm.size());

// filter 

for(int ord=0; ord < ggm.size(); ord++){

gg_vec_t ggv;
	for(int group=0; group<ggm[ord].size(); group++){


	graph_group gg; 

	for(int num=0; num<ggm[ord][group].graph_vec.size(); num++){

	if(!has_fock_insertion(ggm[ord][group].graph_vec[num])){

	gg.graph_vec.push_back(ggm[ord][group].graph_vec[num]);
		
	}

		
	}

		if(gg.graph_vec.size()!=0){
		ggv.push_back(gg);	
		}

	}

			if(ggv.size()!=0){
			ggmc[ord]=ggv;
			}

		}

ggm=ggmc; 	
	
}




void AmiGraph::ggm_to_ggamim(gg_matrix_t &ggm, NewAmiCalc::gg_solution_set_matrix_t &GG_AMI_MATRIX, int min){

GG_AMI_MATRIX.clear();

GG_AMI_MATRIX.resize(ggm.size());

for(int ord=min; ord< ggm.size(); ord++){
	GG_AMI_MATRIX[ord].resize(ggm[ord].size());
	// std::cout<<"On ord "<<ord<<" GGM size is "<<ggm[ord].size();
	for(int group=0; group<ggm[ord].size(); group++){
	
	GG_AMI_MATRIX[ord][group]=ggm[ord][group].ss_vec;
	
	}
}
	
	
}


void AmiGraph::ggm_to_sp_ggamim(gg_matrix_t &ggm, AmiSpec::sp_gg_solution_set_matrix_t &sp_GG_AMI_MATRIX, int min){

sp_GG_AMI_MATRIX.clear();

sp_GG_AMI_MATRIX.resize(ggm.size());

for(int ord=min; ord< ggm.size(); ord++){
	sp_GG_AMI_MATRIX[ord].resize(ggm[ord].size());
	// std::cout<<"On ord "<<ord<<" GGM size is "<<ggm[ord].size();
	for(int group=0; group<ggm[ord].size(); group++){
	
	sp_GG_AMI_MATRIX[ord][group]=ggm[ord][group].sss_vec;
	
	}
}
	
	
}


void AmiGraph::merge_ggm(gg_matrix_t &ggm1, gg_matrix_t &ggm2){

if(ggm1.size()!= ggm2.size()){throw std::runtime_error("Cannot merge ggms of different size - aborting");}

for(int ord=0; ord<ggm1.size(); ord++){

for(int group=0; group< ggm2[ord].size(); group++){
	if(ggm2[ord][group].graph_vec.size()!=0){ // don't merge empty groups 
ggm1[ord].push_back(ggm2[ord][group]);	
	}
}
	
}
	
	
	
}

void AmiGraph::ct_construct_ami_sol(graph_t &g, NewAmiCalc::solution_set &ss, double ereg){



graph_to_R0(g, ss.R0_);

extract_bose_alphas(g, ss.bose_alphas_);

ss.loops_=count_fermi_loops(g)-g[boost::graph_bundle].ct_count;


int graph_ord=graph_order(g)-g[boost::graph_bundle].ct_count;

ss.ct_count_=g[boost::graph_bundle].ct_count;
ss.sigma_ct_count_=g[boost::graph_bundle].sigma_ct_count;
ss.order_shift=0;

ss.prefactor_=get_prefactor(g, graph_ord);

AmiBase::ami_parms amiparms(graph_ord, ereg, ami_parameters.TYPE_);

amiparms.N_INT_=alpha_size(g)-1;

// if((graph_type == AmiBase::Pi_phuu) || (graph_type == AmiBase::Pi_phud) || graph_type== AmiBase::density || graph_type==AmiBase::doubleocc || ami_parameters.TYPE_==AmiBase::Pi_ppud || ami_parameters.TYPE_==AmiBase::Pi_ppuu ){
// amiparms.N_INT_++; // increase the integration order by one for bubble graphs.
// also increase it by one if it is a closed density graph.  
// }

ss.ami_parms_=amiparms;

// std::cout<<"construct"<<std::endl;
// ami.print_g_prod_info( gg.ss_vec[m].R0_);
ami.amibase.construct(ss.ami_parms_, ss.R0_, ss.R_, ss.P_, ss.S_);
// std::cout<<"constructed"<<std::endl;
	
	
	
return; 	
	
}

int AmiGraph::alpha_size(graph_t &g){

boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
int output=-1;
for( boost::tie(ei,ei_end)=edges(g); ei!=ei_end; ++ei){

output=g[*ei].g_struct_.alpha_.size();

}	
	
if(output==-1){ throw std::runtime_error("Couldn't count alpha_size correctly");}	
	
return output;	
	
}

void AmiGraph::construct_ami_sol(graph_t &g, NewAmiCalc::solution_set &ss, double ereg){



graph_to_R0(g, ss.R0_);

extract_bose_alphas(g, ss.bose_alphas_);
ss.loops_=count_fermi_loops(g)-g[boost::graph_bundle].ct_count;

int graph_ord=graph_order(g)-g[boost::graph_bundle].ct_count;
ss.ct_count_=g[boost::graph_bundle].ct_count;
ss.sigma_ct_count_=g[boost::graph_bundle].sigma_ct_count;
ss.order_shift=0;
ss.prefactor_=get_prefactor(g, graph_ord);


AmiBase::ami_parms amiparms(graph_ord, ereg, ami_parameters.TYPE_);


// amiparms.N_INT_=alpha_size(g)-amiparms.N_EXT_;


if((graph_type == AmiBase::Pi_phuu) || (graph_type == AmiBase::Pi_phud) || graph_type== AmiBase::density|| graph_type== AmiBase::ENERGY || graph_type==AmiBase::doubleocc || ami_parameters.TYPE_==AmiBase::Pi_ppud || ami_parameters.TYPE_==AmiBase::Pi_ppuu || ami_parameters.TYPE_==AmiBase::FORCE ){
amiparms.N_INT_++; // increase the integration order by one for bubble graphs.
// also increase it by one if it is a closed density graph.  
}


// std::cout<<"Assigning solution with alpha_size "<<alpha_size(g)<<" and N_INT_="<<amiparms.N_INT_<<std::endl;
ss.ami_parms_=amiparms;

// std::cout<<"construct"<<std::endl;
// ami.print_g_prod_info( gg.ss_vec[m].R0_);
ami.amibase.construct(ss.ami_parms_, ss.R0_, ss.R_, ss.P_, ss.S_);
// std::cout<<"constructed"<<std::endl;
	
	
	
return; 	
	
}


void AmiGraph::construct_ami_terms_sol(graph_t &g, NewAmiCalc::solution_set &ss, double ereg){

// std::cout<<"Constructing ami terms sol in function"<<std::endl;

graph_to_R0(g, ss.R0_);

extract_bose_alphas(g, ss.bose_alphas_);
ss.loops_=count_fermi_loops(g)-g[boost::graph_bundle].ct_count;

int graph_ord=graph_order(g)-g[boost::graph_bundle].ct_count;
ss.ct_count_=g[boost::graph_bundle].ct_count;
ss.sigma_ct_count_=g[boost::graph_bundle].sigma_ct_count;
ss.order_shift=0;
ss.prefactor_=get_prefactor(g, graph_ord);


AmiBase::ami_parms amiparms(graph_ord, ereg, ami_parameters.TYPE_);


// amiparms.N_INT_=alpha_size(g)-amiparms.N_EXT_;


if((graph_type == AmiBase::Pi_phuu) || (graph_type == AmiBase::Pi_phud) || graph_type== AmiBase::density|| graph_type== AmiBase::ENERGY || graph_type==AmiBase::doubleocc || ami_parameters.TYPE_==AmiBase::Pi_ppud || ami_parameters.TYPE_==AmiBase::Pi_ppuu || ami_parameters.TYPE_==AmiBase::FORCE ){
amiparms.N_INT_++; // increase the integration order by one for bubble graphs.
// also increase it by one if it is a closed density graph.  
}


// std::cout<<"Assigning solution with alpha_size "<<alpha_size(g)<<" and N_INT_="<<amiparms.N_INT_<<std::endl;
ss.ami_parms_=amiparms;

// std::cout<<"construct"<<std::endl;
// ami.print_g_prod_info( gg.ss_vec[m].R0_);
ami.amibase.construct(ss.ami_parms_.N_INT_, ss.R0_, ss.ami_terms_);
// std::cout<<"constructed"<<std::endl;
	
	
	
return; 	
	
}


void AmiGraph::gg_construct_ami_term_sol(graph_group &gg, double ereg){

gg.ss_vec.resize(gg.graph_vec.size());

// std::cout<<"Constructing"<<std::endl;
// std::cout<<gg.graph_vec.size()<<" "<<gg.ss_vec.size()<<std::endl;
	
for(int m=0; m< gg.graph_vec.size(); m++){	

// std::cout<<"Calling R0"<<std::endl;

// std::cout<<"Converting graph m="<<m<<std::endl;
// std::cout<<"Graph order is "<<graph_order(gg.graph_vec[m])<<std::endl;
// print_all_edge_info(gg.graph_vec[m]);

graph_to_R0(gg.graph_vec[m], gg.ss_vec[m].R0_);

ami.print_g_prod_info( gg.ss_vec[m].R0_);

// std::cout<<"R0 done"<<std::endl;
// std::cout<<"Graph order is "<<graph_order(gg.graph_vec[m])<<std::endl;
// stash the bosonic lines in the solution_set_matrix_t

// std::cout<<std::flush;
extract_bose_alphas(gg.graph_vec[m], gg.ss_vec[m].bose_alphas_);
gg.ss_vec[m].loops_=count_fermi_loops(gg.graph_vec[m])-gg.graph_vec[m][boost::graph_bundle].ct_count;
// std::cout<<"Graph order is "<<graph_order(gg.graph_vec[m])<<std::endl;
//
// std::cout<<"Resulting count is "<< gg.ss_vec[m].loops_<<std::endl;

// std::cout<<"Checking "<<std::endl;
// std::cout<<std::flush;
//int graph_ord=graph_order(g)-g[boost::graph_bundle].ct_count;

// TODO: Unsure about the prefactor associated with the counter-term
// std::cout<<"Graph order is "<<graph_order(gg.graph_vec[m])<<std::endl;
int graph_ord=graph_order(gg.graph_vec[m])-gg.graph_vec[m][boost::graph_bundle].ct_count+gg.order_shift;
// std::cout<<"Graph order is "<<graph_order(gg.graph_vec[m])<<std::endl;
gg.ss_vec[m].ct_count_=gg.graph_vec[m][boost::graph_bundle].ct_count;
gg.ss_vec[m].sigma_ct_count_=gg.graph_vec[m][boost::graph_bundle].sigma_ct_count;
gg.ss_vec[m].order_shift=gg.order_shift;

// std::cout<<"Order is "<<graph_ord<<" "<<graph_order(gg.graph_vec[m])<<std::endl;
// std::cout<<"Sol set vector is size "<< gg.ss_vec.size()<<std::endl;

gg.ss_vec[m].prefactor_=get_prefactor(gg.graph_vec[m], graph_ord);

// std::cout<<"ss_vec prefactor "<<gg.ss_vec[m].prefactor_ <<std::endl;

AmiBase::ami_parms amiparms(graph_ord, ereg, ami_parameters.TYPE_);

// std::cout<<"NINT is "<<amiparms.N_INT_<<std::endl;

if((graph_type == AmiBase::Pi_phuu) || (graph_type == AmiBase::Pi_phud) || graph_type== AmiBase::density || graph_type== AmiBase::ENERGY || graph_type==AmiBase::doubleocc || graph_type==AmiBase::Pi_ppud || graph_type==AmiBase::Pi_ppuu || graph_type==AmiBase::FORCE){
amiparms.N_INT_++; // increase the integration order by one for bubble graphs.
// also increase it by one if it is a closed density graph.  
}
// std::cout<<"NINT is "<<amiparms.N_INT_<<std::endl;

// For Doubleocc, the N_int_ ends up counting the external line already, so it only needs to be one larger, not two
// if(graph_type==AmiBase::doubleocc){
	// amiparms.N_INT_++;
	// amiparms.N_INT_++;
// }

// std::cout<<"Assigning solution with alpha_size "<<alpha_size(gg.graph_vec[m])<<" and N_INT_="<<amiparms.N_INT_<<std::endl;
// std::cout<<"Bubble and sigma terms are "<< gg.graph_vec[m][boost::graph_bundle].ct_count<<" and "<<gg.graph_vec[m][boost::graph_bundle].sigma_ct_count<<std::endl;
// std::cout<<"two"<<std::endl;
gg.ss_vec[m].ami_parms_=amiparms;

// std::cout<<"construct"<<std::endl;
// ami.print_g_prod_info( gg.ss_vec[m].R0_);
// ami.amibase.construct(gg.ss_vec[m].ami_parms_, gg.ss_vec[m].R0_, gg.ss_vec[m].R_, gg.ss_vec[m].P_, gg.ss_vec[m].S_);
// std::cout<<"constructed"<<std::endl;


ami.amibase.construct(gg.ss_vec[m].ami_parms_.N_INT_, gg.ss_vec[m].R0_, gg.ss_vec[m].ami_terms_);
// std::cout<<"constructed"<<std::endl;
// std::cout<<"Number of terms is "<<gg.ss_vec[m].ami_terms_.size()<<std::endl;

}	

// std::cout<<"Done group moving to pairs"<<std::endl;

// construct solutions for pairs 
// std::cout<<"pair size is "<< gg.gp_vec.size()<<std::endl;
for(int m=0; m< gg.gp_vec.size(); m++){

// g1 
graph_to_R0(gg.gp_vec[m].g1_, gg.gp_vec[m].s1_.R0_);
int graph_ord=graph_order(gg.gp_vec[m].g1_);
gg.gp_vec[m].s1_.prefactor_=get_prefactor(gg.gp_vec[m].g1_, graph_ord);

AmiBase::ami_parms amiparms(graph_ord, ereg, ami_parameters.TYPE_);
if((graph_type == AmiBase::Pi_phuu) || (graph_type == AmiBase::Pi_phud) || graph_type== AmiBase::density || graph_type== AmiBase::ENERGY || graph_type==AmiBase::doubleocc || graph_type==AmiBase::Pi_ppud || graph_type==AmiBase::Pi_ppuu || graph_type==AmiBase::FORCE){
amiparms.N_INT_++; // increase the integration order by one for bubble graphs. or closed density graphs 
}
// increment twice for double-occ 
// if(graph_type==AmiBase::doubleocc){
	// amiparms.N_INT_++;
	// amiparms.N_INT_++;
// }


gg.gp_vec[m].s1_.ami_parms_=amiparms;
ami.amibase.construct(gg.gp_vec[m].s1_.ami_parms_, gg.gp_vec[m].s1_.R0_, gg.gp_vec[m].s1_.R_, gg.gp_vec[m].s1_.P_, gg.gp_vec[m].s1_.S_);	
// g1 end

// g2 
graph_to_R0(gg.gp_vec[m].g2_, gg.gp_vec[m].s2_.R0_);
graph_ord=graph_order(gg.gp_vec[m].g2_);
gg.gp_vec[m].s2_.prefactor_=get_prefactor(gg.gp_vec[m].g2_, graph_ord);

// amiparms the same for the pair 
gg.gp_vec[m].s2_.ami_parms_=amiparms;
ami.amibase.construct(gg.gp_vec[m].s2_.ami_parms_, gg.gp_vec[m].s2_.R0_, gg.gp_vec[m].s2_.R_, gg.gp_vec[m].s2_.P_, gg.gp_vec[m].s2_.S_);
//g2 end 

	
}




	
}


void AmiGraph::gg_construct_ami_sol(graph_group &gg, double ereg){

gg.ss_vec.resize(gg.graph_vec.size());

// std::cout<<"Constructing"<<std::endl;
// std::cout<<gg.graph_vec.size()<<" "<<gg.ss_vec.size()<<std::endl;
	
for(int m=0; m< gg.graph_vec.size(); m++){	

// std::cout<<"Calling R0"<<std::endl;

// std::cout<<"Converting graph m="<<m<<std::endl;
// print_all_edge_info(gg.graph_vec[m]);

graph_to_R0(gg.graph_vec[m], gg.ss_vec[m].R0_);

// ami.print_g_prod_info( gg.ss_vec[m].R0_);

// std::cout<<"R0 done"<<std::endl;

// stash the bosonic lines in the solution_set_matrix_t

// std::cout<<std::flush;
extract_bose_alphas(gg.graph_vec[m], gg.ss_vec[m].bose_alphas_);
gg.ss_vec[m].loops_=count_fermi_loops(gg.graph_vec[m])-gg.graph_vec[m][boost::graph_bundle].ct_count;
gg.ss_vec[m].loops_=count_fermi_loops(gg.graph_vec[m])-gg.graph_vec[m][boost::graph_bundle].ct_count;

//
// std::cout<<"Resulting count is "<< gg.ss_vec[m].loops_<<std::endl;

// std::cout<<"Checking "<<std::endl;
// std::cout<<std::flush;
//int graph_ord=graph_order(g)-g[boost::graph_bundle].ct_count;

// TODO: Unsure about the prefactor associated with the counter-term
int graph_ord=graph_order(gg.graph_vec[m])-gg.graph_vec[m][boost::graph_bundle].ct_count+gg.order_shift;
gg.ss_vec[m].ct_count_=gg.graph_vec[m][boost::graph_bundle].ct_count;
gg.ss_vec[m].sigma_ct_count_=gg.graph_vec[m][boost::graph_bundle].sigma_ct_count;
gg.ss_vec[m].order_shift=gg.order_shift;
// std::cout<<"Order is "<<graph_ord<<" "<<graph_order(gg.graph_vec[m])<<std::endl;
// std::cout<<"Sol set vector is size "<< gg.ss_vec.size()<<std::endl;

gg.ss_vec[m].prefactor_=get_prefactor(gg.graph_vec[m], graph_ord);

// std::cout<<"ss_vec prefactor "<<gg.ss_vec[m].prefactor_ <<std::endl;

AmiBase::ami_parms amiparms(graph_ord, ereg, ami_parameters.TYPE_);

// std::cout<<"Line 1634: NINT is "<<amiparms.N_INT_<<std::endl;
// exit(0);
if((graph_type == AmiBase::Pi_phuu) || (graph_type == AmiBase::Pi_phud) || graph_type== AmiBase::density || graph_type== AmiBase::ENERGY || graph_type==AmiBase::doubleocc || graph_type==AmiBase::Pi_ppud || graph_type==AmiBase::Pi_ppuu || graph_type==AmiBase::FORCE){
amiparms.N_INT_++; // increase the integration order by one for bubble graphs.
// also increase it by one if it is a closed density graph.  
}


// std::cout<<"NINT is "<<amiparms.N_INT_<<std::endl;

// For Doubleocc, the N_int_ ends up counting the external line already, so it only needs to be one larger, not two
// if(graph_type==AmiBase::doubleocc){
	// amiparms.N_INT_++;
	// amiparms.N_INT_++;
// }

// std::cout<<"Assigning solution with alpha_size "<<alpha_size(gg.graph_vec[m])<<" and N_INT_="<<amiparms.N_INT_<<std::endl;
// std::cout<<"Bubble and sigma terms are "<< gg.graph_vec[m][boost::graph_bundle].ct_count<<" and "<<gg.graph_vec[m][boost::graph_bundle].sigma_ct_count<<std::endl;
// std::cout<<"two"<<std::endl;
gg.ss_vec[m].ami_parms_=amiparms;

// std::cout<<"construct"<<std::endl;
// ami.print_g_prod_info( gg.ss_vec[m].R0_);
ami.amibase.construct(gg.ss_vec[m].ami_parms_, gg.ss_vec[m].R0_, gg.ss_vec[m].R_, gg.ss_vec[m].P_, gg.ss_vec[m].S_);
// std::cout<<"constructed"<<std::endl;


}	

// std::cout<<"Done group moving to pairs"<<std::endl;

// construct solutions for pairs 
// std::cout<<"pair size is "<< gg.gp_vec.size()<<std::endl;
for(int m=0; m< gg.gp_vec.size(); m++){

// g1 
graph_to_R0(gg.gp_vec[m].g1_, gg.gp_vec[m].s1_.R0_);
int graph_ord=graph_order(gg.gp_vec[m].g1_);
gg.gp_vec[m].s1_.prefactor_=get_prefactor(gg.gp_vec[m].g1_, graph_ord);

AmiBase::ami_parms amiparms(graph_ord, ereg, ami_parameters.TYPE_);
if((graph_type == AmiBase::Pi_phuu) || (graph_type == AmiBase::Pi_phud) || graph_type== AmiBase::density || graph_type== AmiBase::ENERGY || graph_type==AmiBase::doubleocc || graph_type==AmiBase::Pi_ppud || graph_type==AmiBase::Pi_ppuu || graph_type==AmiBase::FORCE){
amiparms.N_INT_++; // increase the integration order by one for bubble graphs. or closed density graphs 
}
// increment twice for double-occ 
// if(graph_type==AmiBase::doubleocc){
	// amiparms.N_INT_++;
	// amiparms.N_INT_++;
// }


gg.gp_vec[m].s1_.ami_parms_=amiparms;
ami.amibase.construct(gg.gp_vec[m].s1_.ami_parms_, gg.gp_vec[m].s1_.R0_, gg.gp_vec[m].s1_.R_, gg.gp_vec[m].s1_.P_, gg.gp_vec[m].s1_.S_);	
// g1 end

// g2 
graph_to_R0(gg.gp_vec[m].g2_, gg.gp_vec[m].s2_.R0_);
graph_ord=graph_order(gg.gp_vec[m].g2_);
gg.gp_vec[m].s2_.prefactor_=get_prefactor(gg.gp_vec[m].g2_, graph_ord);

// amiparms the same for the pair 
gg.gp_vec[m].s2_.ami_parms_=amiparms;
ami.amibase.construct(gg.gp_vec[m].s2_.ami_parms_, gg.gp_vec[m].s2_.R0_, gg.gp_vec[m].s2_.R_, gg.gp_vec[m].s2_.P_, gg.gp_vec[m].s2_.S_);
//g2 end 

	
}




	
}

void AmiGraph::ggm_remove_pairs(gg_matrix_t &ggm, int min){
	
for(int ord=min; ord< ggm.size(); ord++){
	for(int group=0; group<ggm[ord].size(); group++){
	  for(int pair=0; pair<ggm[ord][group].gp_vec.size(); pair++){

		ggm[ord][group].graph_vec.push_back(ggm[ord][group].gp_vec[pair].g1_);
		ggm[ord][group].graph_vec.push_back(ggm[ord][group].gp_vec[pair].g2_);
		
		  
	  }
	// after all pairs pushing into graph_vec. we delete the gp_vecs
	ggm[ord][group].gp_vec.clear();

	}
}	
	
	
	
}

void AmiGraph::ggm_clear_labels(gg_matrix_t &ggm){
	for(int ord=1; ord< ggm.size(); ord++){
	for(int group=0; group<ggm[ord].size(); group++){
		
		ggm[ord][group].labels.clear();
		
		
	}
	}
	
}



// void AmiGraph::gg_eval_ami_gg(graph_group &gg){
	
	
	
	
	
// }



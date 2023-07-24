#include "amigraph.hpp"



// returns true if passes filter (non-zero).  false if it is zero. 
bool AmiGraph::residue_filter(graph_t &g){
	// TODO should add a check that the graph is actually labelled 
int ord=graph_order(g);	
AmiBase::g_prod_t R0=graph_to_R0(g);

//bool result=true;

AmiBase::pole_array_t poles, multipoles;

for(int ind=0; ind<ord; ind++){
	poles.clear();
	multipoles.clear();

poles=ami.amibase.find_poles(ind,R0);

for(int i=0; i< poles.size(); i++){
if(poles[i].multiplicity_>1){ multipoles.push_back(poles[i]);}	
}

for(int i=0; i< multipoles.size(); i++){
AmiBase::g_prod_t Wi=ami.amibase.reduce_gprod(R0,multipoles[i]);
bool found=false;
for(int m=0; m< Wi.size(); m++){
if(Wi[m].alpha_[ind]!=0){found=true; continue;}
}

if(found==false){
	// std::cout<<"Residue returning false for index "<<ind<<" with multiplicity "<< multipoles[i].multiplicity_<<  std::endl;


	return false;
}	
	
}
	
	
	
	
}
	
	
	
	
	
return true;	
}


#include "amigraph.hpp"

void AmiGraph::ggm_check_pair_momenta(gg_matrix_t &ggm, int mpi_rank){


std::cout<<ggm.size()<<std::endl;	
for(int ord=0; ord< ggm.size(); ord++){
	if(mpi_rank==0){
std::cout<<"Checking pair momenta on ord "<< ord<<" with size "<<ggm[ord].size()<<std::endl;
	}	
for(int i=0; i< ggm[ord].size(); i++){
	// std::cout<<".";
	
		 
	 
	 for(int m=0; m< ggm[ord][i].gp_vec.size(); m++){
	 bool result=true;
	 bool result2=true;
	 check_momentum_conservation(ggm[ord][i].gp_vec[m].g1_, result);
	 check_momentum_conservation(ggm[ord][i].gp_vec[m].g2_, result2);
	 
	 // std::cout<<"On graph "<<ord<<" "<<i<<" pair number "<<m<<std::endl;
	 // std::cout<<std::endl;
	 
	 // std::cout<<"momentum conservation on pairs is "<<result<<" "<<result2<<std::endl;
	 
	 // std::cout<<std::endl;
	 
	 
	 
	 }
	 
	 

}
// std::cout<<std::endl;
}	
	
return;		
	
	
}

void AmiGraph::check_label_sets(gg_matrix_t &ggm,  int mpi_rank, int min){

std::cout<<ggm.size()<<std::endl;	
for(int ord=min; ord< ggm.size(); ord++){
	if(mpi_rank==0){
std::cout<<"Double checking labels: Working on order "<< ord<<" with size "<<ggm[ord].size()<<std::endl;
	}	
for(int i=0; i< ggm[ord].size(); i++){
	// std::cout<<".";
	 for(int j=0; j< ggm[ord][i].labels.size(); j++){
	// std::cout<<ggm[ord][i].labels.size()<<" "<<ggm[ord][i].labels[j].size()<<" "<<ggm[ord].size()<<" "<<j<<std::endl;
	
	 if(ggm[ord][i].labels[j].size()==0){
		 std::cout<<"Label error "<<ord<<" "<<i<< " "<<j<<std::endl;
		 throw std::runtime_error("Found an empty label - must exit");
		 	 }
	 
	 for(int m=0; m< ggm[ord][i].labels[j].size(); m++){
	 bool result=true;
	 check_momentum_conservation(ggm[ord][i].labels[j][m], result);
	 
	 // std::cout<<"On graph "<<ord<<" "<<i<<" term "<<j<<" label number "<<m<<std::endl;
	 // std::cout<<std::endl;
	 // print_all_edge_info(ggm[ord][i].labels[j][m]);
	 // std::cout<<std::endl;
	 
	 if(!result){
		 throw std::runtime_error("Some label is not conserving momentum ");
	 }
	 
	 }
	 
	 
} 
}
// std::cout<<std::endl;
}	
	
return;	
}

void AmiGraph::construct_label_sets(gg_matrix_t &ggm, int mpi_rank, int min, int size, int trim){
	if(mpi_rank==0){
	std::cout<<"Preconstructing label sets "<<std::endl;}
	
	
for(int ord=min; ord< ggm.size(); ord++){
	if(mpi_rank==0){
std::cout<<"Working on order "<< ord<<" with size "<<ggm[ord].size()<<std::endl;
	}	
for(int i=0; i< ggm[ord].size(); i++){
	// std::cout<<".";
	ggm[ord][i].labels.resize(ggm[ord][i].graph_vec.size());
	for(int j=0; j< ggm[ord][i].graph_vec.size(); j++){
	
construct_label_set(ggm[ord][i].graph_vec[j], ggm[ord][i].labels[j], size);

if(ggm[ord][i].labels.size()>trim){
trim_label_set(ggm[ord][i].labels[j], trim);
}

}
}
// std::cout<<std::endl;
}
}	

// TODO: Label sets does not function at all 
void AmiGraph::construct_label_set(graph_t &g, labels_t &L, int size){

L.clear();
graph_t g_temp=g;
bool result=false;
bool added=false;

// int initial,fin;
int check=0;
// std::cout<<"Next one"<<std::endl;

sys_label_sets(L,g, result, size);

if(!result){
	// std::cerr<<"Warning: Failed to produce label sets systematically "<<std::endl;
// std::cerr<<"Falling back to repeated half-random labelling - does this graph have a tadpole?"<<std::endl;

do{
	
repeated_labelling(g_temp, result);	

if(result){
// initial=L.size();
append_to_labels(L, g_temp, added);
	//fin=L.size();
// std::cout<<"Added is "<< added<<std::endl;
if(!added){check++;}else{check=0;}	

// std::cout<<"Checking check="<<check<<std::endl;
}
}while(check<5 && L.size()<size); 


if(!result){throw std::runtime_error("Random label sets failed. exiting ");}
// if(result){std::cerr<<"Warning Resolved: Random labelling succeeded."<<std::endl;}
	
}

// replaced with systematic generation of labels 
/* do{
	
repeated_labelling(g_temp, result);	

if(result){
// initial=L.size();
append_to_labels(L, g_temp, added);
	//fin=L.size();
// std::cout<<"Added is "<< added<<std::endl;
if(!added){check++;}else{check=0;}	

// std::cout<<"Checking check="<<check<<std::endl;
}
}while(check<5 && L.size()<size); */



}

void AmiGraph::append_to_labels(labels_t &L, graph_t &g, bool &added){

bool in_list=false;
added=false;
	
for(int i=0; i< L.size(); i++){

if(equal_graph_labels(L[i],g)){in_list=true; break;}

}	
	
if(!in_list){
L.push_back(g);
added=true;
}//else{added=false;}	
	
}

bool AmiGraph::equal_graph_labels(graph_t &g1, graph_t &g2){
	
// This will assume that the order in which edges are found are equivalent
//These are safetly checks
if( num_vertices(g1)!= num_vertices(g2)){ 
return false;}
if( num_edges(g1)!= num_edges(g2)){ 
return false;}

bool result=true;

edge_vector_t e1, e2;

get_edges(e1,g1);
get_edges(e2,g2);

for(int i=0; i< e1.size(); i++){

if(	!edge_labels_are_equal(e1[i], g1, e2[i],g2)){ result=false; return result;}
	
}

return result;
	
}

void AmiGraph::get_edges(edge_vector_t &e, graph_t &g){

boost::graph_traits<graph_t>::edge_iterator ei, ei_end;

for(boost::tie(ei,ei_end)=edges(g); ei!=ei_end; ++ei){

e.push_back(*ei);

}	
	
	
	
}


void AmiGraph::reset_g(graph_t &g){

// reset labelled counters	
g[boost::graph_bundle].n_indep=0;
g[boost::graph_bundle].n_labelled=0;	

edge_vector_t bose_edges;	
edge_vector_t fermi_edges;

// edge_t bose_edge, fermi_edge;
// find_one_fermi_bose_edge(g, bose_edge, fermi_edge);
find_bose_fermi_edges(g, bose_edges, fermi_edges);
// question is this internal only?


// find_bosonic_edges(g, bose_edges);
// find_fermionic_edges(g, fermi_edges);

//std::vector<int> test1;
	
	
	
boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
// Looking through all edges to find source and targets : this could be improved?
for (boost::tie(ei,ei_end) = edges(g); ei != ei_end; ++ei){




//if( g[*ei].g_struct_.stat_==AmiBase::Fermi){  /// DONT need THIS LOGIC? Bosonic lines get same reset of labels.
// std::vector<int> v{0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
 
// std::fill(v.begin(), v.end(), -1);
//test1.resize(bose_edges.size()+1);
//std::fill(test1.begin, test1.end,0);

if(graph_type == AmiBase::Pi_phuu || graph_type == AmiBase::Pi_phud || (graph_type== AmiBase::doubleocc) || ami_parameters.TYPE_==AmiBase::Pi_ppud || ami_parameters.TYPE_==AmiBase::Pi_ppuu){
	
		// std::cout<<"Bose edges are "<<bose_edges.size()<<std::endl;
	// std::cout<<"Fermi edges are "<<fermi_edges.size()<<std::endl;

g[*ei].g_struct_.alpha_.resize(bose_edges.size());
std::fill(g[*ei].g_struct_.alpha_.begin(),g[*ei].g_struct_.alpha_.end(), 0); 
// TODO: figure out ACTUALLY how big eps vectors are because this is not correct
g[*ei].g_struct_.eps_.resize(fermi_edges.size());
std::fill(g[*ei].g_struct_.eps_.begin(), g[*ei].g_struct_.eps_.end(), 0);

g[*ei].label=unlabelled;	
	
	

}else if(graph_type==AmiBase::density || graph_type==AmiBase::Greens || graph_type==AmiBase::DOS || graph_type==AmiBase::ENERGY){


if( num_edges(g)==1){

g[*ei].g_struct_.alpha_.resize(1);
g[*ei].g_struct_.eps_.resize(1);
}else {
	
	g[*ei].g_struct_.alpha_.resize(bose_edges.size()+1);
	g[*ei].g_struct_.eps_.resize(fermi_edges.size()-1);
	
}

std::fill(g[*ei].g_struct_.alpha_.begin(),g[*ei].g_struct_.alpha_.end(), 0); 
std::fill(g[*ei].g_struct_.eps_.begin(), g[*ei].g_struct_.eps_.end(), 0);

g[*ei].label=unlabelled;	

}else if(graph_type==AmiBase::FORCE){

int this_ord=graph_order(g)+2;
g[*ei].g_struct_.alpha_.resize(this_ord);
std::fill(g[*ei].g_struct_.alpha_.begin(),g[*ei].g_struct_.alpha_.end(), 0); 
// TODO: figure out ACTUALLY how big eps vectors are because this is not correct
g[*ei].g_struct_.eps_.resize(fermi_edges.size());
std::fill(g[*ei].g_struct_.eps_.begin(), g[*ei].g_struct_.eps_.end(), 0);

g[*ei].label=unlabelled;	
	
}
// else
// if(graph_type== AmiCalc::density){

	
// g[*ei].g_struct_.alpha_.resize(bose_edges.size()+1);
// std::fill(g[*ei].g_struct_.alpha_.begin(),g[*ei].g_struct_.alpha_.end(), 0); 
////TODO: figure out ACTUALLY how big eps vectors are because this is not correct
// g[*ei].g_struct_.eps_.resize(fermi_edges.size()-3);
// std::fill(g[*ei].g_struct_.eps_.begin(), g[*ei].g_struct_.eps_.end(), 0);

// g[*ei].label=unlabelled;
	
	
// }
else{
	
	
g[*ei].g_struct_.alpha_.resize(bose_edges.size()+1);
std::fill(g[*ei].g_struct_.alpha_.begin(),g[*ei].g_struct_.alpha_.end(), 0); 
// TODO: figure out ACTUALLY how big eps vectors are because this is not correct
g[*ei].g_struct_.eps_.resize(fermi_edges.size()-2);
std::fill(g[*ei].g_struct_.eps_.begin(), g[*ei].g_struct_.eps_.end(), 0);

g[*ei].label=unlabelled;
	
	
}



}

}

void AmiGraph::label_half_random(graph_t &g){
	
reset_g(g);	
// initialize some counters on the interior of the graph
g[boost::graph_bundle].n_indep=0;
g[boost::graph_bundle].n_labelled=0;

//number_vertices(g);  // assume it is already numbered. - moved to repeated labelling at the start
	
	
// label the external legs
vertex_vector_t extern_vect_list;
edge_vector_t extern_edge_list;
find_external_vertices(g, extern_vect_list, extern_edge_list);

// if(ami_parameters.TYPE_!=AmiBase::density && ami_parameters.TYPE_!= AmiBase::doubleocc){
if(extern_vect_list.size()<1){  throw std::runtime_error("Can't find external leg.");}
// }

label_extern_legs(extern_edge_list,g);

// Need function to find and label tadpoles
vertex_vector_t tp_vert, tp_conn_vert; 
edge_vector_t tp_bose_edges;
label_and_find_tadpoles_ami(g,tp_vert, tp_conn_vert, tp_bose_edges);

// now get all the internal vertices 

vertex_vector_t internal;
find_internal_vertices(g, internal);

// fermionic loop stuff vs label along line stuff
vertex_vector_t labelled_adj_vertices, unlabelled_adj_vertices;
edge_vector_t labelled_edges, unlabelled_edges;


do{	
int rand_choice=random_int(0,internal.size()-1);
sort_labelled_unlabelled_adjacent(g, internal[rand_choice],labelled_adj_vertices, unlabelled_adj_vertices, labelled_edges, unlabelled_edges);
//
// std::cout<<"Currently at "<< g[internal[rand_choice]].index_<<std::endl;
// std::cout<<"Of list with length "<< internal.size()<<" chose random entry "<< rand_choice<<std::endl;
// std::cout<<"Labelled are "<< labelled_edges.size()<<"unlabelled are "<< unlabelled_edges.size()<<std::endl;
bool junk=label_consv_momentum(g, internal[rand_choice], labelled_adj_vertices, unlabelled_adj_vertices, labelled_edges, unlabelled_edges);	


internal.erase(internal.begin()+rand_choice);
	
}while(internal.size()>0)	;
	
	
}


void AmiGraph::label_systematic_redone(graph_t &g){
reset_g(g);

reset_visited(g);

// initialize some counters on the interior of the graph
g[boost::graph_bundle].n_indep=0;
g[boost::graph_bundle].n_labelled=0;

number_vertices(g);
vertex_vector_t vert_list;

collect_unvisited_vertices(g, vert_list);

vertex_vector_t current_list;
vertex_vector_t next_list;

current_list.push_back(vert_list[0]);



// fermionic loop stuff vs label along line stuff
vertex_vector_t labelled_adj_vertices, unlabelled_adj_vertices;
edge_vector_t labelled_edges, unlabelled_edges;

// label the external legs
vertex_vector_t extern_vect_list;
edge_vector_t extern_edge_list;
find_external_vertices(g, extern_vect_list, extern_edge_list);

// if(ami_parameters.TYPE_!=AmiBase::density && ami_parameters.TYPE_ != AmiBase::doubleocc){
if(extern_vect_list.size()<1){  throw std::runtime_error("Can't find external leg.");}
// }

label_extern_legs(extern_edge_list,g);

// Need function to find and label tadpoles
vertex_vector_t tp_vert, tp_conn_vert; 
edge_vector_t tp_bose_edges;
//std::cout<<std::endl;
//std::cout<<"Found some tadpoles, starting to label them labelling "<< tp_conn_vert.size()<<std::endl;
label_and_find_tadpoles_ami(g,tp_vert, tp_conn_vert, tp_bose_edges);

int num_tadpoles=tp_vert.size();
////


// BFS to label everything 
bool exit_do=false;
do{
get_adjacent_lists(g, current_list, next_list);
label_visited(g, current_list);


do{
int rand_choice=random_int(0,current_list.size()-1);
sort_labelled_unlabelled_adjacent(g, current_list[rand_choice],labelled_adj_vertices, unlabelled_adj_vertices, labelled_edges, unlabelled_edges);
//
std::cout<<"Currently at "<< g[current_list[rand_choice]].index_<<std::endl;
std::cout<<"Of list with length "<< current_list.size()<<" chose random entry "<< rand_choice<<std::endl;
//std::cout<<"Labelled are "<< labelled_edges.size()<<"unlabelled are "<< unlabelled_edges.size()<<std::endl;
bool junk=label_consv_momentum(g, current_list[rand_choice], labelled_adj_vertices, unlabelled_adj_vertices, labelled_edges, unlabelled_edges);	


current_list.erase(current_list.begin()+rand_choice);
}while(current_list.size()>0);


if(next_list.size()==0){
exit_do=true;	
}else{
current_list=next_list;
}	
	
}while(!exit_do);



	
	

	
	
	
	
}

int AmiGraph::count_tadpoles(graph_t &g){
	

// Need function to find and label tadpoles
vertex_vector_t tp_vert, tp_conn_vert; 
edge_vector_t tp_bose_edges;
//std::cout<<std::endl;
//std::cout<<"Found some tadpoles, starting to label them labelling "<< tp_conn_vert.size()<<std::endl;
//label_and_find_tadpoles_ami(g,tp_vert, tp_conn_vert, tp_bose_edges);

find_tadpoles(g,tp_vert, tp_conn_vert, tp_bose_edges);

return tp_vert.size();	
	
	
}

// void AmiGraph::reset_epsilons(graph_t &g){
	
	
	
	
	
	
// }
// June 07 2022 
// Modified to resize epsilon to the number of internal edges
// Question: Is there any case where this is not what is wanted? 
void AmiGraph::reset_epsilons(graph_t &g){
	
edge_vector_t internal;
find_internal_fermionic_edges(g, internal);
int this_size=internal.size();//g[internal[i]].g_struct_.eps_.size();
for(int i=0; i< internal.size(); i++){
	
	// this resets the epsilons 
	// int this_size=g[internal[i]].g_struct_.eps_.size();
g[internal[i]].g_struct_.eps_.clear();
g[internal[i]].g_struct_.eps_.resize(this_size,0);
g[internal[i]].g_struct_.eps_[i]=1;
	
	}	

return;	
	
	
}


void AmiGraph::reset_epsilons(AmiBase::g_prod_t &R0){

int this_size=R0.size();//g[internal[i]].g_struct_.eps_.size();
for(int i=0; i< R0.size(); i++){
	
	// this resets the epsilons 
	// int this_size=g[internal[i]].g_struct_.eps_.size();
R0[i].eps_.clear();
R0[i].eps_.resize(this_size,0);
R0[i].eps_[i]=1;
	
	}	

return;	
	
	
}


void AmiGraph::syk_epsilons(graph_t &g){

edge_vector_t internal_fermi;

find_internal_fermionic_edges(g,internal_fermi);	

for(int i=1; i< internal_fermi.size(); i++){


	g[internal_fermi[i]].g_struct_.eps_=g[internal_fermi[0]].g_struct_.eps_;

}

	
	
}



void AmiGraph::fix_epsilons(graph_t &g){

edge_vector_t internal_fermi;

find_internal_fermionic_edges(g,internal_fermi);	

for(int i=0; i< internal_fermi.size(); i++){
for(int j=i; j< internal_fermi.size(); j++){


if( i!=j){

bool check=edge_alphas_are_equal(internal_fermi[i], internal_fermi[j],g);
bool check2=edge_alphas_are_negative(internal_fermi[i],internal_fermi[j],g);
if(check || check2){
	g[internal_fermi[j]].g_struct_.eps_=g[internal_fermi[i]].g_struct_.eps_;
}
	
	
}




}	
	
	
}

	
	
}



// order of labelling should be - external legs - then tadpoles - then things connected to tadpoles - then easy labels - then random
void AmiGraph::label_systematic(graph_t &g){
	
	
// std::cout<<"Labelling!"<<std::endl;	
// reset any existing labels - sets alpha and eps sizes appropriately
reset_g(g);
// initialize some counters on the interior of the graph
g[boost::graph_bundle].n_indep=0;
g[boost::graph_bundle].n_labelled=0;

// label the external legs
vertex_vector_t extern_vect_list;
edge_vector_t extern_edge_list;
find_external_vertices(g, extern_vect_list, extern_edge_list);

// if(ami_parameters.TYPE_!=AmiBase::density && ami_parameters.TYPE_ != AmiBase::doubleocc){
if(extern_vect_list.size()<1){  throw std::runtime_error("Can't find external leg.");}
// }

label_extern_legs(extern_edge_list,g);

// Need function to find and label tadpoles
vertex_vector_t tp_vert, tp_conn_vert; 
edge_vector_t tp_bose_edges;
//std::cout<<std::endl;
// std::cout<<"Found some tadpoles, starting to label them labelling "<< tp_conn_vert.size()<<std::endl;
//label_and_find_tadpoles_ami(g,tp_vert, tp_conn_vert, tp_bose_edges);
find_and_label_tadpole_strings(g,tp_vert, tp_conn_vert, tp_bose_edges);

int num_tadpoles=tp_vert.size();
// std::cout<<"Finished labelling tadpoles "<< tp_conn_vert.size()<<std::endl;
//std::cout<<std::endl;
//


////std::cout<<"so far have labelled "<< g[boost::graph_bundle].n_labelled<< " edges"<<std::endl;	
//print_all_edge_info();


vertex_t current;
vertex_t next;

// find inward external leg
if(source(extern_edge_list[0],g)==extern_vect_list[0]){current=target(extern_edge_list[0],g);}else{current=target(extern_edge_list[1],g);}
////std::cout<<"Starting with vertex "<< g[current].index_<<std::endl;


// fermionic loop stuff vs label along line stuff
vertex_vector_t labelled_adj_vertices, unlabelled_adj_vertices;
edge_vector_t labelled_edges, unlabelled_edges;

edge_vector_t unlab_fermi;
vertex_vector_t to_label;

bool keep_going;
bool skip_do;
int counter=0;
int outer_counter=0;
int marker=0;

keep_going=true;
skip_do=true;

do{
outer_counter++;
counter=0;
do{
counter++;

// first we are going to pick which vertex to go to next 
get_adjacent_vertex_stat(g, current, next, AmiBase::Fermi);
////std::cout<<"Proposed next vertex to be "<< g[next].index_<<std::endl;

////std::cout<<"Label surrounding vertex "<< g[current].index_<<std::endl;
	
// first get all the information around current vertex - labelled and unlabelled edges

sort_labelled_unlabelled_adjacent(g, current,labelled_adj_vertices, unlabelled_adj_vertices, labelled_edges, unlabelled_edges);
//
//std::cout<<"Currently at "<< g[current].index_<<std::endl;
//std::cout<<"Labelled are "<< labelled_edges.size()<<"unlabelled are "<< unlabelled_edges.size()<<std::endl;
keep_going=label_consv_momentum(g, current, labelled_adj_vertices, unlabelled_adj_vertices, labelled_edges, unlabelled_edges);	
if(!keep_going){marker++;}

	
//
////std::cout<<"so far have labelled "<< g[boost::graph_bundle].n_labelled<< " edges"<<std::endl;	
//print_all_edge_info();
//	

// sanity check and extra exit condition
if(	g[boost::graph_bundle].n_labelled==num_edges(g)){ keep_going=false;}
	
current=next;	

// Added on june 22 to try to fix some labelling issues

to_label.clear();
find_unlabelled_fermionic_edges(g, unlab_fermi);
collect_vertices_to_label(g, to_label);	
for(int i=0; i< to_label.size(); i++){
sort_labelled_unlabelled_adjacent(g, to_label[i],labelled_adj_vertices, unlabelled_adj_vertices, labelled_edges, unlabelled_edges);
bool junk=label_consv_momentum(g, to_label[i], labelled_adj_vertices, unlabelled_adj_vertices, labelled_edges, unlabelled_edges);
}
//

	
}while(marker<10);//keep_going);// && counter<10);
marker=0;
// ask are there any unlabelled fermionic edges?
to_label.clear();
find_unlabelled_fermionic_edges(g, unlab_fermi);
collect_vertices_to_label(g, to_label);	

//std::cout<<"remaining labels "<<unlab_fermi.size() << std::endl;	

if(unlab_fermi.size()!=0){


// This became redundant because we label ALL tadpoles at the start now.	
// if(tp_conn_vert.size()!=0){
// //pick current to start at labelled tadpole end 
// int tint=random_int(0, tp_conn_vert.size()-1);
// current=tp_conn_vert[tint];
 // std::cout<<"Moving to tadpole connect "<< g[tp_conn_vert[tint]].index_<<std::endl;
// //std::cout<<"bosonic line is "<< print_edge(tp_bose_edges[tint],g)<<" which has label "<< g[tp_bose_edges[tint]].label<<std::endl;
// tp_conn_vert.erase(tp_conn_vert.begin()+tint);

// }else
if ( to_label.size()!=0 ){
int tint=random_int(0, to_label.size()-1);
current=to_label[tint];
//std::cout<<"Moving to easy label "<< g[to_label[tint]].index_ <<" from list of length " << to_label.size() <<std::endl;
}
else{	
// pick edge from unlabelled fermi edges, let current =source of that edge 
int rint=random_int(0, unlab_fermi.size()-1);
//std::cout<<"Proposed moving to unlabelled source of "<<rint<<" "<< print_edge(unlab_fermi[rint],g) <<std::endl;
current=source(unlab_fermi[rint] ,g);	
}
	
	
}else{skip_do=false;}
	
	
}while(skip_do);// && outer_counter<10);	
	
// conceptually everything is labelled at this point. if not 

	
	
	
}







void AmiGraph::get_edges_fermi_bose(graph_t &g, vertex_t v, edge_vector_t fermi_edges, edge_vector_t bose_edges){
	// for safety, clear the edge lists first
	fermi_edges.clear();
	bose_edges.clear();
	boost::graph_traits<graph_t>::in_edge_iterator ei, ei_end;
	boost::graph_traits<graph_t>::out_edge_iterator eo, eo_end;
	
for (boost::tie(ei,ei_end) = in_edges(v,g); ei != ei_end; ++ei){

if( g[*ei].g_struct_.stat_==AmiBase::Fermi){
fermi_edges.push_back(*ei);
}else{
	bose_edges.push_back(*ei);
}
}

for (boost::tie(eo,eo_end) = out_edges(v,g); eo != eo_end; ++eo){

if( g[*eo].g_struct_.stat_==AmiBase::Fermi){
fermi_edges.push_back(*eo);
}else{
bose_edges.push_back(*eo);	
	
}	
	
}
		
}





void AmiGraph::find_unlabelled_fermionic_edges(graph_t &g, edge_vector_t &vector){

vector.clear();

////std::cout<< "Finding Fermionic edges" <<std::endl;
boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
// Looking through all edges to find source and targets : this could be improved?
for (boost::tie(ei,ei_end) = edges(g); ei != ei_end; ++ei){

if( g[*ei].g_struct_.stat_==AmiBase::Fermi && g[*ei].label==unlabelled){
vector.push_back(*ei);
}
////std::cout << "Edge " << *ei << " of stat_type "<< g[*ei].g_struct_.stat_ << std::endl;

////std::cout<<"edge "<<jg[*ei].edge_number_<<" loop Total is "<< total<<std::endl;
}


}

void AmiGraph::find_unlabelled_edges(graph_t &g, edge_vector_t &vector){

vector.clear();

////std::cout<< "Finding Fermionic edges" <<std::endl;
boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
// Looking through all edges to find source and targets : this could be improved?
for (boost::tie(ei,ei_end) = edges(g); ei != ei_end; ++ei){

if( g[*ei].label==unlabelled){
vector.push_back(*ei);
}
////std::cout << "Edge " << *ei << " of stat_type "<< g[*ei].g_struct_.stat_ << std::endl;

////std::cout<<"edge "<<jg[*ei].edge_number_<<" loop Total is "<< total<<std::endl;
}


}

/* void AmiGraph::find_entry_vertices(graph_t &g, vertex_vector_t &v){
v.clear();

boost::graph_traits<graph_t>::vertex_iterator vi, vi_end;
for (boost::tie(vi,vi_end) = vertices(g); vi != vi_end; ++vi){
	
	if (degree(*vi,g)==1 && in_degree(*vi,g)==0){ v=*vi;}

	}

}	
 */
 /// This finds the vertices that are actually external???
void AmiGraph::find_external_vertices(graph_t &g, vertex_vector_t &v, edge_vector_t &edges){

v.clear();
edges.clear();
	
boost::graph_traits<graph_t>::in_edge_iterator iei, iei_end;
boost::graph_traits<graph_t>::out_edge_iterator oei, oei_end;	
	
boost::graph_traits<graph_t>::vertex_iterator ei, ei_end;
for (boost::tie(ei,ei_end) = vertices(g); ei != ei_end; ++ei){
	
	if (degree(*ei,g)==1  ){ v.push_back(*ei);
	for (boost::tie(iei,iei_end) = in_edges(*ei,g); iei != iei_end; ++iei){
		edges.push_back(*iei);
		
	}
	
	for (boost::tie(oei,oei_end) = out_edges(*ei,g); oei != oei_end; ++oei){
		edges.push_back(*oei);
	}
	
	
	}
	
}
	
	
}

// this function is redundant, since it needs the same operations as the next part of the step...
// the order matters
void AmiGraph::get_adjacent_vertex_stat(graph_t &g, vertex_t &vin, vertex_t &next,  AmiBase::stat_type stat){
	
edge_vector_t edges;	
	
boost::graph_traits<graph_t>::adjacency_iterator ai, ai_end;

for (boost::tie(ai,ai_end) = adjacent_vertices(vin, g); ai != ai_end; ++ai){

if(edge(vin,*ai,g).second){	
if(g[edge(vin,*ai,g).first].g_struct_.stat_==stat){
	
next=*ai;	
	
}

}
}

}

// this function is redundant, since it needs the same operations as the next part of the step...
// the order matters
void AmiGraph::get_adjacent_vertex_stat(graph_t &g, vertex_t &vin, vertex_t &next, edge_t &eout,  AmiBase::stat_type stat){
	
edge_vector_t edges;	
	
boost::graph_traits<graph_t>::adjacency_iterator ai, ai_end;

for (boost::tie(ai,ai_end) = adjacent_vertices(vin, g); ai != ai_end; ++ai){

if(edge(vin,*ai,g).second){	
if(g[edge(vin,*ai,g).first].g_struct_.stat_==stat){
	
next=*ai;	
eout=edge(vin,*ai,g).first;
	
}

}
}

}


void AmiGraph::get_adjacent_vertex_stat(graph_t &g, vertex_t &vin, vertex_t &next,  AmiBase::stat_type stat, bool &success){
	
edge_vector_t edges;	
success=false;
	
boost::graph_traits<graph_t>::adjacency_iterator ai, ai_end;

for (boost::tie(ai,ai_end) = adjacent_vertices(vin, g); ai != ai_end; ++ai){

if(edge(vin,*ai,g).second){	
if(g[edge(vin,*ai,g).first].g_struct_.stat_==stat){
	
next=*ai;
success=true;
break;	
	
}

}
}

}





bool AmiGraph::take_a_step(vertex_t &current, vertex_t &next, edge_t &connect, graph_t& g){
	
//label_edge(current, next, connect,g);
	
	
	
return false;	
	
}
//get_ab_edges(g, query, edge_a, edge_b);
void AmiGraph::get_ab_edges(graph_t &g, vertex_t &v, edge_t &edge_a, edge_t &edge_b){
	
boost::graph_traits<graph_t>::edge_iterator ei, ei_end;

for(boost::tie(ei,ei_end)=edges(g); ei!=ei_end; ++ei){

if(target(*ei,g)==v && g[*ei].g_struct_.stat_==AmiBase::Fermi){ edge_a=*ei;}
if(source(*ei,g)==v && g[*ei].g_struct_.stat_==AmiBase::Fermi){ edge_b=*ei;}



}	
	
	
	
}

void AmiGraph::treat_ab_edges(graph_t &g, edge_vector_t &edges_a,edge_vector_t &edges_b, edge_vector_t &new_edges_a,edge_vector_t &new_edges_b){
new_edges_a.clear();
new_edges_b.clear();


for(int i=0; i< edges_a.size(); i++){

if(source(edges_a[i],g)==target(edges_b[i],g)){

label_indep_edge(edges_a[i],g);

edge_t bose;
find_bosonic_edge_on_three_pointvert(g,source(edges_a[i],g), bose);

vertex_t v=source(edges_a[i],g);
assign_cons_label( g, v ,bose, edges_a[i], edges_b[i]);

vertex_t query;

if(target(bose,g)==v){query=source(bose,g);}
if(source(bose,g)==v){query=target(bose,g);}

edge_t edge_a, edge_b;
get_ab_edges(g, query, edge_a, edge_b);
new_edges_a.push_back(edge_a);
new_edges_b.push_back(edge_b);

//void find_bosonic_edge_on_three_pointvert(graph_t &g, vertex_t vert, edge_t &bose_edge);


}	
	
	
}
	
}

void AmiGraph::find_and_label_tadpole_strings(graph_t &g, vertex_vector_t &tp_vec, vertex_vector_t &tp_conn_vec, edge_vector_t &tp_bose_edges){
	
//	vertex_vector_t tp_vert_vec;
//	edge_vector_t tp_bose_line_vec;
//	vertex_vector_t tp_connect_vert_vec;

edge_vector_t edges_a, edges_b;
find_tadpoles(g, tp_vec, tp_conn_vec, tp_bose_edges, edges_a, edges_b);

// std::cout<<"Found "<< tp_vec.size()<<" tadpoles"<<std::endl;

if(edges_a.size()!= tp_conn_vec.size() || edges_b.size()!=tp_conn_vec.size()){throw std::runtime_error("Did not find edges_a and edges_b correctly");}	
	
edge_t tp_edge;	
for (int i=0; i< tp_vec.size(); i++){
	
	// std::cout<<"Found some tadpoles "<<i<<std::endl;

tp_edge=edge(tp_vec[i],tp_vec[i],g).first;

label_indep_edge(tp_edge,g);

// print_all_edge_info(g);

g[tp_bose_edges[i]].label=labelled;



if(g[edges_a[i]].label==unlabelled && g[edges_b[i]].label==unlabelled){
//std::cout<<"Logic b triggered"<<std::endl;
label_indep_edge(edges_a[i],g);
// assign_cons_label( g, tp_conn_vec[i], tp_bose_edges[i], edges_b[i], edges_a[i]);
		
}


if(g[edges_a[i]].label==unlabelled && g[edges_b[i]].label==labelled){
//std::cout<<"Logic b triggered"<<std::endl;
//label_indep_edge(edges_a[i],g);
assign_cons_label( g, tp_conn_vec[i], tp_bose_edges[i], edges_b[i], edges_a[i]);
		
}

if(g[edges_a[i]].label==labelled && g[edges_b[i]].label==unlabelled){
//std::cout<<"Logic c triggered"<<std::endl;
//label_indep_edge(edges_a[i],g);
assign_cons_label( g, tp_conn_vec[i], tp_bose_edges[i], edges_a[i], edges_b[i]);
		
}



}


// check for string - input is edge_a, edge_b,  
// if edge_a and b are parts of a tp_string. then label edge_a and edge_b
// Then label the bose line they attach to
// then return new edge_a and new edge_b.
// do while edge_a.size()==0

// 2022_03-25 This part not working I think 
// if(edges_a.size()>0){
// edge_vector_t new_edges_a, new_edges_b;
// do{
// treat_ab_edges(g, edges_a, edges_b, new_edges_a, new_edges_b);
// edges_a=new_edges_a;
// edges_b=new_edges_b;
// }while(new_edges_a.size()!=0);
// }
	
	
}

void AmiGraph::label_and_find_tadpoles_ami(graph_t &g, vertex_vector_t &tp_vec, vertex_vector_t &tp_conn_vec, edge_vector_t &tp_bose_edges){

//	vertex_vector_t tp_vert_vec;
//	edge_vector_t tp_bose_line_vec;
//	vertex_vector_t tp_connect_vert_vec;

edge_vector_t edges_a, edges_b;

find_tadpoles(g, tp_vec, tp_conn_vec, tp_bose_edges, edges_a, edges_b);

// std::cout<<"Tadpoles and edges found "<<tp_conn_vec.size()<<" "<<edges_a.size()<<" "<<edges_b.size()<<std::endl;

if(edges_a.size()!= tp_conn_vec.size() || edges_b.size()!=tp_conn_vec.size()){throw std::runtime_error("Did not find edges_a and edges_b correctly");}	
	
edge_t tp_edge;	
for (int i=0; i< tp_vec.size(); i++){

tp_edge=edge(tp_vec[i],tp_vec[i],g).first;

label_indep_edge(tp_edge,g);

g[tp_bose_edges[i]].label=labelled;

// logic to label edges connected to tadpoles 

// if(g[edges_a[i]].label==unlabelled && g[edges_b[i]].label==unlabelled){
// std::cout<<"Logic a triggered"<<std::endl;
// label_indep_edge(edges_a[i],g);
// assign_cons_label( g, tp_conn_vec[i], tp_bose_edges[i], edges_a[i], edges_b[i]);
		
// }

if(g[edges_a[i]].label==unlabelled && g[edges_b[i]].label==labelled){
//std::cout<<"Logic b triggered"<<std::endl;
//label_indep_edge(edges_a[i],g);
assign_cons_label( g, tp_conn_vec[i], tp_bose_edges[i], edges_b[i], edges_a[i]);
		
}

if(g[edges_a[i]].label==labelled && g[edges_b[i]].label==unlabelled){
//std::cout<<"Logic c triggered"<<std::endl;
//label_indep_edge(edges_a[i],g);
assign_cons_label( g, tp_conn_vec[i], tp_bose_edges[i], edges_a[i], edges_b[i]);
		
}



}	
}

void AmiGraph::label_tadpoles_ami(graph_t &g){

	vertex_vector_t tp_vert_vec;
	edge_vector_t tp_bose_line_vec;
	vertex_vector_t tp_connect_vert_vec;
	find_tadpoles(g, tp_vert_vec, tp_connect_vert_vec, tp_bose_line_vec);

edge_t tp_edge;	
for (int i=0; i< tp_vert_vec.size(); i++){

tp_edge=edge(tp_vert_vec[i],tp_vert_vec[i],g).first;

label_indep_edge(tp_edge,g);

}	

for (int i=0;i<tp_bose_line_vec.size(); i++){
	
	g[tp_bose_line_vec[i]].label=labelled;
	
}

	
	
}

void AmiGraph::label_indep_edge(edge_t &e, graph_t &g){
	
	//std::cout<<"Attempting to add to alpha, and independent ends in element "<< g[boost::graph_bundle].n_indep<<std::endl;
	//std::cout<<"Alpha size is "<<g[e].g_struct_.alpha_.size()<<std::endl;
	
		//std::cout<<"Attempting to add to eps, in element "<< g[boost::graph_bundle].n_labelled<<std::endl;
	//std::cout<<"eps size is "<<g[e].g_struct_.eps_.size()<<std::endl;
	
	if(g[e].g_struct_.alpha_.size()<=g[boost::graph_bundle].n_indep || g[e].g_struct_.eps_.size()<=g[boost::graph_bundle].n_labelled){
	return;	
	}
	g[e].g_struct_.alpha_[g[boost::graph_bundle].n_indep]=1;
	g[e].g_struct_.eps_[g[boost::graph_bundle].n_labelled]=1;
	
	
	// update tracking integers
	g[e].label=labelled;
	g[boost::graph_bundle].n_indep++;
	g[boost::graph_bundle].n_labelled++; 
	
	
	
	
}

bool AmiGraph::is_labelled_bose_zero(edge_t &e,graph_t &g){

if(g[e].g_struct_.stat_ != AmiBase::Bose){return false;}

if (g[e].label==unlabelled){return false;}

if (g[e].label== labelled){

for(int i=0; i< g[e].g_struct_.alpha_.size();i++){
//std::cout<<i<<std::endl;	
if(g[e].g_struct_.alpha_[i]!=0){return false;}
	

}

}	
	
	
return true;	
	
}



void AmiGraph::label_indep_epsilon(edge_t &e, graph_t &g){
	
	
	//std::cout<<"Attempting to add to eps, in element "<< g[boost::graph_bundle].n_labelled<<std::endl;
	//std::cout<<"eps size is "<<g[e].g_struct_.eps_.size()<<std::endl;	
	
//	g[e].g_struct_.alpha_[g[boost::graph_bundle].n_indep]=1;
	g[e].g_struct_.eps_[g[boost::graph_bundle].n_labelled]=1;
	
	
	// update tracking integers
	g[e].label=labelled;
//	g[boost::graph_bundle].n_indep++;
	g[boost::graph_bundle].n_labelled++;
	
	
	
	
}

void AmiGraph::mark_labelled(edge_t &e, graph_t &g){
	
	
	// update tracking integers
	g[e].label=labelled;
	
}

void AmiGraph::sort_labelled_unlabelled_adjacent(graph_t &g, vertex_t &vin,  vertex_vector_t &labelled_adj_vertices,vertex_vector_t &unlabelled_adj_vertices, edge_vector_t &labelled_edges, edge_vector_t &unlabelled_edges){
	
	
////std::cout<<"sorting labelled/unlabelled"<<std::endl;	
	
// reset vectors of labelled/unlabelled sorting
labelled_adj_vertices.clear();
unlabelled_adj_vertices.clear();
labelled_edges.clear();
unlabelled_edges.clear();

// in principle this should always be true
// if(g[edge(last,vin,g).first].label==labelled){
// labelled_adj_vertices.push_back(last);
// labelled_edges.push_back(edge(last,vin,g).first);
// }else{
// unlabelled_adj_vertices.push_back(last);
// unlabelled_edges.push_back(edge(last,vin,g).first);
// }	
	
boost::graph_traits<graph_t>::in_edge_iterator iei,iei_end;
boost::graph_traits<graph_t>::out_edge_iterator oei,oei_end;	
	
for (boost::tie(iei,iei_end) = in_edges(vin, g); iei != iei_end; ++iei){	
	
if (g[*iei].label==labelled){
labelled_adj_vertices.push_back(source(*iei,g));
labelled_edges.push_back(*iei);
}else{
unlabelled_adj_vertices.push_back(source(*iei,g));
unlabelled_edges.push_back(*iei);
}

}

for (boost::tie(oei,oei_end) = out_edges(vin, g); oei != oei_end; ++oei){	
	
if (g[*oei].label==labelled){
labelled_adj_vertices.push_back(target(*oei,g));
labelled_edges.push_back(*oei);
}else{
unlabelled_adj_vertices.push_back(target(*oei,g));
unlabelled_edges.push_back(*oei);
}


}



	
// boost::graph_traits<graph_t>::adjacency_iterator ai, ai_end;

// for (boost::tie(ai,ai_end) = adjacent_vertices(vin, g); ai != ai_end; ++ai){
	
// if( edge(vin, *ai, g).second){

// if(g[edge(vin,*ai,g).first].label==labelled){
// labelled_adj_vertices.push_back(target(edge(vin,*ai,g).first,g));
// labelled_edges.push_back(edge(vin,*ai,g).first);
	
// }else{
// unlabelled_adj_vertices.push_back(target(edge(vin,*ai,g).first,g));
// unlabelled_edges.push_back(edge(vin,*ai,g).first);
	
// }


// }else
	// if (edge(*ai, vin, g).second){
		// if(g[edge(*ai, vin,g).first].label==labelled){
	// labelled_adj_vertices.push_back(source(edge(*ai, vin,g).first,g));
	// labelled_edges.push_back(edge(*ai, vin,g).first);
		// }else {
	// unlabelled_adj_vertices.push_back(source(edge(*ai, vin,g).first,g));
	// unlabelled_edges.push_back(edge(*ai, vin,g).first);
			
		// }
	
	// }
// }
	
	
	
	
//std::cout<<"done sorting labelled/unlabelled"<<std::endl;
//std::cout<<"Labelled is "<<  labelled_adj_vertices.size()<<"  unLabelled is "<<  unlabelled_adj_vertices.size()<<std::endl;		
	
	
}

bool AmiGraph::label_consv_momentum(graph_t &g, vertex_t &vin, vertex_vector_t &labelled_adj_vertices, vertex_vector_t &unlabelled_adj_vertices, edge_vector_t &labelled_edges, edge_vector_t &unlabelled_edges){

edge_t bose;	
bool output=true;	
if( unlabelled_edges.size()==3){
	//std::cout<<"found 3 unlabelled edges"<<std::endl;
	output=false; return output;} // immediately break function since we can't label.
if( unlabelled_edges.size()==0){
	//std::cout<<"no unlabelled edges"<<std::endl;
	output=false; return output;}  // immediately break function since nothing to label.
	
if (unlabelled_edges.size()==2){
		// sanity check
		if (labelled_edges.size()!=1){  
		std::cout<<"Labelled are "<< labelled_edges.size()<<"unlabelled are "<< unlabelled_edges.size()<<std::endl;
		print_all_edge_info(g);
		throw std::runtime_error("In Labelling: Number of labelled plus unlabelled not equal to total number? Can't be true :) ");} 	

// both fermionic lines
	if (g[unlabelled_edges[0]].g_struct_.stat_==AmiBase::Fermi && g[unlabelled_edges[1]].g_struct_.stat_==AmiBase::Fermi){ 	

	// label one of the two fermi edges independently and set the other according to consv of momentum 

	label_indep_edge(unlabelled_edges[0],g);
	assign_cons_label(g,vin,labelled_edges[0], unlabelled_edges[0], unlabelled_edges[1]); // 0 is now labelled
	}


// first fermi second bose
	if (g[unlabelled_edges[0]].g_struct_.stat_==AmiBase::Fermi && g[unlabelled_edges[1]].g_struct_.stat_==AmiBase::Bose){

	// bose=unlabelled_edges[1];
	// if(is_labelled_bose_zero(bose,g)){ std::cout<<"Found false logic"<<std::endl;
 	// add_labels(unlabelled_edges[0],labelled_edges[0],g);
	// label_indep_epsilon(unlabelled_edges[0],g);
	// g[unlabelled_edges[0]].label=labelled;
    // }else{
	
	label_indep_edge(unlabelled_edges[0],g);
	assign_cons_label(g,vin,labelled_edges[0], unlabelled_edges[0], unlabelled_edges[1]); 
	// }
	 }	

// first bose second fermi 
	if (g[unlabelled_edges[0]].g_struct_.stat_==AmiBase::Bose && g[unlabelled_edges[1]].g_struct_.stat_==AmiBase::Fermi){

	// bose=unlabelled_edges[0];
	// if(is_labelled_bose_zero(bose,g)){std::cout<<"Found false logic"<<std::endl;
 	// add_labels(unlabelled_edges[1],labelled_edges[0],g);
	// label_indep_epsilon(unlabelled_edges[1],g);
	// g[unlabelled_edges[0]].label=labelled;
    // }else{	
	label_indep_edge(unlabelled_edges[1],g);
	assign_cons_label(g,vin,labelled_edges[0], unlabelled_edges[1], unlabelled_edges[0]); 
	// }
	 }	

}


if (unlabelled_edges.size()==1){

		if (labelled_edges.size()!=2){  
		//std::cout<<"Labelled are "<< labelled_edges.size()<<"unlabelled are "<< unlabelled_edges.size()<<std::endl;
		throw std::runtime_error("In Label_consv: Number of labelled plus unlabelled not equal to total number? Can't be true :) ");} 

assign_cons_label( g, vin, labelled_edges[0], labelled_edges[1], unlabelled_edges[0]);


}	
	
	
	
	
return output;	
}	





// This function takes two labelled fermi lines at a three point vertex and assigns a label to the bosonic line 
// this has assigned only the momentum conserving alpha. the epsilon is just like the others.
void AmiGraph::assign_cons_label(graph_t &g, vertex_t &vin, edge_t &fermi_one, edge_t &fermi_two, edge_t &bose){
	
	
	// std::cout<<print_edge(fermi_one,g)<<" "<<print_edge(fermi_two,g)<<" "<<print_edge(bose,g)<<" "<<std::endl;
	
	if(fermi_one==fermi_two){
	
mark_labelled(bose,g);
return;	
		
	}
	
	
	if(source(bose,g)==vin){
		
		if(target(fermi_one,g)==vin ){
			// std::cout<<"Adding"<<std::endl;
			add_labels(bose,fermi_one,g);
		}
		else{
			subtract_labels(bose,fermi_one, g);
		}
		
		if (target(fermi_two,g)==vin ){
			add_labels(bose,fermi_two,g);
			// std::cout<<"Adding"<<std::endl;
			}else {
				subtract_labels(bose,fermi_two,g);
			}
		
		
	}else  // don't need the following if but for debugging i'm going to put it.
		if (target(bose,g)==vin){
			
		if(target(fermi_one,g)==vin){
			subtract_labels(bose,fermi_one,g);
		}
		else{
			add_labels(bose,fermi_one, g);
			// std::cout<<"Adding"<<std::endl;
		}
		
		if (target(fermi_two,g)==vin){
			subtract_labels(bose,fermi_two,g);
			}else {
				add_labels(bose,fermi_two,g);
				// std::cout<<"Adding"<<std::endl;
			}
			
			
			
			
		}
	
	if(g[bose].g_struct_.stat_==AmiBase::Fermi){
	label_indep_epsilon(bose,g);}else{mark_labelled(bose,g);}
	
	
	
}

// this only checks that each label is in the other graph. it won't fail on multiple poles, but won't pick up the difference between two graphs with different numbers of G's
// should only be used on two graphs which are in fact Equal in graph structure - not just isomorphic
bool AmiGraph::edge_labels_are_equal(edge_t &one,graph_t &g1,  edge_t &two, graph_t &g2){

for (int i=0; i< g1[one].g_struct_.alpha_.size(); i++){
if(g1[one].g_struct_.alpha_[i]!=g2[two].g_struct_.alpha_[i]){return false;}
}	
for (int i=0; i< g1[one].g_struct_.eps_.size(); i++){
if(g1[one].g_struct_.eps_[i]!=g2[two].g_struct_.eps_[i]){return false;}
}


return true;	
}



bool AmiGraph::edge_labels_are_equal(edge_t &one, edge_t &two, graph_t &g){

for (int i=0; i< g[one].g_struct_.alpha_.size(); i++){
if(g[one].g_struct_.alpha_[i]!=g[two].g_struct_.alpha_[i]){return false;}
}	
for (int i=0; i< g[one].g_struct_.eps_.size(); i++){
if(g[one].g_struct_.eps_[i]!=g[two].g_struct_.eps_[i]){return false;}
}


return true;	
}

bool AmiGraph::edge_alphas_are_equal(edge_t &one, edge_t &two, graph_t &g){

for (int i=0; i< g[one].g_struct_.alpha_.size(); i++){
if(g[one].g_struct_.alpha_[i]!=g[two].g_struct_.alpha_[i]){return false;}
}	

return true;	
}

bool AmiGraph::edge_alphas_are_negative(edge_t &one, edge_t &two, graph_t &g){

for (int i=0; i< g[one].g_struct_.alpha_.size(); i++){
if(g[one].g_struct_.alpha_[i]!=-g[two].g_struct_.alpha_[i]){return false;}
}	

return true;	
}


void AmiGraph::add_labels(edge_t &one, edge_t &two, graph_t &g){
	
for (int i=0; i< g[one].g_struct_.alpha_.size(); i++){
g[one].g_struct_.alpha_[i]+=g[two].g_struct_.alpha_[i];
}	
	
}

void AmiGraph::subtract_labels(edge_t &one, edge_t &two, graph_t &g){
	
for (int i=0; i< g[one].g_struct_.alpha_.size(); i++){
g[one].g_struct_.alpha_[i]-=g[two].g_struct_.alpha_[i];
}	
	
}


// TODO: Issue. what if you give this a graph where the alphas and epsilons are empty? This won't work then.
void AmiGraph::label_extern_legs(edge_vector_t &extern_vect_list,graph_t &g){
	
	// if(graph_type != AmiBase::Bose){ 
	edge_t e1;
	
	// std::cout<<"Supposed to be labelling legs here"<<std::endl;
	
	for (int i=0; i< extern_vect_list.size(); i++){
		
		// std::cout<<"Supposed to be labelling legs here at i="<<i<<std::endl;
		e1=extern_vect_list[i];
		g[e1].g_struct_.alpha_[g[e1].g_struct_.alpha_.size()-1]=1;
	    
		// TODO these lines are hardcoded for single external lines
		if(graph_type==AmiBase::density || graph_type==AmiBase::DOS || graph_type==AmiBase::Greens|| graph_type==AmiBase::ENERGY){
		g[e1].g_struct_.eps_[g[e1].g_struct_.eps_.size()-1]=1; 
		}
		
		// external legs have no epsilon value??? we set the last epsilon for the external legs 
		
		g[e1].label=labelled;
		
		
	}
	
	// }
	// else{
		
		// edge_t e1;
	
	// std::cout<<"Supposed to be labelling bose legs here"<<std::endl;
		
	// e1=extern_vect_list[i];
	// g[e1].g_struct_.alpha_[g[e1].g_struct_.alpha_.size()-1]=1;	
		
		
		
		
	// }
	
	//g[boost::graph_bundle].n_labelled++;// we don't increment for external legs
	
	 // TODO this is not general - since they all get the same label if there are not two external legs this won't be right.	
	
	
}

void AmiGraph::random_labelling(graph_t &g, bool &result){
reset_g(g);
// initialize some counters on the interior of the graph
g[boost::graph_bundle].n_indep=0;
g[boost::graph_bundle].n_labelled=0;

// label the external legs
vertex_vector_t extern_vect_list;
edge_vector_t extern_edge_list;
find_external_vertices(g, extern_vect_list, extern_edge_list);


// if(ami_parameters.TYPE_!=AmiBase::density && ami_parameters.TYPE_ != AmiBase::doubleocc){
if(extern_vect_list.size()<1){  throw std::runtime_error("Can't find external leg.");}
// }

label_extern_legs(extern_edge_list,g);

vertex_vector_t to_label;
edge_vector_t unlab_fermi;	
find_unlabelled_fermionic_edges(g, unlab_fermi);

// fermionic loop stuff vs label along line stuff
vertex_vector_t labelled_adj_vertices, unlabelled_adj_vertices;
edge_vector_t labelled_edges, unlabelled_edges;	

int ord=graph_order(g);

for(int m=0; m<ord; m++){

int ran=random_int(0,unlab_fermi.size()-1);
if(ran>=0){
label_indep_edge(unlab_fermi[ran],g);	
unlab_fermi.erase(unlab_fermi.begin(), unlab_fermi.begin()+ran);	
	
}

}
// now we have assigned at random 'ord' number of independent fermionic edges.

// find edges with 2 labelled, and conserve momentum.  Repeat until no more with 2.  Check unlabelled edges=0, else fail.  check momentum conservation else fail 

// ask are there any unlabelled fermionic edges?

do{
to_label.clear();
unlab_fermi.clear();
find_unlabelled_fermionic_edges(g, unlab_fermi);
collect_vertices_to_label(g, to_label);		

// label the to_label list 
for(int i=0; i< to_label.size(); i++){
sort_labelled_unlabelled_adjacent(g, to_label[i],labelled_adj_vertices, unlabelled_adj_vertices, labelled_edges, unlabelled_edges);
bool junk=label_consv_momentum(g, to_label[i], labelled_adj_vertices, unlabelled_adj_vertices, labelled_edges, unlabelled_edges);
}	
}while(to_label.size()!=0);

// finally - check for unlabelled edges
unlab_fermi.clear();
find_unlabelled_fermionic_edges(g, unlab_fermi);

if(unlab_fermi.size()!=0){result=false;}

if(unlab_fermi.size()==0){check_momentum_conservation(g, result);}	
	
	
}

void AmiGraph::print_label_counts(std::vector<int> in){
	
std::cout<<"( ";

for(int i=0; i< in.size(); i++){

std::cout<<in[i]<<" ";
}	

std::cout<<")"<<std::endl;
	
	
}

void  AmiGraph::label_counts(graph_t &g, std::vector<int> &out){
	out.clear();
edge_vector_t internal;
find_internal_fermionic_edges(g, internal);

out.resize(g[internal[0]].g_struct_.alpha_.size());

for(int i=0; i< internal.size(); i++){

for(int m=0; m< g[internal[i]].g_struct_.alpha_.size(); m++){

out[m]+=std::abs(g[internal[i]].g_struct_.alpha_[m]);
	
	
}

}	
	
	
	
}

void AmiGraph::repeated_labelling(graph_t &g, bool &result){
int iter=0;
result=false;

number_vertices(g);
// std::cout<<"Num vertices is "<<num_vertices(g)<<std::endl;
// std::cout<<"Num edges is "<<num_edges(g)<<std::endl;

/* 
for(int i=0; i<10; i++){
	i++;
label_systematic(g);	
check_momentum_conservation(g, result);
if(result==true){
	iter+=i;
	break;
	}
} */

if(result==false){
	// std::cout<<"Starting half-random labelling"<<std::endl;
// for(int i=0; i<2000;i++){
do{	
//label_systematic_redone(g);

// std::cout<<"Next attempt"<<std::endl;

label_half_random(g);
check_momentum_conservation(g, result);
edge_vector_t unlabelled;
find_unlabelled_fermionic_edges(g, unlabelled);
if(unlabelled.size()!=0){result=false;}
}while(result==false);
}
/* 
if(result==false){
	// std::cout<<"Starting half-random labelling"<<std::endl;
for(int i=0; i<20000;i++){
	
//label_systematic_redone(g);
random_labelling(g,result);
// check_momentum_conservation(g, result);
if(result==true){
	iter+=i;
	break;
	}
}
} */

fix_epsilons(g);

if(result==false){std::cerr<<"Warning: Labelling Failed."<<std::endl;}
//std::cout<<"Labelling took iter="<<iter<<std::endl;
}

void AmiGraph::check_momentum_conservation(graph_t &g, bool &result){
	
bool cons=true;
bool vert_cons=true;
	
boost::graph_traits<graph_t>::vertex_iterator vi,vi_end;

boost::graph_traits<graph_t>::in_edge_iterator iei,iei_end;
boost::graph_traits<graph_t>::out_edge_iterator oei,oei_end;	


std::vector<int> sum_vec(g[random_edge(g,rand_gen)].g_struct_.alpha_.size(),0);

//std::cout<<std::endl;
//std::cout<<"Checking momentum conservation"<<std::endl;

for (boost::tie(vi,vi_end) = vertices(g); vi != vi_end; ++vi){
vert_cons=true;
//std::cout<<"On vertex "<<g[*vi].index_ <<std::endl;
std::fill(sum_vec.begin(), sum_vec.end(), 0);


	if(degree(*vi,g)!=1){
		//std::cout<<"On vertex "<<g[*vi].index_ << " checking conserv"<<std::endl;

		for (boost::tie(iei,iei_end) = in_edges(*vi,g); iei != iei_end; ++iei){
		for (int i=0; i< g[*iei].g_struct_.alpha_.size(); i++){
		//	//std::cout<<"in edges i= "<<i<<std::endl;
		sum_vec[i]+=g[*iei].g_struct_.alpha_[i];
		
		}
		}

		for (boost::tie(oei,oei_end) = out_edges(*vi,g); oei != oei_end; ++oei){
		
		for (int i=0; i< g[*oei].g_struct_.alpha_.size(); i++){
		//	//std::cout<<"in edges i= "<<i<<std::endl;
		sum_vec[i]-=g[*oei].g_struct_.alpha_[i];
		
		}
		
		

			
		}	
			
	}

for (int i=0; i< sum_vec.size(); i++){

if (sum_vec[i]!=0){
//std::cout<<"Momentum not conserved at vertex "<< g[*vi].index_ <<std::endl;	

cons=false;
vert_cons=false;
}

}	
// if(vert_cons==false){	
// print_edgeson_vert(g[*vi].index_);	
// }
}


//if(cons==true){std::cout<<"Momentum was conserved everywhere "<<std::endl;}
//if(cons==false){std::cout<<"Momentum was NOT conserved  "<<std::endl;}	


result=cons;	
}

int AmiGraph::count_labelled_edges_of_vert(vertex_t &v, graph_t &g){
	
int output=0;	
	
boost::graph_traits<graph_t>::in_edge_iterator iei,iei_end;
boost::graph_traits<graph_t>::out_edge_iterator oei,oei_end;
	
	for (boost::tie(iei,iei_end) = in_edges(v,g); iei != iei_end; ++iei){

		if(g[*iei].label==labelled){output++;}
	
	}
	
	for (boost::tie(oei,oei_end) = out_edges(v,g); oei != oei_end; ++oei){
	if(g[*oei].label==labelled){output++;}
	
	}
return output;	
}

void AmiGraph::collect_vertices_to_label(graph_t &g, vertex_vector_t &vert_list){
vert_list.clear();	
int num=0;
boost::graph_traits<graph_t>::vertex_iterator vi, vi_end;

for (boost::tie(vi,vi_end)= vertices(g); vi  != vi_end; ++vi){

num=count_labelled_edges_of_vert(*vi,g);

if(num==2){ vert_list.push_back(*vi);}

}	
//std::cout<<"Found some easy vertices "<<vert_list.size()<<std::endl;

}



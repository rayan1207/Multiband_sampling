#include "amigraph.hpp"


// need function to systematically generate labels until we have some specific number

void AmiGraph::trim_label_set(labels_t &L, int trim){

if(L.size()>trim){	
do{

int rand=random_int(0, L.size()-1);
L.erase(L.begin()+rand);
}while(L.size()>trim);	

	
L.shrink_to_fit();	
}
	
}

void AmiGraph::sys_label_sets(labels_t &L, graph_t g, bool &result, int max){

int ord=graph_order(g);

ord=graph_order(g);	
// ord is one bigger
if( (graph_type == AmiBase::Pi_phuu) || (graph_type == AmiBase::Pi_phud) || (graph_type==AmiBase::doubleocc)|| graph_type ==AmiBase::Pi_ppud || graph_type ==AmiBase::Pi_ppuu || graph_type ==AmiBase::FORCE){
ord=graph_order(g)+1;
}
// if( graph_type==AmiBase::density){
	// ord=graph_order(g)+1;
// }
// if(graph_type==AmiBase::doubleocc){
	// ord=graph_order(g)+2;
// }


bool has_external_legs=true;
// if a force graph remove the external legs completely.

vertex_vector_t in_vv,out_vv;
if(graph_type==AmiBase::FORCE){
	
	// std::cout<<"Looking for LR vertices"<<std::endl;
	// print_all_edge_info(g);

find_force_LR_vertices(g, in_vv, out_vv);
	
vertex_vector_t extern_vect_list;// warning there is a naming duplication below 
edge_vector_t extern_edge_list;	
find_external_vertices(g, extern_vect_list, extern_edge_list);
delete_legs(g,extern_vect_list, extern_edge_list);
has_external_legs=false;
// std::cout<<"Legs removed check"<<std::endl;
}





std::vector<int> v;

result=false;

	
reset_g(g);	

// std::cout<<"Removed labels"<<std::endl;
// print_all_edge_info(g);
// initialize some counters on the interior of the graph
g[boost::graph_bundle].n_indep=0;
g[boost::graph_bundle].n_labelled=0;

// label the external legs



vertex_vector_t extern_vect_list;
edge_vector_t extern_edge_list;

if(has_external_legs){
find_external_vertices(g, extern_vect_list, extern_edge_list);
// if(ami_parameters.TYPE_!=AmiBase::density && ami_parameters.TYPE_!= AmiBase::doubleocc){
if(extern_vect_list.size()<1){  throw std::runtime_error("Can't find external leg.");}
// }
// if(extern_vect_list.size()<1){  throw std::runtime_error("Can't find external leg.");}
label_extern_legs(extern_edge_list,g);
}


edge_vector_t internal_F_edges;
find_internal_fermionic_edges(g,internal_F_edges);

vertex_vector_t internal_v;
find_internal_vertices(g, internal_v);

// fermionic loop stuff vs label along line stuff
vertex_vector_t labelled_adj_vertices, unlabelled_adj_vertices;
edge_vector_t labelled_edges, unlabelled_edges;

int count=0;

// take entire list of F edges. 
for(int i=0; i< internal_F_edges.size(); i++){
v.push_back(i);	
}


// if(scramble){
// random_shuffle(v.begin(), v.end());
// }

// std::cout<<"Using list "<<std::endl;
// for(int i=0; i<v.size(); i++){
// std::cout<<v[i]<<" ";	
// }
// std::cout<<std::endl;


do{
	// std::cout<<"Start of do "<<count<<std::endl;
	count++;
reset_g(g);	

// std::cout<<"G was reset and is "<<std::endl;
// print_all_edge_info(g);

if(has_external_legs){
label_extern_legs(extern_edge_list,g);
}
// std::cout<<"-----------"<<std::endl;
// print_all_edge_info(g);
// std::cout<<"-----------"<<std::endl;
// initialize some counters on the interior of the graph
g[boost::graph_bundle].n_indep=0;
g[boost::graph_bundle].n_labelled=0;	

// std::cout<<"On permutation list "<<std::endl;
// std::cout<<"About to assign independent edges of size "<<ord<<std::endl;
for(int j=0; j< ord; j++){
// std::cout<<v[j]<<" ";	
label_indep_edge(internal_F_edges[v[j]],g);
// std::cout<<"Labelling independent edge "<< v[j]<<std::endl;
}
// std::cout<<std::endl;
// now we have the independent edges labelled. So now we iterate through vertices, looking for a vertex with two edges labelled 

// print_all_edge_info(g);

bool cont=true;
edge_vector_t in, out;
int first,second;

do{
	find_unlabelled_edges(g, in);
	first=in.size();

for(int i=0; i<internal_v.size();i++){
	
// get labelled and unlabelled edges attached to the i'th vertex 
sort_labelled_unlabelled_adjacent(g, internal_v[i],labelled_adj_vertices, unlabelled_adj_vertices, labelled_edges, unlabelled_edges);

if( unlabelled_edges.size()!=1){continue;} // immediately break function since we can't label.

// Can only now do something if 2 of three edges are labelled.
// if (unlabelled_edges.size()==1){

		if (labelled_edges.size()!=2){  
		// std::cout<<"Labelled are "<< labelled_edges.size()<<"unlabelled are "<< unlabelled_edges.size()<<std::endl;
		throw std::runtime_error("In Syslabel: Number of labelled plus unlabelled not equal to total number? Can't be true :) ");} 

// std::cout<<"Assigning conserving label to edge "<<print_edge(unlabelled_edges[0],g)<<  std::endl;
		assign_cons_label( g, internal_v[i], labelled_edges[0], labelled_edges[1], unlabelled_edges[0]);

// }

}

find_unlabelled_edges(g, out);
second=out.size();	

// std::cout<<first<<" "<<second<<std::endl;

if(first==second){cont=false;}


}while(cont);

// std::cout<<"Out of do loop"<<std::endl;
// print_all_edge_info(g);

// at this point the graph cannot be labelled further. check two things

// 1) Has every line been labelled? aka, out.size()==0
// 2) is momentum conserved. 

// std::cout<<"out.size()="<<out.size()<<std::endl;
bool append=false;
if(out.size()==0){

check_momentum_conservation(g, append);
// std::cout<<"Momentum conserved? "<<append<<std::endl;	

if(append){

fix_epsilons(g);	

if(graph_type==AmiBase::FORCE){
	// std::cout<<"Fixing force labels with and and out sizes "<<in_vv.size()<<" "<< out_vv.size()<<std::endl;
fix_force_labels(g,in_vv, out_vv);
// std::cout<<"Done"<<std::endl;
// put_back_legs(g,in_vv, out_vv);// don't need to put the legs back here 
}	
	// std::cout<<"Pushing into L on graph order "<< graph_order(g) <<std::endl;
L.push_back(g);	
result=true;

if(L.size()>max){return;}

}	

}

// break out of do loop if
// 1) size of L is > max 
// std::cout<<"-----------"<<std::endl;
// print_all_edge_info(g);
// std::cout<<"-----------"<<std::endl;

	 std::reverse(v.begin()+ord, v.end());
}while(std::next_permutation(v.begin(),v.end()));	
	
// std::cout<<"Ran for count="<<count<<std::endl;	
	
}


void AmiGraph::sys_label(graph_t &g, bool &result){

// std::cout<<"Systematic bubble labelling for graph"<<std::endl;
// print_all_edge_info(g);

graph_t gc;
gc=g;
int ord;

result=false;

ord=graph_order(gc);
// std::cout<<"Ord is "<<ord<<std::endl;
// std::cout<<"Graph type is "<<graph_type<<std::endl;
// ord is one bigger

if( (graph_type == AmiBase::Pi_phuu) || (graph_type == AmiBase::Pi_phud) || (graph_type==AmiBase::doubleocc)|| graph_type==AmiBase::Pi_ppud || graph_type==AmiBase::Pi_ppuu || graph_type==AmiBase::FORCE){
	// std::cout<<"Increasing order by 1"<<std::endl;
ord=graph_order(gc)+1;
}

bool has_external_legs=true;
// if a force graph remove the external legs completely.

vertex_vector_t in_vv,out_vv;
if(graph_type==AmiBase::FORCE){

find_force_LR_vertices(gc, in_vv, out_vv);
	
vertex_vector_t extern_vect_list;
edge_vector_t extern_edge_list;	
find_external_vertices(gc, extern_vect_list, extern_edge_list);
delete_legs(gc,extern_vect_list, extern_edge_list);
has_external_legs=false;
}



// if( graph_type==AmiBase::density){
	// ord=graph_order(gc)+1;
// }
// if(graph_type==AmiBase::doubleocc){
	// ord=graph_order(gc)+2;
// }


// std::cout<<"Ord is "<<ord<<std::endl;



std::vector<int> v;

reset_g(gc);


// print_all_edge_info(gc);

// if(ord==1){
	
// }


	
// label_indep_edge(internal_F_edges[v[j]],gc);	
	
// }

// std::cout<<"Graph was reset and is now"<<std::endl;
// print_all_edge_info(gc);




// label the external legs
vertex_vector_t extern_vect_list;
edge_vector_t extern_edge_list;


if(has_external_legs){
// std::cout<<"Looking for external "<<std::endl;
find_external_vertices(gc, extern_vect_list, extern_edge_list);
// std::cout<<"exited"<<std::endl;
// std::cout<<"Extern edge list of size "<<extern_vect_list.size()<<std::endl;
// if(ami_parameters.TYPE_!=AmiBase::density && ami_parameters.TYPE_!= AmiBase::doubleocc){
if(extern_vect_list.size()<1){  throw std::runtime_error("Can't find external leg.");}

}
// }
// if(extern_vect_list.size()<1){  throw std::runtime_error("Can't find external leg.");}

// TODO: I think this is redundant 
// label_extern_legs(extern_edge_list,gc);


// std::cout<<std::endl<<std::endl<<"-----external labelled-------"<<std::endl;
// print_all_edge_info(gc);

///



edge_vector_t internal_F_edges;
find_internal_fermionic_edges(gc,internal_F_edges);


// std::cout<<"-----Prelabelled------"<<std::endl;
// print_all_edge_info(gc);
// std::cout<<"-----------"<<std::endl;

// std::cout<<"Number of internal edges is "<< internal_F_edges.size()<<std::endl;


vertex_vector_t internal_v;
find_internal_vertices(gc, internal_v);

// std::cout<<"Number of internal vertices is "<< internal_v.size()<<std::endl;

// fermionic loop stuff vs label along line stuff
vertex_vector_t labelled_adj_vertices, unlabelled_adj_vertices;
edge_vector_t labelled_edges, unlabelled_edges;

int count=0;

// take entire list of F edges. 
for(int i=0; i< internal_F_edges.size(); i++){
v.push_back(i);	
}


///
do{
	count++;
reset_g(gc);

// std::cout<<"count is "<<count<<std::endl;


// std::cout<<"-----Prelabelled------"<<std::endl;
// print_all_edge_info(gc);
// std::cout<<"-----------"<<std::endl;

if(has_external_legs){
label_extern_legs(extern_edge_list,gc);
}
// std::cout<<"-----------"<<std::endl;
// print_all_edge_info(gc);
// std::cout<<"-----------"<<std::endl;
// initialize some counters on the interior of the graph
gc[boost::graph_bundle].n_indep=0;
gc[boost::graph_bundle].n_labelled=0;	

// 2022-03-25 : Added tadpole labelling 
// Need function to find and label tadpoles
// vertex_vector_t tp_vert, tp_conn_vert; 
// edge_vector_t tp_bose_edges;
// label_and_find_tadpoles_ami(gc,tp_vert, tp_conn_vert, tp_bose_edges);

	// vertex_vector_t tp_vec;
	// edge_vector_t tp_bose_edges;
	// vertex_vector_t tp_conn_vec;

// edge_vector_t edges_a, edges_b;

// find_tadpoles(gc, tp_vec, tp_conn_vec, tp_bose_edges, edges_a, edges_b);
// std::cout<<"Found number of tadpoles="<<tp_bose_edges.size()<<std::endl;
// for(int m=0; m<tp_bose_edges.size(); m++){
	// mark_labelled(tp_bose_edges[m],gc);
// }
// vertex_vector_t tp_vert, tp_conn_vert; 
// edge_vector_t tp_bose_edges;

// print_all_edge_info(gc);

// find_and_label_tadpole_strings(gc,tp_vert, tp_conn_vert, tp_bose_edges);

// std::cout<<"After tadpole labelling "<<std::endl;
// print_all_edge_info(gc);
// exit(0);

// how is this even systematic?
for(int j=0; j< ord; j++){

label_indep_edge(internal_F_edges[v[j]],gc);
// std::cout<<"Labelling independent edge "<< v[j]<<std::endl;
// print_all_edge_info(gc);
}

// std::cout<<"Indep edges labelled"<<std::endl;
// print_all_edge_info(gc);
// now we have the independent edges labelled. So now we iterate through vertices, looking for a vertex with two edges labelled 

bool cont=true;
edge_vector_t in, out;
int first,second;

// print_all_edge_info(gc);

do{
	find_unlabelled_edges(gc, in);
	first=in.size();

for(int i=0; i<internal_v.size();i++){
	
// get labelled and unlabelled edges attached to the i'th vertex 
sort_labelled_unlabelled_adjacent(gc, internal_v[i],labelled_adj_vertices, unlabelled_adj_vertices, labelled_edges, unlabelled_edges);

if( unlabelled_edges.size()!=1){continue;} // immediately break function since we can't label.

// Can only now do something if 2 of three edges are labelled.
// if (unlabelled_edges.size()==1){

		if (labelled_edges.size()!=2){

			continue; // immediatly break instead of error?
		// std::cout<<"Labelled are "<< labelled_edges.size()<<"unlabelled are "<< unlabelled_edges.size()<<std::endl;
		// print_all_edge_info(gc);
		// std::cout<<"In and out vertices are"<<std::endl;
		// std::cout<<gc[in_vv[0]].index_<<" "<<gc[in_vv[1]].index_<<" "<< gc[out_vv[0]].index_<<" "<<gc[out_vv[1]].index_<<std::endl;
		throw std::runtime_error("Syslabel: Number of labelled plus unlabelled not equal to total number? Can't be true :) ");} 

// std::cout<<"Assigning conserving label to edge "<<print_edge(unlabelled_edges[0],g)<<  std::endl;
// std::cout<<"Labelled edges are "<<print_edge(labelled_edges[0],g)<<" "<<print_edge(labelled_edges[1],g)<<std::endl;
		assign_cons_label( gc, internal_v[i], labelled_edges[0], labelled_edges[1], unlabelled_edges[0]);
// std::cout<<"Done"<<std::endl;
// }

}

find_unlabelled_edges(gc, out);
second=out.size();	

// std::cout<<first<<" "<<second<<std::endl;

if(first==second){cont=false;


}


}while(cont);

// at this point the graph cannot be labelled further. check two things

// 1) Has every line been labelled? aka, out.size()==0
// 2) is momentum conserved. 

// std::cout<<"out.size()="<<out.size()<<std::endl;
bool append=false;
find_unlabelled_edges(gc, out);
if(out.size()==0){

check_momentum_conservation(gc, append);
// std::cout<<"Momentum conserved? "<<append<<std::endl;

// print_all_edge_info(gc);	

if(append){

fix_epsilons(gc);	
	// std::cout<<"Pushing into L"<<std::endl;
	
	
if(graph_type==AmiBase::FORCE){
  std::cout<<"Fixing force labels"<<std::endl;
fix_force_labels(gc,in_vv, out_vv);
put_back_legs(gc,in_vv, out_vv);
}	

// std::cout<<"Final print before copy"<<std::endl;
// print_all_edge_info(g);
// std::cout<<std::endl;
// print_all_edge_info(gc);
	
g=gc;
result=true;



// std::cout<<"-----Postlabelled------"<<std::endl;
// print_all_edge_info(gc);
// std::cout<<"-----------"<<std::endl;




return;
}
}
// std::cout<<"one..";
	 std::reverse(v.begin()+ord, v.end());
// std::cout<<"two"<<std::endl;	 
}while(std::next_permutation(v.begin(),v.end()));	





}


/* 
// systematically moves through labelling until it finds one label and then exits. This will often return the same label over and over.
void AmiGraph::sys_label(graph_t &g, bool &result){

// std::cout<<"Sys label start"<<std::endl;

graph_t gc;
gc=g;	

int ord=graph_order(gc);
std::vector<int> v;

	
reset_g(gc);	

std::cout<<"Removed labels"<<std::endl;
print_all_edge_info(g);
// initialize some counters on the interior of the graph
gc[boost::graph_bundle].n_indep=0;
gc[boost::graph_bundle].n_labelled=0;

// label the external legs
vertex_vector_t extern_vect_list;
edge_vector_t extern_edge_list;
find_external_vertices(gc, extern_vect_list, extern_edge_list);
if(extern_vect_list.size()<1){  throw std::runtime_error("Can't find external leg.");}
label_extern_legs(extern_edge_list,gc);


edge_vector_t internal_F_edges;
find_internal_fermionic_edges(gc,internal_F_edges);

vertex_vector_t internal_v;
find_internal_vertices(gc, internal_v);

// fermionic loop stuff vs label along line stuff
vertex_vector_t labelled_adj_vertices, unlabelled_adj_vertices;
edge_vector_t labelled_edges, unlabelled_edges;

int count=0;

// take entire list of F edges. 
for(int i=0; i< internal_F_edges.size(); i++){
v.push_back(i);	
}

do{
	count++;
reset_g(gc);	
label_extern_legs(extern_edge_list,gc);
std::cout<<"-----------"<<std::endl;
print_all_edge_info(gc);
std::cout<<"-----------"<<std::endl;
// initialize some counters on the interior of the graph
gc[boost::graph_bundle].n_indep=0;
gc[boost::graph_bundle].n_labelled=0;	

for(int j=0; j< ord; j++){

label_indep_edge(internal_F_edges[v[j]],gc);
std::cout<<"Labelling independent edge "<< v[j]<<std::endl;
}
// now we have the independent edges labelled. So now we iterate through vertices, looking for a vertex with two edges labelled 

bool cont=true;
edge_vector_t in, out;
int first,second;

do{
	find_unlabelled_edges(gc, in);
	first=in.size();

for(int i=0; i<internal_v.size();i++){
	
// get labelled and unlabelled edges attached to the i'th vertex 
sort_labelled_unlabelled_adjacent(gc, internal_v[i],labelled_adj_vertices, unlabelled_adj_vertices, labelled_edges, unlabelled_edges);

if( unlabelled_edges.size()!=1){continue;} // immediately break function since we can't label.

// Can only now do something if 2 of three edges are labelled.
// if (unlabelled_edges.size()==1){

		if (labelled_edges.size()!=2){  
		std::cout<<"Labelled are "<< labelled_edges.size()<<"unlabelled are "<< unlabelled_edges.size()<<std::endl;
		throw std::runtime_error("Number of labelled plus unlabelled not equal to total number? Can't be true :) ");} 

std::cout<<"Assigning conserving label to edge "<<print_edge(unlabelled_edges[0],g)<<  std::endl;
std::cout<<"Labelled edges are "<<print_edge(labelled_edges[0],g)<<" "<<print_edge(labelled_edges[1],g)<<std::endl;
		assign_cons_label( gc, internal_v[i], labelled_edges[0], labelled_edges[1], unlabelled_edges[0]);
std::cout<<"Done"<<std::endl;
// }

}

find_unlabelled_edges(gc, out);
second=out.size();	

// std::cout<<first<<" "<<second<<std::endl;

if(first==second){cont=false;}


}while(cont);

// at this point the graph cannot be labelled further. check two things

// 1) Has every line been labelled? aka, out.size()==0
// 2) is momentum conserved. 

// std::cout<<"out.size()="<<out.size()<<std::endl;
bool append=false;
if(out.size()==0){

check_momentum_conservation(gc, append);
// std::cout<<"Momentum conserved? "<<append<<std::endl;	

if(append){

fix_epsilons(gc);	
	// std::cout<<"Pushing into L"<<std::endl;
g=gc;
return;
}
}
// std::cout<<"one..";
	 std::reverse(v.begin()+ord, v.end());
// std::cout<<"two"<<std::endl;	 
}while(std::next_permutation(v.begin(),v.end()));	
	
std::cout<<"Ran for count="<<count<<std::endl;		


print_all_edge_info(g);
	
}


 */


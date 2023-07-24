#include "amigraph.hpp"

// TODO - check if there is needless memory duplication in is_isomorphic AND the sub functions


bool AmiGraph::is_isomorphic(graph_t g1, graph_t g2){
	
//	number_vertices(g1);
//	number_vertices(g2);
	
	// print_all_edge_info(g1);
	// std::cout<<"--In isomorphism check--"<< std::endl;
	// print_all_edge_info(g2);
	
bool isbose=false;	
// these can only occur prior to symmetrization of bosonic lines

if(is_pi_graph(g1)!= is_pi_graph(g2)){
	// std::cout<<"One pi one not"<<std::endl;
	return false;
}
if(ami_parameters.TYPE_!=AmiBase::FORCE){
if(!is_pi_graph(g1) && !is_pi_graph(g2)){
	// std::cout<<"Both non-pi"<<std::endl;
if( get_length_principle_line(g1) != get_length_principle_line(g2)){ return false;}
}else{ isbose=true;}
}	
///////////	
// std::cout<<"Made it to here"<<std::endl;
		
		
// print_all_edge_info(g1);		

 // symmetrize_bosonic_lines(g1);
 // symmetrize_bosonic_lines(g2);

// std::cout<<"Symmetrized "<<std::endl;
// print_all_edge_info(g1);		
	
	
// CHANGE: jan 23rd - swapped to only symmetrizing internal lines
// aug 3rd edit - try not symmetrizing if pp bubble
// if(graph_type!=AmiBase::Pi_ppuu && graph_type !=AmiBase::Pi_ppud){
	
/* if(graph_type==AmiBase::Pi_ppuu || graph_type==AmiBase::Pi_ppud){
number_vertices(g1);
number_vertices(g2);



// std::cout<<"pp barcode test resulted in "<< pptest<<std::endl;

bool spincheck=spin_iso(g1,g2);	
// std::cout<<"exited spin iso "<<spincheck<<std::endl;

// if(spincheck){return true;}// if this shows false it could be a false false but true is robust. 	
 
if(!spincheck){

if(full_iso(g1,g2)){
AmiGraph::pp_barcode_t ppb1, ppb2;

pp_barcode(g1, ppb1);

pp_barcode(g2, ppb2);

bool pptest=compare_pp_barcodes(ppb1,ppb2);

return pptest;	
}


return false;
	
 
 
}	 */
	
	
symmetrize_internal_bosonic_lines(g1);
symmetrize_internal_bosonic_lines(g2);	
// }

// std::cout<<"Symmetrized g1"<<std::endl;
// print_all_edge_info(g1);	

// std::cout<<"Symmetrized g2"<< std::endl;
// print_all_edge_info(g2);
	
		
//num nodes1=nodes 2?
// std::cout<<"Checking vertices"<<std::endl;
if( num_vertices(g1)!= num_vertices(g2)){ 
// std::cout<<"Non isomorphic due to vertices"<<std::endl;
return false;}

//num edges 1 = edges 2?
// std::cout<<"Checking edges"<<std::endl;
if( num_edges(g1)!= num_edges(g2)){ 
// std::cout<<"Non isomorphic due to num_edges"<<std::endl;
return false;}




//num F_edges1 = num F_edges 2
// std::cout<<"Checking Fermi edges"<<std::endl;

edge_vector_t fermi_edges1, fermi_edges2;

find_fermionic_edges(g1, fermi_edges1);
find_fermionic_edges(g2, fermi_edges2);

if( fermi_edges1.size() != fermi_edges2.size()){ 
// std::cout<<"Non isomorphic due to num fermi lines"<<std::endl;
return false;}

if(count_tadpoles(g1)!= count_tadpoles(g2)){ return false;}



// num B_edges1=num B_edges 2
// std::cout<<"Checking Bose edges"<<std::endl;
edge_vector_t bose_edges1, bose_edges2;

find_bosonic_edges(g1, bose_edges1);
find_bosonic_edges(g2, bose_edges2);

// std::cout<<"Found edges"<<std::endl;

if( bose_edges1.size() != bose_edges2.size()){ 
// std::cout<<"Non isomorphic due to num bose lines"<<std::endl;
return false;}



if(graph_type!=AmiBase::Pi_ppuu && graph_type !=AmiBase::Pi_ppud && graph_type!=AmiBase::FORCE){
// function to extract the F or B array for each graph 
// std::cout<<"Constructing ForB arrays"<<std::endl;
systematic_vertex_label(g1);
systematic_vertex_label(g2);

// std::cout<<"Systematic complete"<<std::endl;

 std::vector< std::vector<int>> e1, e2;

get_labelled_bose_lists(e1,g1);
get_labelled_bose_lists(e2, g2);


bool first=bose_lists_equal(e1,e2);


if(!first){
	// std::cout<<"Non isomorphic due to first"<<std::endl;
	return false;}



 barcode_t B1, B2;
tree_t T1, T2;
// std::cout<<"Constructing barcodes"<<std::endl;
construct_graph_barcode(g1, B1, T1);
construct_graph_barcode(g2, B2, T2);

// std::cout<<"Checking equal"<<std::endl;
// bool second=true;
bool second=barcodes_are_equal(B1, B2);

// TODO: The question is do you actually WANT the different spin graphs or do you sum over them since they are the same?
if(ami_parameters.int_type_==AmiBase::coulomb){
	second=true;
}


		if(!second){
			// std::cout<<"Non isomorphic due to second "<<std::endl;
			return false;}  
	// std::cout<<"Full iso"<<std::endl;
	bool final=full_iso(g1,g2);
	if(!final){
		// std::cout<<"Non isomorphic in full iso"<<std::endl;
		return false;}

}

if(graph_type==AmiBase::Pi_ppuu || graph_type==AmiBase::Pi_ppud|| graph_type==AmiBase::FORCE){
	
// std::cout<<"In pp stuff"<<std::endl;	
	
symmetrize_spin_up(g1);
symmetrize_spin_up(g2);

// print_all_edge_info(g1);
// print_all_edge_info(g2);	
	
number_vertices(g1);
number_vertices(g2);



// std::cout<<"pp barcode test resulted in "<< pptest<<std::endl;

bool spincheck=full_iso(g1,g2); //spin_iso(g1,g2);	
// std::cout<<"exited spin iso "<<spincheck<<std::endl;

// if(spincheck){return true;}// if this shows false it could be a false false but true is robust. 	
 
if(!spincheck){

// if(full_iso(g1,g2)){
// AmiGraph::pp_barcode_t ppb1, ppb2;

// pp_barcode(g1, ppb1);

// pp_barcode(g2, ppb2);

// bool pptest=compare_pp_barcodes(ppb1,ppb2);

// return pptest;	
// }


return false;
	
 
 
}	
// if(!spincheck){return false;}	

// as a final 


	
	// bool final=full_iso(g1,g2);
	// std::cout<<"exited full iso "<< final<<std::endl;
	// bool spincheck=spin_iso(g1,g2);	
	// std::cout<<"exited spin iso "<<spincheck<<std::endl;
	
	
	// if(!final){
		// std::cout<<"Non isomorphic in full iso "<<std::endl;
		// return false;}
	
	
}


// std::cout<<"End iso"<<std::endl;
return true;
	
}

void AmiGraph::symmetrize_spin_up(graph_t &g){

// std::cout<<"Symetrizing graph before"<<std::endl;

// print_all_edge_info(g);

boost::graph_traits<graph_t>::edge_iterator ei, ei_end;

edge_vector_t add_list;

for(boost::tie(ei,ei_end)=edges(g); ei!=ei_end; ++ei){

if(g[*ei].g_struct_.stat_==AmiBase::Fermi && g[*ei].spin==0){

add_list.push_back(*ei);

// std::cout<<"Adding edge"<<std::endl;
// add_edge( source(*ei,g), target(*ei,g), edge_info(AmiBase::Fermi,g[*ei].fermi_loop_id, g[*ei].spin), g); 	
	
}
	
	
}

for(int i=0; i<add_list.size(); i++){

add_edge( source(add_list[i],g), target(add_list[i],g), edge_info(AmiBase::Fermi,g[add_list[i]].fermi_loop_id, g[add_list[i]].spin), g);	
	
}


// std::cout<<"after"<<std::endl;
// print_all_edge_info(g);

return; 
	
}



bool AmiGraph::is_pi_graph(graph_t &g){

boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
bool isbose=false;
for (boost::tie(ei,ei_end)=edges(g); ei!=ei_end; ++ei){
	
	// std::cout<<"Source and targets are "<< g[source(*ei,g)].index_<<" "<< g[target(*ei,g)].index_<<std::endl;
	// std::cout<<"Degree is "<<degree(source(*ei,g),g)<<std::endl;
	
if(degree(source(*ei,g),g)==1 && (source(*ei,g) != target(*ei,g) )){
	
	
	if(g[*ei].g_struct_.stat_==AmiBase::Bose){ isbose=true;}
	
	 break;}
	
}

return isbose;
	
	
}


bool AmiGraph::bose_lists_equal(std::vector< std::vector<int>> b1, std::vector< std::vector<int>> b2){

if(b1.size()!= b2.size()){return false;}
	
for(int i=0; i< b1.size(); i++){

for(int j=0; j< b2.size(); j++){	

// std::cout<<"Comparing "<< b1[i][0]<<","<<b1[i][1]<<std::endl;
// std::cout<<"with "<< b2[j][0]<<","<<b2[j][1]<<std::endl;

if( b1[i]==b2[j]){
	// std::cout<<"Found equal"<<std::endl;
//b1.erase(b1.begin()+i);
b2.erase(b2.begin()+j);
continue;

}	
	
}
}
	
if( b2.size()==0){ 
//std::cout<<"Found isomorphic"<<std::endl;
return true;}	
	
return false;	
}

int AmiGraph::get_length_principle_line(graph_t &g){
reset_visited(g);
vertex_t current;//=vert_list[0]; // this will not work if vertex zero is no longer the external vertex 
vertex_t next;
vertex_vector_t fermi_vert_list;

boost::graph_traits<graph_t>::edge_iterator ei, ei_end;

for (boost::tie(ei,ei_end)=edges(g); ei!=ei_end; ++ei){
	
if(degree(source(*ei,g),g)==1){current= source(*ei,g); break;}
	
	
}

bool instop=false;
do{
 // std::cout<<"Source is "<<graph[current].index_<<std::endl;	

g[current].visited=1;
//g[current].index_=counter;
//counter++;
fermi_vert_list.push_back(current);
get_adjacent_vertex_stat(g, current, next, AmiBase::Fermi); 
// std::cout<<"Source is "<< graph[current].index_<<" next is "<<graph[next].index_<<std::endl;
 
 if( g[next].visited==1){instop=true; g[next].visited=1;};
 current=next;

} while(!instop);


// std::cout<<"Found principle line length "<< fermi_vert_list.size()<<std::endl;	
return fermi_vert_list.size();	
	
}

// THIS FUNCTION FOR SURE DID NOT WORK PRIOR TO Feb 11 2020
void AmiGraph::get_principle_line(graph_t &g, edge_vector_t &ev){

vertex_t current;//=vert_list[0]; // this will not work if vertex zero is no longer the external vertex 
vertex_t next;
vertex_vector_t fermi_vert_list;
// edge_vector_t ev;

boost::graph_traits<graph_t>::edge_iterator ei, ei_end;

for (boost::tie(ei,ei_end)=edges(g); ei!=ei_end; ++ei){
	
if(degree(source(*ei,g),g)==1){current= source(*ei,g); break;}
	
	
}

bool instop=false;
do{
 // std::cout<<"Source is "<<graph[current].index_<<std::endl;	

g[current].visited=1;
//g[current].index_=counter;
//counter++;
fermi_vert_list.push_back(current);
get_adjacent_vertex_stat(g, current, next, AmiBase::Fermi); 
// std::cout<<"Source is "<< graph[current].index_<<" next is "<<graph[next].index_<<std::endl;
 if(edge(current,next,g).second){ev.push_back(edge(current,next,g).first);}
 if( g[next].visited==1){instop=true; g[next].visited=1;};
 current=next;

} while(!instop);

	
//return fermi_vert_list.size();	
	
}

void AmiGraph::get_labelled_bose_lists(std::vector< std::vector<int>> &bose_list,  graph_t &g){
	

boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
std::vector<int> entry;
for( boost::tie(ei, ei_end)=edges(g); ei!=ei_end; ++ei){
if( g[*ei].g_struct_.stat_==AmiBase::Bose){

// std::cout<<"Bose line from "<<	g[source(*ei,g)].index_ << " to "<<g[target(*ei,g)].index_<<std::endl;
entry.push_back(g[source(*ei,g)].index_);
entry.push_back(g[target(*ei,g)].index_);
bose_list.push_back(entry);
entry.clear();	
}

}
}	
	
	
	

	
	
void AmiGraph::get_labelled_edge_list(std::vector< std::vector<int>> &bose_list,  graph_t &g){
	

boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
std::vector<int> entry;
for( boost::tie(ei, ei_end)=edges(g); ei!=ei_end; ++ei){
entry.push_back(g[source(*ei,g)].index_);
entry.push_back(g[target(*ei,g)].index_);
bose_list.push_back(entry);
entry.clear();	


}	
	
	
	
}




void AmiGraph::symmetrize_internal_bosonic_lines(graph_t &g){
	
edge_vector_t edge_list;	

find_internal_edges_stat(g, edge_list, AmiBase::Bose);

for(int i=0; i< edge_list.size(); i++){

//add_edge(source(pickone,g), new_vert, edge_info(AmiBase::Fermi,g[pickone].fermi_loop_id, g[pickone].spin), g);
add_edge(target(edge_list[i],g), source(edge_list[i],g), edge_info(g[edge_list[i]].g_struct_.stat_,g[edge_list[i]].fermi_loop_id, g[edge_list[i]].spin),g); 	
	
}

	
	
}


void AmiGraph::symmetrize_bosonic_lines(graph_t &g){
	
edge_vector_t edge_list;	

find_bosonic_edges(g, edge_list);

for(int i=0; i< edge_list.size(); i++){

//add_edge(source(pickone,g), new_vert, edge_info(AmiBase::Fermi,g[pickone].fermi_loop_id, g[pickone].spin), g);
add_edge(target(edge_list[i],g), source(edge_list[i],g), edge_info(g[edge_list[i]].g_struct_.stat_,g[edge_list[i]].fermi_loop_id, g[edge_list[i]].spin),g); 	
	
}

	
	
}

void AmiGraph::symmetrize_non_principle_lines(graph_t &g){
	
edge_vector_t edge_list, e_temp;	
edge_vector_t pl_list;
find_fermionic_edges(g, e_temp);
get_principle_line(g,pl_list);

bool add;
for(int i=0; i< e_temp.size(); i++){
	add=true;
	for(int j=0; i< pl_list.size(); j++){
	
	if(e_temp[i]==pl_list[j]){add=false;}
	
	}
	if(add){
		edge_list.push_back(e_temp[i]);
	}
}



for(int i=0; i< edge_list.size(); i++){

//add_edge(source(pickone,g), new_vert, edge_info(AmiBase::Fermi,g[pickone].fermi_loop_id, g[pickone].spin), g);
add_edge(target(edge_list[i],g), source(edge_list[i],g), edge_info(g[edge_list[i]].g_struct_.stat_,g[edge_list[i]].fermi_loop_id, g[edge_list[i]].spin),g); 	
	
}

	
	
}

void AmiGraph::symmetrize_fermionic_lines(graph_t &g){
	
edge_vector_t edge_list;	

find_fermionic_edges(g, edge_list);

for(int i=0; i< edge_list.size(); i++){

//add_edge(source(pickone,g), new_vert, edge_info(AmiBase::Fermi,g[pickone].fermi_loop_id, g[pickone].spin), g);
add_edge(target(edge_list[i],g), source(edge_list[i],g), edge_info(g[edge_list[i]].g_struct_.stat_,g[edge_list[i]].fermi_loop_id, g[edge_list[i]].spin),g); 	
	
}

	
	
}

bool AmiGraph::is_AID(graph_t g1, graph_t g2){

 symmetrize_bosonic_lines(g1);
 symmetrize_bosonic_lines(g2);

  symmetrize_fermionic_lines(g1);
  symmetrize_fermionic_lines(g2);
// symmetrize_non_principle_lines(g1);
// symmetrize_non_principle_lines(g2);  
	
	
bool final=full_iso(g1,g2);
if(!final){return false;}



return true;
	
}

int AmiGraph::connected_components(graph_t &g){
int out	;
number_vertices(g);
	
cc_graph_t ccg(num_vertices(g));

std::vector<boost::graph_traits<cc_graph_t>::vertex_descriptor> v1(num_vertices(g));

boost::graph_traits<cc_graph_t>::vertex_iterator i, end;
  int id = 0;
  for (boost::tie(i, end) = vertices(ccg); i != end; ++i, ++id) {
        v1[id] = *i;
  }
  
boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
  
  for(boost::tie(ei,ei_end)=edges(g); ei!=ei_end; ++ei){
	
	int sour=g[source(*ei,g)].index_;
	int targ=g[target(*ei,g)].index_;
	
	add_edge(v1[sour], v1[targ] , ccg);	  
	  
  }
  
// conceptually now have a copy of the graph. in cc_graph_t 

std::vector<int> component(boost::num_vertices (ccg));
size_t num_components = boost::connected_components (ccg, &component[0]);  
  
out=(int)num_components;
	
	
return out;	
}

int AmiGraph::fermi_connected_components(graph_t &g){

int out	;
graph_t gc=g;

vertex_vector_t extern_vect_list;
edge_vector_t extern_edge_list;
find_external_vertices(gc, extern_vect_list, extern_edge_list);
// std::cout<<"Deleting legs"<<std::endl;
// print_all_edge_info(gc);
delete_legs(gc, extern_vect_list, extern_edge_list);
// std::cout<<"After Deleting legs"<<std::endl;
// print_all_edge_info(gc);

number_vertices(gc);
	
cc_graph_t ccg(num_vertices(gc));

std::vector<boost::graph_traits<cc_graph_t>::vertex_descriptor> v1(num_vertices(gc));

boost::graph_traits<cc_graph_t>::vertex_iterator i, end;
  int id = 0;
  for (boost::tie(i, end) = vertices(ccg); i != end; ++i, ++id) {
        v1[id] = *i;
  }
  
boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
  
  for(boost::tie(ei,ei_end)=edges(gc); ei!=ei_end; ++ei){
	
  if(gc[*ei].g_struct_.stat_==AmiBase::Fermi){  
  
	int sour=gc[source(*ei,gc)].index_;
	int targ=gc[target(*ei,gc)].index_;
	
	add_edge(v1[sour], v1[targ] , ccg);	  
  }
	  
  }
  
// conceptually now have a copy of the graph. in cc_graph_t 

std::vector<int> component(boost::num_vertices (ccg));
size_t num_components = boost::connected_components (ccg, &component[0]);  
  
out=(int)num_components;
	
	
return out;	
}



// this function removes all the bosonic edges before it checks isomorphism - since those don't have spin labels they just get in the way.
bool AmiGraph::spin_iso(graph_t g1, graph_t g2){
	
// std::cout<<"In spin iso printing graphs"<<std::endl;

// print_all_edge_info(g1);
// std::cout<<"-----"<<std::endl;
// print_all_edge_info(g2);	
	

spin_iso_graph_t 	gc1(num_vertices(g1)), gc2(num_vertices(g2));	
int n= num_vertices(g1);
 std::vector<boost::graph_traits<iso_graph_t>::vertex_descriptor> v1(n), v2(n);

// boost::property_map < spin_iso_graph_t, boost::edge_weight_t >::type EdgeWeightMap1 = get(boost::edge_weight, gc1),
// EdgeWeightMap2 = get(boost::edge_weight, gc2);

// typedef boost::property_map<spin_iso_graph_t, boost::edge_weight_t>::type edge_weight_map_t;

// typedef boost::property_map_equivalent<edge_weight_map_t, edge_weight_map_t> edge_comp_t;

  // edge_comp_t edge_comp = boost::make_property_map_equivalent(
    // boost::get(boost::edge_weight, gc1), boost::get(boost::edge_weight, gc2));

 boost::property_map<spin_iso_graph_t, boost::vertex_index_t>::type 
    v1_index_map = boost::get(boost::vertex_index, gc1),
    v2_index_map = boost::get(boost::vertex_index, gc2);

  boost::graph_traits<spin_iso_graph_t>::vertex_iterator i, end;
  int id = 0;
  for (boost::tie(i, end) = vertices(gc1); i != end; ++i, ++id) {
    put(v1_index_map, *i, id);
    v1[id] = *i;
  }
  id = 0;
  for (boost::tie(i, end) = vertices(gc2); i != end; ++i, ++id) {
    put(v2_index_map, *i, id);
    v2[id] = *i;
  }


  boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
  
  for(boost::tie(ei,ei_end)=edges(g1); ei!=ei_end; ++ei){
	
	// if( g1[*ei].g_struct_.stat_==AmiBase::Fermi){
	int sour=g1[source(*ei,g1)].index_;
	int targ=g1[target(*ei,g1)].index_;
	
	if( g1[*ei].g_struct_.stat_==AmiBase::Fermi){
	add_edge(v1[sour], v1[targ] , EdgeWeightProperty(g1[*ei].spin), gc1);	  
	}
	if( g1[*ei].g_struct_.stat_==AmiBase::Bose){
	add_edge(v1[sour], v1[targ] , EdgeWeightProperty(2), gc1);		
	}
	// }
	
	
  }
  
  for(boost::tie(ei,ei_end)=edges(g2); ei!=ei_end; ++ei){
	  // if( g2[*ei].g_struct_.stat_==AmiBase::Fermi){
	 int sour=g2[source(*ei,g2)].index_;
	int targ=g2[target(*ei,g2)].index_;
	
	if( g2[*ei].g_struct_.stat_==AmiBase::Fermi){
	add_edge(v2[sour], v2[targ] ,EdgeWeightProperty(g2[*ei].spin), gc2);	  
	}
	if( g2[*ei].g_struct_.stat_==AmiBase::Bose){
		add_edge(v2[sour], v2[targ] ,EdgeWeightProperty(2), gc2);	  
	}
	  // } 
  }
  
  
  std::vector<boost::graph_traits<iso_graph_t>::vertex_descriptor> f(n);
  
#if defined(BOOST_MSVC) && BOOST_MSVC <= 1300
  bool ret = isomorphism
    (gc1, gc2, boost::make_iterator_property_map(f.begin(), v1_index_map, f[0]),
     boost::degree_vertex_invariant(), boost::get(boost::vertex_index, gc1), boost::get(boost::vertex_index, gc2));
#else
  bool ret = isomorphism
    (gc1, gc2, isomorphism_map
     (boost::make_iterator_property_map(f.begin(), v1_index_map, f[0])));
#endif	

if(!ret){
	
	return ret;
}

// std::cout << "isomorphic without spins " << ret << std::endl;

// up to here is equivalent to the non-spin variant. Now lets just check that every source and target from one graph is in the other graph



// this part does not work because if there are two edges with same source and target only one is considered





for (int i=0; i< f.size(); i++){
	for(int j=0; j< f.size(); j++){
	
	edge_vector_t ev1, ev2;
	//get(v1_index_map,v1[i]) -> these just return 'i' in the case of v[i]
	ev1=get_edges(i,j, g1);
	
	int one, two;
	
	one=get(v2_index_map, f[i]);
	two=get(v2_index_map, f[j]);
	
	ev2=get_edges(one,two, g2);
	
	// now ev1 and ev2 contain all edges from each graph - check if equivalent - that is that one of the edges in each set matches the other
	
	bool ee=edges_equivalent(ev1, ev2, g1,g2);
	
	if(!ee){ 
	// std::cout<<"Found spin isomorphism false "<<std::endl;
	return false;}
	
	}
	
}


/* 
for (int i=0; i< f.size(); i++){
	for(int j=0; j< f.size(); j++){
	
	if(edge( v1[i],v1[j], gc1).second ){
		
		size_t one, two;
		
		one=get(v2_index_map, f[i]);
		two=get(v2_index_map, f[j]);
	
		// std::cout<<"Edge weight is "<<get(get(boost::edge_weight, gc1) ,edge( v1[i],v1[j], gc1).first) <<std::endl;
		// std::cout<<"Does this exist in other graph: "<< edge(v2[one], v2[two],gc2).second<<std::endl;
		// std::cout<<"Compared to weight "<<get(get(boost::edge_weight, gc1) ,edge( v2[one],v2[two], gc2).first);
		
		int spin_one, spin_two;
		spin_one=get(get(boost::edge_weight, gc1) ,edge( v1[i],v1[j], gc1).first) ;
		spin_two=get(get(boost::edge_weight, gc1) ,edge( v2[one],v2[two], gc2).first);
		
		if(spin_one != spin_two){ 
		std::cout<<"Eliminated graph due to spin mismatch"<<std::endl;
		return false;}
		
		// if( get(get(boost::edge_weight, gc1) ,edge( v1[i],v1[j], gc1).first)
		
		
	}
	
	}
} */

// std::cout<<"Found spin isomorphism true "<< ret<<std::endl;

return ret;	
	
	
}

bool AmiGraph::edges_equivalent( edge_vector_t ev1, edge_vector_t ev2, graph_t &g1, graph_t &g2){

// std::cout<<"Comparing edge lists of sizes "<< ev1.size()<<" "<< ev2.size()<<std::endl;

if(ev1.size()!= ev2.size()){return false;}

if(ev1.size()==0 && ev2.size()==0){ return true;}

// print_all_edge_info(g1);
// std::cout<<"---"<<std::endl;
// print_all_edge_info(g2);

// check vectors - when found in ev2, remove from vector list and decriment iterator 

for( int i = 0 ; i< ev1.size(); i++){
	bool found=false;
	
	int spin1=g1[ev1[i]].spin;
	int stat1=g1[ev1[i]].g_struct_.stat_;
	
	// std::cout<<"On line from "<< g1[source(ev1[i],g1)].index_ <<" to "<< g1[target(ev1[i],g1)].index_<<std::endl;
	// std::cout<<"Searching for spin1 and stat1 "<< spin1<<" "<<stat1<<std::endl;
	// std::cout<<"ev2 size="<<ev2.size()<<std::endl;

	for(int j=0; j< ev2.size(); j++){

	
	
	int spin2=g2[ev2[j]].spin;
	int stat2=g2[ev2[j]].g_struct_.stat_;
	
	// std::cout<<"Comparing to spin2 stat2 "<< spin2<<" "<< stat2<<std::endl;
	// std::cout<<"From vertex "<< g2[source(ev2[j],g2)].index_<<" to "<<g2[target(ev2[j],g2)].index_<<std::endl;
	
	if( (spin1 == spin2) && stat1==stat2){
		// found a match 
		// std::cout<<"Found match"<<std::endl;
		found=true;
		ev2.erase(ev2.begin()+j);
		break;// break out of the j loop 
	}
	}

	if(!found){return false;} // edge vectors are not equivalent because an entry was not found 

}	
	
	// std::cout<<"----"<<std::endl<<std::endl;
// if it gets here then the edge sets are equivalent 
return true;	
	
	
}

AmiGraph::edge_vector_t AmiGraph::get_edges(int si, int ti, graph_t &g){

edge_vector_t output;

boost::graph_traits<graph_t>::edge_iterator ei, ei_end;

for( boost::tie(ei,ei_end)=edges(g); ei!= ei_end; ++ei){	
	
if( g[source(*ei,g)].index_==si && g[target(*ei,g)].index_==ti){

// std::cout<<"Found edge from "<<si <<" to "<< ti<<std::endl;
output.push_back(*ei);

}	
	
	
}

return output;
	
}


bool AmiGraph::full_iso(graph_t g1, graph_t g2){
	
iso_graph_t gc1(num_vertices(g1)), gc2(num_vertices(g2));	


// assign labels to gc1 and gc2 
int n= num_vertices(g1);
 std::vector<boost::graph_traits<iso_graph_t>::vertex_descriptor> v1(n), v2(n);

  boost::property_map<iso_graph_t, boost::vertex_index_t>::type 
    v1_index_map = boost::get(boost::vertex_index, gc1),
    v2_index_map = boost::get(boost::vertex_index, gc2);

  boost::graph_traits<iso_graph_t>::vertex_iterator i, end;
  int id = 0;
  for (boost::tie(i, end) = vertices(gc1); i != end; ++i, ++id) {
    put(v1_index_map, *i, id);
    v1[id] = *i;
  }
  id = 0;
  for (boost::tie(i, end) = vertices(gc2); i != end; ++i, ++id) {
    put(v2_index_map, *i, id);
    v2[id] = *i;
  }
/////

  boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
  
  for(boost::tie(ei,ei_end)=edges(g1); ei!=ei_end; ++ei){
	
	int sour=g1[source(*ei,g1)].index_;
	int targ=g1[target(*ei,g1)].index_;
	
	add_edge(v1[sour], v1[targ] , gc1);	  
	  
  }
  
  for(boost::tie(ei,ei_end)=edges(g2); ei!=ei_end; ++ei){
	 int sour=g2[source(*ei,g2)].index_;
	int targ=g2[target(*ei,g2)].index_;
	
	add_edge(v2[sour], v2[targ] , gc2);	  
	  
  }
  
  // now should have copy of original
  
  std::vector<boost::graph_traits<iso_graph_t>::vertex_descriptor> f(n);
  
#if defined(BOOST_MSVC) && BOOST_MSVC <= 1300
  bool ret = isomorphism
    (gc1, gc2, boost::make_iterator_property_map(f.begin(), v1_index_map, f[0]),
     boost::degree_vertex_invariant(), boost::get(boost::vertex_index, gc1), boost::get(boost::vertex_index, gc2));
#else
  bool ret = isomorphism
    (gc1, gc2, isomorphism_map
     (boost::make_iterator_property_map(f.begin(), v1_index_map, f[0])));
#endif	

//std::cout << "isomorphic? " << ret << std::endl;

return ret;	
	
	
}


void AmiGraph::remove_childless(graph_t &g, vertex_vector_t &vertex_list){
	std::cout<<"Removing childless"<<std::endl;

	for (int i=0; i< vertex_list.size(); i++){
		
	clear_vertex(vertex_list[i],g);
	//remove_vertex(vertex_list[i],g);
		
		
	}
	std::cout<<"Removing childless done"<<std::endl;
	
	
}

void AmiGraph::childless_parent(tree_t &tree, vertex_vector_t &vertex_list){
	
	for(int i=0; i< tree.size(); i++){
		std::cout<<"Tree geometry "<< tree.size()<<" long with "<<tree[i].size()<<" in row "<<i<<std::endl;
		for(int j=0; j< tree[i].size(); j++){
		std::cout<<tree[i][j].child.size()<<std::endl;
		if(tree[i][j].child.size()==0){ 
		 bool add=true;
		for(int j=0; j< vertex_list.size(); j++){
			
			std::cout<<tree[i][j].parent<<" "<< vertex_list[j]<<std::endl;
			if(tree[i][j].parent == vertex_list[j]){ add=false;}
		}
		if(add){
		std::cout<<"Pushing into length of "<< vertex_list.size()<<std::endl;
		vertex_list.push_back(tree[i][j].parent);
		//vertex_list.push_back(tree[i][j].grandpa);
		}
		}
		}	
}

	}
	
void AmiGraph::label_via_tree(tree_t &tree, graph_t &g){
	
int counter=0;	
for( int i=0; i< tree.size(); i++){
	
for(int j=0; j< tree[i].size(); j++){

g[tree[i][j].parent].index_=counter;
counter++;

}
	
}
	
}

void AmiGraph::compare_trees(tree_t &T1, graph_t &g1, tree_t &T2, graph_t &g2){

for( int i=0; i< T1.size(); i++){
	
for(int j=0; j< T1[i].size(); j++){

std::cout<< "Parents "<<g1[T1[i][j].parent].index_<<" "<<g2[T2[i][j].parent].index_<<std::endl;


}}	
	
	
}

void AmiGraph::pp_barcode(graph_t &g, pp_barcode_t &barcode){

	reset_visited(g);
	number_vertices(g);	
	
	vertex_t start_v;
	//int index=find_start_index(g, start_v);
	
	find_entry_vertex(g,start_v);
	
	vertex_vector_t vv, next_vv;
	vv.push_back(start_v);
	
	// std::cout<<"Start index "<< index<<std::endl;
	
	g[start_v].visited=1;
	// get_pp_line(vv, next_vv, graph, barcode);
	
vertex_vector_t vert_list;	
	
	do{
		
	// for(int i=0; i< vv.size(); i++){
		// g[vv[i]].visited=1;	
	// }
		
	get_pp_line(vv, next_vv, g, barcode);
	
	
	for(int i=0; i< next_vv.size(); i++){
		g[next_vv[i]].visited=1;	
	}
	
	// std::cout<<"Next vertices are "<<std::endl;
	// for(int i=0; i<next_vv.size(); i++){

		// std::cout<<g[next_vv[i]].index_<<" "<<g[next_vv[i]].visited<<std::endl;
	
	// }

	vv=next_vv;
	
		
	vert_list.clear();	
	collect_unvisited_vertices(g, vert_list);

// std::cout<<"Remaining unvisited vertices is "<<vert_list.size()<<std::endl;	

// for(int i=0; i< vert_list.size(); i++){
// std::cout<<"Vertex: "<<g[vert_list[i]].index_<<std::endl;	
	
// }
		
	}while(vert_list.size()>0);
	
// collect_unvisited_vertices(g, vert_list);
	
	
}

bool AmiGraph::compare_pp_barcodes(pp_barcode_t &B1, pp_barcode_t &B2){

if(B1.size()!= B2.size()){return false;}

for(int i=0; i< B1.size(); i++){
// std::cout<<"b1:";
// print_bc_line(B1[i]);

// std::cout<<"b2:";
// print_bc_line(B2[i]);

bool line_equal=compare_bc_line(B1[i],B2[i]);

if(!line_equal){
	// std::cout<<"Lines not equal on i="<<i<<std::endl;
	return false;}

}	
	

return true;	
}

void AmiGraph::print_bc_line(pp_barcode_line_t b){

for(int i=0; i<b.size(); i++){

print_pc(b[i]);



}	
	
}

bool AmiGraph::compare_bc_line(pp_barcode_line_t bl1, pp_barcode_line_t bl2){

if(bl1.size()!= bl2.size()){ return false;}

for(int i=0; i< bl1.size()	; i++){

// std::cout<<"Comparing "<<std::endl;
// print_pc(bl1[i]);

// look in bl2 for a matching parent and child and kick it out 
bool found=false;
for(int j=0; j< bl2.size(); j++){

// std::cout<<"to"<<std::endl;
// print_pc(bl2[j]);
// std::cout<<std::endl;
	
	
bool result=pc_equal(bl1[i], bl2[j]);	

// std::cout<<"resulted in "<<result<<std::endl;

if( result ){

bl2.erase(bl2.begin()+j);
found=true;
break; // found a match so break out of 'j' loop 

}	
	
	
}

if(!found){
	std::cout<<"Parent in bl1 at entry "<<i<<std::endl;
	return false;}
	
	
	
}

if(bl2.size()!=0){ return false;}

return true;// if it got here then it is true 
	
	
}

void AmiGraph::print_pc(pp_pc pc){
	
std::cout<<"P[";
for(int j=0; j< pc.parent.size(); j++){
	std::cout<<"("<<pc.parent[j].first <<","<<pc.parent[j].second <<") ";
}
std::cout<<"]"<<std::endl;

std::cout<<"C[";
for(int j=0; j< pc.child.size(); j++){
	std::cout<<"("<<pc.child[j].first <<","<<pc.child[j].second <<") ";
}
std::cout<<"]"<<std::endl;
	
	
}

bool AmiGraph::pc_equal(pp_pc pc1, pp_pc pc2){

if(pc1.parent.size()!= pc2.parent.size()){return false;}

if(pc1.child.size()!= pc2.child.size()){return false;}

for(int i=0; i< pc1.parent.size(); i++){

bool found=false;
for(int j=0; j< pc2.parent.size(); j++){
	
// std::cout<<"Comparing "<<i<<" to "<<j<<std::endl;	
	
// bool result=	(pc1.parent[i]==pc2.parent[j]);
	
// std::cout<<"comparison returned "<< result<<std::endl;	
	
if(pc1.parent[i]==pc2.parent[j]){

pc2.parent.erase(pc2.parent.begin()+j);
found=true;
break;

}


	
}
if(!found){
// std::cout<<"Did not find item "<<i<<std::endl;	
	return false;}	


}

for( int i=0; i< pc1.child.size(); i++){
bool found=false;
for(int j=0; j< pc2.child.size(); j++){	
	
if(pc1.child[i]== pc2.child[j]){

pc2.child.erase(pc2.child.begin()+j);
found=true;
break;

}	
	
}

if(!found){return false;}

}


if(pc2.parent.size()!=0 || pc2.child.size()!=0){

return false;
	
}

return true;	
	
	
	
}

void AmiGraph::get_pp_line(vertex_vector_t &vv, vertex_vector_t &next_vv, graph_t &g, pp_barcode_t &barcode){
	
	next_vv.clear();
	
	pp_barcode_line_t pcline;
	// pp_barcode_line_t bc_line;
	
	// for each vertex find all edges in graph that point to vertex - make these parents
	// for each vertex find all edges that leave the vertex - make these childs 
	// combine parents and childs for each vertex
	for(int i=0; i< vv.size(); i++){
	pp_pc pcstruct;
	
		
	boost::graph_traits<graph_t>::edge_iterator ei, ei_end;

for(boost::tie(ei,ei_end)=edges(g); ei!=ei_end; ++ei){	
	
// child	
if(source(*ei,g)==vv[i]){

std::pair<int,int> c;

c.first=g[*ei].g_struct_.stat_;
c.second=g[*ei].spin;

// std::cout<<"assigned spin "<< c.second<<std::endl;

pcstruct.child.push_back(c);	

if(g[target(*ei,g)].visited==0){
next_vv.push_back(target(*ei,g));	
}
	
	
}

//parent 
if(target(*ei,g)==vv[i]){
	
std::pair<int,int> p;

p.first=g[*ei].g_struct_.stat_;
p.second=g[*ei].spin;

pcstruct.parent.push_back(p);	
	
	
}

	
		
		
}

pcline.push_back(pcstruct);
pcstruct.child.clear();
pcstruct.parent.clear();

	}
	
barcode.push_back(pcline);	
	
	
}

void AmiGraph::construct_graph_barcode(graph_t &g, barcode_t &barcode, tree_t &tree){
// first set all vertices in both not visited 	
	reset_visited(g);
	number_vertices(g);
    vertex_vector_t vert_list;
	tree.clear();
	barcode.clear();
	//tree_t tree;

collect_unvisited_vertices(g, vert_list);
//
vertex_vector_t current_list;
vertex_vector_t next_list;

boost::graph_traits<graph_t>::edge_iterator ei, ei_end;

for (boost::tie(ei,ei_end)=edges(g); ei!=ei_end; ++ei){
	
if(degree(source(*ei,g),g)==1){current_list.push_back(source(*ei,g)); break;}
	
	
}

 // std::cout<<"Starting at "<<g[current_list[0]].index_<<std::endl;

//current_list.push_back(vert_list[0]);
//if(g[vert_list[0]].index_!=0){std::cout<<"Possible error in starting point"<<std::endl;}

//

// NEED function get_ordered_adjacency lists - only difference is that it needs to be in the F, B order

bool exit_do=false;
do{
	// std::cout<<"Getting line "<<std::endl;
//get_ForB_line(g, current_list, next_list, single_line);

// std::cout<<"Current list is size "<< current_list.size()<<std::endl;
// std::cout<<"Labelling"<<std::endl;

get_tree_leg(g, current_list, next_list, tree);
label_visited(g, current_list);
//FBarray.push_back(single_line);
//single_line.clear();
// std::cout<<"Finished "<<std::endl;

// set current_list to next _list if next_list.size()!=0

if(next_list.size()==0){
exit_do=true;	
}else{
current_list=next_list;
}	
	
}while(!exit_do);


// function to convert tree to barcode 
// std::cout<<"Convert to barcode"<<std::endl;

convert_tree_to_barcode(tree, barcode, g );
// std::cout<<"Convert to barcode finished"<<std::endl;

}


void AmiGraph::reconstruct_graph_barcode(graph_t &g, barcode_t &barcode, tree_t &tree){
// first set all vertices in both not visited 	
	reset_visited(g);
	number_vertices(g);
    vertex_vector_t vert_list;
	tree.clear();
	barcode.clear();
	//tree_t tree;

collect_unvisited_vertices(g, vert_list);
//
vertex_vector_t current_list;
vertex_vector_t next_list;

boost::graph_traits<graph_t>::edge_iterator ei, ei_end;

for (boost::tie(ei,ei_end)=edges(g); ei!=ei_end; ++ei){
	
if(degree(source(*ei,g),g)==1){current_list.push_back(source(*ei,g)); break;}
	
	
}


bool exit_do=false;
do{


get_tree_leg(g, current_list, next_list, tree);
label_visited(g, current_list);

if(next_list.size()==0){
exit_do=true;	
}else{
current_list=next_list;
}	
	
}while(!exit_do);



reconvert_tree_to_barcode(tree, barcode, g );
// std::cout<<"Convert to barcode finished"<<std::endl;

}

void AmiGraph::reconvert_tree_to_barcode(tree_t &tree, barcode_t &barcode, graph_t &g){
	
	std::cout<<"Summarizing tree parent structure"<<std::endl;
	for(int i=1; i< tree.size(); i++){
std::cout<<"On tree level"<<i<<std::endl;
for(int j=0; j<tree[i].size(); j++){
	std::cout<<"Line from  "<< g[tree[i][j].grandpa].index_<<" to "<< g[tree[i][j].parent].index_<< std::endl;
	
	}}

barcode_line_t line;
pc_struct entry;

if(tree.size()>0){	
for(int i=1; i< tree.size(); i++){

for(int j=0; j<tree[i].size(); j++){

// std::cout<<"GP "<< g[tree[i][j].grandpa].index_<<std::endl;
	
entry.parent.push_back(g[edge(tree[i][j].grandpa, tree[i][j].parent, g).first].g_struct_.stat_);
entry.parent.push_back(g[edge(tree[i][j].grandpa, tree[i][j].parent, g).first].spin);	

for(int k=0; k< tree[i][j].child.size(); k++){
child_t kid;

	
	kid.push_back(	g[edge(tree[i][j].parent,tree[i][j].child[k],g).first].g_struct_.stat_);
	kid.push_back(	g[edge(tree[i][j].parent,tree[i][j].child[k],g).first].spin);
	entry.child.push_back(kid);
	kid.clear();	
	
	
}


line.push_back(entry);
entry.parent.clear();
entry.child.clear();	
	
	// parent_list[i].push_back(	g[edge(current_list[i],*ai,g).first].g_struct_.stat_);
	// parent_list[i].push_back(	g[edge(current_list[i],*ai,g).first].spin);	
}

barcode.push_back(line);
line.clear();

}
}	

}

void AmiGraph::get_tree_leg(graph_t &g, vertex_vector_t &current_list, vertex_vector_t &next_list, tree_t &tree){

next_list.clear();	

tree_leg_t this_leg;

pc_vert_struct this_entry;
 // std::cout<<"Tree size is "<<tree.size()<<std::endl;

if(tree.size()==0){
	for (int i=0; i< current_list.size(); i++){
	 // std::cout<<"On element "<<i<<std::endl;
	this_entry.parent=current_list[i];
	
	boost::graph_traits<graph_t>::adjacency_iterator ai, ai_end;
	for(boost::tie(ai,ai_end)=adjacent_vertices(current_list[i],g); ai!=ai_end; ++ai){

	if (g[*ai].visited==0){
	next_list.push_back(*ai);
	
	this_entry.child.push_back(*ai);
	}
	// this_entry.child.push_back(*ai);
	
	}
	
	this_leg.push_back(this_entry);
	this_entry.child.clear();
	
	}
	
	
}

if(tree.size()>0){
	// std::cout<<"Previous tree line has j entries = "<< tree[tree.size()-1].size()<<std::endl;
	
	// get parents list from previous children list 
for(int j=0; j<tree[tree.size()-1].size(); j++){

// std::cout<<"Child size from last leg is "<<	tree[tree.size()-1][j].child.size() <<std::endl;
for ( int k=0; k<tree[tree.size()-1][j].child.size(); k++){
 // std::cout<<"k is "<<k<<std::endl;
 // if (g[tree[tree.size()-1][j].child[k]].visited==0){
this_entry.parent= tree[tree.size()-1][j].child[k];
this_leg.push_back(this_entry);
// }
}
}


// this leg now has a list of parents so for each parent in this leg, find it's adjacent 

for(int j=0; j< this_leg.size(); j++){
	
	boost::graph_traits<graph_t>::adjacency_iterator ai, ai_end;
	for(boost::tie(ai,ai_end)=adjacent_vertices(this_leg[j].parent,g); ai!=ai_end; ++ai){

	if (g[*ai].visited==0){
	next_list.push_back(*ai);
	
	}
	if(g[edge(this_leg[j].parent, *ai, g).first].g_struct_.stat_==AmiBase::Fermi){
	this_leg[j].child.push_back(*ai);
	}
		
	}	
	
// do it again for bosonic edges 

	for(boost::tie(ai,ai_end)=adjacent_vertices(this_leg[j].parent,g); ai!=ai_end; ++ai){

	if(g[edge(this_leg[j].parent, *ai, g).first].g_struct_.stat_==AmiBase::Bose){
	this_leg[j].child.push_back(*ai);
	}
		
	}	
	
	
	
	// find grandparent for this_leg[j].parent 
	
	for(int l=0; l<tree[tree.size()-1].size(); l++){

// std::cout<<"Child size is "<<	tree[tree.size()-1][j].child.size() <<std::endl;
for ( int k=0; k<tree[tree.size()-1][l].child.size(); k++){

if(this_leg[j].parent== tree[tree.size()-1][l].child[k]){
	// std::cout<<"Found grandparent"<<std::endl;
this_leg[j].grandpa=tree[tree.size()-1][l].parent	;
break;
}

}	
	
	
}
	
	
	
}








}
	
	
// std::cout<<"This leg has length "<< this_leg.size()<<std::endl;	
if(this_leg.size()>0){
tree.push_back(this_leg);
}	
	
	
}
/* 
void AmiGraph::get_tree_leg(graph_t &g, vertex_vector_t &current_list, vertex_vector_t &next_list, tree_t &tree){
	
next_list.clear();	

tree_leg_t this_leg;

pc_vert_struct this_entry;

for (int i=0; i< current_list.size(); i++){
	std::cout<<"On element "<<i<<std::endl;
	this_entry.parent=current_list[i];
	
	boost::graph_traits<graph_t>::adjacency_iterator ai, ai_end;
	for(boost::tie(ai,ai_end)=adjacent_vertices(current_list[i],g); ai!=ai_end; ++ai){

	if (g[*ai].visited==0){
	next_list.push_back(*ai);
	
	
	}
	this_entry.child.push_back(*ai);
	
	}

// find this parent's grandparent 

// look through the children of tree[i-1] to get vertex=current_list[i], then grandparent = that vertex.parent
std::cout<<"Tree size is now "<<tree.size()<< " and i is "<< i<<std::endl;

if(tree.size()>0){
	// std::cout<<"Previous tree line has j entries = "<< tree[tree.size()-1].size()<<std::endl;
for(int j=0; j<tree[tree.size()-1].size(); j++){

std::cout<<"Child size is "<<	tree[tree.size()-1][j].child.size() <<std::endl;
for ( int k=0; k<tree[tree.size()-1][j].child.size(); k++){

if(this_entry.parent== tree[tree.size()-1][j].child[k]){
	std::cout<<"Found grandparent"<<std::endl;
this_entry.grandpa=tree[tree.size()-1][j].parent	;
break;
}

}	
	
	
}

	
}

this_leg.push_back(this_entry);
	//this_entry.parent.clear();
	this_entry.child.clear();

}

if(this_leg.size()>0){
tree.push_back(this_leg);
}
} */

void AmiGraph::convert_tree_to_barcode(tree_t &tree, barcode_t &barcode, graph_t &g){

barcode_line_t line;
pc_struct entry;

if(tree.size()>0){	
for(int i=1; i< tree.size(); i++){

for(int j=0; j<tree[i].size(); j++){

	
entry.parent.push_back(g[edge(tree[i][j].grandpa, tree[i][j].parent, g).first].g_struct_.stat_);
entry.parent.push_back(g[edge(tree[i][j].grandpa, tree[i][j].parent, g).first].spin);	

for(int k=0; k< tree[i][j].child.size(); k++){
child_t kid;

	
	kid.push_back(	g[edge(tree[i][j].parent,tree[i][j].child[k],g).first].g_struct_.stat_);
	kid.push_back(	g[edge(tree[i][j].parent,tree[i][j].child[k],g).first].spin);
	entry.child.push_back(kid);
	kid.clear();	
	
	
}


line.push_back(entry);
entry.parent.clear();
entry.child.clear();	
	
	// parent_list[i].push_back(	g[edge(current_list[i],*ai,g).first].g_struct_.stat_);
	// parent_list[i].push_back(	g[edge(current_list[i],*ai,g).first].spin);	
}

barcode.push_back(line);
line.clear();

}
}	

}
	
	


bool AmiGraph::barcodes_are_equal(barcode_t &B1, barcode_t &B2){
	
if(B1.size()!= B2.size()){
	//std::cout<<"Barcode sizes are different"<<std::endl;
	return false;}

for (int i=0; i< B1.size(); i++){

// std::cout<<"Comparing lines in barcode "<<i<<" out of "<< B1.size()<<std::endl;

if (!barcode_lines_equal(B1[i],B2[i])){ return false;}



}	
	
return true;	
}

void AmiGraph::print_barcodes(barcode_t B1, barcode_t B2){
	
for(int c=0; c<B1.size(); c++){

std::cout<<"B1_line size is "<<B1[c].size()<<std::endl;
for(int i=0; i< B1[c].size(); i++){
	std::cout<<"B1 parent "<<B1[c][i].parent[0]<<" "<<B1[c][i].parent[1]<< std::endl;
	for(int j=0; j<B1[c][i].child.size(); j++){
		std::cout<<"B1 child "<<B1[c][i].child[j][0]	<<" "<<B1[c][i].child[j][1]<<" | ";
	}	
	std::cout<<std::endl;
}
}	


for(int c=0; c<B2.size(); c++){

std::cout<<"B2_line size is "<<B2[c].size()<<std::endl;
for(int i=0; i< B2[c].size(); i++){
	std::cout<<"B2 parent "<<B2[c][i].parent[0]<<" "<<B2[c][i].parent[1]<< std::endl;
	for(int j=0; j<B2[c][i].child.size(); j++){
		std::cout<<"B2 child "<<B2[c][i].child[j][0]	<<" "<<B2[c][i].child[j][1]<<" | ";
	}	
	std::cout<<std::endl;
}
}	
	
	
	
	
}

// in case destructive, B1_line are not pointers
bool AmiGraph::barcode_lines_equal(barcode_line_t B1_line, barcode_line_t B2_line){
	
if(B1_line.size() != B2_line.size()){return false;}

for  (int i=0; i< B1_line.size(); i++){

// for each element in B1_line, find a matching element in B2_line and remove it from B2.  If at the end B2_line.size()==0 then entire line is the same.
// std::cout<<"B1_line size is "<<B1_line.size()<<std::endl;
// std::cout<<"B2_line size is "<<B2_line.size()<<std::endl;
for( int j=0; j< B2_line.size(); j++){

// std::cout<<"Comparing element of B1 to B2 "<<i<<" "<<j << " of totals "<< B1_line.size()<<" "<< B2_line.size()<<std::endl;
// std::cout<<"Parent size is "<< B1_line[i].parent.size()<<" and "<< B2_line[j].parent.size()<<std::endl;

// std::cout<<"B1 parent "<<B1_line[i].parent[0]<<" "<<B1_line[i].parent[1]<< std::endl;
// std::cout<<"B2 parent "<<B2_line[j].parent[0]<<" "<<B2_line[j].parent[1]<< std::endl;

// for(int c=0; c< B1_line[i].child.size(); c++){

// std::cout<<"B1 child "<<B1_line[i].child[c][0]	<<" "<<B1_line[i].child[c][1]<<" | ";
	
// }
// std::cout<<std::endl;

// for(int c=0; c< B2_line[j].child.size(); c++){

// std::cout<<"B2 child "<<B2_line[j].child[c][0]	<<" "<<B2_line[j].child[c][1]<<" | ";
	
// }
// std::cout<<std::endl; 

//if( std::memcmp( &B1_line[i], &B2_line[j], sizeof(&B1_line[i]))==0){

if ( pc_are_equal(B1_line[i], B2_line[j])){

// the jth element of B2 is the i'th of B1.
// std::cout<<"Elements are equal "<<i<<" "<<j<<std::endl;
// B1_line.erase(B1_line.begin()+(i));	
B2_line.erase(B2_line.begin()+(j));	
continue; // break j loop since we reduced the size of B2 we need to do this anyways 
	
}
	
	
}
	

// here is the tricky part.  I have a set of parent-child sets in one vector.  They need to be exactly in the other vector but the order does not matter.
//  std::memcmp( pointer1, pointer2, sizeof(pointer1)) = 0 if they are the same 
// can use to compare pc_structs

	
	
}	
	
if(B2_line.size()==0){ return true;}
else{ return false;}
	
	
}

bool AmiGraph::pc_are_equal( pc_struct pc1, pc_struct pc2){
	

// std::cout<<"Parent size "<< pc1.parent.size()<<" "<< pc2.parent.size()<<std::endl;
// std::cout<<"Childred size "<< pc1.child.size()<<" "<< pc2.child.size()<<std::endl;


if( ! parents_are_equal(pc1.parent, pc2.parent)){ return false;}

if (! children_are_equal(pc1.child, pc2.child)){ return false;}
	
// std::cout<<"Found that PC struct was equal"<<std::endl;	
return true;	
}

bool AmiGraph::parents_are_equal( parent_t p1, parent_t p2){
if( p1.size()!= p2.size()){ return false;}
for(int i=0; i< p1.size(); i++){

if (p1[i] != p2[i]){return false;}

} 	
	
return true;	
}

bool AmiGraph::children_are_equal( std::vector<child_t> c1, std::vector<child_t> c2){

if (c1.size()!= c2.size()){ return false;}

if(c1.size()==0){ return true;}

for(int i=0; i< c1.size(); i++){

for(int j=0; j< c2.size(); j++){
	
if( c1[i]==c2[j]){ // if child vectors are the same, kick out of list then continue

c2.erase(c2.begin()+j);
continue;

}	
	
	
	
}



}	
	
if( c2.size()==0){return true;}
else{ return false;}	
	
}


void AmiGraph::get_barcode_line(graph_t &g, vertex_vector_t &current_list, vertex_vector_t &next_list, barcode_t &barcode){
	
next_list.clear();	


std::vector< pc_struct> pc_vector;
std::vector< parent_t > parent_list;

for (int i=0; i< current_list.size(); i++){
	pc_struct pc;
	
	vertex_vector_t child_list;
	child_list.clear();
	
	bool parent_needed=true;
	boost::graph_traits<graph_t>::adjacency_iterator ai, ai_end;
	for(boost::tie(ai,ai_end)=adjacent_vertices(current_list[i],g); ai!=ai_end; ++ai){

	if (g[*ai].visited==0){
	next_list.push_back(*ai);
	child_list.push_back(*ai);
	
	
	parent_list[i].push_back(	g[edge(current_list[i],*ai,g).first].g_struct_.stat_);
	parent_list[i].push_back(	g[edge(current_list[i],*ai,g).first].spin);
	
	// std::cout<<"Parent is edge "<< g[source(edge(current_list[i],*ai,g).first, g)].index_<<" to "<< g[target(edge(current_list[i],*ai,g).first, g)].index_<<std::endl;
	}
	
	}
	
	// std::cout<<"Next list has size "<<next_list.size()<<std::endl;
	// std::cout<<"Child list has size "<<child_list.size()<<std::endl;
		
	pc.parent=parent_list[i];
	
	for (int j=0; j< child_list.size(); j++){
	// std::cout<<"i and j are "<<i<<" "<<j<<std::endl;	
	boost::graph_traits<graph_t>::adjacency_iterator ai2, ai_end2;
	for(boost::tie(ai2,ai_end2)=adjacent_vertices(next_list[j],g); ai2!=ai_end2; ++ai2){
		// std::cout<<"In tie loop"<<std::endl;
    child_t kid;
	if (g[*ai2].visited==0){
	kid.push_back(	g[edge(next_list[j],*ai2,g).first].g_struct_.stat_);
	kid.push_back(	g[edge(next_list[j],*ai2,g).first].spin);
	pc.child.push_back(kid);
	
	kid.clear();
	// std::cout<<"Added child from "<<g[next_list[j]].index_<<" to "<< g[*ai2].index_<<std::endl;
	}
	
		
		
		
	}
	}
	
if(pc.parent.size()!=0){	
pc_vector.push_back(pc);
}
pc.parent.clear();
pc.child.clear();

	
}

barcode.push_back(pc_vector);	
	
}

void AmiGraph::get_ForB_line(graph_t &g, vertex_vector_t &current_list, vertex_vector_t &next_list, std::vector< std::vector<int>> &single_line ){
	
next_list.clear();	

std::vector<int> fermi_info, bose_info, combined;
	
for (int i=0 ; i< current_list.size();i++){
fermi_info.clear();
bose_info.clear();

boost::graph_traits<graph_t>::adjacency_iterator ai, ai_end;

	for(boost::tie(ai,ai_end)=adjacent_vertices(current_list[i],g); ai!=ai_end; ++ai){

	if (g[*ai].visited==0){
	next_list.push_back(*ai);

	if(g[edge(current_list[i],*ai,g).first].g_struct_.stat_==AmiBase::Fermi)
	{
	combined.push_back(	g[edge(current_list[i],*ai,g).first].g_struct_.stat_);
	combined.push_back(	g[edge(current_list[i],*ai,g).first].spin);
	single_line.push_back(combined);
	combined.clear();
		
	}

		
	}	
		
		
	}

	for(boost::tie(ai,ai_end)=adjacent_vertices(current_list[i],g); ai!=ai_end; ++ai){

	if (g[*ai].visited==0){

	if(g[edge(current_list[i],*ai,g).first].g_struct_.stat_==AmiBase::Bose)
	{
		combined.push_back(g[edge(current_list[i],*ai,g).first].g_struct_.stat_);
		single_line.push_back(combined);
	    combined.clear();
	}


		
	}	
		
		
	}




}	
// if(combined.size()!=0){
// FBarray.push_back(combined);
// combined.clear();
// }
}
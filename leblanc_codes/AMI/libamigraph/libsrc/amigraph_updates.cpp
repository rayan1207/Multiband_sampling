
#include "amigraph.hpp"


void AmiGraph::AIL_defined(graph_t &g, edge_t &one, edge_t &two){

// We need two cases.  One where the edges are the same.  and if not, they must be different.
if (one == two){



vertex_t new_left = add_vertex(vertex_info(three_leg), g);
vertex_t new_right = add_vertex(vertex_info(three_leg), g);

//add_edge(source(one,g), new_onevert, edge_info(AmiBase::Fermi,g[one].fermi_loop_id), g);

add_edge(source(one,g),new_left,edge_info(AmiBase::Fermi,g[one].fermi_loop_id, g[one].spin),g);
add_edge(new_left,new_right,edge_info(AmiBase::Fermi,g[one].fermi_loop_id, g[one].spin),g);
add_edge(new_right,target(one,g),edge_info(AmiBase::Fermi,g[one].fermi_loop_id, g[one].spin),g);

add_edge(new_left,new_right,edge_info(AmiBase::Bose),g);

remove_edge(one,g);
}
else{




vertex_t new_onevert = add_vertex(vertex_info(three_leg), g);
vertex_t new_twovert = add_vertex(vertex_info(three_leg), g);

// TODO: This will make lots of cases where we get non-hubbard diagrams
// Add edges
add_edge(source(one,g), new_onevert, edge_info(AmiBase::Fermi,g[one].fermi_loop_id, g[one].spin), g);
add_edge(new_onevert, target(one,g), edge_info(AmiBase::Fermi,g[one].fermi_loop_id, g[one].spin), g);


add_edge(source(two,g), new_twovert, edge_info(AmiBase::Fermi,g[two].fermi_loop_id, g[two].spin), g);
add_edge(new_twovert, target(two,g), edge_info(AmiBase::Fermi,g[two].fermi_loop_id, g[two].spin), g);


// Add bosonic edge
add_edge(new_onevert,new_twovert, edge_info(AmiBase::Bose),g);



remove_edge(one,g);
remove_edge(two,g);
}// end else line



}


// TODO: It might be better to separate these two updates. 
void AmiGraph::AIL_rand(graph_t &g){

edge_t one, two;

std::cout<<"Selecting two edges"<<std::endl;

// need function - pick two fermionic lines

find_two_edges(g, AmiBase::Fermi , one,  two);

//std::cout<<"two edges selected"<<std::endl;

//std::cout<< one<< " "<< two << std::endl;


// We need two cases.  One where the edges are the same.  and if not, they must be different.
if (one == two){


vertex_t new_left = add_vertex(vertex_info(three_leg), g);
vertex_t new_right = add_vertex(vertex_info(three_leg), g);

//add_edge(source(one,g), new_onevert, edge_info(AmiBase::Fermi,g[one].fermi_loop_id), g);

add_edge(source(one,g),new_left,edge_info(AmiBase::Fermi,g[one].fermi_loop_id, g[one].spin),g);
add_edge(new_left,new_right,edge_info(AmiBase::Fermi,g[one].fermi_loop_id, g[one].spin),g);
add_edge(new_right,target(one,g),edge_info(AmiBase::Fermi,g[one].fermi_loop_id, g[one].spin),g);

add_edge(new_left,new_right,edge_info(AmiBase::Bose),g);

remove_edge(one,g);
}
else{


vertex_t new_onevert = add_vertex(vertex_info(three_leg), g);
vertex_t new_twovert = add_vertex(vertex_info(three_leg), g);

// TODO: This will make lots of cases where we get non-hubbard diagrams
// Add edges
add_edge(source(one,g), new_onevert, edge_info(AmiBase::Fermi,g[one].fermi_loop_id, g[one].spin), g);
add_edge(new_onevert, target(one,g), edge_info(AmiBase::Fermi,g[one].fermi_loop_id, g[one].spin), g);


add_edge(source(two,g), new_twovert, edge_info(AmiBase::Fermi,g[two].fermi_loop_id, g[two].spin), g);
add_edge(new_twovert, target(two,g), edge_info(AmiBase::Fermi,g[two].fermi_loop_id, g[two].spin), g);


// Add bosonic edge
add_edge(new_onevert,new_twovert, edge_info(AmiBase::Bose),g);


remove_edge(one,g);
remove_edge(two,g);

}// end else line



}


/////////////////////////////////////////////////



void AmiGraph::DIL_rand(graph_t &g){

// get list of bosonic edges


AmiGraph::edge_vector_t edge_list;
AmiGraph::find_bosonic_edges(g, edge_list);

// for(int i=0;i< edge_list.size(); i++){

// std::cout<<edge_list[i] <<std::endl;
// }


//std::cout<<"List of bosonic edges is "<<std::endl;

//std::uniform_int_distribution<int> int_dist(0,edge_list.size()-1);

int randint= random_int(0, edge_list.size()-1);
edge_t pickone=edge_list[randint];

//std::cout <<"Picked edge "<< pickone << " to remove" << std::endl;

remove_edge(pickone,g);

// that logic can occur later
// need logic to do this. properly -> to know if this is a sigma or pi diagram. not reduce sigma below Nb=2, or remove the external legs


}

void AmiGraph::AT_defined(graph_t &g, edge_t &pickone, AmiGraph::spin_type spin){
	
		
g[boost::graph_bundle].fermi_loop_id++;	

// insert vertex - attach source and target appropriately - remove edge

vertex_t new_vert=add_vertex(g);

add_edge(source(pickone,g), new_vert, edge_info(AmiBase::Fermi,g[pickone].fermi_loop_id, g[pickone].spin), g);
add_edge(new_vert, target(pickone,g), edge_info(AmiBase::Fermi,g[pickone].fermi_loop_id, g[pickone].spin), g);
remove_edge(pickone,g);

// insert bubble vertex - attach bosonic line to iter_swap
// TODO: question.  Does the direction of hartree bosonic lines matter? later when we label in and out bosonic lines we will reun into problems.  We don't correctly find bubbles if bosonic lines added later have the wrong direction.
vertex_t top=add_vertex(g);
add_edge(new_vert, top, edge_info(AmiBase::Bose), g);


// insert fermionic loop with source and target the same	

	
add_edge(top, top, edge_info(AmiBase::Fermi,g[boost::graph_bundle].fermi_loop_id, spin), g);	



}


void AmiGraph::AT_rand(graph_t &g){
	
	
g[boost::graph_bundle].fermi_loop_id++;	
// pick random fermionic edge - exclude 'external edges'
edge_t pickone;

find_rand_internal_edge(g, AmiBase::Fermi, pickone);


// insert vertex - attach source and target appropriately - remove edge

vertex_t new_vert=add_vertex(g);

add_edge(source(pickone,g), new_vert, edge_info(AmiBase::Fermi,g[pickone].fermi_loop_id, g[pickone].spin), g);
add_edge(new_vert, target(pickone,g), edge_info(AmiBase::Fermi,g[pickone].fermi_loop_id, g[pickone].spin), g);
remove_edge(pickone,g);

// insert bubble vertex - attach bosonic line to iter_swap
// TODO: question.  Does the direction of hartree bosonic lines matter? later when we label in and out bosonic lines we will reun into problems.  We don't correctly find bubbles if bosonic lines added later have the wrong direction.
vertex_t top=add_vertex(g);
add_edge(new_vert, top, edge_info(AmiBase::Bose), g);


// insert fermionic loop with source and target the same	

if(g[pickone].spin==up){	
add_edge(top, top, edge_info(AmiBase::Fermi,g[boost::graph_bundle].fermi_loop_id, dn), g);	
}else{
add_edge(top, top, edge_info(AmiBase::Fermi,g[boost::graph_bundle].fermi_loop_id, up), g);	
}
	
	
}

void AmiGraph::DT_rand(graph_t &g){
	
	vertex_vector_t tp_vert_vec;
	edge_vector_t tp_bose_line_vec;
	vertex_vector_t tp_connect_vert_vec;
	find_tadpoles(g, tp_vert_vec, tp_connect_vert_vec, tp_bose_line_vec);
	
//	std::cout<<"Number of tadpoles found is "<< tp_vert_vec.size()<<std::endl;
	
	edge_t V1_in, V1_out;
	
	if(tp_vert_vec.size()!=0){
		
		int randint= random_int(0, tp_vert_vec.size()-1);
		find_in_out_edges_oftype(g, tp_connect_vert_vec[randint], V1_in,V1_out, AmiBase::Fermi); 
		add_edge(source(V1_in,g),target(V1_out,g), edge_info(AmiBase::Fermi,g[V1_in].fermi_loop_id, g[V1_in].spin), g);	
		
		
		
		clear_vertex(tp_connect_vert_vec[randint],g);
		remove_vertex(tp_connect_vert_vec[randint],g);
		clear_vertex(tp_vert_vec[randint],g);
		remove_vertex(tp_vert_vec[randint],g);
			
		
	}
	
	
}

void AmiGraph::RBI_rand(graph_t &g){
	
AmiGraph::vertex_vector_list_t bubble_vertex_list;
AmiGraph::edge_vector_list_t bubble_edge_list;
AmiGraph::edge_vector_list_t legs_edge_list;
int bubble_choice;

vertex_t left, right;


bubble_finder(g, bubble_vertex_list,  bubble_edge_list, legs_edge_list);

//std::cout<<"Found nb = "<< bubble_vertex_list.size()<<std::endl;

if( bubble_vertex_list.size() != 0){
	
	bubble_choice=random_int(0, bubble_vertex_list.size()-1);
	
	if(target(legs_edge_list[bubble_choice][0], g)== bubble_vertex_list[bubble_choice][0] || target(legs_edge_list[bubble_choice][0], g)== bubble_vertex_list[bubble_choice][1]){ left=source(legs_edge_list[bubble_choice][0], g);
	}else{
left=target(legs_edge_list[bubble_choice][0], g);
	}

	if(target(legs_edge_list[bubble_choice][1], g)== bubble_vertex_list[bubble_choice][0] || target(legs_edge_list[bubble_choice][1], g)== bubble_vertex_list[bubble_choice][1]){ right=source(legs_edge_list[bubble_choice][1], g);
	}else{
right=target(legs_edge_list[bubble_choice][1], g);
	}	
	
// now want to deal with left and right 	
	// add new interaction line
	add_edge(left,right,edge_info(AmiBase::Bose),g); 
	
	// clear the bubble vertices and remove them.
	clear_vertex(bubble_vertex_list[bubble_choice][0],g);
	clear_vertex(bubble_vertex_list[bubble_choice][1],g);
	remove_vertex(bubble_vertex_list[bubble_choice][0],g);
	remove_vertex(bubble_vertex_list[bubble_choice][1],g);
	
	
}
	
}

/* 
void AmiGraph::RB_RPA_defined(graph_t &g, vertex_vector_t &bubble_vertices, edge_vector_t &bubble_edges,edge_vector_t legs_edges){
	
// AmiGraph::vertex_vector_list_t bubble_vertex_list;
// AmiGraph::edge_vector_list_t bubble_edge_list;
// AmiGraph::edge_vector_list_t legs_edge_list;
// int bubble_choice;

vertex_t left, right;


	if(target(legs_edges[0], g)== bubble_vertices[0] || target(legs_edges[0], g)== bubble_vertices[1]){ left=target(legs_edges[0], g);
	}else{
left=source(legs_edges[0], g);
	}

	if(source(legs_edges[1], g)== bubble_vertices[0] || source(legs_edges[1], g)== bubble_vertices[1]){ right=source(legs_edges[1], g);
	}else{
right=target(legs_edges[1], g);
	}	
	
	// create two new vertices 
	
// now want to deal with left and right 	
	// add new interaction line
	add_edge(left,right,edge_info(AmiBase::Bose),g); 
	
	// clear the bubble vertices and remove them.
	clear_vertex(bubble_vertex_list[bubble_choice][0],g);
	clear_vertex(bubble_vertex_list[bubble_choice][1],g);
	remove_vertex(bubble_vertex_list[bubble_choice][0],g);
	remove_vertex(bubble_vertex_list[bubble_choice][1],g);
	
	

	
} */


void AmiGraph::RIB_defined(graph_t &g, edge_t &one, spin_type spin){

	g[boost::graph_bundle].fermi_loop_id++;
	// edge_vector_t bose_edges;
	
	// find_bosonic_edges(g,bose_edges);
//	std::cout<<"Found m interaction lines "<< bose_edges.size()<<std::endl;
	
	// int int_choice;
	// int_choice=random_int(0,bose_edges.size()-1);
	
	// have an interaction line now. so do something with it.
	
	// creat left and right bubble vertices
	vertex_t left=add_vertex(vertex_info(three_leg), g);
	vertex_t right=add_vertex(vertex_info(three_leg), g);
	// attach source of chosen int to left
	add_edge(source(one,g), left, edge_info(AmiBase::Bose), g);
	// attach right to target chosen int
	add_edge(right, target(one,g), edge_info(AmiBase::Bose), g);

	// attach left to right and right to left as fermionic
	// TODO: Same problem as adding interaction lines or bubbles - can't prevent moving to bad configuration
	// int rnum=random_int(0,1);
	add_edge(left, right, edge_info(AmiBase::Fermi, g[boost::graph_bundle].fermi_loop_id, spin), g);
	add_edge(right, left, edge_info(AmiBase::Fermi, g[boost::graph_bundle].fermi_loop_id, spin), g);
	// delete int edge
	remove_edge(one,g);
	
		
}

void AmiGraph::RIB_rand(graph_t &g){
	
	g[boost::graph_bundle].fermi_loop_id++;
	edge_vector_t bose_edges;
	
	find_bosonic_edges(g,bose_edges);
//	std::cout<<"Found m interaction lines "<< bose_edges.size()<<std::endl;
	
	int int_choice;
	int_choice=random_int(0,bose_edges.size()-1);
	
	// have an interaction line now. so do something with it.
	
	// creat left and right bubble vertices
	vertex_t left=add_vertex(vertex_info(three_leg), g);
	vertex_t right=add_vertex(vertex_info(three_leg), g);
	// attach source of chosen int to left
	add_edge(source(bose_edges[int_choice],g), left, edge_info(AmiBase::Bose), g);
	// attach right to target chosen int
	add_edge(right, target(bose_edges[int_choice],g), edge_info(AmiBase::Bose), g);

	// attach left to right and right to left as fermionic
	// TODO: Same problem as adding interaction lines or bubbles - can't prevent moving to bad configuration
	int rnum=random_int(0,1);
	add_edge(left, right, edge_info(AmiBase::Fermi, g[boost::graph_bundle].fermi_loop_id, AmiGraph::spin_type(rnum)), g);
	add_edge(right, left, edge_info(AmiBase::Fermi, g[boost::graph_bundle].fermi_loop_id, AmiGraph::spin_type(rnum)), g);
	// delete int edge
	remove_edge(bose_edges[int_choice],g);
	
	
	
	
}

// This is a completely optional update.
void AmiGraph::SN_rand(graph_t &g){
	
}






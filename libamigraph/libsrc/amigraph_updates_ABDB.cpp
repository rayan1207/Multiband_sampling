
#include "amigraph.hpp"

void AmiGraph::AB_defined(graph_t &g, edge_t &one, edge_t &two, spin_type spin){
	
	int global_loop=g[boost::graph_bundle].fermi_loop_id;
// We need two cases.  One where the edges are the same.  and if not, they must be different.
if (one == two){
//	std::cout<<"Add bubble case 1"<< std::endl;



vertex_t new_leftb = add_vertex(vertex_info(three_leg), g);
vertex_t new_rightb = add_vertex(vertex_info(three_leg), g);

vertex_t new_leftt = add_vertex(vertex_info(three_leg), g);
vertex_t new_rightt = add_vertex(vertex_info(three_leg), g);

add_edge(source(one,g),new_leftb,edge_info(AmiBase::Fermi,g[one].fermi_loop_id, g[one].spin),g);
add_edge(new_leftb,new_rightb,edge_info(AmiBase::Fermi, g[one].fermi_loop_id, g[one].spin),g);
add_edge(new_rightb,target(one,g),edge_info(AmiBase::Fermi,g[one].fermi_loop_id, g[one].spin),g);

add_edge(new_leftb, new_leftt,edge_info(AmiBase::Bose),g);
add_edge(new_rightt,new_rightb,edge_info(AmiBase::Bose),g);

// add correct spin to bubble - opposite spin to connection 
// TODO: Is this what we actually want?  In general no. In general - want random spin...
// TODO: Check, does this mess with detailed balance at all?

//int rnum=random_int(0,1);
add_edge(new_leftt, new_rightt, edge_info(AmiBase::Fermi, global_loop+1, spin), g);
add_edge(new_rightt, new_leftt, edge_info(AmiBase::Fermi, global_loop+1, spin), g);

// This is commented out because it doesn't use random spin which I now want.
// if(g[one].spin==up){
// add_edge(new_leftt, new_rightt, edge_info(AmiBase::Fermi, global_loop+1, dn), g);
// add_edge(new_rightt, new_leftt, edge_info(AmiBase::Fermi, global_loop+1, dn), g); 
// }else{
	
	// add_edge(new_leftt, new_rightt, edge_info(AmiBase::Fermi, global_loop+1, up), g);
	// add_edge(new_rightt, new_leftt, edge_info(AmiBase::Fermi, global_loop+1, up), g); 
// }

remove_edge(one,g);


g[boost::graph_bundle].fermi_loop_id ++;


}
else{
//	std::cout<<"Add bubble case 2"<< std::endl;
// TODO: when these are different - we don't know how to label the spin of the new bubble to ensure that it is Hubbard like. So pick a random spin - 0 or 1
// TODO: must filter to check for 'Hubbard-like' diagrams




vertex_t new_onevert = add_vertex(vertex_info(three_leg), g);
vertex_t new_twovert = add_vertex(vertex_info(three_leg), g);

vertex_t new_bubblel = add_vertex(vertex_info(three_leg), g);
vertex_t new_bubbler = add_vertex(vertex_info(three_leg), g);


// Add edges
add_edge(source(one,g), new_onevert, edge_info(AmiBase::Fermi,g[one].fermi_loop_id, g[one].spin), g);
add_edge(new_onevert, target(one,g), edge_info(AmiBase::Fermi,g[one].fermi_loop_id,g[one].spin), g);


add_edge(source(two,g), new_twovert, edge_info(AmiBase::Fermi,g[two].fermi_loop_id, g[two].spin), g);
add_edge(new_twovert, target(two,g), edge_info(AmiBase::Fermi,g[two].fermi_loop_id, g[two].spin), g);


// Add bosonic edge
add_edge(new_onevert,new_bubblel, edge_info(AmiBase::Bose),g);
add_edge(new_bubbler, new_twovert, edge_info(AmiBase::Bose),g);

if(g[one].spin==g[two].spin){
if(g[one].spin==up){	
add_edge(new_bubblel,new_bubbler, edge_info(AmiBase::Fermi, global_loop+1, dn), g);
add_edge(new_bubbler,new_bubblel, edge_info(AmiBase::Fermi, global_loop+1,dn), g);
}else
{
add_edge(new_bubblel,new_bubbler, edge_info(AmiBase::Fermi, global_loop+1, up), g);
add_edge(new_bubbler,new_bubblel, edge_info(AmiBase::Fermi, global_loop+1,up), g);
}	


}else{ // THIS ELSE CAUSES NON-PHYSICAL DIAGRAMS WE WILL REJECT LATER?
int rnum=random_int(0,1);
add_edge(new_bubblel,new_bubbler, edge_info(AmiBase::Fermi, global_loop+1, AmiGraph::spin_type(rnum)), g);
add_edge(new_bubbler,new_bubblel, edge_info(AmiBase::Fermi, global_loop+1,AmiGraph::spin_type(rnum)), g);	
}

remove_edge(one,g);
remove_edge(two,g);

g[boost::graph_bundle].fermi_loop_id ++;

}// end else line
	
	
	
	
	

	
	
	
}


void AmiGraph::AB_rand(graph_t &g){
	
edge_t one, two;

int global_loop=g[boost::graph_bundle].fermi_loop_id;

//std::cout<<"Selecting two edges"<<std::endl;

// need function - pick two fermionic lines

find_two_edges(g, AmiBase::Fermi , one,  two);

// TODO find_two_edges with SAME spin, would fix this, but mess with detailed balance

// std::cout<<"two edges selected"<<std::endl;

// std::cout<< one<< " "<< two << std::endl;


// We need two cases.  One where the edges are the same.  and if not, they must be different.
if (one == two){
//	std::cout<<"Add bubble case 1"<< std::endl;



vertex_t new_leftb = add_vertex(vertex_info(three_leg), g);
vertex_t new_rightb = add_vertex(vertex_info(three_leg), g);

vertex_t new_leftt = add_vertex(vertex_info(three_leg), g);
vertex_t new_rightt = add_vertex(vertex_info(three_leg), g);

add_edge(source(one,g),new_leftb,edge_info(AmiBase::Fermi,g[one].fermi_loop_id, g[one].spin),g);
add_edge(new_leftb,new_rightb,edge_info(AmiBase::Fermi, g[one].fermi_loop_id, g[one].spin),g);
add_edge(new_rightb,target(one,g),edge_info(AmiBase::Fermi,g[one].fermi_loop_id, g[one].spin),g);

add_edge(new_leftb, new_leftt,edge_info(AmiBase::Bose),g);
add_edge(new_rightt,new_rightb,edge_info(AmiBase::Bose),g);

// add correct spin to bubble - opposite spin to connection 
// TODO: Is this what we actually want?  In general no. In general - want random spin...
// TODO: Check, does this mess with detailed balance at all?

int rnum=random_int(0,1);
add_edge(new_leftt, new_rightt, edge_info(AmiBase::Fermi, global_loop+1, AmiGraph::spin_type(rnum)), g);
add_edge(new_rightt, new_leftt, edge_info(AmiBase::Fermi, global_loop+1,AmiGraph::spin_type(rnum)), g);

// This is commented out because it doesn't use random spin which I now want.
// if(g[one].spin==up){
// add_edge(new_leftt, new_rightt, edge_info(AmiBase::Fermi, global_loop+1, dn), g);
// add_edge(new_rightt, new_leftt, edge_info(AmiBase::Fermi, global_loop+1, dn), g); 
// }else{
	
	// add_edge(new_leftt, new_rightt, edge_info(AmiBase::Fermi, global_loop+1, up), g);
	// add_edge(new_rightt, new_leftt, edge_info(AmiBase::Fermi, global_loop+1, up), g); 
// }

remove_edge(one,g);


g[boost::graph_bundle].fermi_loop_id ++;


}
else{
//	std::cout<<"Add bubble case 2"<< std::endl;
// TODO: when these are different - we don't know how to label the spin of the new bubble to ensure that it is Hubbard like. So pick a random spin - 0 or 1
// TODO: must filter to check for 'Hubbard-like' diagrams




vertex_t new_onevert = add_vertex(vertex_info(three_leg), g);
vertex_t new_twovert = add_vertex(vertex_info(three_leg), g);

vertex_t new_bubblel = add_vertex(vertex_info(three_leg), g);
vertex_t new_bubbler = add_vertex(vertex_info(three_leg), g);


// Add edges
add_edge(source(one,g), new_onevert, edge_info(AmiBase::Fermi,g[one].fermi_loop_id, g[one].spin), g);
add_edge(new_onevert, target(one,g), edge_info(AmiBase::Fermi,g[one].fermi_loop_id,g[one].spin), g);


add_edge(source(two,g), new_twovert, edge_info(AmiBase::Fermi,g[two].fermi_loop_id, g[two].spin), g);
add_edge(new_twovert, target(two,g), edge_info(AmiBase::Fermi,g[two].fermi_loop_id, g[two].spin), g);


// Add bosonic edge
add_edge(new_onevert,new_bubblel, edge_info(AmiBase::Bose),g);
add_edge(new_bubbler, new_twovert, edge_info(AmiBase::Bose),g);

if(g[one].spin==g[two].spin){
if(g[one].spin==up){	
add_edge(new_bubblel,new_bubbler, edge_info(AmiBase::Fermi, global_loop+1, dn), g);
add_edge(new_bubbler,new_bubblel, edge_info(AmiBase::Fermi, global_loop+1,dn), g);
}else
{
add_edge(new_bubblel,new_bubbler, edge_info(AmiBase::Fermi, global_loop+1, up), g);
add_edge(new_bubbler,new_bubblel, edge_info(AmiBase::Fermi, global_loop+1,up), g);
}	


}else{ // THIS ELSE CAUSES NON-PHYSICAL DIAGRAMS WE WILL REJECT LATER?
int rnum=random_int(0,1);
add_edge(new_bubblel,new_bubbler, edge_info(AmiBase::Fermi, global_loop+1, AmiGraph::spin_type(rnum)), g);
add_edge(new_bubbler,new_bubblel, edge_info(AmiBase::Fermi, global_loop+1,AmiGraph::spin_type(rnum)), g);	
}

remove_edge(one,g);
remove_edge(two,g);

g[boost::graph_bundle].fermi_loop_id ++;

}// end else line
	
	
	
	
	
}

void AmiGraph::DB_rand(graph_t &g){


bubble_type b_type;
	
AmiGraph::vertex_vector_list_t bubble_vertex_list;
AmiGraph::edge_vector_list_t bubble_edge_list;
AmiGraph::edge_vector_list_t legs_edge_list;
int bubble_choice;

bubble_finder(g, bubble_vertex_list,  bubble_edge_list, legs_edge_list);

vertex_vector_t bubble_vect;

// std::cout << "Legs edge list stuff" <<std::endl;
// std::cout<< legs_edge_list.size()<<std::endl;
// std::cout<< legs_edge_list[0].size()<<std::endl;

// for(int i=0;i<legs_edge_list.size(); i++){
	// std::cout<< legs_edge_list[i].size()<<std::endl;
	
// }


//std::cout<< "Found bubbles Nb = "<< bubble_vertex_list.size()<< " expecting to find "<<g[boost::graph_bundle].fermi_loop_id<<std::endl;

if( bubble_vertex_list.size() != 0){
	
	bubble_choice=random_int(0, bubble_vertex_list.size()-1);
	//std::cout<<"Bubble choice is "<< bubble_choice << std::endl;
	
	bubble_vect=bubble_vertex_list[bubble_choice];
	get_bubble_type(g, b_type, bubble_vertex_list[bubble_choice], legs_edge_list[bubble_choice]);
	
	vertex_t B1;
	vertex_t B2;
	
	vertex_t V1, V2;
	edge_t e1,e2;
	
	//std::cout<<"Bubble vect "<<std::endl;
	//std::cout<< bubble_vect[0]<<" "<<bubble_vect[1]<<std::endl;
	//std::cout<<g[bubble_vect[0]].index_<<" "<< g[bubble_vect[1]].index_<<std::endl;
	
	
	// for (int i=0; i<bubble_edge_list[bubble_choice].size(); i++){
	// std::cout<<"Bubble edges are "<<std::endl;
	// std::cout<< bubble_edge_list[bubble_choice][i]<<std::endl;
	// std::cout<<AmiGraph::print_edge(bubble_edge_list[bubble_choice][i],g)<<std::endl;
	
	
	
// }

// for (int i=0; i<bubble_vertex_list[bubble_choice].size(); i++){
	// std::cout<<"Bubble vertices are "<<std::endl;
	// std::cout<<g[bubble_vertex_list[bubble_choice][i]].index_<<std::endl;
	
// }

// for (int i=0; i<legs_edge_list[bubble_choice].size(); i++){
	// std::cout<<"Bubble bosonic legs are"<<std::endl;
	// std::cout<<AmiGraph::print_edge(legs_edge_list[bubble_choice][i],g)<<std::endl;
	
// }
	
	
	// std::cout<<"Made it here"<<std::endl;
	
	// logic to artificially label directionality of bosonic lines- please check
	// Want to define V1, B1, B2, V2 according to picture
	// logic sometimes fails
	//
	// V1-bose-B1-loop-B2-V2 
if(b_type==directed){	
	if( target(legs_edge_list[bubble_choice][0],g)== bubble_vect[0] || target(legs_edge_list[bubble_choice][0],g)== bubble_vect[1] ){
	// std::cout<<"logic 1"<<std::endl;
		V1=source(legs_edge_list[bubble_choice][0],g);
		B1=target(legs_edge_list[bubble_choice][0],g);
		
		B2=source(legs_edge_list[bubble_choice][1],g);
		V2=target(legs_edge_list[bubble_choice][1],g);
		
		e1=legs_edge_list[bubble_choice][0];
		e2=legs_edge_list[bubble_choice][1];
		
				
	}else 	
	if(source(legs_edge_list[bubble_choice][0],g)== bubble_vect[0] || source(legs_edge_list[bubble_choice][0],g)== bubble_vect[1]){
	// std::cout<<"logic 2"<<std::endl;
		V1=source(legs_edge_list[bubble_choice][1],g);
		B1=target(legs_edge_list[bubble_choice][1],g);
		
		B2=source(legs_edge_list[bubble_choice][0],g);
		V2=target(legs_edge_list[bubble_choice][0],g);
		
		e1=legs_edge_list[bubble_choice][0];
		e2=legs_edge_list[bubble_choice][1];
				
	}
}

if(b_type==both_in){
	
	V1=source(legs_edge_list[bubble_choice][0],g);
	V2=source(legs_edge_list[bubble_choice][1],g);
	
	B1=target(legs_edge_list[bubble_choice][0],g);
	B2=target(legs_edge_list[bubble_choice][1],g);	
	e1=legs_edge_list[bubble_choice][0];
	e2=legs_edge_list[bubble_choice][1];
	
	
}

if(b_type==both_out){
	
	V1=target(legs_edge_list[bubble_choice][0],g);
	V2=target(legs_edge_list[bubble_choice][1],g);
	
	B1=source(legs_edge_list[bubble_choice][0],g);
	B2=source(legs_edge_list[bubble_choice][1],g);	
	e1=legs_edge_list[bubble_choice][0];
	e2=legs_edge_list[bubble_choice][1];
	
	
}




	
		
	// next, check if V1 or V2 share an in or out edge 	
		
		
edge_t V1_in, V1_out;
edge_t V2_in, V2_out;

edge_t shared;

// std::cout<<"Finding Edges for V1 "<< g[V1].index_ <<std::endl;
	find_in_out_edges_oftype(g, V1, V1_in,V1_out, AmiBase::Fermi); 
	// std::cout<< V1_in<<" "<<V1_out<<std::endl;
	 // std::cout<< AmiGraph::print_edge(V1_in,g)<<" "<<AmiGraph::print_edge(V1_out,g) <<std::endl;
	
	 // std::cout<<"Finding Edges for V2 "<< g[V2].index_<<std::endl;
	find_in_out_edges_oftype(g, V2, V2_in,V2_out, AmiBase::Fermi); 
	// std::cout<< V2_in<<" "<< V2_out<<std::endl;
     // std::cout<< AmiGraph::print_edge(V2_in,g) <<" "<< AmiGraph::print_edge(V2_out,g) <<std::endl;
// std::cout<<"Made it here 2"<<std::endl;
	
	// add edges around V1
	if(V1_out==V2_in || V1_in==V2_out){
		
		if(V1_out==V2_in){shared=V1_out;  
// std::cout<<"shared line is "<<AmiGraph::print_edge(shared,g)<<std::endl;
// std::cout<<"Adding edge from "<<g[source(V1_in,g)].index_<<" to "<<g[target(V2_out,g)].index_<<std::endl;		
		add_edge(source(V1_in,g),target(V2_out,g),edge_info(AmiBase::Fermi,g[V1_in].fermi_loop_id, g[V1_in].spin),g);
		}
		if(V1_in==V2_out){shared=V2_out;
		add_edge(source(V2_in,g),target(V1_out,g),edge_info(AmiBase::Fermi,g[V2_in].fermi_loop_id, g[V2_in].spin),g);
		}
		
		
		// std::cout<<"Made it here 3"<<std::endl;
		
	// std::cout<<"Removing vertices "<< V1<<" "<<V2 <<" "<<B1<<" "<<B2 <<std::endl;
  	
	

	clear_vertex(B1,g);
	clear_vertex(B2,g);
	remove_edge(shared,g);
	clear_vertex(V1,g);
	clear_vertex(V2,g);
	
	remove_vertex(B1,g);
	remove_vertex(B2,g);
	remove_vertex(V1,g);
	remove_vertex(V2,g);

	//g[boost::graph_bundle].fermi_loop_id --;
	
	
		
		// std::cout<<"Made it here 4"<<std::endl;
	}else{
		// std::cout<<"Made it here 5"<<std::endl;
	// if(source(V1_in,g)==target(V2_out,g)){
    // add_edge(source(V1_in,g),target(V2_out,g),edge_info(AmiBase::Fermi,g[V1_in].fermi_loop_id),g);
	// }else{
	add_edge(source(V1_in,g),target(V1_out,g),edge_info(AmiBase::Fermi,g[V1_in].fermi_loop_id, g[V1_in].spin),g);
	add_edge(source(V2_in,g),target(V2_out,g),edge_info(AmiBase::Fermi,g[V2_in].fermi_loop_id, g[V2_in].spin),g);
	// }
	// std::cout<<"added line from "<<g[source(V1_in,g)].index_<<" to "<< g[target(V1_out,g)].index_<<std::endl;
	// std::cout<<"added line from "<<g[source(V2_in,g)].index_<<" to "<< g[target(V2_out,g)].index_<<std::endl;
	// std::cout<<"Removing vertices "<< g[V1].index_<<" "<<g[V2].index_ <<" "<<g[B1].index_<<" "<<g[B2].index_ <<std::endl;
	

	clear_vertex(V1,g);
	clear_vertex(V2,g);
	clear_vertex(B1,g);
	clear_vertex(B2,g);
	remove_vertex(B1,g);
	remove_vertex(B2,g);
	remove_vertex(V1,g);
	remove_vertex(V2,g);

	//g[boost::graph_bundle].fermi_loop_id --;
	
	// std::cout<<"Made it here 5.1"<<std::endl;
	
	// std::cout<<"Made it here 5.2"<<std::endl;
		
	
	// std::cout<<"Made it here 6"<<std::endl;
	
	   
		
		
	}	
		



	
	


	
		
}

	
}







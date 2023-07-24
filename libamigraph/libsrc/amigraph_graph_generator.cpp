#include "amigraph.hpp"


void AmiGraph::generate_force_graphs(std::vector< std::vector< AmiGraph::graph_t>> &graph_matrix, int max, int mpi_rank){
	
	
// resize the order container 
graph_matrix.resize(max+1);

std::cout<<"Generating force Graphs Starting with graph order "<< graph_order(current_graph)<<std::endl;


graph_matrix[graph_order(current_graph)].push_back(current_graph);
graph_matrix[graph_order(f2)].push_back(f2);
graph_matrix[graph_order(f3)].push_back(f3);
//int max=6;

for (int order=graph_order(current_graph); order<max; order++){
if(mpi_rank==0){
std::cout<<"Working on order "<< order <<std::endl;}
if(	graph_matrix[order].size() >0){
	for (int i=0; i< graph_matrix[order].size(); i++){
		
		std::cout<<i<<std::endl;

	// g.current_graph=graph_matrix[order][i];
	
	
	// AmiGraph::edge_vector_t edge_list;
	// g.find_internal_fermionic_edges(g.current_graph,edge_list);
	// std::cout<<"There are "<< edge_list.size()<<" fermionic edges"<<std::endl;
	
	// g.proposed_graph=g.current_graph;
	// std::cout<<source(edge_list[0], g.current_graph)<<" "<< source(edge_list[0], g.proposed_graph)<<std::endl;
	// Now have edge list for a specific graph.  First, do tadpoles. then do interaction lines
	// std::cout<<"On graph "<<i+1<< " out of "<< graph_matrix[order].size() <<" for order "<< order <<std::endl;
	
	current_graph=graph_matrix[order][i];
	
	// std::cout<<"Attempting systematic vertex numbering"<<std::endl;
	systematic_vertex_label(current_graph);
	// std::cout<<"Completed"<<std::endl;
	
	// print_all_edge_info(current_graph);
	
	AmiGraph::edge_vector_t edge_list, translate_list;
	find_fermionic_edges(current_graph,edge_list);
	
	// only internal edges for tadpoles 
	// std::cout<<"Inserting tadpoles"<<std::endl;
	// std::cout<<"Edge list size is "<< edge_list.size()<<std::endl;
	for (int spin=0; spin<2; spin++){
	for(int tint=0; tint<edge_list.size();tint++){
		 // std::cout<<"tint and spin are "<<tint <<" "<< spin<<std::endl;
	
	proposed_graph=current_graph;
	//g.number_vertices(g.proposed_graph);
	find_fermionic_edges(proposed_graph,translate_list); // TODO: this is pretty dangerous... but it works so....
	
	// std::cout<<"Attempting to add tadpole between source/target "<< proposed_graph[source(translate_list[tint],proposed_graph)].index_<<" "<< proposed_graph[target(translate_list[tint],proposed_graph)].index_<<std::endl;
	
	AT_defined(proposed_graph,translate_list[tint], AmiGraph::spin_type(spin));
	
	// std::cout<<"Did it"<<std::endl;
	number_vertices(proposed_graph);
	// print_all_edge_info(proposed_graph);

	//g.number_vertices(g.proposed_graph);	
	//g.print_all_edge_info(g.proposed_graph);

// std::cout<<"Updating proposed"<<std::endl;
	graph_matrix_update(proposed_graph, graph_matrix); // this function checks isomorphisms and inserts the graph if it is new. 
// std::cout<<"Complete"<<std::endl;
  // std::cout<<"After attempting to add the graph list size is "<<graph_matrix[order+1].size()<<std::endl;

	// std::cout<<"Exit tint loop"<<std::endl;		
	}
	}
	
	// now do all interaction line insertions, this includes the external degree 1 vertices 
	
	find_fermionic_edges(current_graph, edge_list);
	
	// for all combinations of edge choices, propose a change and check isomorphism 
	// std::cout<<"Inserting AIL"<<std::endl;
	// std::cout<<"Edge list size is "<< edge_list.size()<<std::endl;
	for(int li=0; li<edge_list.size();li++){
		
		for(int lj=li; lj<edge_list.size(); lj++){
			// std::cout<<li<<" "<< lj<<std::endl;
		proposed_graph=current_graph;
	    //g.number_vertices(g.proposed_graph);
	    find_fermionic_edges(proposed_graph,translate_list);
		
		AIL_defined(proposed_graph, translate_list[li], translate_list[lj]);
	// std::cout<<"----"<<std::endl;
// std::cout<<"Proposed is "<<std::endl;	
	// g.number_vertices(g.proposed_graph);	
	// g.print_all_edge_info(g.proposed_graph);
// std::cout<<"----"<<std::endl;
	graph_matrix_update(proposed_graph, graph_matrix);

		
		}		
	}
	
	// now replace each interaction line with a bubble - probably only need to do the last line, but for safety, lets check them all
	// each bubble could have spin up or spin down
	
	// TODO: I think this is not necessary in general
	// Definitely causes problems for force graph type
	find_non_external_bosonic_edges(current_graph,edge_list);
	// std::cout<<"Found "<<edge_list.size()<<" non-entry bosonic edges"<<std::endl;
	for (int spin=0; spin<2; spin++){
	for(int li=0; li< edge_list.size(); li++){
	proposed_graph=current_graph;
	// find_bosonic_edges(proposed_graph,translate_list);
	find_non_external_bosonic_edges(proposed_graph,translate_list);
	
	RIB_defined(proposed_graph,translate_list[li], spin_type(spin));
	
	graph_matrix_update(proposed_graph, graph_matrix);
		
	}
	} 
	
	
	

	}	
		
}
if(mpi_rank==0){
// std::cout<<"Finished order "<< order <<std::endl;
std::cout<<"Total number of diagrams of order "<<order+1<<" is "<< graph_matrix[order+1].size()<<std::endl;	
}
} 	
	
	
	
	
	
}


void AmiGraph::generate_bubble_graphs(std::vector< std::vector< AmiGraph::graph_t>> &graph_matrix, int max, int mpi_rank){
	
	
// resize the order container 
graph_matrix.resize(max+1);

std::cout<<"Generating Bubble Graphs Starting with graph order "<< graph_order(current_graph)<<std::endl;
number_vertices(current_graph);
// print_all_edge_info(current_graph);

graph_matrix[graph_order(current_graph)].push_back(current_graph);
//int max=6;

for (int order=graph_order(current_graph); order<max; order++){
if(mpi_rank==0){
std::cout<<"Working on order "<< order <<std::endl;}
if(	graph_matrix[order].size() >0){
	for (int i=0; i< graph_matrix[order].size(); i++){
		
		// std::cout<<i<<std::endl;

	// g.current_graph=graph_matrix[order][i];
	
	
	// AmiGraph::edge_vector_t edge_list;
	// g.find_internal_fermionic_edges(g.current_graph,edge_list);
	// std::cout<<"There are "<< edge_list.size()<<" fermionic edges"<<std::endl;
	
	// g.proposed_graph=g.current_graph;
	// std::cout<<source(edge_list[0], g.current_graph)<<" "<< source(edge_list[0], g.proposed_graph)<<std::endl;
	// Now have edge list for a specific graph.  First, do tadpoles. then do interaction lines
	// std::cout<<"On graph "<<i+1<< " out of "<< graph_matrix[order].size() <<" for order "<< order <<std::endl;
	
	current_graph=graph_matrix[order][i];
	
	systematic_vertex_label(current_graph);
	// print_all_edge_info(current_graph);
	
	AmiGraph::edge_vector_t edge_list, translate_list;
	find_fermionic_edges(current_graph,edge_list);
	
	// only internal edges for tadpoles 
	// std::cout<<"Inserting tadpoles"<<std::endl;
	// std::cout<<"Edge list size is "<< edge_list.size()<<std::endl;
	for (int spin=0; spin<2; spin++){
	for(int tint=0; tint<edge_list.size();tint++){
		 // std::cout<<"tint and spin are "<<tint <<" "<< spin<<std::endl;
	
	proposed_graph=current_graph;
	//g.number_vertices(g.proposed_graph);
	find_fermionic_edges(proposed_graph,translate_list); // TODO: this is pretty dangerous... but it works so....
	
	// std::cout<<"Attempting to remove edge from source/target "<< proposed_graph[source(translate_list[tint],proposed_graph)].index_<<" "<< proposed_graph[target(translate_list[tint],proposed_graph)].index_<<std::endl;
	
	AT_defined(proposed_graph,translate_list[tint], AmiGraph::spin_type(spin));
	
	// std::cout<<"Did it"<<std::endl;
	// number_vertices(proposed_graph);
	// print_all_edge_info(proposed_graph);

	//g.number_vertices(g.proposed_graph);	
	//g.print_all_edge_info(g.proposed_graph);

// std::cout<<"Updating proposed"<<std::endl;
	graph_matrix_update(proposed_graph, graph_matrix); // this function checks isomorphisms and inserts the graph if it is new. 
// std::cout<<"Complete"<<std::endl;

	// std::cout<<"Exit tint loop"<<std::endl;		
	}
	}
	
	// now do all interaction line insertions, this includes the external degree 1 vertices 
	
	find_fermionic_edges(current_graph, edge_list);
	
	// for all combinations of edge choices, propose a change and check isomorphism 
	// std::cout<<"Inserting AIL"<<std::endl;
	// std::cout<<"Edge list size is "<< edge_list.size()<<std::endl;
	for(int li=0; li<edge_list.size();li++){
		
		for(int lj=li; lj<edge_list.size(); lj++){
			// std::cout<<li<<" "<< lj<<std::endl;
		proposed_graph=current_graph;
	    //g.number_vertices(g.proposed_graph);
	    find_fermionic_edges(proposed_graph,translate_list);
		
		AIL_defined(proposed_graph, translate_list[li], translate_list[lj]);
	// std::cout<<"----"<<std::endl;
// std::cout<<"Proposed is "<<std::endl;	
	// g.number_vertices(g.proposed_graph);	
	// g.print_all_edge_info(g.proposed_graph);
// std::cout<<"----"<<std::endl;
	graph_matrix_update(proposed_graph, graph_matrix);

		
		}		
	}
	
	// now replace each interaction line with a bubble - probably only need to do the last line, but for safety, lets check them all
	// each bubble could have spin up or spin down
	
	find_non_entry_bosonic_edges(current_graph,edge_list);
	
	for (int spin=0; spin<2; spin++){
	for(int li=0; li< edge_list.size(); li++){
	proposed_graph=current_graph;
	// find_bosonic_edges(proposed_graph,translate_list);
	find_non_entry_bosonic_edges(proposed_graph,translate_list);
	
	RIB_defined(proposed_graph,translate_list[li], spin_type(spin));
	
	graph_matrix_update(proposed_graph, graph_matrix);
		
	}
	}
	
	

	}	
		
}
if(mpi_rank==0){
// std::cout<<"Finished order "<< order <<std::endl;
std::cout<<"Total number of diagrams of order "<<order+1<<" is "<< graph_matrix[order+1].size()<<std::endl;	
}
} 	
	
	
	
	
	
}


void AmiGraph::generate_sigma_graphs(std::vector< std::vector< AmiGraph::graph_t>> &graph_matrix, int max, int mpi_rank){
	
	
// resize the order container 
graph_matrix.resize(max+1);

graph_matrix[graph_order(current_graph)].push_back(current_graph);
//int max=6;

for (int order=graph_order(current_graph); order<max; order++){
if(mpi_rank==0){
std::cout<<"Working on order "<< order <<std::endl;}
if(	graph_matrix[order].size() >0){
	for (int i=0; i< graph_matrix[order].size(); i++){

	// g.current_graph=graph_matrix[order][i];
	
	
	// AmiGraph::edge_vector_t edge_list;
	// g.find_internal_fermionic_edges(g.current_graph,edge_list);
	// std::cout<<"There are "<< edge_list.size()<<" fermionic edges"<<std::endl;
	
	// g.proposed_graph=g.current_graph;
	// std::cout<<source(edge_list[0], g.current_graph)<<" "<< source(edge_list[0], g.proposed_graph)<<std::endl;
	// Now have edge list for a specific graph.  First, do tadpoles. then do interaction lines
	// std::cout<<"On graph "<<i+1<< " out of "<< graph_matrix[order].size() <<" for order "<< order <<std::endl;
	
	current_graph=graph_matrix[order][i];
	
	AmiGraph::edge_vector_t edge_list, translate_list;
	find_fermionic_edges(current_graph,edge_list);
	
	// only internal edges for tadpoles 
	// std::cout<<"Inserting tadpoles"<<std::endl;
	// std::cout<<"Edge list size is "<< edge_list.size()<<std::endl;
	for (int spin=0; spin<2; spin++){
	for(int tint=0; tint<edge_list.size();tint++){
		 // std::cout<<tint <<" "<< spin<<std::endl;
	
	proposed_graph=current_graph;
	//g.number_vertices(g.proposed_graph);
	find_fermionic_edges(proposed_graph,translate_list); // TODO: this is pretty dangerous... but it works so....
	
	//std::cout<<"Attempting to remove edge from source/target "<< g.proposed_graph[source(translate_list[tint],g.proposed_graph)].index_<<" "<< g.proposed_graph[target(translate_list[tint],g.proposed_graph)].index_<<std::endl;
	
	AT_defined(proposed_graph,translate_list[tint], AmiGraph::spin_type(spin));

	//g.number_vertices(g.proposed_graph);	
	//g.print_all_edge_info(g.proposed_graph);



	graph_matrix_update(proposed_graph, graph_matrix); // this function checks isomorphisms and inserts the graph if it is new. 
			
	}
	}
	
	// now do all interaction line insertions, this includes the external degree 1 vertices 
	
	find_fermionic_edges(current_graph, edge_list);
	
	// for all combinations of edge choices, propose a change and check isomorphism 
	// std::cout<<"Inserting AIL"<<std::endl;
	// std::cout<<"Edge list size is "<< edge_list.size()<<std::endl;
	for(int li=0; li<edge_list.size();li++){
		
		for(int lj=li; lj<edge_list.size(); lj++){
			// std::cout<<li<<" "<< lj<<std::endl;
		proposed_graph=current_graph;
	    //g.number_vertices(g.proposed_graph);
	    find_fermionic_edges(proposed_graph,translate_list);
		
		AIL_defined(proposed_graph, translate_list[li], translate_list[lj]);
	// std::cout<<"----"<<std::endl;
// std::cout<<"Proposed is "<<std::endl;	
	// g.number_vertices(g.proposed_graph);	
	// g.print_all_edge_info(g.proposed_graph);
// std::cout<<"----"<<std::endl;
	graph_matrix_update(proposed_graph, graph_matrix);

		
		}		
	}
	
	
	


	}	
		
}
if(mpi_rank==0){
std::cout<<"Finished order "<< order <<std::endl;
std::cout<<"Total number of diagrams of order "<<order+1<<" is "<< graph_matrix[order+1].size()<<std::endl;	
}
} 	
	
	
	
	
	
}



// Greens function graphs look almost identical to sigma graphs - but allow for 1p reducible graphs 
// No separate function required. Can just use sigma - and flag for Greens type to keep 1p reducible and to close the loops
void AmiGraph::generate_greens_graphs(std::vector< std::vector< AmiGraph::graph_t>> &graph_matrix, int max, int mpi_rank){
	
	
// resize the order container 
graph_matrix.resize(max+1);

graph_matrix[graph_order(current_graph)].push_back(current_graph);

// std::cout<<"Generating graphs starting from "<<std::endl;
// number_vertices(current_graph);
// print_all_edge_info(current_graph);

//int max=6;

for (int order=graph_order(current_graph); order<max; order++){
if(mpi_rank==0){
std::cout<<"Working on order "<< order <<std::endl;}
if(	graph_matrix[order].size() >0){
	for (int i=0; i< graph_matrix[order].size(); i++){

		
	current_graph=graph_matrix[order][i];
	
	AmiGraph::edge_vector_t edge_list, translate_list;
	find_fermionic_edges(current_graph,edge_list);
	
	// only internal edges for tadpoles 
	// std::cout<<"Inserting tadpoles"<<std::endl;
	// std::cout<<"Edge list size is "<< edge_list.size()<<std::endl;
	for (int spin=0; spin<2; spin++){
	for(int tint=0; tint<edge_list.size();tint++){
		 // std::cout<<tint <<" "<< spin<<std::endl;
	
	proposed_graph=current_graph;
	//g.number_vertices(g.proposed_graph);
	find_fermionic_edges(proposed_graph,translate_list); // TODO: this is pretty dangerous... but it works so....
	
	
	AT_defined(proposed_graph,translate_list[tint], AmiGraph::spin_type(spin));



	graph_matrix_update(proposed_graph, graph_matrix); // this function checks isomorphisms and inserts the graph if it is new. 
			
	}
	}
	
	// now do all interaction line insertions, this includes the external degree 1 vertices 
	
	find_fermionic_edges(current_graph, edge_list);
	
	// for all combinations of edge choices, propose a change and check isomorphism 
	// std::cout<<"Inserting AIL"<<std::endl;
	// std::cout<<"Edge list size is "<< edge_list.size()<<std::endl;
	for(int li=0; li<edge_list.size();li++){
		
		for(int lj=li; lj<edge_list.size(); lj++){
			// std::cout<<li<<" "<< lj<<std::endl;
		proposed_graph=current_graph;
	    //g.number_vertices(g.proposed_graph);
	    find_fermionic_edges(proposed_graph,translate_list);
		
		AIL_defined(proposed_graph, translate_list[li], translate_list[lj]);
	// std::cout<<"----"<<std::endl;
// std::cout<<"Proposed is "<<std::endl;	
	// g.number_vertices(g.proposed_graph);	
	// g.print_all_edge_info(g.proposed_graph);
// std::cout<<"----"<<std::endl;
	graph_matrix_update(proposed_graph, graph_matrix);

		
		}		
	}
	
	
	


	}	
		
}
if(mpi_rank==0){
std::cout<<"Finished order "<< order <<std::endl;
std::cout<<"Total number of diagrams of order "<<order+1<<" is "<< graph_matrix[order+1].size()<<std::endl;	
}
} 	
	

// if(ami_parameters.TYPE_==AmiBase::density){	
// if(mpi_rank==0){
// std::cout<<"Finished Graph Generation - closing Greens function loops " <<std::endl;
// }

// gm_close_loops(graph_matrix);	
// }	
	
}

void AmiGraph::gm_close_loops(std::vector< std::vector< AmiGraph::graph_t>> &gm){
	
std::cout<<"Closing loops"<<std::endl;
for (int ord=0; ord<gm.size(); ord++){
for(int m=0; m< gm[ord].size(); m++){	
	
close_loop(gm[ord][m]);	
	
}
}
	
std::cout<<"Closing loops complete"<<std::endl;	
	
	
}



void AmiGraph::gm_attach_probe_lines(std::vector< std::vector< AmiGraph::graph_t>> &gm,std::vector< std::vector< AmiGraph::graph_t>> &gm_out){
gm_out.clear();
gm_out.resize(gm.size());
	
std::cout<<"Attaching probe lines"<<std::endl;
for (int ord=0; ord<gm.size(); ord++){
for(int m=0; m< gm[ord].size(); m++){	

std::cout<<ord <<" "<<m<<std::endl;
	
std::vector< AmiGraph::graph_t> with_probes;	
	
attach_probe_lines(gm[ord][m],with_probes);	

gm_out[ord].insert(gm_out[ord].end(),with_probes.begin(), with_probes.end());
	
}
}
	
std::cout<<"Closing loops complete"<<std::endl;	
	
	
}

void AmiGraph::attach_probe_lines(AmiGraph::graph_t &g, std::vector< AmiGraph::graph_t> &gm_out){
	
// get list of bosonic lines 
edge_vector_t bose_edges;

find_bosonic_edges(g, bose_edges);

if(bose_edges.size()<2){
	throw std::runtime_error("Can't probe if less than two bosonic edges");
}

// given length of bosonic list 

for(int i=0; i< bose_edges.size()-1; i++){
	for(int j=i+1; j< bose_edges.size(); j++){
		
		std::cout<<"On "<<i<<" "<< j<<std::endl;
	
		AmiGraph::graph_t gc=g;
		find_bosonic_edges(gc, bose_edges);
		
		probe_lines(gc, bose_edges[i],bose_edges[j]);
	
		bool add=true;
		for(int m=0; m< gm_out.size(); m++){
			
		if(is_isomorphic(gm_out[m],gc)){
		add=false;
		break;
		}			
			
		}
		
		if(add){
			gm_out.push_back(gc);
		}
		
	
	}
}

// if length<2 throw error 

	
	
	
}

void AmiGraph::probe_lines(graph_t &g, edge_t &e1, edge_t &e2){
	
	// make two vertices
	// add bosonic lines to the source and target of e1 
	// repeat for e2 
vertex_t in_vert1 = add_vertex(g);
vertex_t in_vert2 = add_vertex(g);

vertex_t out_vert1 = add_vertex(g);
vertex_t out_vert2 = add_vertex(g);

// Add interaction lines that are Bose not Fermi
add_edge(in_vert1,source(e1,g), edge_info(AmiBase::Bose),g);
add_edge(in_vert2,target(e1,g), edge_info(AmiBase::Bose),g);

// Add interaction lines that are Bose not Fermi
add_edge(source(e2,g),out_vert1, edge_info(AmiBase::Bose),g);
add_edge(target(e2,g),out_vert2, edge_info(AmiBase::Bose),g);		
	

	
}



void AmiGraph::close_ggm(gg_matrix_t &ggm){
	
for(int ord=0; ord<ggm.size(); ord++){
for(int group=0; group<ggm[ord].size(); group++){

for(int m=0; m< ggm[ord][group].graph_vec.size(); m++){

close_loop(ggm[ord][group].graph_vec[m]);

}

for(int p=0; p< ggm[ord][group].gp_vec.size(); p++){
	
close_loop(ggm[ord][group].gp_vec[p].g1_);
close_loop(ggm[ord][group].gp_vec[p].g2_);	
}


}
}	
	
	
}

// TODO: the labelling system may not be able to label closed loops 

// TODO: this does not properly close a bare G line.

// THE close_loop function is designed to take a labelled or unlabelled graph and close the loop.  But the labelling functions do not know how to label closed diagrams.  So in practice, the diagrams should be labelled and then closed, then converted to ami. 

// TODO: this can only be called with double-occ graphs 
void AmiGraph::close_loop(AmiGraph::graph_t &g){

AmiBase::stat_type this_stat;

if(graph_type==AmiBase::density){
this_stat=AmiBase::Fermi;}
else if (graph_type==AmiBase::doubleocc){
	this_stat=AmiBase::Bose;
}
else{
	throw std::runtime_error("Can't close graphs that are not density or doubleocc types");
}

	
// std::cout<<"Closing graph with properties"<<std::endl;
// print_all_edge_info(g);	
	
if(num_edges(g)>1){
// find external legs 

// label the external legs
vertex_vector_t extern_vert_list;
edge_vector_t extern_edge_list;
find_external_vertices(g, extern_vert_list, extern_edge_list);

// std::cout<<"Checking exit"<<std::endl;
// exit condition
assert(extern_vert_list.size() == 2 && "More than 2 external vertices - can't close loop");
//

// decide source and target - track degree==1 vertices
// std::cout<<"Deciding source and target"<<std::endl;

vertex_t vin, vout;

if(in_degree(extern_vert_list[0],g)==0){
vin=extern_vert_list[0];
vout=extern_vert_list[1];	
}
else{
vin=extern_vert_list[1];
vout=extern_vert_list[0];	
}

// std::cout<<"Vin and Vout are "<<g[vin].index_<<" "<<g[vout].index_<<std::endl;

vertex_t tar, sour;


edge_t ext_edge1, ext_edge2;

// boost::graph_traits<graph_t>::vertex_iterator vi, vi_end;
// for (boost::tie(vi,vi_end) = vertices(g); vi != vi_end; ++vi){

// if( edge(vin,*vi,g).second && ){
// tar=*vi;
// ext_edge1=edge(vin,*vi,g).first;
// }

// if(edge(*vi,vout,g).second){
// sour=*vi;
// ext_edge2=edge(*vi,vout,g).first;
// }	
	
// }

boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
for(boost::tie(ei,ei_end)=edges(g); ei!=ei_end; ++ei){

// std::cout<<"Comparing "<< g[source(*ei,g)].index_<<" "<< g[target(*ei,g)].index_<<std::endl;
// std::cout<<"To "<<g[vin].index_ <<" and "<< g[vout].index_<<std::endl;

if(source(*ei,g)==vin){ // && g[*ei].g_struct_.stat_==this_stat){
	// std::cout<<"Assigning edge1"<<std::endl;
	ext_edge1=*ei;
	tar=target(*ei,g);
}

if (target(*ei,g)==vout){ // && g[*ei].g_struct_.stat_==this_stat ){
	// std::cout<<"Assigning edge2"<<std::endl;
	ext_edge2=*ei;
	sour=source(*ei,g);
}

	
	
}




// std::cout<<"Sour and tar are "<< g[sour].index_<<" "<< g[tar].index_<<std::endl;


// std::cout<<"Closing loop edge"<<std::endl;
// this edge closes the loop 
edge_t added_edge=add_edge(sour,tar,edge_info(g[ext_edge1].g_struct_.stat_,g[ext_edge1].fermi_loop_id, g[ext_edge1].spin),g).first;
// now adopt the properties of ext_edge1 

// std::cout<<"Added edge"<<std::endl;
// number_vertices(g);
// print_all_edge_info(g);


g[added_edge].g_struct_=g[ext_edge1].g_struct_;

// std::cout<<"Set gstruct"<<std::endl;

// print_all_edge_info(g);

// TODO: So now there is an issue - the external leg had no specified external epsilon attached to it, and now we need to add that.
// for now we presume that it is the last entry - but this should be user specified.
// only do this for this_stat==AmiBase::Fermi 
if(this_stat==AmiBase::Fermi){
int last=g[added_edge].g_struct_.eps_.size()-1;
g[added_edge].g_struct_.eps_[last]=1;
}
///
// for hubbard interaction in doubleocc that line gets dropped so it does not matter
// for density that line does NOT get dropped - so it might matter 

// clear degree==1 vertices 
// get properties of the external line. add new edge with those properties 


// std::cout<<"Before cleanup graph is "<<std::endl;
// number_vertices(g);
// print_all_edge_info(g);

// std::cout<<"Cleaning up"<<std::endl;
remove_edge(ext_edge1,g);
remove_edge(ext_edge2,g);
clear_vertex(vin,g);
clear_vertex(vout,g);
remove_vertex(vin,g);
remove_vertex(vout,g);

// std::cout<<"loop closed"<<std::endl;
// ALL DONE!

}
else{

edge_t e1=random_edge(g,rand_gen);
vertex_t tar, sour;
tar=target(e1,g);
sour=source(e1,g);

// std::cout<<"Triggered random closure "<<std::endl;
edge_t added_edge=add_edge(sour,sour,edge_info(g[e1].g_struct_.stat_,g[e1].fermi_loop_id, g[e1].spin),g).first;
// now adopt the properties of ext_edge1 

// std::cout<<"Added edge"<<std::endl;
// number_vertices(g);
// print_all_edge_info(g);


g[added_edge].g_struct_=g[e1].g_struct_;

//TODO: this is hardcoded and probably should be in the labelling 
// g[added_edge].g_struct_.eps_[0]=1;

// add_edge(sour,sour,edge_info(g[e1].g_struct_.stat_,g[e1].fermi_loop_id, g[e1].spin),g);

remove_edge(e1,g);
remove_vertex(tar,g);


}	
	
	
// std::cout<<"On output graph is"<<std::endl;
// number_vertices(g);
// print_all_edge_info(g);	
	
}


void AmiGraph::generate_graphs(std::vector< std::vector< AmiGraph::graph_t>> &graph_matrix, int max, int mpi_rank){

if((ami_parameters.TYPE_==AmiBase::Sigma) || ami_parameters.TYPE_==AmiBase::Bare){
	generate_sigma_graphs(graph_matrix, max, mpi_rank);
}

if(ami_parameters.TYPE_==AmiBase::Pi_phuu || ami_parameters.TYPE_==AmiBase::Pi_phud || ami_parameters.TYPE_==AmiBase::doubleocc || ami_parameters.TYPE_==AmiBase::Pi_ppud || ami_parameters.TYPE_==AmiBase::Pi_ppuu){
	generate_bubble_graphs(graph_matrix, max, mpi_rank);
}

if(ami_parameters.TYPE_==AmiBase::Greens || ami_parameters.TYPE_==AmiBase::density || ami_parameters.TYPE_==AmiBase::ENERGY || ami_parameters.TYPE_==AmiBase::DOS ){
	generate_greens_graphs(graph_matrix, max, mpi_rank);
}

if(ami_parameters.TYPE_==AmiBase::FORCE){
	generate_force_graphs(graph_matrix,max,mpi_rank);
}
	
}

/* 
void AmiGraph::bubble_label_graphs_sys(std::vector< std::vector< AmiGraph::graph_t>> &graph_matrix, int min, int max){
	
int top=max+1;

//TODO: min, max maybe?

for(int i=min; i<top; i++){
	std::cout<<"Labelling order "<<i<<std::endl;

for(int m=0; m< graph_matrix[i].size(); m++){
number_vertices(graph_matrix[i][m]);
// systematic_vertex_label(graph_matrix[i][m]);
// if(m%10==0){	
std::cout<<"Labelling graph ord"<<i<<" num "<<m<<std::endl;
// }
bool success=true;
bubble_sys_label(graph_matrix[i][m], success);
// std::cout<<"Exited sys label function"<<std::endl;
if(!success){std::cerr<<"Failed to label graph ord"<<i<<" num "<<m<<std::endl;	
	
}


}	
		
}
}
 */


void AmiGraph::label_graphs_sys(std::vector< std::vector< AmiGraph::graph_t>> &graph_matrix, int min, int max){
	
int top=max+1;

//TODO: min, max maybe?
// std::cout<<"Graph matrix "<< graph_matrix.size()<<std::endl;
for(int i=min; i<top; i++){
	std::cout<<"Labelling order "<<i<<std::endl;
	// std::cout<<graph_matrix[i].size()<<std::endl;

for(int m=0; m< graph_matrix[i].size(); m++){
number_vertices(graph_matrix[i][m]);
// std::cout<<"Entering systematic_ver..."<<std::endl;


// this no longer needs to be enforced - since no labelling of closed graphs 
// if(ami_parameters.TYPE_!= AmiBase::density && ami_parameters.TYPE_ != AmiBase::doubleocc){
// systematic_vertex_label(graph_matrix[i][m]);
// }


if(m%10==0){	
std::cout<<"Labelling graph ord"<<i<<" num "<<m<<std::endl;
}
bool success=true;
// std::cout<<"Entering sys_label function"<<std::endl;
sys_label(graph_matrix[i][m], success);
// std::cout<<"Exited sys label function"<<std::endl;
if(!success){std::cerr<<"Warning: Failed to label systematically graph ord"<<i<<" num "<<m<<std::endl;
std::cerr<<"Falling back to repeated half-random labelling - does this graph have a tadpole?"<<std::endl;
systematic_vertex_label(graph_matrix[i][m]);	
repeated_labelling(graph_matrix[i][m], success);

if(!success){std::cerr<<"Warning: Random labelling failed for graph ord"<<i<<" num "<<m<<std::endl;}
if(success){std::cerr<<"Warning Resolved: Random labelling succeeded for graph ord"<<i<<" num "<<m<<std::endl;}
	
}

// print_all_edge_info(graph_matrix[i][m]);


}	
		
}


}

void AmiGraph::label_graphs(std::vector< std::vector< AmiGraph::graph_t>> &graph_matrix, int min, int max){
	
int top=max+1;
std::cout<<"Graph matrix "<< graph_matrix.size()<<std::endl;
for(int i=min; i<top; i++){
	// std::cout<<"Labelling order "<<i<<std::endl;
	// std::cout<<graph_matrix[i].size()<<std::endl;

for(int m=0; m< graph_matrix[i].size(); m++){
	
	number_vertices(graph_matrix[i][m]);



if(ami_parameters.TYPE_!= AmiBase::density && ami_parameters.TYPE_!= AmiBase::ENERGY && ami_parameters.TYPE_ != AmiBase::doubleocc){
// std::cout<<"Entering systematic_ver..."<<std::endl;
systematic_vertex_label(graph_matrix[i][m]);
}
// systematic_vertex_label(graph_matrix[i][m]);
// std::cout<<"Labelling graph ord"<<i<<" num "<<m<<std::endl;
// if(m%10==0){	
std::cout<<"Labelling graph ord"<<i<<" num "<<m<<std::endl;
std::cout<<"Prelabel"<<std::endl;
print_all_edge_info(graph_matrix[i][m]);
// }
bool success=true;

if(graph_type==AmiBase::density || graph_type==AmiBase::Greens|| graph_type==AmiBase::DOS|| graph_type==AmiBase::ENERGY || graph_type==AmiBase::FORCE){

// NOT SURE why we wouldn't just try the systematic case anyways. 
	// std::cout<<"In Sys label"<<std::endl;
sys_label(graph_matrix[i][m], success);
// std::cout<<"Exit Sys label"<<std::endl;
if(!success){ repeated_labelling(graph_matrix[i][m], success);}

}else{
repeated_labelling(graph_matrix[i][m], success);
}
if(!success){std::cerr<<"Failed to label graph ord"<<i<<" num "<<m<<std::endl;	
	
}


}	
		
}
}



void AmiGraph::reduce_gm_oneleg(std::vector< std::vector< AmiGraph::graph_t>> &graph_matrix, int mpi_rank){
	if(mpi_rank==0){
	std::cout<<"---------"<<std::endl<<"Filtering one-leg diagrams out"<<std::endl<<"---------"<<std::endl;	}	
  
  // label_graphs_sys(graph_matrix,0, graph_matrix.size());
  
for(int i=0; i<graph_matrix.size(); i++){

for(int m=0; m< graph_matrix[i].size(); m++){
	
// std::cout<<"Principle line length is "<<	get_length_principle_line(graph_matrix[i][m])<<std::endl;
if(is_oneleg(graph_matrix[i][m])){
graph_matrix[i].erase( graph_matrix[i].begin()+m);	
m--;
}


}	
if(mpi_rank==0){
	std::cout<<"Finished order "<< i <<std::endl;
std::cout<<"Total number of diagrams of order "<<i<<" is "<< graph_matrix[i].size()<<std::endl;	
}
}		
	
}

void AmiGraph::reduce_gm_rf(std::vector< std::vector< AmiGraph::graph_t>> &graph_matrix, int mpi_rank){
	if(mpi_rank==0){
	std::cout<<"---------"<<std::endl<<"Applying RF filter"<<std::endl<<"---------"<<std::endl;	}
for(int i=2; i<graph_matrix.size(); i++){

for(int m=0; m< graph_matrix[i].size(); m++){
	
	
bool success=true;
success=residue_filter(graph_matrix[i][m]);
if(success==false){
// std::cout<<"Rf returned false for graph ord "<<i<<" num "<<m<<std::endl;	
graph_matrix[i].erase( graph_matrix[i].begin()+m);	
m--;
}


}	
if(mpi_rank==0){
	std::cout<<"Finished order "<< i <<std::endl;
std::cout<<"Total number of diagrams of order "<<i<<" is "<< graph_matrix[i].size()<<std::endl;	}
	
}	



	
	
}

void AmiGraph::reduce_gm_1PBose(std::vector< std::vector< AmiGraph::graph_t>> &graph_matrix, int mpi_rank, int min){
	if(mpi_rank==0){
	std::cout<<"---------"<<std::endl<<"Applying 1P Bose filter"<<std::endl<<"---------"<<std::endl;	}
for(int i=min; i<graph_matrix.size(); i++){

for(int m=0; m< graph_matrix[i].size(); m++){
	
	
bool success=true;
success=!is_1P_bose_reducible_graph(graph_matrix[i][m]);
if(success==false){
std::cout<<"1p bose returned false for graph ord "<<i<<" num "<<m<<std::endl;	
graph_matrix[i].erase( graph_matrix[i].begin()+m);	
m--;
}


}	
if(mpi_rank==0){
	std::cout<<"Finished order "<< i <<std::endl;
std::cout<<"Total number of diagrams of order "<<i<<" is "<< graph_matrix[i].size()<<std::endl;	}
	
}	



	
	
}


void AmiGraph::reduce_gm_hubbard(std::vector< std::vector< AmiGraph::graph_t>> &graph_matrix, int mpi_rank, int min){
	if(mpi_rank==0){
	std::cout<<"---------"<<std::endl<<"Applying Hubbard filter"<<std::endl<<"---------"<<std::endl;	}
for(int i=min; i<graph_matrix.size(); i++){

for(int m=0; m< graph_matrix[i].size(); m++){
	
	
bool success=true;
success=is_hubbard(graph_matrix[i][m]);
if(success==false){
// std::cout<<"is_hubbard returned false for graph ord "<<i<<" num "<<m<<std::endl;	
graph_matrix[i].erase( graph_matrix[i].begin()+m);	
m--;
}


}	
if(mpi_rank==0){
	std::cout<<"Finished order "<< i <<std::endl;
std::cout<<"Total number of diagrams of order "<<i<<" is "<< graph_matrix[i].size()<<std::endl;	}
	
}	



	
	
}



void AmiGraph::reduce_gm_1PFermi(std::vector< std::vector< AmiGraph::graph_t>> &graph_matrix, int mpi_rank, int min){
	if(mpi_rank==0){
	std::cout<<"---------"<<std::endl<<"Applying 1P Fermi filter"<<std::endl<<"---------"<<std::endl;	}
for(int i=min; i<graph_matrix.size(); i++){

for(int m=0; m< graph_matrix[i].size(); m++){
	
	
bool success=true;
success=!is_1P_reducible_graph(graph_matrix[i][m]);
if(success==false){
std::cout<<"1p Fermi returned false for graph ord "<<i<<" num "<<m<<std::endl;	
graph_matrix[i].erase( graph_matrix[i].begin()+m);	
m--;
}


}	
if(mpi_rank==0){
	std::cout<<"Finished order "<< i <<std::endl;
std::cout<<"Total number of diagrams of order "<<i<<" is "<< graph_matrix[i].size()<<std::endl;	}
	
}	



	
	
}




void AmiGraph::reduce_gm_ppvertex(std::vector< std::vector< AmiGraph::graph_t>> &graph_matrix, int mpi_rank, int min){
	if(mpi_rank==0){
	std::cout<<"---------"<<std::endl<<"Applying vertex filter"<<std::endl<<"---------"<<std::endl;	}
for(int i=min; i<graph_matrix.size(); i++){

for(int m=0; m< graph_matrix[i].size(); m++){
	
std::cout<<"On graph "<<i<<" "<<m<<std::endl;	
bool success=true;
success=is_ppvertex_graph(graph_matrix[i][m]);
if(!success){
std::cout<<"vertex filter returned false for graph ord "<<i<<" num "<<m<<std::endl;	
graph_matrix[i].erase( graph_matrix[i].begin()+m);	
m--;
}


}	
if(mpi_rank==0){
	std::cout<<"Finished order "<< i <<std::endl;
std::cout<<"Total number of diagrams of order "<<i<<" is "<< graph_matrix[i].size()<<std::endl;	}
	
}	



	
	
}




void AmiGraph::reduce_gm_RPA_chains(std::vector< std::vector< AmiGraph::graph_t>> &graph_matrix, int mpi_rank, int min){
	if(mpi_rank==0){
	std::cout<<"---------"<<std::endl<<"Applying 1P Bose filter"<<std::endl<<"---------"<<std::endl;	}
for(int i=min; i<graph_matrix.size(); i++){

for(int m=0; m< graph_matrix[i].size(); m++){
	
int ord=graph_order(graph_matrix[i][m]);
int num_bubbles=count_bubbles(	graph_matrix[i][m]);


bool success=true;

if(ord>0 && num_bubbles==ord+1){success=false;}

// success=!is_RPA_reducible_graph(graph_matrix[i][m]);
if(success==false){
std::cout<<"Removing RPA chain for graph ord "<<i<<" num "<<m<<std::endl;	
graph_matrix[i].erase( graph_matrix[i].begin()+m);	
m--;
}


}	
if(mpi_rank==0){
	std::cout<<"Finished order "<< i <<std::endl;
std::cout<<"Total number of diagrams of order "<<i<<" is "<< graph_matrix[i].size()<<std::endl;	}
	
}	



	
	
}




void AmiGraph::reduce_gm_fock(std::vector< std::vector< AmiGraph::graph_t>> &graph_matrix, int mpi_rank, int min){
	if(mpi_rank==0){
	std::cout<<"---------"<<std::endl<<"Applying Fock filter"<<std::endl<<"---------"<<std::endl;	}
for(int i=min; i<graph_matrix.size(); i++){

for(int m=0; m< graph_matrix[i].size(); m++){
	
	
// bool success=true;
// success=is_1P_bose_reducible_graph(graph_matrix[i][m]);
if(has_fock_insertion(graph_matrix[i][m])){
// std::cout<<"Rf returned false for graph ord "<<i<<" num "<<m<<std::endl;	
graph_matrix[i].erase( graph_matrix[i].begin()+m);	
m--;
}


}	
if(mpi_rank==0){
	std::cout<<"Finished order "<< i <<std::endl;
std::cout<<"Total number of diagrams of order "<<i<<" is "<< graph_matrix[i].size()<<std::endl;	}
	
}	



	
	
}


// Requires labelled diagram 
void AmiGraph::reduce_gm_nested(std::vector< std::vector< AmiGraph::graph_t>> &graph_matrix, int mpi_rank, int min){
	if(mpi_rank==0){
	std::cout<<"---------"<<std::endl<<"Applying Nested filter"<<std::endl<<"---------"<<std::endl;	}
for(int i=min; i<graph_matrix.size(); i++){

for(int m=0; m< graph_matrix[i].size(); m++){
	
// if(i==4 && m==6){	
bool success=true;
success=is_nested(graph_matrix[i][m]);//is_1P_bose_reducible_graph(graph_matrix[i][m]);
// std::cout<<"Is nested function returned "<<success<<std::endl;
if(success){
// std::cout<<"Rf returned false for graph ord "<<i<<" num "<<m<<std::endl;	
graph_matrix[i].erase( graph_matrix[i].begin()+m);	
m--;
}
}

	
if(mpi_rank==0){
	std::cout<<"Finished order "<< i <<std::endl;
std::cout<<"Total number of diagrams of order "<<i<<" is "<< graph_matrix[i].size()<<std::endl;	}
	
}	



	
	
}



void AmiGraph::reduce_gm_ladder(std::vector< std::vector< AmiGraph::graph_t>> &graph_matrix, int mpi_rank, int min){
	if(mpi_rank==0){
	std::cout<<"---------"<<std::endl<<"Applying Ladder filter"<<std::endl<<"---------"<<std::endl;	}
for(int i=min; i<graph_matrix.size(); i++){

for(int m=0; m< graph_matrix[i].size(); m++){
	
	
bool success=true;
success=!is_ladder_graph(graph_matrix[i][m]);
if(success==false){
// std::cout<<"Rf returned false for graph ord "<<i<<" num "<<m<<std::endl;	
graph_matrix[i].erase( graph_matrix[i].begin()+m);	
m--;
}


}	
if(mpi_rank==0){
	std::cout<<"Finished order "<< i <<std::endl;
std::cout<<"Total number of diagrams of order "<<i<<" is "<< graph_matrix[i].size()<<std::endl;	}
	
}	



	
	
}



void AmiGraph::reduce_gm_skeleton(std::vector< std::vector< AmiGraph::graph_t>> &graph_matrix, int mpi_rank){
	if(mpi_rank==0){
	std::cout<<"---------"<<std::endl<<"Applying Skeleton filter"<<std::endl<<"---------"<<std::endl;	}
for(int i=2; i<graph_matrix.size(); i++){

for(int m=0; m< graph_matrix[i].size(); m++){
	
	
bool success=true;
success=is_skeleton(graph_matrix[i][m]);
if(success==false){
// std::cout<<"Skeleton returned false for graph ord "<<i<<" num "<<m<<std::endl;	
graph_matrix[i].erase( graph_matrix[i].begin()+m);	
m--;
}


}	
if(mpi_rank==0){
	std::cout<<"Finished order "<< i <<std::endl;
std::cout<<"Total number of diagrams of order "<<i<<" is "<< graph_matrix[i].size()<<std::endl;	}
	
}	



	
	
}






void AmiGraph::reduce_gm_tp(std::vector< std::vector< AmiGraph::graph_t>> &graph_matrix,  int mpi_rank, int min){
	if(mpi_rank==0){
	std::cout<<"---------"<<std::endl<<"Filtering tadpoles out"<<std::endl<<"---------"<<std::endl;	}
for(int i=min; i<graph_matrix.size(); i++){

for(int m=0; m< graph_matrix[i].size(); m++){
	
// std::cout<<"Found "<<count_tadpoles(graph_matrix[i][m])<<" on "<<i<<" "<<m<<std::endl;

if(count_tadpoles(graph_matrix[i][m])>0){
// std::cout<<"Rf returned false for graph ord "<<i<<" num "<<m<<std::endl;	
graph_matrix[i].erase( graph_matrix[i].begin()+m);	
m--;
}


}	
if(mpi_rank==0){
	std::cout<<"Finished tp order "<< i <<std::endl;
std::cout<<"Total number of diagrams with no tadpoles of order "<<i<<" is "<< graph_matrix[i].size()<<std::endl;}	
	
}	

	
}


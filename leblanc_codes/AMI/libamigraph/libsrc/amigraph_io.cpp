
#include "amigraph.hpp"

void AmiGraph::print_edgeson_vert(int index){
	
graph_t g=current_graph;
vertex_t vert;

boost::graph_traits<graph_t>::vertex_iterator ei, edge_end;

for (boost::tie(ei,edge_end) = vertices(g); ei != edge_end; ++ei){
if( g[*ei].index_ == index){vert=*ei;}
   }

boost::graph_traits<graph_t>::in_edge_iterator iei, iedge_end;
boost::graph_traits<graph_t>::out_edge_iterator oei, oedge_end;	
std::cout<<"In on "<< index<<std::endl;
for (boost::tie(iei,iedge_end) = in_edges(vert, g); iei != iedge_end; ++iei){	
	
for( int i=0; i<g[*iei].g_struct_.alpha_.size();i++){
	
std::cout<<	g[*iei].g_struct_.alpha_[i]<<" ";
	
}

std::cout<<std::endl;	
}
std::cout<<"out on "<< index<<std::endl;
for (boost::tie(oei,oedge_end) = out_edges(vert, g); oei != oedge_end; ++oei){	
for( int i=0; i<g[*oei].g_struct_.alpha_.size();i++){
	
std::cout<<	g[*oei].g_struct_.alpha_[i]<<" ";
	
}

std::cout<<std::endl;	
	
}



}

void AmiGraph::print_all_vertex_info(){

graph_t g=current_graph;

boost::graph_traits<graph_t>::vertex_iterator ei, edge_end;
//set_all_vertex_freq(g);

//
std::cout<<"Current graph Contains the following vertices and info"<<std::endl;
std::cout<<num_vertices(g)<<std::endl;
for (boost::tie(ei,edge_end) = vertices(g); ei != edge_end; ++ei){
std::cout << "Vertex " << g[*ei].index_ << " of vertex type "<< g[*ei].type_<<" with degree "<< in_degree(*ei,g)+out_degree(*ei,g) << std::endl;
   }



}


	
	
	
}

std::string AmiGraph::print_edge(edge_t &e, graph_t &g){
	
	std::string output;
	output="("+std::to_string(g[source(e,g)].index_)+","+std::to_string(g[target(e,g)].index_)+")";
	//output="("<<g[source(e,g)]<<","<<g[target(e,g)]<<")";
	return output;
}

void AmiGraph::print_edge_info(edge_t &e, graph_t &g){
	

std::cout << "Edge (" << g[source(e,g)].index_ << ","<< g[target(e,g)].index_<<") of stat_type "<< g[e].g_struct_.stat_ <<" with loop id "<< g[e].fermi_loop_id<<" with spin "<<g[e].spin<<" with label "<< g[e].label<< " label=[ ";// << std::endl;

for( int i=0; i<g[e].g_struct_.alpha_.size();i++){
	
std::cout<<	g[e].g_struct_.alpha_[i]<<" ";
	
}
std::cout<<"] and eps=[";

for(int i=0; i<g[e].g_struct_.eps_.size();i++){
	
std::cout<<g[e].g_struct_.eps_[i]<<" ";	
	
}

std::cout<<"]"<<std::endl;
	
	
	
}


void AmiGraph::print_all_edge_info(graph_t &g){

boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
// Looking through all edges to find source and targets : this could be improved?
for (boost::tie(ei,ei_end) = edges(g); ei != ei_end; ++ei){


std::cout << "Edge (" << g[source(*ei,g)].index_ << ","<< g[target(*ei,g)].index_<<") with bm "<<g[*ei].band_element.first<<","<<g[*ei].band_element.second<<" of stat_type "<< g[*ei].g_struct_.stat_ <<" with loop id "<< g[*ei].fermi_loop_id<<" with spin "<<g[*ei].spin<<" with label "<< g[*ei].label<< " label=[ ";// << std::endl;

for( int i=0; i<g[*ei].g_struct_.alpha_.size();i++){
	
std::cout<<	g[*ei].g_struct_.alpha_[i]<<" ";
	
}
std::cout<<"] and eps=[";

for(int i=0; i<g[*ei].g_struct_.eps_.size();i++){
	
std::cout<<g[*ei].g_struct_.eps_[i]<<" ";	
	
}

std::cout<<"]"<<std::endl;

//std::cout<<"edge "<<jg[*ei].edge_number_<<" loop Total is "<< total<<std::endl;
}

}


void AmiGraph::print_all_edge_info(){

graph_t g=current_graph;

boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
// Looking through all edges to find source and targets : this could be improved?
for (boost::tie(ei,ei_end) = edges(g); ei != ei_end; ++ei){


std::cout << "Edge (" << g[source(*ei,g)].index_ << ","<< g[target(*ei,g)].index_<<") of stat_type "<< g[*ei].g_struct_.stat_ <<" with loop id "<< g[*ei].fermi_loop_id<<" with spin "<<g[*ei].spin<<" with label "<< g[*ei].label<< " label=[ ";// << std::endl;

for( int i=0; i<g[*ei].g_struct_.alpha_.size();i++){
	
std::cout<<	g[*ei].g_struct_.alpha_[i]<<" ";
	
}
std::cout<<"] and eps=[";

for(int i=0; i<g[*ei].g_struct_.eps_.size();i++){
	
std::cout<<g[*ei].g_struct_.eps_[i]<<" ";	
	
}

std::cout<<"]"<<std::endl;


//std::cout<<"edge "<<jg[*ei].edge_number_<<" loop Total is "<< total<<std::endl;
}


/////////////////                         ////////////////////////////
////////////////          RAYAN's Code   /////////////////////////////



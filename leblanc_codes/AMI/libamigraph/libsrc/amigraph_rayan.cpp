#include "amigraph.hpp"

void AmiGraph::Rfind_bosonic_edges(graph_t &g, edge_vector_t &vector){

vector.clear();

//std::cout<< "Finding bosonic edges" <<std::endl;
boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
// Looking through all edges to find source and targets : this could be improved?
for (boost::tie(ei,ei_end) = edges(g); ei != ei_end; ++ei){

if( g[*ei].g_struct_.stat_==AmiBase::Bose){
vector.push_back(*ei);
}
//std::cout << "Edge " << *ei << " of stat_type "<< g[*ei].g_struct_.stat_ << std::endl;

//std::cout<<"edyge "<<jg[*ei].edge_number_<<" loop Total is "<< total<<std::endl;
}



}

void AmiGraph::Rfind_fermionic_edges(graph_t &g, edge_vector_t &vector){

vector.clear();

//std::cout<< "Finding Fermionic edges" <<std::endl;
boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
// Looking through all edges to find source and targets : this could be improved?
for (boost::tie(ei,ei_end) = edges(g); ei != ei_end; ++ei){

if( g[*ei].g_struct_.stat_==AmiBase::Fermi){
vector.push_back(*ei);
}
//std::cout << "Edge " << *ei << " of stat_type "<< g[*ei].g_struct_.stat_ << std::endl;

//std::cout<<"edge "<<jg[*ei].edge_number_<<" loop Total is "<< total<<std::endl;
}




}


#include "amigraph.hpp"

// This is all for multi-band stuff 



// this will assign each g_struct a species entry 
void AmiGraph::mb_setup(graph_t &g, std::vector< std::vector< int >> &V_prod){
  V_prod.clear();
  if(graph_type != AmiBase::Sigma){throw std::runtime_error("Multi-band only set up for sigma graphs - exiting.");}
  
  edge_vector_t ev, external_ev;
  vertex_vector_t vv;
  
  find_internal_fermionic_edges(g,ev);
  find_external_vertices(g, vv, external_ev);
  
  for(int i=0; i< ev.size(); i++){
    g[ev[i]].g_struct_.species_=i;
  }
  
  for(int i=0; i< external_ev.size(); i++){
    g[ev[i]].g_struct_.species_=i+int(ev.size());
  }
  
  edge_vector_t bev;
  find_bosonic_edges(g,bev);
  
  get_vprod(g, bev,V_prod);
  
  
  
}


void AmiGraph::get_vprod(graph_t &g, edge_vector_t &bev,std::vector< std::vector< int >> &V_prod ){
  
 for(int i=0; i< bev.size(); i++){
   V_prod.push_back(get_v(g,bev[i]));
 }
  
  
}

std::vector< int > AmiGraph::get_v(graph_t &g, edge_t &be){
  
  std::vector<int> output;
  output.resize(4,-1);
  
  vertex_t v1,v2;
  
  v1=source(be,g);
  v2=target(be,g);
  
boost::graph_traits<graph_t>::edge_iterator ei, ei_end;

int found=0;
for(boost::tie(ei,ei_end)=edges(g); ei!=ei_end; ++ei){
if(source(*ei,g)==v1){
output[0]=g[*ei].g_struct_.species_;
found++;
}
if(target(*ei,g)==v1){
output[1]=g[*ei].g_struct_.species_;
found++;
}
if(source(*ei,g)==v2){
output[2]=g[*ei].g_struct_.species_;
found++;
}
if(target(*ei,g)==v2){
output[3]=g[*ei].g_struct_.species_;
found++;
}


}	

if(found!=4){throw std::runtime_error("Didn't find one of the V species indexes - exiting.");}
  
return output;  
  
  
}

// this is the most likely place to cause a segfault.
std::complex< double > AmiGraph::eval_symbolic_Vprod(std::vector< int > &V_values, std::vector< std::vector< int >> &V_indices,AmiGraph::U_matrix_type &UM){
  
  std::complex<double> output(1,0);
  
  for(int i=0; i< V_indices.size(); i++){
    
    output=output*UM[V_values[V_indices[i][0]]][V_values[V_indices[i][1]]][V_values[V_indices[i][2]]][V_values[V_indices[i][3]]];
        
  }
  
  return output;
  
  
  
}

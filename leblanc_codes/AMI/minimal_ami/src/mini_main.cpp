#include "mini_ami.hpp"


int main(int argc, char** argv)
{ 


mband::params_param params;
params_loader("../loader/params.txt", params);
    
     std::vector<std::vector<int>> interaction = readFile("../loader/Hubbard_U.txt");
	 std::vector<double> interaction_value = readFile1("../loader/Hubbard_U.txt",5);
	 std::vector<double> band_energy = readFile1("../loader/ccpdvz_h.txt",2);
	 std::cout<<interaction.size();
    

	
	
////////////////////////////James loading  stuff/////////////////////////////////
        AmiBase ami;
		int seed=3;	
		std::cout<< "Constructing AmiGraph object using seed: "<<seed <<" " <<std::endl;

	
	AmiGraph gp(AmiBase::Pi_ppud, seed);	
	std::string infile("ext_vars.dat");
	
		std::cout<<"Reading external parameters from ext_vars.dat"<<std::endl;
		NewAmiCalc::external_variable_list extern_list;

	// read column data into a file - you can get external variables however you would like	
	gp.ami.read_external(infile, extern_list);	
	



	AmiGraph::gg_matrix_t ggm; 
	std::cout<<"Attempting to load self-energy graphs from example_graphs"<<std::endl;
	int max=4;
	//g.read_ggmp("../example_graphs/ggm_sigma_nofock_notp/",ggm, max);
	gp.read_ggmp("../graphs/ggm_ppvertex/",ggm, max);
	std::cout<<"Completed read"<<std::endl;
	std::cout<<std::endl;
	
 
	
	gp.ggm_label(ggm,0);
	AmiGraph::graph_t graph = ggm[4][0].graph_vec[0];

	
/////////////////////////////////////loading and labelling ends////////////////////
//std::vector<mband::sampler_collector> sigma_ToSum;
//std::vector<AmiGraph::graph_t>sigma_FromGraph;

std::cout<<"External parameters read are"<<std::endl;
for(int i=0; i<extern_list.size();i++){
		std::cout<<extern_list[i].BETA_<<" "<<extern_list[i].MU_<<" "<< extern_list[i].H_<<" "<<extern_list[i].KDIM_<<" "<<extern_list[i].external_k_list_[0][0]<<" "<<extern_list[i].external_k_list_[0][1] <<" "<<extern_list[i].external_freq_[0]<<std::endl;
	
		}	



///////////////////////////////constructing mband ///////////////////////////////////////
	
bool hf = true; ///set to true if you are expecting a hatree or fock type interaction. True works for all type of graph so its set to default. Setting to false 
	//speeds up the process
	//construction 
mband mb(interaction,interaction_value,band_energy,hf);
int  min_ord = params.min_ord;		
int  max_ord = params.max_ord;

/*
AmiGraph::edge_vector_t fermionic_edge;


gp.find_non_external_bosonic_edges(graph, bosonic_edge);

std::cout <<"Printing the internal bosonic_edges \n"; 
for (auto edge  : bosonic_edge){

gp.print_edge_info(edge,graph);
}

AmiGraph::edge_vector_t bosonic_ext_edge;
AmiGraph::vertex_vector_t bosonic_vert;

gp.find_external_vertices(graph, bosonic_vert,bosonic_ext_edge);
std::cout << " extenral bosonic edges are " <<std::endl;
for (auto edge  : bosonic_ext_edge){
gp.print_edge_info(edge,graph);
}



std::cout <<" All Fermionic edges are \n" ;
boost::graph_traits<AmiGraph::graph_t>::edge_iterator  ei,ei_end;
for (boost::tie(ei,ei_end) = edges(graph); ei !=ei_end; ++ei){
	if(graph[*ei].g_struct_.stat_ == AmiBase::Fermi){
	AmiGraph::edge_t e = *ei;
		gp.print_edge_info(e, graph);
		
	}

}

*/

/*
AmiGraph::edge_vector_t pair1;
AmiGraph::edge_vector_t pair2;
std::vector<std::vector<int>> alpha_p1;
std::vector<std::vector<int>> alpha_p2;

std::vector<int> bandindex = {1,2};
mb.find_connected_pp_edge(graph,pair1,pair2);
mb.assign_bandindex_fourline_pp(graph,pair1,pair2, bandindex);

std::cout << "Assigning species are done, lets look at what we got \n";
AmiGraph::edge_vector_t fermionic_edges;

gp.find_fermionic_edges(graph,fermionic_edges);

for (auto e : fermionic_edges){
	gp.print_edge_info(e, graph);
	std::cout <<"species assigned are "  << graph[e].g_struct_.species_ <<std::endl;	
}




AmiGraph::edge_vector_t fermionic_edge1;
std::vector<std::vector<int>> fermionic_species;
std::vector<std::vector<std::vector<int>>> interaction_species;
std::vector<std::vector<int>> bosonic_Alpha;
std::vector<std::vector<int>> gkkp_Alpha;

std::cout << "Sampling starts here \n\n\nn";
mb.solve_pp_ord(graph,fermionic_edge1,fermionic_species,interaction_species,bosonic_Alpha,gkkp_Alpha,bandindex);



std::cout << "printing all edge info \n\n" ;
for (auto e: fermionic_edge1){
	gp.print_edge_info(e,graph);
		
}
std::cout << " Now corresponding fermionic edge species are \n\n" ;
print2d(fermionic_species);*/
/*
std::cout << "Next outgoing vertex is ";

start_vertex = find_next_ppv_outedge(start_vertex, graph);
start_vertex = find_next_ppv_outedge(start_vertex, graph);

std::cout << graph[start_vertex].index_ << std::endl;
*/
	
	
std::vector<int> bandindex = {1,2};	
std::vector<mband::sampler_collector> pp_ToSum;
std::vector<AmiGraph::graph_t> pp_FromGraph;
for (int i = min_ord; i < max_ord+1; ++i) {
    for (int j = 0; j < ggm[i].size(); ++j) {
        for (int k = 0; k < ggm[i][j].graph_vec.size(); ++k) {
            mband::sampler_collector pp_collector;
            mb.pp_sampler(ggm[i][j].graph_vec[k], pp_collector,bandindex);
			if (!pp_collector.fermionic_edge_species.empty()){
				if ( params.mfreq_indp == 0){
			pp_ToSum.push_back(pp_collector);
			pp_FromGraph.push_back(ggm[i][j].graph_vec[k]);
				}
			}
		}
	}
}












//mb.solve_multiband(graph,fermionic_edge,fermionic_species,interaction_species,bosonic_Alpha,ext_legs);



}
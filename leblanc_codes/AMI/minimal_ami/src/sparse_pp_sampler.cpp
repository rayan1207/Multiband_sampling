#include "mini_ami.hpp"

AmiGraph gp(AmiBase::Pi_ppud, 0);	

void mband::find_four_pp_lines(AmiGraph::graph_t &graph, AmiGraph::edge_vector_t &outgoing_fedge,
 AmiGraph::edge_vector_t &incoming_fedge){
	AmiGraph::vertex_t start;
	AmiGraph::vertex_t end;
	gp.find_start_end_vertices( graph, start,end);
	//std::cout << " Start vertices are [" <<graph[start].index_ << "]"<< std::endl;
	//std::cout << " End vertices are [" <<graph[end].index_ << "]"<< std::endl;
	//std::cout<< " Outcoming edges of first are " <<std::endl;    
	auto range = boost::out_edges(start, graph);
		for (auto it = range.first; it != range.second; ++it) {
			if(graph[*it].g_struct_.stat_ == AmiBase::Fermi){
			AmiGraph::edge_t e = *it;
			outgoing_fedge.push_back(e);
			//gp.print_edge_info(e, graph);	
		}	 
	}


	//std::cout<< " Incoming edges of first are " <<std::endl;
	auto range1 = boost::in_edges(end,graph);
		for (auto it = range1.first; it != range1.second; ++it) {
			if(graph[*it].g_struct_.stat_ == AmiBase::Fermi){
			AmiGraph::edge_t e = *it;
			incoming_fedge.push_back(e);
			//gp.print_edge_info(e, graph);	
		}	 
	}

}

void mband::check_iso_pp(AmiGraph::graph_t g, AmiGraph::gg_matrix_t ggm, int min_ord, int max_ord) {
    bool match_found = false; // Add a flag to track if a match is found

    for (int i = min_ord; i < max_ord + 1; ++i) {
        for (int j = 0; j < ggm[i].size(); ++j) {
            for (int k = 0; k < ggm[i][j].graph_vec.size(); ++k) {
                if (gp.spin_iso(g, ggm[i][j].graph_vec[k])) {
                    std::cout << "matches-> o" << std::to_string(i) << "_g" << std::to_string(j)
                              << "_n" << std::to_string(k) << std::endl;
                    match_found = true;
                    break; 
                }
            }
            if (match_found) {
                break; 
            }
        }
        if (match_found) {
            break; 
        }
    }

    if (!match_found) {
        std::cout << "no match found" << std::endl;
    }
}


AmiGraph::vertex_t mband::find_next_ppv_outedge(AmiGraph::vertex_t &start, AmiGraph::graph_t &graph){ 
AmiGraph::edge_t e;
auto range = boost::out_edges(start, graph);
	for (auto it = range.first; it != range.second; ++it) {
		if(graph[*it].g_struct_.stat_ == AmiBase::Fermi){
		 e = *it;
	}	 
}
AmiGraph::vertex_t v = boost::target(e,graph);
return v;

}


void mband::find_connected_pp_edge(AmiGraph::graph_t &graph, AmiGraph::edge_vector_t &pair1, AmiGraph::edge_vector_t &pair2){
AmiGraph::edge_vector_t outgoing_fedge;
AmiGraph::edge_vector_t incoming_fedge;
mband::find_four_pp_lines(graph,outgoing_fedge,incoming_fedge);

AmiGraph::vertex_t start;
AmiGraph::vertex_t end;
gp.find_start_end_vertices( graph, start,end);

bool matchlook = true;
AmiGraph::vertex_t start_vertex = target(outgoing_fedge[0], graph);
std::cout << "The start vertex is " << graph[start_vertex].index_ << std::endl;

AmiGraph::vertex_t target_vertex = source(incoming_fedge[0], graph);
std::cout << "The target vertex is " << graph[source(incoming_fedge[0], graph)].index_ << std::endl;	

if (outgoing_fedge[0] == incoming_fedge[0]){
		std::cout << "the first edges are same, no vertices in between " << std::endl;
		pair1 = {outgoing_fedge[0],incoming_fedge[0]};
		pair2 = {outgoing_fedge[1],incoming_fedge[1]};
		matchlook = false;	
}

else if (outgoing_fedge[0] == incoming_fedge[1]){
		std::cout << "the first edges are same, no vertices in between " << std::endl;
		pair1 = {outgoing_fedge[0],incoming_fedge[1]};
		pair2 = {outgoing_fedge[1],incoming_fedge[0]};
		matchlook = false;	
	}
	


else {	
do {

	if (start_vertex == target_vertex){
		
		std::cout << "destination reached, match found" << std::endl;
		pair1 = {outgoing_fedge[0],incoming_fedge[0]};
		pair2 = {outgoing_fedge[1],incoming_fedge[1]};
		matchlook = false;
		
	}
	if (start_vertex == end){
		std::cout << "no match found";
		pair1 = {outgoing_fedge[0],incoming_fedge[1]};
		pair2 = {outgoing_fedge[1],incoming_fedge[0]};
		matchlook = false;
		
	}
	else {
		
		start_vertex= mband::find_next_ppv_outedge(start_vertex, graph);
		
	}
		
		
	}while(matchlook);
}	
/*
	for (auto i : pair1){
		alpha_p1.push_back(graph[i].g_struct_.alpha_);
	}
	for (auto i : pair2){
		alpha_p2.push_back(graph[i].g_struct_.alpha_);
	}
*/	
	std::cout << " set of pair1 are " <<std::endl;
	
	for  (auto vec: pair1){
		
 	gp.print_edge_info(vec, graph);
		
	}
	//std::cout << " set of pair1 alpha  are " <<std::endl;
	//print2d(alpha_p1);
	
	std::cout << " set of pair2 are " <<std::endl;
	
	for  (auto vec: pair2){
		
 	gp.print_edge_info(vec, graph);
			
	}
	
	//std::cout << " set of pair2 alpha  are " <<std::endl;
	//print2d(alpha_p2);
		
}	




void mband::assign_bandindex_fourline_pp(AmiGraph::graph_t &graph,AmiGraph::edge_vector_t &pair1,AmiGraph::edge_vector_t &pair2, std::vector<int> &bandindex){
std::vector<int> band1;
std::vector<int> band2;
/*
if (bandindex[0] <10 && bandindex[1]<10) {
	for (auto edge : pair1){
		int spin = graph[edge].spin;
		graph[edge].g_struct_.species_ = bandindex[spin];	
	}
	for (auto edge : pair2){
		int spin = graph[edge].spin; 
		graph[edge].g_struct_.species_ = bandindex[spin];	
	}
}
*/
	int a1 = bandindex[0]/10;
	int a2 = bandindex[0]%10;
	int b1 = bandindex[1]/10;
	int b2 = bandindex[1]%10;
	band1 ={a1,b1};
	band2 ={a2,b2};
	print1d(band1);
	print1d(band2);
	
   if (graph[pair1[0]].spin==0){
	   graph[pair1[0]].g_struct_.species_ = band1[0];
	   graph[pair1[1]].g_struct_.species_ = band1[1];
	   graph[pair2[0]].g_struct_.species_ = band2[0];
	   graph[pair2[1]].g_struct_.species_ = band2[1];  
   }
   
   else {
	   graph[pair1[0]].g_struct_.species_ = band2[0];
	   graph[pair1[1]].g_struct_.species_ = band2[1];
	   graph[pair2[0]].g_struct_.species_ = band1[0];
	   graph[pair2[1]].g_struct_.species_ = band1[1];  
	   
		}
	
	}



void mband::solve_pp_ord1(AmiGraph::graph_t &graph,AmiGraph::edge_vector_t &fermionic_edge,
std::vector<std::vector<int>> &fermionic_species,std::vector<std::vector<std::vector<int>>> &interaction_species,std::vector<std::vector<int>> &bosonic_Alpha,std::vector<std::vector<int>> &gkkp_Alpha,std::vector<int> &bandindex){
    AmiGraph::edge_vector_t pair1;
	AmiGraph::edge_vector_t pair2;
	mband::find_connected_pp_edge(graph,pair1,pair2);
	mband::assign_bandindex_fourline_pp(graph,pair1,pair2, bandindex);
	gp.find_internal_fermionic_edges(graph,fermionic_edge);
	std::vector<int> orginal_species = mband::generate_edge_species(graph, fermionic_edge);
	

	AmiGraph::edge_vector_t bvector;
	std::vector<AmiGraph::edge_vector_t> int_vector;
	gp.find_non_external_bosonic_edges(graph,bvector);
	 std::cout <<"bosonic edges are \n";
	 for (auto b: bvector){
		gp.print_edge_info(b,graph);
		std::cout<<std::endl;
		bosonic_Alpha.push_back(graph[b].g_struct_.alpha_);
		
	}
	print2d(bosonic_Alpha);
  for (auto b: pair1){
		gkkp_Alpha.push_back(graph[b].g_struct_.alpha_);
	}	
	for (int i=0;i <bvector.size();i++){
		//g.print_edge_info(bvector[i],gself2);
		std::cout << "Edge = ("  << graph[source(bvector[i],graph)].index_<<"," << graph[target(bvector[i],graph)].index_ <<")" 
		<<std::endl;
	}
	for (auto b: bvector){
		bosonic_Alpha.push_back(graph[b].g_struct_.alpha_);
	}	
	mband::find_interaction(graph,bvector,int_vector);
	mband::print_interactions(graph,bvector, int_vector);
	std::vector<int> initial_species_1 = mband::generate_edge_species(graph, int_vector[0]);
	std::vector<std::vector<int>> possible_species_1;
	if (mband::Hartee_fock){
		possible_species_1 = findmatch(initial_species_1,int_vector[0]);
	}
	else{
		possible_species_1 = mband::interaction_legs;
	}
	std::cout<< " Possible external species in curly bracket and internal species in () is shown below \n";
	if (!possible_species_1.empty()) {
		for (int i = 0; i<possible_species_1.size();i++){

		mband::assign_label(graph,int_vector[0],possible_species_1[i]);
        interaction_species.push_back({possible_species_1[i]});
		std::vector<int> v;
		std::cout<<"(";
		for (int x =0; x < fermionic_edge.size(); x++){
			v.push_back(graph[fermionic_edge[x]].g_struct_.species_);
			std::cout<<graph[fermionic_edge[x]].g_struct_.species_;				
		}
		fermionic_species.push_back(v);
		v.clear();
		std::cout<<")"<<std::endl;			    
		}
		
		//mband::reset_species(graph,int_vector);
		mband::assign_label(graph,fermionic_edge,orginal_species);	
		
	}
	else {std::cout <<"No match found" <<std::endl;}

mband::print_assigned_species(interaction_species);

  
}


void mband::solve_pp_ord2(AmiGraph::graph_t &graph,AmiGraph::edge_vector_t &fermionic_edge,
std::vector<std::vector<int>> &fermionic_species,std::vector<std::vector<std::vector<int>>> &interaction_species,std::vector<std::vector<int>> &bosonic_Alpha,std::vector<std::vector<int>> &gkkp_Alpha,std::vector<int> &bandindex){
    ///do pp stuff here first///
	AmiGraph::edge_vector_t pair1;
	AmiGraph::edge_vector_t pair2;
	mband::find_connected_pp_edge(graph,pair1,pair2);
	mband::assign_bandindex_fourline_pp(graph,pair1,pair2, bandindex);
	gp.find_internal_fermionic_edges(graph,fermionic_edge);
	std::vector<int> orginal_species = mband::generate_edge_species(graph, fermionic_edge);
	

	AmiGraph::edge_vector_t bvector;
	std::vector<AmiGraph::edge_vector_t> int_vector;
	gp.find_non_external_bosonic_edges(graph,bvector);
	
	 std::cout <<"bosonic edges are \n";
	 for (auto b: bvector){
		gp.print_edge_info(b,graph);
		std::cout<<std::endl;
		bosonic_Alpha.push_back(graph[b].g_struct_.alpha_);
	}
	print2d(bosonic_Alpha);
  for (auto b: pair1){
		gkkp_Alpha.push_back(graph[b].g_struct_.alpha_);
	}		
	///pp stuff ends here///
	for (int i=0;i <bvector.size();i++){
		//g.print_edge_info(bvector[i],gself2);
		std::cout << "Edge = ("  << graph[source(bvector[i],graph)].index_<<"," << graph[target(bvector[i],graph)].index_ <<")" 
		<<std::endl;
	}
	for (auto b: bvector){
		print1d(graph[b].g_struct_.alpha_);
	}
	mband::find_interaction(graph,bvector,int_vector);
	mband::print_interactions(graph,bvector, int_vector);
	std::vector<int> initial_species_1 = mband::generate_edge_species(graph, int_vector[0]);
	std::vector<std::vector<int>> possible_species_1;
	if (mband::Hartee_fock){
		possible_species_1 = findmatch(initial_species_1,int_vector[0]);
	}
	else{
		possible_species_1 = mband::interaction_legs;
	}
	std::cout<< " Possible external species in curly bracket and internal species in () is shown below \n";
	if (!possible_species_1.empty()) {
		for (int i = 0; i<possible_species_1.size();i++){

			mband::assign_label(graph,int_vector[0],possible_species_1[i]);

			std::vector<int> initial_species_2 = mband::generate_edge_species(graph, int_vector[1]);

			std::vector<std::vector<int>> possible_species_2 = mband::findmatch(initial_species_2,int_vector[1]);
			//end_code
			if (!possible_species_2.empty()){
				for (int j = 0; j < possible_species_2.size();j++){
				mband::assign_label(graph,int_vector[1],initial_species_2);				
				mband::assign_label(graph,int_vector[1],possible_species_2[j]);
				interaction_species.push_back({possible_species_1[i],possible_species_2[j]});
				std::cout<< fermionic_edge.size();
				std::vector<int> v;
				std::cout<<"(";
				for (int x =0; x < fermionic_edge.size(); x++){
					v.push_back(graph[fermionic_edge[x]].g_struct_.species_);
					std::cout<<graph[fermionic_edge[x]].g_struct_.species_;				
				}
				fermionic_species.push_back(v);
				v.clear();
				std::cout<<")"<<std::endl;			    
				}
				mband::assign_label(graph,int_vector[1],initial_species_2);
				
			}
			else {mband::assign_label(graph,int_vector[1],initial_species_2);}
			
		//mband::reset_species(graph,int_vector);
		mband::assign_label(graph,fermionic_edge,orginal_species);		
		}
	}
	else {std::cout <<"No match found" <<std::endl;}

mband::print_assigned_species(interaction_species);
}




void mband::solve_pp_ord3(AmiGraph::graph_t &graph,AmiGraph::edge_vector_t &fermionic_edge,
std::vector<std::vector<int>> &fermionic_species,std::vector<std::vector<std::vector<int>>> &interaction_species,std::vector<std::vector<int>> &bosonic_Alpha,std::vector<std::vector<int>> &gkkp_Alpha,std::vector<int> &bandindex){
      ///do pp stuff here first///
	AmiGraph::edge_vector_t pair1;
	AmiGraph::edge_vector_t pair2;
	mband::find_connected_pp_edge(graph,pair1,pair2);
	mband::assign_bandindex_fourline_pp(graph,pair1,pair2, bandindex);
	gp.find_internal_fermionic_edges(graph,fermionic_edge);
	std::vector<int> orginal_species = mband::generate_edge_species(graph, fermionic_edge);
	

	AmiGraph::edge_vector_t bvector;
	std::vector<AmiGraph::edge_vector_t> int_vector;
	gp.find_non_external_bosonic_edges(graph,bvector);
	 std::cout <<"bosonic edges are \n";
	 for (auto b: bvector){
		gp.print_edge_info(b,graph);
		std::cout<<std::endl;
		bosonic_Alpha.push_back(graph[b].g_struct_.alpha_);	
	}
	print2d(bosonic_Alpha);
	for (auto b: pair1){
		gkkp_Alpha.push_back(graph[b].g_struct_.alpha_);
	}		
	///pp stuff ends here///
	for (int i=0;i <bvector.size();i++){
		//g.print_edge_info(bvector[i],gself2);
		std::cout << "Edge = ("  << graph[source(bvector[i],graph)].index_<<"," << graph[target(bvector[i],graph)].index_ <<")" 
		<<std::endl;
	}
	mband::find_interaction(graph,bvector,int_vector);
	mband::print_interactions(graph,bvector, int_vector);
	std::vector<int> initial_species_1 = mband::generate_edge_species(graph, int_vector[0]);
	std::vector<std::vector<int>> possible_species_1;
	if (mband::Hartee_fock){
		possible_species_1 = findmatch(initial_species_1,int_vector[0]);
	}
	else{
		possible_species_1 = mband::interaction_legs;
	}
	std::cout<< " Possible external species in curly bracket and internal species in () is shown below \n";
	if (!possible_species_1.empty()) {
		for (int i = 0; i<possible_species_1.size();i++){
			std::vector<int> initial_species_1 = mband::generate_edge_species(graph, int_vector[0]);	
			mband::assign_label(graph,int_vector[0],possible_species_1[i]);

			std::vector<int> initial_species_2 = mband::generate_edge_species(graph, int_vector[1]);

			std::vector<std::vector<int>> possible_species_2 = mband::findmatch(initial_species_2,int_vector[1]);
			if (!possible_species_2.empty()){
				for (int j = 0; j < possible_species_2.size();j++){
				mband::assign_label(graph,int_vector[1],initial_species_2);				
				mband::assign_label(graph,int_vector[1],possible_species_2[j]);
				std::vector<int> initial_species_3 = mband::generate_edge_species(graph, int_vector[2]);
				std::vector<std::vector<int>> possible_species_3 = mband::findmatch(initial_species_3, int_vector[2]);
					if (!possible_species_3.empty()){
						for (int k = 0; k < possible_species_3.size();k++){
							mband::assign_label(graph,int_vector[2],initial_species_3);
							mband::assign_label(graph,int_vector[2],possible_species_3[k]);											
							interaction_species.push_back({possible_species_1[i],possible_species_2[j],possible_species_3[k]});
							
							std::cout<< fermionic_edge.size();
							std::cout<<j << " printing internal fermionic band index" << std::endl;
							std::cout<<"(";
							std::vector<int> v;
				for (int x =0; x < fermionic_edge.size(); x++){
					v.push_back(graph[fermionic_edge[x]].g_struct_.species_);
					std::cout<<graph[fermionic_edge[x]].g_struct_.species_ ;
				}
				
				fermionic_species.push_back(v);
				v.clear();
				std::cout<<")"<<std::endl;
															
													
						}mband::assign_label(graph,int_vector[2],initial_species_3);
										
					}
					else{mband::assign_label(graph,int_vector[2],initial_species_3);}
				}
				mband::assign_label(graph,int_vector[1],initial_species_2);
				
			}
			else {mband::assign_label(graph,int_vector[1],initial_species_2);}
			
		//mband::reset_species(graph,int_vector);
		mband::assign_label(graph,fermionic_edge,orginal_species);		
			}
    }
	else {std::cout << " No match found \n";}
mband::print_assigned_species(interaction_species);  
}




void mband::solve_pp_ord4(AmiGraph::graph_t &graph,AmiGraph::edge_vector_t &fermionic_edge,
std::vector<std::vector<int>> &fermionic_species,std::vector<std::vector<std::vector<int>>> &interaction_species,std::vector<std::vector<int>> &bosonic_Alpha,std::vector<std::vector<int>> &gkkp_Alpha,std::vector<int> &bandindex){
         ///do pp stuff here first///
	AmiGraph::edge_vector_t pair1;
	AmiGraph::edge_vector_t pair2;
	mband::find_connected_pp_edge(graph,pair1,pair2);
	mband::assign_bandindex_fourline_pp(graph,pair1,pair2, bandindex);
	gp.find_internal_fermionic_edges(graph,fermionic_edge);
	std::vector<int> orginal_species = mband::generate_edge_species(graph, fermionic_edge);
	

	AmiGraph::edge_vector_t bvector;
	std::vector<AmiGraph::edge_vector_t> int_vector;
	gp.find_non_external_bosonic_edges(graph,bvector);
	
     std::cout <<"bosonic edges are \n";
	 for (auto b: bvector){
		gp.print_edge_info(b,graph);
		std::cout<<std::endl;
		bosonic_Alpha.push_back(graph[b].g_struct_.alpha_);
	}
	print2d(bosonic_Alpha);
  for (auto b: pair1){
		gkkp_Alpha.push_back(graph[b].g_struct_.alpha_);
	}	
	///pp stuff ends here///
	for (int i=0;i <bvector.size();i++){
		//g.print_edge_info(bvector[i],gself2);
		std::cout << "Edge = ("  << graph[source(bvector[i],graph)].index_<<"," << graph[target(bvector[i],graph)].index_ <<")" 
		<<std::endl;
	}
	 for (auto b: bvector){
		bosonic_Alpha.push_back(graph[b].g_struct_.alpha_);
	}
	mband::find_interaction(graph,bvector,int_vector);
	mband::print_interactions(graph,bvector, int_vector);
	std::vector<std::vector<int>> possible_species_1;
	std::vector<int> initial_species_1 = mband::generate_edge_species(graph, int_vector[0]);
	if (mband::Hartee_fock){
		possible_species_1 = findmatch(initial_species_1,int_vector[0]);
	}
	else{
		possible_species_1 = mband::interaction_legs;
	}
	
	std::cout<< " Possible external species o4 in curly bracket and internal species in () is shown below \n";
	if (!possible_species_1.empty()){
		for (int i = 0; i<possible_species_1.size();i++){
			
			std::vector<int> initial_species_1 = mband::generate_edge_species(graph, int_vector[0]);	
			mband::assign_label(graph,int_vector[0],possible_species_1[i]);
			std::vector<int> initial_species_2 = mband::generate_edge_species(graph, int_vector[1]);
			std::vector<std::vector<int>> possible_species_2 = mband::findmatch(initial_species_2,int_vector[1]);
		

			if (!possible_species_2.empty()){
				for (int j = 0; j < possible_species_2.size();j++){
				//mband::print_match(possible_species_2[j],2);
				mband::assign_label(graph,int_vector[1],initial_species_2);				
				mband::assign_label(graph,int_vector[1],possible_species_2[j]);
				std::vector<int> initial_species_3 = mband::generate_edge_species(graph, int_vector[2]);
				std::vector<std::vector<int>> possible_species_3 = mband::findmatch(initial_species_3,int_vector[2]);
				
					if (!possible_species_3.empty()){
						for (int k = 0; k < possible_species_3.size();k++){
							//mband::print_match(possible_species_3[k],3);
							mband::assign_label(graph,int_vector[2],initial_species_3);
							mband::assign_label(graph,int_vector[2],possible_species_3[k]);											
						//interaction_species.push_back({mband::interaction_legs[i],possible_species_2[j],possible_species_3[k]});
											
							std::vector<int> initial_species_4 = mband::generate_edge_species(graph, int_vector[3]);
							std::vector<std::vector<int>> possible_species_4 = mband::findmatch(initial_species_4,int_vector[3]);					
							if (!possible_species_4.empty()){
								for (int l = 0; l < possible_species_4.size();l++){
									//num++
								//mband::print_match(possible_species_4[l],4);
								mband::assign_label(graph,int_vector[3],initial_species_4);							
								mband::assign_label(graph,int_vector[3],possible_species_4[l]);
								interaction_species.push_back({possible_species_1[i],possible_species_2[j],possible_species_3[k],possible_species_4[l]});
								
								std::cout<< fermionic_edge.size();
								std::cout<<j << " printing internal fermionic band index" << std::endl;
								std::cout<<"(";		
								std::vector<int> v;
								for (int x =0; x < fermionic_edge.size(); x++){
									v.push_back(graph[fermionic_edge[x]].g_struct_.species_);
									std::cout<<graph[fermionic_edge[x]].g_struct_.species_;
					
				}
				fermionic_species.push_back(v);
				v.clear();
								std::cout<<")"<<std::endl;
								}mband::assign_label(graph,int_vector[3],initial_species_4);
								
							}
							else { mband::assign_label(graph,int_vector[3],initial_species_4);}				
						}mband::assign_label(graph,int_vector[2],initial_species_3);
										
					}
					else{mband::assign_label(graph,int_vector[2],initial_species_3);}
				}
				mband::assign_label(graph,int_vector[1],initial_species_2);			
			}
			else {mband::assign_label(graph,int_vector[1],initial_species_2);}
			
		//mband::reset_species(graph,int_vector);
		mband::assign_label(graph,fermionic_edge,orginal_species);	
		}
	}
	else{std::cout<< "NO match found \n";}
mband::print_assigned_species(interaction_species);
}



void mband::solve_pp_ord(AmiGraph::graph_t &graph,AmiGraph::edge_vector_t &fermionic_edge,std::vector<std::vector<int>> &fermionic_species,std::vector<std::vector<std::vector<int>>> &interaction_species,std::vector<std::vector<int>> &bosonic_Alpha,std::vector<std::vector<int>> &gkkp_Alpha,std::vector<int> &bandindex)
{
	int ord =  gp.graph_order(graph);
	if (ord ==0)
	std::cerr << "0-th order not yet possible \n";
	else if (ord==1){
	 mband::solve_pp_ord1(graph,fermionic_edge,fermionic_species,interaction_species,bosonic_Alpha,gkkp_Alpha,bandindex);	
	}
	else if (ord==2){	
	 mband::solve_pp_ord2(graph,fermionic_edge,fermionic_species,interaction_species,bosonic_Alpha,gkkp_Alpha,bandindex);	
	}
	else if (ord==3){	
	 mband::solve_pp_ord3(graph,fermionic_edge,fermionic_species,interaction_species,bosonic_Alpha,gkkp_Alpha,bandindex);
	}
	else if (ord==4){
	 mband::solve_pp_ord4(graph,fermionic_edge,fermionic_species,interaction_species,bosonic_Alpha,gkkp_Alpha,bandindex);	
	}
	else{
		std::cout << "this order not possible";}
	
}


void mband::pp_sampler( AmiGraph::graph_t &graph, mband::sampler_collector& collector,std::vector<int> &bandindex){ 
    mband::solve_pp_ord(graph,collector.fermionic_edge,collector.fermionic_edge_species,collector.interaction_species,collector.bosonic_Alpha,
	collector.gkkp_Alpha,bandindex);	
	mband::generate_eps_alpha(graph,collector.fermionic_edge,collector.Epsilon,collector.Alpha);
	collector.graph = graph;
	for (auto interac: collector.interaction_species){
	collector.Uindex.push_back(interaction_index(interac));
	}
    std::vector<std::vector<int>> ext;
	for (int i =0; i <collector.fermionic_edge_species.size();i++){	
		ext.push_back(bandindex);
	}
	
	std::cout << "gkkp_alpha are" <<std::endl;
	print2d(collector.gkkp_Alpha);
	collector.external_line=ext;	
}



std::tuple<std::complex<double>, std::complex<double>, int> mband::lcalc_sampled_pp(AmiGraph::graph_t &gself, std::vector<AmiBase::epsilon_t>& Epsilon, std::vector<AmiBase::alpha_t>& Alpha,std::vector<std::vector<int>> &bosonic_Alpha,std::vector<std::vector<int>> &gkkp_Alpha,std::vector<int> &Utype,
 std::vector<int>& Species,NewAmiCalc::ext_vars& ext_params,int MC_num,params_param& param) {
	int cutoff_num = 0;
    AmiBase ami;
	std::vector<int> gg_type = {param.in%10, param.out%10};
    int ord = gp.graph_order(gself);
	int loop = gp.count_fermi_loops(gself);
	
	double power= (double) (ord+loop);
    //double prefactor = std::pow(-1,power);
	double prefactor = gp.get_prefactor(gself,ord);
	//std::cout << "1prefactor used is" << prefactor << "with loop "<< loop<<"ord is "<< ord<<"bubble no is " <<gp.count_bubbles(gself) <<std::endl;
	//std::cout << "2prefactor used is" << prefactor1 <<std::endl; 
	
    int n = 2 * ord +2; // number of fermionic lines
    AmiBase::g_struct gs[n];
    std::vector<AmiBase::g_struct> gs_vec;
	

    for (int i = 0; i < n; i++) {
        gs[i] = {Epsilon[i], Alpha[i]};
        gs_vec.push_back(gs[i]);
    }

    AmiBase::g_prod_t R0 = gs_vec;
    AmiBase::S_t S_array;
    AmiBase::P_t P_array;
    AmiBase::R_t R_array;
    
    double E_REG = param.E_reg;
    int N_INT = ord+1;
    ami.precision_cutoff=param.set_precision;
    AmiBase::graph_type bose=AmiBase::Pi_ppud;
	AmiBase::ami_parms test_amiparms(N_INT, E_REG,bose);
    ami.construct(test_amiparms, R0, R_array, P_array, S_array);
    AmiBase::frequency_t frequency;
	
    for (int i = 0; i < N_INT; i++) {
        frequency.push_back(std::complex<double>(0, 0));
    }
    frequency.push_back(ext_params.external_freq_[0]);

    std::random_device rd;
    std::default_random_engine engine(rd());
    std::uniform_real_distribution<double> distribution(0, 1.0);
    int kspace = Alpha[0].size() - 1;
    std::complex<double> localSum(0.0, 0.0);
    std::complex<double> localSumOfSquares(0.0, 0.0);
	
	
	// Storage Structures
	AmiBase::g_prod_t unique;
	AmiBase::R_ref_t rref;
	AmiBase::ref_eval_t eval_list;

	// Take existing solution from first part and factorize it 
	ami.factorize_Rn(R_array.back(), unique, rref, eval_list);

	
	
	
    // Generate random samples and calculate the sum
    for (int i = 0; i < MC_num; i++) {
        std::vector<std::vector<double>> momenta;
        momenta.reserve(kspace);
        for (int j = 0; j < kspace; j++) {
		if (param.lattice_type != 5 && param.lattice_type !=6){
            double momentum1 = 2 * M_PI * distribution(engine);
            double momentum2 = 2 * M_PI * distribution(engine);
            momenta.push_back({momentum1, momentum2});
		}
		else {
			//std::cout<<"Generating momenta for Hexagon \n";
			std::pair<double,double> ks =  mband::generate_hex_bz();
			double momentum1 = ks.first;
            double momentum2 = ks.second ;
            momenta.push_back({momentum1, momentum2});
			
			}
        }
        momenta.push_back(ext_params.external_k_list_[0]);

        std::vector<std::vector<double>> summed_momenta;
        summed_momenta.reserve(Alpha.size());
        for (const auto& alpha : Alpha) {
            double qx = 0;
            double qy = 0;
            for (int j = 0; j < alpha.size(); j++) {
                qx += static_cast<double>(alpha[j]) * momenta[j][0];
                qy += static_cast<double>(alpha[j]) * momenta[j][1];
            }
            summed_momenta.push_back({qx, qy});
        }
	double form_factor = 1;
	

    

	if (param.lattice_type ==2 ||param.lattice_type ==3 || param.lattice_type ==4 || param.lattice_type ==6 ||param.lattice_type ==7 || param.lattice_type ==8  ){		
		std::vector<std::vector<double>> V_momenta;
		V_momenta.reserve(bosonic_Alpha.size());

		for (const auto& b_alpha : bosonic_Alpha) {
			double qx = 0;
			double qy = 0;
			for (int j = 0; j < b_alpha.size(); j++) {
				qx += static_cast<double>(b_alpha[j]) * momenta[j][0];
				qy += static_cast<double>(b_alpha[j]) * momenta[j][1];
			}
		   V_momenta.push_back({ qx, qy });
		}
		if (param.lattice_type==4){
                        for (int i =0;i<Utype.size();i++){
                                if (Utype[i]>21 && Utype[i]<30){
                                form_factor= form_factor*(1.0+  2*param.V*(std::cos(V_momenta[i][0])+std::cos(V_momenta[i][1])));
                                }
                        }
                }
		if (param.lattice_type==8){
			//std::cout<< "triggering pair hopping for quad layer\n";
                        for (int i =0;i<Utype.size();i++){
                                if (Utype[i]>31 && Utype[i]<44){
                                form_factor= form_factor*(1.0+  2*param.V*(std::cos(V_momenta[i][0])+std::cos(V_momenta[i][1])));
                                }
                        }
                }

  
		if (param.lattice_type ==3){
			for (int i = 0; i<Utype.size();i++){
				if (Utype[i]==0 || Utype[i]==1){
				form_factor= form_factor*2*ext_params.MU_.imag()*(std::cos(V_momenta[i][0])+std::cos(V_momenta[i][1]));
						}
				else {
				form_factor= form_factor*(1.0+  2*ext_params.H_*(std::cos(V_momenta[i][0])+std::cos(V_momenta[i][1])));	
	
				}	
			}
		}
		if (param.lattice_type ==6){
			for (int i = 0; i<Utype.size();i++){
				if (Utype[i]==0 || Utype[i]==1){
				form_factor= form_factor*2*param.V*(std::cos(V_momenta[i][1]) + 2*std::cos(V_momenta[i][1]/2)*std::cos(std::sqrt(3)*V_momenta[i][0]/2) );
						}
				else {
				form_factor= form_factor*(1.0+  2*param.V*(std::cos(V_momenta[i][1]) + 2*std::cos(V_momenta[i][1]/2)*std::cos(std::sqrt(3)*V_momenta[i][0]/2)));	
	
				}	
			}
		}
	    if (param.lattice_type==2){
			for (int i =0;i<Utype.size();i++){
				if (Utype[i]>11 && Utype[i]<16){
				form_factor= form_factor*(1.0+  2*param.V*(std::cos(V_momenta[i][0])+std::cos(V_momenta[i][1])));
				}
				
			}
			
		}
		if (param.lattice_type==7){
			for (int i = 0; i<Utype.size();i++){
				if (Utype[i] <= 5){
			    //std::cout << " triggering extended result  with Utype " << Utype[i] << "with value " << ext_params.MU_.imag() << "and " <<ext_params.H_ <<"respectively" << std::endl;
				form_factor= form_factor*( 2*ext_params.MU_.imag()*(std::cos(V_momenta[i][0])+std::cos(V_momenta[i][1])) +
				4*ext_params.H_*(std::cos(V_momenta[i][0])*std::cos(V_momenta[i][1])));
						}
				else if(Utype[i] >=6 && Utype[i]<=11)  {
				//std::cout << " triggering extended result  with Utype " << Utype[i] << "with value " << ext_params.MU_.imag() << "and " <<ext_params.H_ <<"respectively" << std::endl;
				form_factor= form_factor*(1.0+  2*ext_params.MU_.imag()*(std::cos(V_momenta[i][0])+std::cos(V_momenta[i][1])) +
				4*ext_params.H_*(std::cos(V_momenta[i][0])*std::cos(V_momenta[i][1])));	
				}
				else if(Utype[i] >=12 && Utype[i]<=35)  {
				double D = 1-2*ext_params.MU_.real();
				//std::cout << " triggering extended result  with Utype " << Utype[i] << "with value " << ext_params.MU_.imag() << "and " <<ext_params.H_ <<"respectively" << std::endl;
				form_factor= form_factor*(D+  2*ext_params.MU_.imag()*(std::cos(V_momenta[i][0])+std::cos(V_momenta[i][1])) +
				4*ext_params.H_*(std::cos(V_momenta[i][0])*std::cos(V_momenta[i][1])));	
				}	
			}
		}
		
	} 
	std::complex<double> gk(1,0);
	
	if (param.G_FUNC > 0){		
		std::vector<std::vector<double>> V1_momenta;
		V1_momenta.reserve(gkkp_Alpha.size());
		for (const auto& gk_alpha : gkkp_Alpha) {
			double qx = 0;
			double qy = 0;
			for (int j = 0; j < gk_alpha.size(); j++) {
				qx += static_cast<double>(gk_alpha[j]) * momenta[j][0];
				qy += static_cast<double>(gk_alpha[j]) * momenta[j][1];
			}
		   V1_momenta.push_back({ qx, qy });
		}
		/*
		std::cout << "print summed momenta\n" ;
		print2d(momenta);
		std::cout <<"\n V1_momenta are" <<std::endl;
		print2d(V1_momenta);
		*/
		
		for (int i = 0; i<gkkp_Alpha.size();i++){
			std::complex<double> g1= mband::gfunc_pp(V1_momenta[i],param,gg_type[i]);
			
		
			if (i ==1){
				g1 = std::conj(g1);
			}
			
			gk = gk*g1;
			//std::cout<<gk;
	
        }
	}
		

	/////////////////////////////////energy/////////////////////////////	
        std::vector<double> energy;
        energy.reserve(summed_momenta.size());
        for (int i = 0; i < summed_momenta.size(); i++) {
            if (param.lattice_type == 1 || param.lattice_type == 3 ) {
                energy.push_back(mband::Hubbard_Energy(ext_params, summed_momenta[i], Species[i],param));
            }
			if (param.lattice_type == 2 ) {
                energy.push_back(mband::Bilayer_Hubbard_Energy(ext_params, summed_momenta[i], Species[i],param));
            }
			if (param.lattice_type == 4 ) {
                energy.push_back(mband::Trilayer_Hubbard_Energy(ext_params, summed_momenta[i], Species[i],param));
            }
			if (param.lattice_type == 8 ) {
				//std::cout<<"triggering quadlayer energy"<<std::endl;
                energy.push_back(mband::Quadlayer_Hubbard_Energy(ext_params, summed_momenta[i], Species[i],param));
            }
			if (param.lattice_type == 5 || param.lattice_type == 6 ) {
                energy.push_back(mband::Triangular_Hubbard_Energy(ext_params, summed_momenta[i], Species[i],param));
            }
			if (param.lattice_type == 7 ) {
				//std::cout << "Using SRO energy";
                energy.push_back(mband::SRO_Hubbard_Energy(ext_params, summed_momenta[i], Species[i],param));
            }
        }
        std::vector<std::complex<double>> energy_t = mband::generate_ept(Epsilon, energy);

        AmiBase::ami_vars external(energy_t, frequency, ext_params.BETA_);
	
        std::complex<double> raw_coeff = ami.evaluate(test_amiparms, R_array, P_array, S_array, external,unique, rref, eval_list);
		std::complex<double> result =form_factor*gk* prefactor *raw_coeff;

	if (std::isnan(raw_coeff.real()) || std::isnan(raw_coeff.imag()) ||
        std::isinf(raw_coeff.real()) || std::isinf(raw_coeff.imag())) {
        std::cout << "Result is NaN or infinity, dropping the result." << std::endl;
        std::cout<<"result is " <<raw_coeff<<std::endl;
		cutoff_num++;    
    }
	else if ((abs(raw_coeff.real()) > param.cutoff_value || abs(raw_coeff.imag()) > param.cutoff_value) && ord>2){
		std::cout << "Large value detected." << std::endl;
        std::cout<<"result is " <<raw_coeff<<std::endl;
		cutoff_num++; 
		
	}	
	else if (ami.overflow_detected) {
			std::cout<<"over flow detected \n ";
			std::cout<<"result is " <<raw_coeff<<std::endl;
			cutoff_num++;     		
		}
	
	else{
        localSum +=  result;
        localSumOfSquares += std::complex<double> (std::pow(result.real(),2),std::pow(result.imag(),2)) ;
		}
    }
    int samples = MC_num-cutoff_num;
    return  std::make_tuple(localSum, localSumOfSquares, samples);
}
	

     
#include "mini_ami.hpp"


AmiGraph g_ph(AmiBase::Pi_phuu, 0);

void mband::findInternalBoson(AmiGraph::graph_t &graph, AmiGraph::edge_vector_t &internal_b_vector){
    AmiGraph::edge_vector_t all_b_vector;
    g_ph.find_bosonic_edges(graph, all_b_vector);
    AmiGraph::vertex_vector_t v;
    AmiGraph::edge_vector_t external_bosonic_edges;
    g_ph.find_external_vertices(graph, v, external_bosonic_edges);
    for (int i = 0; i < all_b_vector.size(); i++){
        bool not_same_vector = true;
        for (int j = 0; j < external_bosonic_edges.size(); j++){
            if ((graph[source(all_b_vector[i], graph)].index_ == graph[source(external_bosonic_edges[j], graph)].index_) &&
                    (graph[target(all_b_vector[i], graph)].index_ == graph[target(external_bosonic_edges[j], graph)].index_)){
                not_same_vector = false;
            }
        }
        if (not_same_vector == true){
            internal_b_vector.push_back(all_b_vector[i]);
        }
    }
}


void mband::findNotAssignedFermion(AmiGraph::graph_t &graph, AmiGraph::edge_vector_t &internal_f_vector){
    AmiGraph::vertex_vector_t v;
    AmiGraph::edge_vector_t external_bosonic_edges;
    g_ph.find_external_vertices(graph, v, external_bosonic_edges);
    AmiGraph::edge_vector_t ext_bosonic_source;
    AmiGraph::edge_vector_t ext_bosonic_target;
    mband::find_interaction_of_source(graph, external_bosonic_edges[1], ext_bosonic_source);
    mband::find_interaction_of_target(graph, external_bosonic_edges[0], ext_bosonic_target);
    AmiGraph::edge_vector_t ext_fermi(ext_bosonic_source.begin(), ext_bosonic_source.end());
    ext_fermi.insert(ext_fermi.end(), ext_bosonic_target.begin(), ext_bosonic_target.end());
    AmiGraph::edge_vector_t all_f_vector;
    g_ph.find_fermionic_edges(graph, all_f_vector);
    for (int i = 0; i < all_f_vector.size(); i++){
        bool not_same_vector = true;
        for(int j = 0; j < ext_fermi.size(); j++){
            if ((graph[source(all_f_vector[i], graph)].index_ == graph[source(ext_fermi[j], graph)].index_) &&
            (graph[target(all_f_vector[i], graph)].index_ == graph[target(ext_fermi[j], graph)].index_)){
                not_same_vector = false;
            }
        }
        if (not_same_vector == true){
            internal_f_vector.push_back(all_f_vector[i]);
        }
    }

}


void mband::generateCombination(std::vector<int> &bands, int &n, int index, std::vector<int> combination, std::vector<std::vector<int>> &all_combinations){
    if (index == n) {
        all_combinations.push_back(combination);
        return;
    }
    for (int i=0; i < bands.size(); i++){
        combination[index] = bands[i];
        mband::generateCombination(bands, n, index+1, combination, all_combinations);
    }
}


void mband::findBandCombination(AmiGraph::graph_t &graph, std::vector<int> &bands, AmiGraph::edge_vector_t &internal_f_vector, std::vector<std::vector<int>> &all_combinations){
    mband::findNotAssignedFermion(graph, internal_f_vector);
    int n = internal_f_vector.size();
    std::vector<int> combination (n, 0);
    mband::generateCombination(bands, n, 0, combination, all_combinations);
}


void mband::findSpeciesForUVector(AmiGraph::graph_t &graph, std::vector<AmiGraph::edge_vector_t> &U_int_vector,
                                  std::vector<std::vector<int>> &U_int_vector_species){
    for (int i = 0; i < U_int_vector.size(); i++){
        std::vector<int> sp_vector;
        for (int j = 0; j < U_int_vector[i].size(); j++){
            sp_vector.push_back(graph[U_int_vector[i][j]].g_struct_.species_);
        }
        U_int_vector_species.push_back(sp_vector);
    }
}


void mband::checkCombination(AmiGraph::graph_t &graph, std::vector<int> &bands,
                             AmiGraph::edge_vector_t &all_f_vector,
                             std::vector<std::vector<int>> &good_combinations,
                             std::vector<std::vector<int>> &bosonic_Alpha,
                             std::vector<std::vector<int>> &eps,
							  std::vector<std::vector<int>> &alp,
                             std::vector<AmiGraph::edge_vector_t> u_int_vector,
                             AmiGraph::edge_vector_t internal_b_vector,
                             std::vector<std::vector<std::vector<int>>> &good_U_combinations){
    std::vector<std::vector<int>> allowed_U_combination = mband::interaction_legs;
    std::vector<std::vector<int>> all_combinations;
    AmiGraph::edge_vector_t internal_f_vector;
    mband::findBandCombination(graph, bands, internal_f_vector, all_combinations);
//    std::vector<AmiGraph::edge_vector_t> u_int_vector;
//    AmiGraph::edge_vector_t internal_b_vector;
//    mband::find4VertexInteractions(graph, internal_b_vector, u_int_vector);
    g_ph.find_fermionic_edges(graph, all_f_vector);
    for (auto b: internal_b_vector){
        bosonic_Alpha.push_back(graph[b].g_struct_.alpha_);
    }
    for (auto f: all_f_vector){
        eps.push_back(graph[f].g_struct_.eps_);
		alp.push_back(graph[f].g_struct_.alpha_);
    }
//    if (internal_f_vector.empty()){
//        std::vector<std::vector<int>> u_int_vector_species;
//        mband::findSpeciesForUVector(graph, u_int_vector, u_int_vector_species);
//        std::vector<int> combination;
//        for (int m = 0; m < all_f_vector.size(); m ++){
//            combination.push_back(graph[all_f_vector[m]].g_struct_.species_);
//        }
//        good_combinations.push_back(combination);
//        good_U_combinations.push_back(u_int_vector_species);
//        return;
//    }
    for (int i = 0; i < all_combinations.size(); i++){
        bool big_check = true;
        for (int j = 0; j < all_combinations[0].size(); j++){
            graph[internal_f_vector[j]].g_struct_.species_ = all_combinations[i][j];
        }
        std::vector<std::vector<int>> u_int_vector_species;
        mband::findSpeciesForUVector(graph, u_int_vector, u_int_vector_species);
        for (int k = 0; k < u_int_vector_species.size(); k++){
            bool check = false;
            for (int l = 0; l < allowed_U_combination.size(); l++){
                if (u_int_vector_species[k] == allowed_U_combination[l]){
                    check = true;
                }
            }
            if (check == false){
                big_check = false;
                goto end_loop;
            }
        }
        end_loop:
        if (big_check == true){
            std::vector<int> combination;
            for (int m = 0; m < all_f_vector.size(); m ++){
                combination.push_back(graph[all_f_vector[m]].g_struct_.species_);
            }
            good_combinations.push_back(combination);
            good_U_combinations.push_back(u_int_vector_species);
        }
    }
}






void mband::find_interaction_of_source(AmiGraph::graph_t &graph, AmiGraph::edge_t &b_vector,
                                AmiGraph::edge_vector_t &v){
    boost::graph_traits<AmiGraph::graph_t>::in_edge_iterator iei, iedge_end;
    boost::graph_traits<AmiGraph::graph_t>::out_edge_iterator oei, oedge_end;
    for (boost::tie(iei, iedge_end) = in_edges(source(b_vector, graph), graph); iei != iedge_end; ++iei){
        if (graph[*iei].g_struct_.stat_==AmiBase::Fermi){
            v.push_back(*iei);
        }
    }
    for (boost::tie(oei, oedge_end) = out_edges(source(b_vector, graph), graph); oei != oedge_end; ++oei) {
        if (graph[*oei].g_struct_.stat_ == AmiBase::Fermi) {
            v.push_back(*oei);
        }
    }
}


void mband::find_interaction_of_target(AmiGraph::graph_t &graph, AmiGraph::edge_t &b_vector,
                                AmiGraph::edge_vector_t &v){

    boost::graph_traits<AmiGraph::graph_t>::in_edge_iterator iei, iedge_end;
    boost::graph_traits<AmiGraph::graph_t>::out_edge_iterator oei, oedge_end;
    for (boost::tie(iei, iedge_end) = in_edges(target(b_vector, graph), graph); iei != iedge_end; ++iei){
        if( graph[*iei].g_struct_.stat_ == AmiBase::Fermi){
            v.push_back(*iei);
        }
    }
    for (boost::tie(oei, oedge_end) = out_edges(target(b_vector, graph), graph); oei != oedge_end; ++oei){
        if( graph[*oei].g_struct_.stat_==AmiBase::Fermi){
            v.push_back(*oei);
        }
    }
}


void mband::assign_initial_species(AmiGraph::graph_t &g,
                            AmiGraph::edge_vector_t ext_bosonic_source,
                            AmiGraph::edge_vector_t ext_bosonic_target,
                            std::vector<int> possible_interaction){
	 
								
    g[ext_bosonic_source[0]].g_struct_.species_ = possible_interaction[1] / 10;
    g[ext_bosonic_source[1]].g_struct_.species_ = possible_interaction[1] % 10;
    g[ext_bosonic_target[1]].g_struct_.species_ = possible_interaction[0] / 10 ;
    g[ext_bosonic_target[0]].g_struct_.species_ = possible_interaction[0] % 10;
}
///////added functionality
bool mband::check_valid_initial_species(AmiGraph::graph_t &g,
                            AmiGraph::edge_vector_t ext_bosonic_source,
                            AmiGraph::edge_vector_t ext_bosonic_target,
                            std::vector<int> possible_interaction){
	 
	if (g[ext_bosonic_source[0]].g_struct_.species_ == possible_interaction[1] / 10 &&
    g[ext_bosonic_source[1]].g_struct_.species_ == possible_interaction[1] % 10 &&
    g[ext_bosonic_target[1]].g_struct_.species_ == possible_interaction[0] / 10 && 
    g[ext_bosonic_target[0]].g_struct_.species_ == possible_interaction[0] % 10 ){
		return true;
		
	}
	else { std::cout <<"Initial assingment band species produced spurious comination\n";


	return false;}
}


void mband::findInitialSpeciesPH(AmiGraph::graph_t &g,
                          std::vector<int> &band_ind,
                          AmiGraph::edge_vector_t &external_bosonic_edges,
                          AmiGraph::edge_vector_t &ext_bosonic_source,
                          AmiGraph::edge_vector_t &ext_bosonic_target){
    AmiGraph::vertex_vector_t v;
    g_ph.find_external_vertices(g, v, external_bosonic_edges);
    mband::find_interaction_of_source(g, external_bosonic_edges[1], ext_bosonic_source);
    mband::find_interaction_of_target(g, external_bosonic_edges[0], ext_bosonic_target);
    mband::assign_initial_species(g, ext_bosonic_source, ext_bosonic_target, band_ind);
}



void mband::find4VertexInteractions(AmiGraph::graph_t &graph, AmiGraph::edge_vector_t &internal_b_vector, std::vector<AmiGraph::edge_vector_t> &f_vector){
//    AmiGraph::edge_vector_t internal_b_vector;
    mband::findInternalBoson(graph, internal_b_vector);
    boost::graph_traits<AmiGraph::graph_t>::in_edge_iterator iei, iedge_end;
    boost::graph_traits<AmiGraph::graph_t>::out_edge_iterator oei, oedge_end;
    for (int i=0; i < internal_b_vector.size(); i++){
        AmiGraph::edge_vector_t v;
        for (boost::tie(iei, iedge_end) = in_edges(source(internal_b_vector[i], graph), graph); iei != iedge_end; ++iei){
            if (graph[*iei].g_struct_.stat_ == AmiBase::Fermi){
                v.push_back(*iei);
            }
        }

        for (boost::tie(oei, oedge_end) = out_edges(source(internal_b_vector[i], graph), graph); oei != oedge_end; ++oei){
            if( graph[*oei].g_struct_.stat_==AmiBase::Fermi){
                v.push_back(*oei);}
        }

        for (boost::tie(iei, iedge_end) = in_edges(target(internal_b_vector[i], graph), graph); iei != iedge_end; ++iei){
            if( graph[*iei].g_struct_.stat_==AmiBase::Fermi){
                v.push_back(*iei);}
        }

        for (boost::tie(oei, oedge_end) = out_edges(target(internal_b_vector[i], graph), graph); oei != oedge_end; ++oei){
            if( graph[*oei].g_struct_.stat_==AmiBase::Fermi){
                v.push_back(*oei);}
        }
        f_vector.push_back(v);
    }
}


void mband::ph_sampler(AmiGraph::graph_t &graph, mband::sampler_collector &collector, std::vector<int> &bandindex, std::vector<int>& possible_bands){

    AmiGraph::edge_vector_t external_bosonic_edges;
    AmiGraph::edge_vector_t ext_bosonic_left;
    AmiGraph::edge_vector_t ext_bosonic_right;
    mband::findInitialSpeciesPH(graph, bandindex, external_bosonic_edges, ext_bosonic_right,  ext_bosonic_left);
    if (mband::check_valid_initial_species(graph,  ext_bosonic_right,  ext_bosonic_left,bandindex)){
    std::cout << "print external bosonic edges" << std::endl;
    for (auto i : external_bosonic_edges){
        g_ph.print_edge_info(i, graph);
    }

    std::cout << "print external interactions source" << std::endl;
    for (int j = 0; j < ext_bosonic_right.size(); j++){
        g_ph.print_edge_info(ext_bosonic_right[j], graph);
    }
    std::cout << std::endl;

    std::cout << "print external interactions target" << std::endl;
    for (int j = 0; j < ext_bosonic_left.size(); j++){
        g_ph.print_edge_info(ext_bosonic_left[j], graph);
    }
    std::cout << std::endl;

    std::vector<AmiGraph::edge_vector_t> u_int_vector;
    AmiGraph::edge_vector_t internal_b_vector;
    mband::find4VertexInteractions(graph, internal_b_vector, u_int_vector);


    
    mband::checkCombination(graph, possible_bands, collector.fermionic_edge, collector.fermionic_edge_species,
                            collector.bosonic_Alpha, collector.Epsilon,collector.Alpha, u_int_vector, internal_b_vector, collector.interaction_species);

    collector.graph = graph;
    for (auto interac: collector.interaction_species){
        collector.Uindex.push_back(interaction_index(interac));
    }
    std::vector<std::vector<int>> ext;
    for (int i =0; i <collector.fermionic_edge_species.size();i++){
        ext.push_back(bandindex);
    }
    collector.external_line=ext;

    for (int k = 0; k < collector.interaction_species.size(); k++){
        std::cout << "Combination " << k+1 << std::endl;
        for (int i=0; i < internal_b_vector.size();i++){
            std::cout << "Bosonic = ("  << graph[source(internal_b_vector[i],graph)].index_<<"," << graph[target(internal_b_vector[i],graph)].index_ <<") "
                      << "alpha = " ;
            print1d(collector.bosonic_Alpha[i]) ;
            std::cout <<std::endl;
            std::cout << "Interaction =" << std::endl;
            for (int j = 0; j < u_int_vector[i].size(); j++){
                std::cout << " (" << graph[source(u_int_vector[i][j],graph)].index_<<"," << graph[target(u_int_vector[i][j],graph)].index_ <<") "
                          << collector.interaction_species[k][i][j] <<std::endl;
            }
            std::cout << "\n";
        }

        for (int j = 0; j < collector.fermionic_edge.size(); j++){
            std::cout << "Edge = ("  << graph[source(collector.fermionic_edge[j],graph)].index_<<"," << graph[target(collector.fermionic_edge[j],graph)].index_ <<") "
                      << collector.fermionic_edge_species[k][j] << " ";
            print1d(collector.Epsilon[j]) ;
            std::cout <<std::endl;
        }
        std::cout << "\n";
    }
	for (auto interac: collector.interaction_species){
	collector.Uindex.push_back(interaction_index(interac));
		}
	}
	else {
		
	std::cout << "Invalid species assigned in the initial band assignment stage. Combinations ignored" <<std::endl;}

}



std::tuple<std::complex<double>, std::complex<double>, int> mband::lcalc_sampled_ph(AmiGraph::graph_t &gself,
                                 std::vector<AmiBase::epsilon_t>& Epsilon,
                                 std::vector<AmiBase::alpha_t>& Alpha,std::vector<std::vector<int>> &bosonic_Alpha,
                                 std::vector<int> &Utype, std::vector<int>& Species,
                                 NewAmiCalc::ext_vars& ext_params,int MC_num,params_param& param) {
									 

	int cutoff_num = 0;
    AmiBase ami;
    int ord = g_ph.graph_order(gself);
        
    //double prefactor = std::pow(-1,power);
    double prefactor = g_ph.get_prefactor(gself,ord);
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
    AmiBase::graph_type bose=AmiBase::Pi_phuu;
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
		std::complex<double> result =form_factor* prefactor *raw_coeff;

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




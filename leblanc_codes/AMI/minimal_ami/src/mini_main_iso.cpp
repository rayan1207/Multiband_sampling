#include "mini_ami.hpp"
/*
AmiGraph gi(AmiBase::Pi_ppud, 3);
void check_iso(AmiGraph::graph_t g, AmiGraph::gg_matrix_t ggm, int min_ord, int max_ord) {
    bool match_found = false; // Add a flag to track if a match is found

    for (int i = min_ord; i < max_ord + 1; ++i) {
        for (int j = 0; j < ggm[i].size(); ++j) {
            for (int k = 0; k < ggm[i][j].graph_vec.size(); ++k) {
                if (gi.full_iso(g, ggm[i][j].graph_vec[k])) {
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
*/

int main(int argc, char** argv)
{ 

std::vector<std::vector<int>> interaction = readFile("../loader/bilayer_interaction.txt");
std::vector<double> interaction_value = readFile1("../loader/Hubbard_U.txt",5);
std::vector<double> band_energy = readFile1("../loader/ccpdvz_h.txt",2);
AmiBase ami;
int seed=0;	

//AmiGraph gp(AmiBase::Pi_ppud, seed);	
AmiGraph gp(AmiBase::Sigma,seed);
std::string infile("ext_vars.dat");
NewAmiCalc::external_variable_list extern_list;
gp.ami.read_external(infile, extern_list);	




AmiGraph::gg_matrix_t ggm_1; 
std::cout<<"Attempting to load self-energy graphs from example_graphs"<<std::endl;
int max=4;
//gp.read_ggmp("../graphs/ggm_ppvertex_ext_hubb/",ggm_1, max);
gp.read_ggmp("../graphs/ggm_test1/",ggm_1, max);
std::cout<<"Completed read"<<std::endl;
std::cout<<std::endl;
	
	
AmiGraph::gg_matrix_t ggm_2; 
std::cout<<"Attempting to load self-energy graphs from example_graphs"<<std::endl;
//gp.read_ggmp("../graphs/ggm_ppvertex_ext_hubb/",ggm_2, max);
gp.read_ggmp("../graphs/ggm_all/",ggm_2, max);
std::cout<<"Completed read"<<std::endl;
std::cout<<std::endl;
	
	


	






///////////////////////////////constructing mband ///////////////////////////////////////
	
bool hf = true; 
mband mb(interaction,interaction_value,band_energy,hf);

int min_ord =4;
int max_ord=4;

for (int i = min_ord; i < max_ord+1; ++i) {
    for (int j = 0; j < ggm_1[i].size(); ++j) {
        for (int k = 0; k < ggm_1[i][j].graph_vec.size(); ++k) {
			std::cout <<"checking o" <<std::to_string(i)<<"_g" <<std::to_string(j)<<"_n"<<std::to_string(k)<<" ";
			mb.check_iso_sigma(ggm_1[i][j].graph_vec[k],ggm_2,min_ord,max_ord);
				
		}
	}
}

			
			
			
			
			
}
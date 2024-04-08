#include "mini_ami.hpp"


int main(int argc, char** argv)
{ 

std::vector<std::vector<int>> interaction = readFile("../loader/bilayer_interaction.txt");
std::vector<double> interaction_value = readFile1("../loader/Hubbard_U.txt",5);
std::vector<double> band_energy = readFile1("../loader/ccpdvz_h.txt",2);
AmiBase ami;
int seed=0;	

AmiGraph g(AmiBase::Pi_ppud, seed);	
//AmiGraph g(AmiBase::Sigma,seed);
std::string infile("ext_vars.dat");
NewAmiCalc::external_variable_list extern_list;
g.ami.read_external(infile, extern_list);	




AmiGraph::gg_matrix_t ggm; 
AmiGraph::gg_matrix_t ggm_new;
std::cout<<"Attempting to load self-energy graphs from example_graphs"<<std::endl;
int max_ord=4;
int min_ord=1;
g.read_ggmp("../graphs/ggm_ppvertex_ext_hubb/",ggm, max_ord);
//g.read_ggmp("../graphs/ggm_sigma_nofock_notp/",ggm, max_ord);
std::cout<<"Completed read"<<std::endl;
std::cout<<std::endl;
g.ggm_label(ggm,seed);


ggm_new.resize(ggm.size());
for (int i = 0; i<ggm_new.size();i++){
	ggm_new[i].resize((int) 1);	
	for (int j =0; j <ggm_new[i].size();j++){
		ggm_new[i][j].graph_vec.resize((int) 2500);
		
	}
}

std::cout << "the size is " << ggm_new.size() <<std::endl;



for (int i = min_ord; i < max_ord+1; ++i) {
	int num = 0;
    for (int j = 0; j < ggm[i].size(); ++j) {
	
        for (int k = 0; k < ggm[i][j].graph_vec.size(); ++k) {
			
			std::cout <<"checking o" <<std::to_string(i)<<"_g" <<std::to_string(j)<<"_n"<<std::to_string(k)<<"";
			if (g.is_oneleg(ggm[i][j].graph_vec[k])){
				std::cout << " -> One legged diagram found - excluding it\n";
			}
			else{
				std::string name = "ggm_ppvertex_ext_hubb_nl/o"+std::to_string(i)+"_g"+std::to_string(0)+"_n"+std::to_string(num)+".graph";
				std::cout <<" -> o" <<std::to_string(i)<<"_g" <<std::to_string(0)<<"_n"<<std::to_string(num)<<" \n";
				
				AmiGraph::graph_t grp = ggm[i][j].graph_vec[k];
				
				ggm_new[i][0].graph_vec.push_back(grp) ;
				g.graph_write(name,grp);
				num++;			
			}
			
			
			
				
		}
	}
}

//g.write_ggmp("graphs",ggm_new,min_ord);			
			
			
			
			
}

#include "mini_ami.hpp"

int main(int argc, char** argv)
{ 
auto   startTime = std::chrono::high_resolution_clock::now();
AmiBase ami;
mband::params_param params;
params_loader("params.txt", params);
int seed =0;
AmiBase::graph_type baseType;
if (params.graph_type == 0) {
    baseType = AmiBase::Sigma;
} else if (params.graph_type ==1) {
    baseType = AmiBase::Pi_ppud;
}
else if (params.graph_type ==2) {
    baseType = AmiBase::Pi_phuu;
}

AmiGraph g(baseType, seed);

AmiGraph::gg_matrix_t ggm;	
NewAmiCalc::external_variable_list extern_list;
std::vector<std::vector<int>> interaction;
std::vector<double> interaction_value;
std::vector<double> band_energy;
std::vector<int> possible_bands;
MPI_Init(&argc, &argv);
int numProcesses, rank;
MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);


if (params.molecular==0&& (params.lattice_type == 1 || params.lattice_type == 5)){
    interaction = readFile("/project/6005092/rfarid/workflow_2023/leblanc_codes/AMI/minimal_ami/loader/Hubbard_U.txt");
	interaction_value = readFile1("/project/6005092/rfarid/workflow_2023/leblanc_codes/AMI/minimal_ami/loader/Hubbard_U.txt",5);
	 band_energy = {0,0};
    std::cout<< "Constructing AmiGraph object using seed: "<<seed <<" " <<std::endl;
	std::string infile("ext_vars.dat");
	std::cout<<"Reading external parameters from ext_vars.dat"<<std::endl;
	g.ami.read_external(infile, extern_list);	
	std::cout<<"Attempting to load self-energy graphs from example_graphs"<<std::endl;
	int max=params.max_ord;
	g.read_ggmp(params.graph,ggm, max);	
	
	
	std::cout<<"Completed read"<<std::endl;
	std::cout<<std::endl;
	g.ggm_label(ggm,0); 
	possible_bands={1,2};

    std::cout<<"External parameters read are"<<std::endl;
for(int i=0; i<extern_list.size();i++){
		std::cout<<extern_list[i].BETA_<<" "<<extern_list[i].MU_<<" "<< extern_list[i].H_<<" "<<extern_list[i].KDIM_<<" "<<extern_list[i].external_k_list_[0][0]<<" "<<extern_list[i].external_k_list_[0][1] <<" "<<extern_list[i].external_freq_[0]<<std::endl;
	}
}


else if (params.molecular==0&& (params.lattice_type ==2 || params.lattice_type ==4 ||params.lattice_type ==7 || params.lattice_type ==8 ) ){
	 if (params.lattice_type==2){
     interaction = readFile("/project/6005092/rfarid/workflow_2023/leblanc_codes/AMI/minimal_ami/loader/bilayer_interaction.txt");
	 interaction_value = readFile1("/project/6005092/rfarid/workflow_2023/leblanc_codes/AMI/minimal_ami/loader/bilayer_interaction.txt",5);
	 possible_bands={1,2,3,4};
     }
	 else if (params.lattice_type==4){
     interaction = readFile("/project/6005092/rfarid/workflow_2023/leblanc_codes/AMI/minimal_ami/loader/trilayer_interaction.txt");
	 interaction_value = readFile1("/project/6005092/rfarid/workflow_2023/leblanc_codes/AMI/minimal_ami/loader/trilayer_interaction.txt",5);
	  possible_bands={1,2,3,4,5,6};
     }
	 else if (params.lattice_type==7){
		 std::cout << " opening interaction file for SRO " <<std::endl;
		 interaction = readFile("/project/6005092/rfarid/workflow_2023/leblanc_codes/AMI/minimal_ami/loader/SRO_interaction.txt");
		 print2d(interaction);
		 interaction_value = readFile1("/project/6005092/rfarid/workflow_2023/leblanc_codes/AMI/minimal_ami/loader/SRO_interaction.txt",5);
		possible_bands={1,2,3,4,5,6};		 
	}
	 else if (params.lattice_type==8){
		 std::cout << " opening interaction file for quadlayer_interaction " <<std::endl;
		 interaction = readFile("/project/6005092/rfarid/workflow_2023/leblanc_codes/AMI/minimal_ami/loader/quadlayer_interaction.txt");
		 print2d(interaction);
		 interaction_value = readFile1("/project/6005092/rfarid/workflow_2023/leblanc_codes/AMI/minimal_ami/loader/quadlayer_interaction.txt",5);	
         possible_bands={1,2,3,4,5,6,7,8};			 
	}
	 
	 band_energy = {0,0};

    std::cout<< "Constructing AmiGraph object using seed: "<<seed <<" " <<std::endl;	
	std::string infile("ext_vars.dat");
	std::cout<<"Reading external parameters from ext_vars.dat"<<std::endl;
	g.ami.read_external(infile, extern_list);	
	std::cout<<"Attempting to load self-energy graphs from example_graphs"<<std::endl;
	int max=params.max_ord;
	

	g.read_ggmp(params.graph,ggm, max);	
	
	g.ggm_label(ggm,0);    
	
	std::cout<<"External parameters read are"<<std::endl;
    for(int i=0; i<extern_list.size();i++){
		std::cout<<extern_list[i].BETA_<<" "<<extern_list[i].MU_<<" "<< extern_list[i].H_<<" "<<extern_list[i].KDIM_<<" "<<extern_list[i].external_k_list_[0][0]<<" "<<extern_list[i].external_k_list_[0][1] <<" "<<extern_list[i].external_freq_[0]<<std::endl;
	
		}	

  }
  
  else if (params.molecular==0&& (params.lattice_type ==3 ||params.lattice_type ==6) ){
    interaction = readFile("/project/6005092/rfarid/workflow_2023/leblanc_codes/AMI/minimal_ami/loader/extended_U.txt");
	interaction_value = readFile1("/project/6005092/rfarid/workflow_2023/leblanc_codes/AMI/minimal_ami/loader/extended_U.txt",5);
    band_energy = {0,0};

    std::cout<< "Constructing AmiGraph object using seed: "<<seed <<" " <<std::endl;	
	std::string infile("ext_vars.dat");
	std::cout<<"Reading external parameters from ext_vars.dat"<<std::endl;
	g.ami.read_external(infile, extern_list);	
	std::cout<<"Attempting to load self-energy graphs from example_graphs"<<std::endl;
	int max=params.max_ord;
	g.read_ggmp(params.graph,ggm,max);
	possible_bands={1,2};
	std::cout<<"Completed read"<<std::endl;
	std::cout<<std::endl;
	g.ggm_label(ggm,0); 
	//g.print_all_edge_info(ggm[4][0].graph_vec[0]);
	
	std::cout<<"External parameters read are"<<std::endl;
    for(int i=0; i<extern_list.size();i++){
		std::cout<<extern_list[i].BETA_<<" "<<extern_list[i].MU_<<" "<< extern_list[i].H_<<" "<<extern_list[i].KDIM_<<" "<<extern_list[i].external_k_list_[0][0]<<" "<<extern_list[i].external_k_list_[0][1] <<" "<<extern_list[i].external_freq_[0]<<std::endl;
	
		}	

  }
  
  else if (params.molecular==1&& params.molecular_type ==1){
     interaction = readFile("/project/6005092/rfarid/workflow_2023/leblanc_codes/AMI/minimal_ami/loader/h2sto_U.txt");
	 interaction_value = readFile1("/project/6005092/rfarid/workflow_2023/leblanc_codes/AMI/minimal_ami/loader/h2sto_U.txt",5);
	 band_energy = readFile1("/project/6005092/rfarid/workflow_2023/leblanc_codes/AMI/minimal_ami/loader/h2sto_e0.txt",2);

	 

    std::cout<< "Constructing AmiGraph object using seed: "<<seed <<" " <<std::endl;		
	std::cout<<"Attempting to load self-energy graphs from example_graphs"<<std::endl;
	int max=params.max_ord;
	g.read_ggmp(params.graph,ggm, max);
	std::cout<<"Completed read"<<std::endl;
	std::cout<<std::endl;
	g.ggm_label(ggm,0);    
  }
  
    else if (params.molecular==1&& params.molecular_type ==2){
     interaction = readFile("/project/6005092/rfarid/workflow_2023/leblanc_codes/AMI/minimal_ami/loader/ccpdvz.txt");
	 interaction_value = readFile1("/project/6005092/rfarid/workflow_2023/leblanc_codes/AMI/minimal_ami/loader/ccpdvz.txt",5);
	 band_energy = readFile1("/project/6005092/rfarid/workflow_2023/leblanc_codes/AMI/minimal_ami/loader/ccpdvz_h.txt",2);
	

	
	std::cout<<"Attempting to load self-energy graphs from example_graphs"<<std::endl;
	int max=params.max_ord;
	g.read_ggmp(params.graph,ggm, max);
	std::cout<<"Completed read"<<std::endl;
	std::cout<<std::endl;
	g.ggm_label(ggm,0);    

  }
  else{
	  std::cout<<"Enter Correct params file";
	  MPI_Finalize();
	  return 0;}






if (params.molecular == 0){
std::cout<<"beginning lattice calculations"<<std::endl;
bool hf = params.hatree_fock; ///set to true if you are expecting a hatree or fock type interaction. True works for all type of graph so its set to default. Setting to false 
	//speeds up the process
	//construction 
mband mb(interaction,interaction_value,band_energy,hf);	
int  min_ord = params.min_ord;		
int  max_ord = params.max_ord;
 
std::vector<int> bandindex = {params.in, params.out};
std::vector<mband::sampler_collector> sigma_ToSum;
std::vector<AmiGraph::graph_t>sigma_FromGraph;
for (int i = min_ord; i < max_ord+1; ++i) {
    for (int j = 0; j < ggm[i].size(); ++j) {
        for (int k = 0; k < ggm[i][j].graph_vec.size(); ++k) {
			std::cout << "sampling graph with " << "o" << i << "_g" << j << "_n" << k <<std::endl;
            mband::sampler_collector sigma_collector;
			if (params.graph_type==0){
			std::cout <<"Sampling sigma with params  " << params.graph_type;
            mb.sigma_sampler(ggm[i][j].graph_vec[k], sigma_collector);}
			else if (params.graph_type ==1){
			        std::cout <<"Sampling particle particle with params  " << params.graph_type;
				mb.pp_sampler(ggm[i][j].graph_vec[k], sigma_collector,bandindex);}
			else if (params.graph_type ==2){
			        std::cout <<"Sampling particle hole with params  " << params.graph_type;
				mb.ph_sampler(ggm[i][j].graph_vec[k], sigma_collector,bandindex,possible_bands);}
			
			if (!sigma_collector.fermionic_edge_species.empty()){
				
				if ( params.mfreq_indp == 0 ){
			sigma_ToSum.push_back(sigma_collector);
			sigma_FromGraph.push_back(ggm[i][j].graph_vec[k]);
				}
				else if (params.mfreq_indp == 1 && !mb.check_mfreq_independent(sigma_collector.Alpha))  {
					sigma_ToSum.push_back(sigma_collector);
					sigma_FromGraph.push_back(ggm[i][j].graph_vec[k]);
				}
			
				
				
			}
		}
	}
}
std::cout <<"graphs have been sampled, beginning calculations" <<std::endl;
int lattice_type = params.lattice_type;
int samplesPerProcess = params.MC_num;
int numIterations = 0;
int totalSampleCount = 0;


std::vector<std::vector<std::vector<std::complex<double>>>> localSums(sigma_ToSum.size(), std::vector<std::vector<std::complex<double>>>(extern_list.size()));
std::vector<std::vector<std::vector<std::complex<double>>>> localSumSquared(sigma_ToSum.size(), std::vector<std::vector<std::complex<double>>>(extern_list.size()));
std::vector<std::vector<std::vector<std::complex<double>>>> globalSums(sigma_ToSum.size(), std::vector<std::vector<std::complex<double>>>(extern_list.size()));
std::vector<std::vector<std::vector<std::complex<double>>>> globalSumSquared(sigma_ToSum.size(), std::vector<std::vector<std::complex<double>>>(extern_list.size()));
std::vector<std::vector<std::vector<int>>> localSample(sigma_ToSum.size(), std::vector<std::vector<int>>(extern_list.size()));
std::vector<std::vector<std::vector<int>>> globalSample(sigma_ToSum.size(), std::vector<std::vector<int>>(extern_list.size()));


std::chrono::time_point<std::chrono::high_resolution_clock> startTime;
std::chrono::time_point<std::chrono::high_resolution_clock> currentTime;

if (rank == 0) {
    startTime = std::chrono::high_resolution_clock::now();
}

// Broadcast start time from rank 0 to all other ranks
MPI_Bcast(&startTime, sizeof(startTime), MPI_BYTE, 0, MPI_COMM_WORLD);

std::cout <<"caculations begins now\n";
while (std::chrono::duration_cast<std::chrono::seconds>(currentTime - startTime).count() < params.time) {
    for (int i = 0; i < sigma_ToSum.size(); i++) {
        for (int j = 0; j < extern_list.size(); j++) {
            int speciesSize = sigma_ToSum[i].fermionic_edge_species.size();
            globalSums[i][j].resize(speciesSize, std::complex<double>(0.0, 0.0));
            localSums[i][j].resize(speciesSize, std::complex<double>(0.0, 0.0));
            globalSumSquared[i][j].resize(speciesSize, std::complex<double>(0.0, 0.0));
            localSumSquared[i][j].resize(speciesSize, std::complex<double>(0.0, 0.0));
			localSample[i][j].resize(speciesSize, 0);
            globalSample[i][j].resize(speciesSize, 0);
            if (params.graph_type==0 && params.in == -1 && params.out == -1) {
                for (int k = 0; k < speciesSize; k++) {
					  std::tuple<std::complex<double>, std::complex<double>, int> result;
                        result = mb.lcalc_sampled_sigma(sigma_ToSum[i].graph, sigma_ToSum[i].Epsilon, sigma_ToSum[i].Alpha, sigma_ToSum[i].bosonic_Alpha, sigma_ToSum[i].Uindex[k], sigma_ToSum[i].fermionic_edge_species[k],
                                                        extern_list[j], samplesPerProcess, params);
                        localSums[i][j][k] += std::get<0>(result);
                        localSumSquared[i][j][k] += std::get<1>(result);
						localSample[i][j][k] += std::get<2>(result);     
                }
            } else if (params.graph_type==0) {
                for (int k = 0; k < speciesSize; k++) {
                    if (sigma_ToSum[i].external_line[k][0] == params.in && sigma_ToSum[i].external_line[k][1] == params.out) {
                        std::tuple<std::complex<double>, std::complex<double>, int> result;
                        result = mb.lcalc_sampled_sigma(sigma_ToSum[i].graph, sigma_ToSum[i].Epsilon, sigma_ToSum[i].Alpha, sigma_ToSum[i].bosonic_Alpha, sigma_ToSum[i].Uindex[k], sigma_ToSum[i].fermionic_edge_species[k],
                                                        extern_list[j], samplesPerProcess, params);
														
												       
                        localSums[i][j][k] += std::get<0>(result);
                        localSumSquared[i][j][k] += std::get<1>(result);
						localSample[i][j][k] += std::get<2>(result);	     
                    }
                }
            }
			
			else if (params.graph_type==1) {
                for (int k = 0; k < speciesSize; k++) {
					    
                        std::tuple<std::complex<double>, std::complex<double>, int> result;
                        result = mb.lcalc_sampled_pp(sigma_ToSum[i].graph, sigma_ToSum[i].Epsilon, sigma_ToSum[i].Alpha, sigma_ToSum[i].bosonic_Alpha,sigma_ToSum[i].gkkp_Alpha, sigma_ToSum[i].Uindex[k], sigma_ToSum[i].fermionic_edge_species[k],
                                                        extern_list[j], samplesPerProcess, params);				       
                        localSums[i][j][k] += std::get<0>(result);
                        localSumSquared[i][j][k] += std::get<1>(result);
						localSample[i][j][k] += std::get<2>(result);
					
                }
            }
			else if (params.graph_type==2) {
                for (int k = 0; k < speciesSize; k++) {
					    
                        std::tuple<std::complex<double>, std::complex<double>, int> result;
                        result = mb.lcalc_sampled_ph(sigma_ToSum[i].graph, sigma_ToSum[i].Epsilon, sigma_ToSum[i].Alpha, sigma_ToSum[i].bosonic_Alpha, sigma_ToSum[i].Uindex[k], sigma_ToSum[i].fermionic_edge_species[k],
                                                        extern_list[j], samplesPerProcess, params);				       
                        localSums[i][j][k] += std::get<0>(result);
                        localSumSquared[i][j][k] += std::get<1>(result);
						localSample[i][j][k] += std::get<2>(result);
					
                }
            }
			
			
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Reduce(localSums[i][j].data(), globalSums[i][j].data(), speciesSize * 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce(localSumSquared[i][j].data(), globalSumSquared[i][j].data(), speciesSize * 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(localSample[i][j].data(), globalSample[i][j].data(), speciesSize, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD);


        }
    }
 
    totalSampleCount += samplesPerProcess;
	MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0) {
        currentTime = std::chrono::high_resolution_clock::now();
    }

    // Broadcast current time from rank 0 to all other ranks
    MPI_Bcast(&currentTime, sizeof(currentTime), MPI_BYTE, 0, MPI_COMM_WORLD);
	std::cout<< "Time elapsed for rank " << rank <<" is " <<  std::chrono::duration_cast<std::chrono::seconds>(currentTime - startTime).count()<<std::endl;
   
}

// Compute the global totalSampleCount on rank 0

int globalTotalSampleCount;
MPI_Reduce(&totalSampleCount, &globalTotalSampleCount, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);


std::cout<<"sample of each process " << totalSampleCount;
// Print the global totalSampleCount on rank 0
if (rank == 0) {
std::cout << "Global Total Sample Count: " << globalTotalSampleCount << std::endl;
double U = 1;

std::ofstream outputFile;
if (params.lattice_type==1){
outputFile.open("output_1.txt");
}
else if (params.lattice_type==2){
outputFile.open("output_2.txt");
	
}
else if (params.lattice_type==4){
outputFile.open("output_4.txt");
	
}
else if (params.lattice_type==3){
outputFile.open("output_3.txt");
	
}
else if (params.lattice_type==5 ){
outputFile.open("output_5.txt");
}
else if (params.lattice_type==6){
outputFile.open("output_6.txt");
}
else if (params.lattice_type==7){
outputFile.open("output_7.txt");
}
else if (params.lattice_type==8){
outputFile.open("output_8.txt");

}

// Open the result file for writing

    if (outputFile.is_open()) {
        std::cout << "Number of processes: " << numProcesses << std::endl;
        outputFile << "#Ord "<<"Beta " << "Mu " <<"H " <<"Kx " <<"Ky " <<"Rfq "<<"Imfq " <<"in " <<"out " <<"Re " << "Re_err " <<"Imag " << "Imag_err " <<"U1 " <<"U2 " <<"U3 "<<"U4 " <<std::endl;
        for (int i = 0; i < sigma_ToSum.size(); i++) {
            for (int j = 0; j < extern_list.size(); j++) {
                int speciesSize = sigma_ToSum[i].fermionic_edge_species.size();
                for (int k = 0; k < speciesSize; k++) {
					double Ufactor = std::pow(U,g.graph_order(sigma_ToSum[i].graph));
					
                    std::complex<double> integral =globalSums[i][j][k] / static_cast<double>(globalSample[i][j][k] );
			        std::complex<double> integral_sq = globalSumSquared[i][j][k] ;
					double integral_err_real =  mb.Stdev(integral_sq.real(), integral.real(),globalSample[i][j][k]  );
					double integral_err_imag =  mb.Stdev(integral_sq.imag(), integral.imag(),globalSample[i][j][k]  );
					
					

                    std::cout << g.graph_order(sigma_ToSum[i].graph) << " ";
                    outputFile << g.graph_order(sigma_ToSum[i].graph) << " ";

                    std::cout << extern_list[j].BETA_ << " " << extern_list[j].MU_.real() << " "<<extern_list[j].MU_.imag() << " " << extern_list[j].H_ << " "
                              << extern_list[j].external_k_list_[0][0] << " " << extern_list[j].external_k_list_[0][1]
                              << " " << extern_list[j].external_freq_[0].real() << " " << extern_list[j].external_freq_[0].imag()<<" ";
                    outputFile << extern_list[j].BETA_ << " " << extern_list[j].MU_.real() << " "<< extern_list[j].MU_.imag() << " "<< extern_list[j].H_ << " "
                               << extern_list[j].external_k_list_[0][0] << " " << extern_list[j].external_k_list_[0][1]
                               << " " << extern_list[j].external_freq_[0].real() << " " << extern_list[j].external_freq_[0].imag()<<" ";

                    std::cout << sigma_ToSum[i].external_line[k][0] << " " << sigma_ToSum[i].external_line[k][1] << " ";
                    outputFile << sigma_ToSum[i].external_line[k][0] << " " << sigma_ToSum[i].external_line[k][1] << " ";

                    std::cout << integral.real() << " "<<integral_err_real<<" "<<integral.imag() << " "<<integral_err_imag<<" " ;
                    outputFile <<integral.real() << " "<<integral_err_real<<" "<< integral.imag() << " "<<integral_err_imag<<" " ;

                    for (int l = 0; l < sigma_ToSum[i].Uindex[k].size(); l++) {
                        std::cout << sigma_ToSum[i].Uindex[k][l] << " ";
                        outputFile << sigma_ToSum[i].Uindex[k][l] << " ";
                    }
					int remainingElements = 4 - sigma_ToSum[i].Uindex[k].size();
						for (int l = 0; l < remainingElements; l++) {
							std::cout << -1 << " ";
							outputFile << -1 << " ";
														}
                    std::cout << i << " "<<globalSample[i][j][k];
                    std::cout << std::endl;
                    outputFile << std::endl;
                }
            }
        }

        outputFile.close(); // Close the result file
	}
	
	else {
        std::cout << "Error: Failed to open the result file." << std::endl;
    }
std::cout << "molecular: " << params.molecular << std::endl;
std::cout << "E_reg: " << params.E_reg << std::endl;
std::cout << "lattice_type: " << params.lattice_type << std::endl;
std::cout << "set_precision: " << params.set_precision << std::endl;
std::cout << "hatree_fock: " << std::boolalpha << params.hatree_fock << std::endl;
std::cout << "tp: " << params.tp << std::endl;
std::cout << "tpp: " << params.tpp << std::endl;
std::cout << "tperp: " << params.tperp << std::endl;
std::cout << "tbs: " << params.tbs << std::endl;
std::cout << "SOC: " << params.SOC << std::endl;
std::cout << "t_orb: " << params.t_orb << std::endl;
std::cout << "Mc_num" << params.MC_num << std::endl;
std::cout << "max: " << params.max_ord << std::endl;
std::cout << "min_ord" << params.min_ord << std::endl;
std::cout << "line1 " << params.in<<" " <<params.out<< std::endl;
std::cout << "V " << params.V<<" "  <<std::endl;
std::cout << "time" << params.time<<" " << std::endl;
std::cout << "G_FUNC " << params.G_FUNC<<" " << std::endl;
std::cout << "Global Total Sample Count: " << globalTotalSampleCount << std::endl;
std::cout << "sample per process: " << samplesPerProcess<< std::endl;
	}
	


}
else if (params.molecular == 1){
std::cout<<"beginning molecular calculations"<<std::endl;
bool hf = params.hatree_fock; ///set to true if you are expecting a hatree or fock type interaction. True works for all type of graph so its set to default. Setting to false 
	//speeds up the process
	//construction 
mband mb(interaction,interaction_value,band_energy,hf);
print1d(band_energy);
print1d(interaction_value);	
std::cout<<params.molec_beta << " " <<params.molec_mfreq<<std::endl;

int  min_ord = params.min_ord;		
int  max_ord = params.max_ord;

std::vector<double> Beta_ext_vec = {params.molec_beta};
std::vector<double> Mfreq_ext_vec;
for (int i =0;i<params.molec_mfreq;i++){Mfreq_ext_vec.push_back(i);}
std::vector<int> line = {params.in,params.out};


 int totalGraphs = 0;
    for (int i = min_ord; i <= max_ord; ++i) {
        for (int j = 0; j < ggm[i].size(); ++j) {
            totalGraphs += ggm[i][j].graph_vec.size();
        }
    }

	// Calculate the workload for each process
    int graphsPerProcess = totalGraphs / numProcesses;
    int remainder = totalGraphs % numProcesses;
    int startGraph = rank * graphsPerProcess;
    int endGraph = startGraph + graphsPerProcess - 1;

    if (rank < remainder) {
        // Assign one extra graph to the remaining processes
        startGraph += rank;
        endGraph += rank + 1;
    } else {
        // Distribute the remaining graphs among the processes
        startGraph += remainder;
        endGraph += remainder;
    }
 

    int currentGraphCount = 0;
   
    for (int i = min_ord; i <= max_ord; ++i) {
        for (int j = 0; j < ggm[i].size(); ++j) {
            for (int k = 0; k < ggm[i][j].graph_vec.size(); ++k) {
                if (currentGraphCount >= startGraph && currentGraphCount <= endGraph) {
                    // Process only the relevant workload
                    mband::output_collector output_collector;
					if (line[0]==-1 && line[1]==-1){
					std::string filename;
					if (params.molecular_type==1){
                    filename = "h2_output_p" + std::to_string(rank)+ "_o" + std::to_string(i) + "_g" + std::to_string(j) +"_n" + std::to_string(k) + ".txt";
					}
					else if (params.molecular_type==2){
                    filename = "h10_output_p" + std::to_string(rank)+ "_o" + std::to_string(i) + "_g" + std::to_string(j) +"_n" + std::to_string(k) + ".txt";
					}
					mb.molecular_solver(ggm[i][j].graph_vec[k], output_collector, Beta_ext_vec, Mfreq_ext_vec,filename);		
                    
					}
					else{
					std::string filename;
                    if (params.molecular_type==1){
                    filename = "h2_output_p" + std::to_string(rank)+ "_o" + std::to_string(i) + "_g" + std::to_string(j) +"_n" + std::to_string(k) + ".txt";
					}
					else if (params.molecular_type==2){
                    filename = "h10_output_p" + std::to_string(rank)+ "_o" + std::to_string(i) + "_g" + std::to_string(j) +"_n" + std::to_string(k) + ".txt";
					}
					mb.molecular_solver_ext(ggm[i][j].graph_vec[k], output_collector, Beta_ext_vec, Mfreq_ext_vec,line,filename);	
					}
                    //mb.write_output(filename.str(), output_collector, Beta_ext_vec, Mfreq_ext_vec);
                }

                ++currentGraphCount;
            }
        }
    }
	
	

	
}


MPI_Finalize();
return 0;

//////////////////////////////////sampling ends////////////////////////


///////////////////////////////using sampler,collecting and evaluating///////////////////////////////////
/*
std::vector<double> Beta_ext_vec = {50};
std::vector<double> Mfreq_ext_vec;
for (int i =0;i<100;i++){Mfreq_ext_vec.push_back(i);}
std::vector<int> line={1,1};
std::vector<int> line1={2,2};
*/
/*
for (int i = min_ord; i < max_ord+1; ++i) {
    for (int j = 0; j < ggm[i].size(); ++j) {
        for (int k = 0; k < ggm[i][j].graph_vec.size(); ++k) {
            mband::sampler_collector sigma_collector;
            mb.sigma_sampler(ggm[i][j].graph_vec[k], sigma_collector);
		    sigma_ToSum.push_back(sigma_collector);
			sigma_FromGraph.push_back(ggm[i][j].graph_vec[k]);
			if(!sigma_collector.fermionic_edge_species.empty()){
				mband::output_collector output_collector;
				mb.calculate_sampled_sigma_ext(ggm[i][j].graph_vec[k], sigma_collector,output_collector,Beta_ext_vec,Mfreq_ext_vec,line);
				std::string filename = "10bandH2_o"+std::to_string(i)+"_g" +std::to_string(j)+"_n" +std::to_string(k)+"11.txt";
				mb.write_output(filename,output_collector,Beta_ext_vec,Mfreq_ext_vec);				
			}		
        }
    }
}

*/

///////////////////////////////////Molecular stuff on hold???//////////////////data  collected already/////////////////////////
/*
for (int i = min_ord; i < max_ord+1; ++i) {
    for (int j = 0; j < ggm[i].size(); ++j) {
        for (int k = 0; k < ggm[i][j].graph_vec.size(); ++k) {
			mband::output_collector output_collector;
        	mb.molecular_solver_ext(ggm[i][j].graph_vec[k], output_collector, Beta_ext_vec, Mfreq_ext_vec,line);
            std::string filename = "2band_o"+std::to_string(i)+"_g" +std::to_string(j)+"_n" +std::to_string(k)+"11.txt";
		    mb.write_output(filename,output_collector,Beta_ext_vec,Mfreq_ext_vec);				
			}		
        }
    }

	
for (int i = min_ord; i < max_ord+1; ++i) {
    for (int j = 0; j < ggm[i].size(); ++j) {
        for (int k = 0; k < ggm[i][j].graph_vec.size(); ++k) {
			mband::output_collector output_collector;
        	mb.molecular_solver_ext(ggm[i][j].graph_vec[k], output_collector, Beta_ext_vec, Mfreq_ext_vec,line1);
            std::string filename = "2band_o"+std::to_string(i)+"_g" +std::to_string(j)+"_n" +std::to_string(k)+"22.txt";
		    mb.write_output(filename,output_collector,Beta_ext_vec,Mfreq_ext_vec);				
			}		
        }
    }
*/
//////////////////////////for a single graph///////////////
/*
mband::output_collector output_collector;
mb.molecular_solver_ext(ggm[3][0].graph_vec[10], output_collector, Beta_ext_vec, Mfreq_ext_vec,line);
std::string filename = "10band_o"+std::to_string(3)+"_g" +std::to_string(0)+"_n" +std::to_string(10)+"11.dat";
mb.write_output(filename,output_collector,Beta_ext_vec,Mfreq_ext_vec);

mband::output_collector output_collector1;
mb.molecular_solver_ext(ggm[3][0].graph_vec[11], output_collector1, Beta_ext_vec, Mfreq_ext_vec,line);
std::string filename1 = "10band_o"+std::to_string(3)+"_g" +std::to_string(0)+"_n" +std::to_string(11)+"11.dat";
mb.write_output(filename1,output_collector1,Beta_ext_vec,Mfreq_ext_vec);
*/
//////molecular stuff ends

}


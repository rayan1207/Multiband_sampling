//=======================================================================
// Copyright 2018 JPF LeBlanc
//=======================================================================
#include "amigraph.hpp"
#include <ctime>
#include <Eigen/Dense>
// #include "cuba.h"
// #include "mpi.h"


// NOte, 3rd order translation is ( 4,2,1,3,5,6)

int main(int argc, char *argv[])
{
  
int mpi_rank=0;

int comm_size=1;

    // MPI_Init(&argc, &argv);
    // MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    // MPI_Comm_size(MPI_COMM_WORLD, &comm_size);	
  
  
  int max=5;
  int min=2;
  int maxeval=10000;
  int seed=0;
  int first=0;
  int second=1;
  int ord=2;
  
  //int seed=0; // but this can be set to anything
auto now =std::chrono::high_resolution_clock::now();
seed = std::chrono::duration_cast<std::chrono::microseconds>(now.time_since_epoch()).count();
//int seed=0;

  
 if(argc==1){ 
 printf("\nNo Extra Command Line Argument Passed Other Than Program Name"); }
    if(argc==2) 
    { 
        seed=atoi(argv[1]);
    } 
if(argc==4){
	
	ord=atoi(argv[1]);	
first=atoi(argv[2]);
second=atoi(argv[3]);

max=ord;

}	

if(argc==5){
ord=atoi(argv[1]);	
first=atoi(argv[2]);
second=atoi(argv[3]);
seed=atoi(argv[4]);
max=ord;
}	
	
  
 
  

std::cout<<"Seed is "<< seed <<" " <<std::endl;
AmiGraph g(AmiCalc::Bare, seed);//g(AmiCalc::Pi,seed);//g(AmiCalc::Hartree,seed);



AmiGraph::vertex_vector_list_t bubble_vertex_list;
AmiGraph::edge_vector_list_t bubble_edge_list;
AmiGraph::edge_vector_list_t legs_edge_list;
g.bubble_finder(g.current_graph, bubble_vertex_list,  bubble_edge_list, legs_edge_list);

AmiGraph::edge_t one, two;




	
		std::cout<<"Reading external parameters from ext_vars.dat"<<std::endl;
		AmiCalc::external_variable_list extern_list;

		std::string infile("ext_vars.dat");
		g.ami.read_external(infile, extern_list);	
	
	
std::vector< std::vector< AmiGraph::graph_t>> graph_matrix;
g.generate_graphs(graph_matrix,max,mpi_rank);

bool success;


g.label_graphs(graph_matrix,1,max);
// g.min_ext_label_counts(graph_matrix);
g.reduce_gm_rf(graph_matrix,mpi_rank);
g.reduce_gm_tp(graph_matrix, mpi_rank);
g.reduce_gm_oneleg(graph_matrix, mpi_rank);
g.reduce_gm_1PBose(graph_matrix, mpi_rank);

AmiGraph::gg_matrix_t ggm;
ggm.resize(graph_matrix.size());




bool result=0;
double prefactor=1;

g.print_all_edge_info(graph_matrix[ord][second]);

// Timing testing 
clock_t total;
total=clock();

AmiGraph::graph_t g2_out;

AmiGraph::git_perm_set_t pst;

std::cout<<"Prefactor Before "<< g.get_prefactor(graph_matrix[ord][second],ord)<<std::endl;

for(int i=0; i<1000; i++){

bool lab=true;
g.repeated_labelling(graph_matrix[ord][second], lab);
// g.repeated_labelling
g.git_compare_v2(graph_matrix[ord][first], graph_matrix[ord][second], result, prefactor, pst, g2_out);

if(result){
std::cout<<"These match with prefactor "<< prefactor <<" and took attempts="<<i<<std::endl;	

// optional
graph_matrix[ord][second]=g2_out;
	break;}

}



std::cout<<"Prefactor after "<< g.get_prefactor(graph_matrix[ord][second],ord)<<std::endl;


g.gm_to_ggm(graph_matrix, ggm);

// g.git_compare(graph_matrix[ord][first], graph_matrix[ord][second], result, prefactor);

total=clock()-total;

 float tim=(float)total/(float)CLOCKS_PER_SEC; ///(float)top;
std::cout<<"It took "<< tim <<" to do the whole loop per attempt clocks="<< (float)CLOCKS_PER_SEC<<std::endl;
std::cout<<"This could be done "<< 1.0/tim*60.0 <<" times per minute :)"<<std::endl;

// g.print_all_edge_info(graph_matrix[ord][first]);
// std::cout<<std::endl;
// g.print_all_edge_info(graph_matrix[ord][second]);






g.ggm_construct_ami_sol(ggm, 1e-8, mpi_rank);

AmiCalc::gg_solution_set_matrix_t GG_AMI_MATRIX;
g.ggm_to_ggamim(ggm, GG_AMI_MATRIX); 




// std::cout<<"Testing git shifts"<<std::endl;

int orde=ord;
int dim=2;

    AmiCalc::internal_state state(orde,dim);
	AmiCalc::internal_state shifted;
	
double s1, s2, si1r, si1i, si2r, si2i;
s1=0;
s2=0;
si1r=0;
si1i=0;
si2i=0;
si2r=0;
double s3=0;
double s4=0;
double s5=0;
double s6=0;
int count=10000;	
	
// for(int i=0;i<count; i++){	
	g.zero_state(state,orde);
	g.randomize_state(state,orde);
	
	AmiCalc::k_vector_t k_shift(dim,M_PI);
	std::cout<<"K-shift vector is "<< k_shift[0]<<" "<<k_shift[1]<<std::endl;
// define some internal state to evaluate
	 // for(int i=0; i< state.internal_k_list_.size(); i++){
		
		// state.internal_k_list_[i][0]=M_PI/20.1+.51*i+M_PI;
		// state.internal_k_list_[i][1]==M_PI/10.1+.51*i;
			
	// }

// Permute the state 
shifted=state;	

g.permute_state(shifted,k_shift, pst);

// g.shift_state(shifted,index,k_shift);	
 // g.shift_state(shifted,0,k_shift);
 
// std::swap(shifted.internal_k_list_[1],shifted.internal_k_list_[2]); 
 
// g.shift_state(shifted,1,k_shift); 
// shifted.internal_k_list_[1][0]=-shifted.internal_k_list_[1][0];
// shifted.internal_k_list_[1][1]=-shifted.internal_k_list_[1][1];
// g.shift_state(shifted,1,k_shift);


// g.print_state(state);
// std::cout<<std::endl;
// g.print_state(shifted);	

double p1, p2;

p1=GG_AMI_MATRIX[ord][first][0].prefactor_;
p2=GG_AMI_MATRIX[ord][second][0].prefactor_;

// std::cout<<"Prefactors of each graph were "<<GG_AMI_MATRIX[ord][first][0].prefactor_<<" and "<<GG_AMI_MATRIX[ord][second][0].prefactor_<<std::endl;

	
    AmiCalc::ami_vars_list vars_list;
	g.ami.construct_ami_vars_list(GG_AMI_MATRIX[ord][first][0].R0_,GG_AMI_MATRIX[ord][first][0].prefactor_, state, extern_list, vars_list);

	std::vector<double> Re_results, Im_results;
	// actual evaluation line - answer is in results 
	//graph.ami.evaluate_solutions(results, AMI_MATRIX[ord][num], vars_list);
	g.ami.evaluate_solutions(Re_results,Im_results, GG_AMI_MATRIX[ord][first][0], vars_list);
	
	
// then evaluate the shifted case 	
	g.ami.construct_ami_vars_list(GG_AMI_MATRIX[ord][second][0].R0_,GG_AMI_MATRIX[ord][second][0].prefactor_, shifted, extern_list, vars_list);

	std::vector<double> Re_shifted, Im_shifted;
	// actual evaluation line - answer is in results 
	//graph.ami.evaluate_solutions(results, AMI_MATRIX[ord][num], vars_list);
	g.ami.evaluate_solutions(Re_shifted, Im_shifted, GG_AMI_MATRIX[ord][second][0], vars_list);
	

// if(sgn(Re_shifted[0]/Re_results[0])>0){ s1++;}else{s3++;}
// if(sgn(Im_shifted[0]/Im_results[0])>0){ s2++;}else{s4++;}
// if(sgn(Re_shifted[0]+Re_results[0])== sgn(Re_results[0])){s5++;}else{s6++;}

// s1+=sgn(Re_shifted[0]/Re_results[0]);
// s2+=sgn(Im_shifted[0]/Im_results[0]);

// s3+=sgn(-Re_shifted[0]/Re_results[0]);
// s4+=sgn(-Im_shifted[0]/Im_results[0]);

// si1r+=sgn(Re_shifted[0]);
// si2r+=sgn(Re_results[0]);

// si1i+=sgn(Im_shifted[0]);
// si2i+=sgn(Im_results[0]);



for(int i=0; i< Re_results.size(); i++){	
std::cout<<"First gave "<<Re_results[i]<<" and "<<Im_results[i]<<std::endl;
std::cout<<"Second gave "<<Re_shifted[i]<<" and "<<Im_shifted[i]<<std::endl;	
}



if(result){
std::cout<<"Was Success with prefactor "<<prefactor*p1*p2<<std::endl;
}

AmiGraph::labels_t L;
bool work=false;

// g.sys_label_sets(L, graph_matrix[ord][first] , work, 10000);

// std::cout<<std::endl;
// g.print_all_edge_info(graph_matrix[ord][first]);
// std::cout<<"Labels sets of size "<< L.size()<<std::endl;

bool scramble=false;
std::cout<<"Testing construct label sets"<<std::endl;
g.construct_label_sets(ggm, mpi_rank, 250, 100);
std::cout<<"Testing label sets"<<std::endl;
g.check_label_sets(ggm, mpi_rank);
std::cout<<"Git reduce"<<std::endl;
g.reduce_ggm_git(ggm,mpi_rank, 15);


return 0;
}


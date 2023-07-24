//=======================================================================
// Copyright 2018 JPF LeBlanc
//=======================================================================
#include "amigraph.hpp"
#include <ctime>
#include <Eigen/Dense>

#include <string> 
#include <sstream> 
#include <fstream>
#include <iostream>

// in amical simple_residue a lot of time is spent allocating memory? worth fixing?


// NOte, 3rd order translation is ( 4,2,1,3,5,6)

int main(int argc, char *argv[])
{
  int max=0;
  int write=0;
  int read=0;
  int filter=0;
  int reduce=0;
  if(argc==1) 
        printf("\nNo Extra Command Line Argument Passed Other Than Program Name"); 
    if(argc>=2) 
    { 
        max=atoi(argv[1]);
		
    }  
	if(argc==4){
		max=atoi(argv[1]);
		write=atoi(argv[2]);
		read=atoi(argv[3]);
	}
	if(argc==5){
		max=atoi(argv[1]);
		write=atoi(argv[2]);
		read=atoi(argv[3]);
		filter=atoi(argv[4]);
	}
	if(argc==6){
		max=atoi(argv[1]);
		write=atoi(argv[2]);
		read=atoi(argv[3]);
		filter=atoi(argv[4]);
		reduce=atoi(argv[5]);
	}
  
  
//int seed=0; // but this can be set to anything
auto now =std::chrono::high_resolution_clock::now();
auto seed = std::chrono::duration_cast<std::chrono::microseconds>(now.time_since_epoch()).count();
//int seed=0;

std::cout<< seed <<" " <<std::endl;
AmiGraph g(AmiCalc::Bare, seed);//g(AmiCalc::Pi,seed);//g(AmiCalc::Hartree,seed);

std::cout<<"Constructor finished"<<std::endl;


AmiGraph::vertex_vector_list_t bubble_vertex_list;
AmiGraph::edge_vector_list_t bubble_edge_list;
AmiGraph::edge_vector_list_t legs_edge_list;
g.bubble_finder(g.current_graph, bubble_vertex_list,  bubble_edge_list, legs_edge_list);

AmiGraph::edge_t one, two;

int mpi_rank=0;

// std::vector<double> test=g.sobol_random_real();
 
// for(int j=0; j< 20; j++){ 
 // test=g.sobol_random_real();
// for (int i=0; i< test.size(); i++){
// std::cout<<test[i]<<" ";
// }	
// std::cout<<std::endl;
// }

// This is all debugging stuff really for generating high order graphs 
AmiGraph::gg_matrix_t preggm;
AmiGraph::gg_matrix_t ggm;

if(max!=0){
std::vector< std::vector< AmiGraph::graph_t>> graph_matrix;
if(read==0){
g.generate_graphs(graph_matrix,max,mpi_rank);
}
if(read==1){
g.read_ggmp("ggm_paired",preggm, max);
g.print_ggm(preggm);
// g.ggm_remove_pairs(preggm);
// g.print_ggm(preggm);
g.ggm_to_gm(preggm,graph_matrix,1);
}


///



if(filter){

std::cout<<"Labelling"<<std::endl;
g.reduce_gm_tp(graph_matrix, mpi_rank,1);
g.label_graphs(graph_matrix,1,max);
// g.min_ext_label_counts(graph_matrix);
// g.reduce_gm_rf(graph_matrix,mpi_rank);
g.reduce_gm_oneleg(graph_matrix, mpi_rank);
g.reduce_gm_1PBose(graph_matrix, mpi_rank,2);


ggm.resize(graph_matrix.size());
}


// 99% sure we only need to number the vertices up to this point. maybe fix this?
// g.number_vertices(Each Graph);
g.label_graphs_sys(graph_matrix,1,max);
g.gm_to_ggm(graph_matrix, ggm);

bool success;
if (write==1){
g.write_ggm("ggm_reduced",ggm);}

if(read && !filter){
ggm=preggm;
}	

// g.reduce_gm_skeleton(graph_matrix,mpi_rank);

// Dec 4 - put this back. unsure if it still works
// std::vector< std::vector< int>> group_lists;
// g.reorder_AID(graph_matrix, group_lists);

// for(int i=0; i<group_lists[max].size();i++){
	// std::cout<<group_lists[max][i]<<std::endl;
// }


// g.gos(graph_matrix, ggm);
//g.gos(graph_matrix, ggm, group_lists);



// std::string testfile="test.graph";
// g.graph_write(testfile, ggm[2][0].graph_vec[0]);
// AmiGraph::graph_t readin;
// g.graph_read(testfile, readin);

// std::cout<<"Original graph was "<< std::endl;
// g.print_all_edge_info(ggm[2][0].graph_vec[0]);
// std::cout<<"Readin is"<<std::endl;
// g.print_all_edge_info(readin);

// g.write_ggm("ggm",ggm);

// Preconstruct labels

// g.construct_label_sets(ggm, mpi_rank);
// g.print_ggm(ggm);

// First round ggm reduction using AID pair info 
// g.ggm_AID_reduce(4, ggm,group_lists);
// g.print_ggm(ggm);

// Continued reduction given existing labels 
// g.reduce_ggm_gos(ggm);
// g.print_ggm(ggm);

std::cout<<"Compare sizes "<< graph_matrix[0].size()<<" "<< ggm[0].size()<<std::endl;

bool scramble=true;

if(reduce){
g.construct_label_sets(ggm, mpi_rank, 1, 10, 10);
// g.trim_label_sets(ggm,mpi_rank,200);
g.check_label_sets(ggm, mpi_rank,1);

// g.reduce_ggm_git_odd_random(ggm,mpi_rank,2);
g.reduce_ggm_git(ggm,mpi_rank, 100);

}

// g.reduce_ggm_gos(ggm,mpi_rank, 80);
// g.reduce_ggm_gos(ggm,mpi_rank, 200);

g.write_ggm("ggm_reduced",ggm);
// g.read_ggm("ggm",ggm, max);
g.print_ggm(ggm);

// g.print_all_edge_info(ggm[2][0].graph_vec[0]);
// std::cout<<"Test systematic"<<std::endl;
// std::cout<<num_vertices(ggm[2][0].graph_vec[0])<<std::endl;
// g.systematic_vertex_label(ggm[2][0].graph_vec[0]);

//TODO: What is the status of the graph_vec graphs? Are they modified by reduce_ggm_git? or no?

g.ggm_label(ggm,1); // don't think this labelling is necessary


g.construct_label_sets(ggm, mpi_rank,1, 100, 100);
g.ggm_populate_git_pairs(ggm,mpi_rank,1);
g.print_ggm(ggm);
g.ggm_construct_ami_sol(ggm, 1e-8, mpi_rank);

g.write_ggmp("ggm_paired", ggm, 1);
// g.read_ggmp("ggm_paired",ggm,max);


AmiCalc::gg_solution_set_matrix_t GG_AMI_MATRIX;
g.ggm_to_ggamim(ggm, GG_AMI_MATRIX,1); 


// g.reduce_ggm_gos(ggm);
//g.print_ggm(ggm);

//AmiGraph::ggss_t ggss; // graph_group solutions set 

// g.ggm_AID_reduce(4, ggm,group_lists);
// g.print_ggm(ggm);
// g.reduce_ggm_gos(ggm);
// g.print_ggm(ggm);

// g.gg_check(ggm);
// g.print_ggm(ggm);
// g.reduce_ggm_gos(ggm);
// g.print_ggm(ggm);

/*  for(int i=0; i< ggm[3].size(); i++){
for(int j=0; j< ggm[3][i].graph_vec.size(); j++){	
	std::cout<<std::endl;
	std::cout<<"Printing graph "<< i<< std::endl;
	std::cout<<std::endl;
// g.print_all_edge_info(ggm[3][i].graph_vec[j]);	

g.ami.print_final(3, ggm[3][i].ss_vec[j].R_, ggm[3][i].ss_vec[j].P_, ggm[3][i].ss_vec[j].S_);
	std::cout<<std::endl;
	std::cout<<ggm[3][i].ss_vec[j].R_[3].size()<<std::endl;
}
}  */

//// Read External Variables////
std::string infile="ext_vars.dat";
std::cout<<"Reading external parameters from "<< infile<<std::endl;
AmiCalc::external_variable_list extern_list;
g.ami.read_external(infile, extern_list);

int orde=4;
int dim=2;

    AmiCalc::internal_state state(orde,dim);
	g.zero_state(state,orde);
	g.randomize_state(state,orde);
	AmiCalc::k_vector_t kx_shift(dim,0);
	AmiCalc::k_vector_t ky_shift(dim,0);
	AmiCalc::k_vector_t minus_k_shift(dim,0);
	int grid=50;
	double length=2*M_PI/double(grid);
	kx_shift[0]=length;	
	kx_shift[1]=0;//length;
	ky_shift[0]=0.0;
	ky_shift[1]=length;
	
	minus_k_shift[0]=-length;	
	minus_k_shift[1]=-length;
	
std::ofstream file;
file.open("kxkyfile.dat",  std::ofstream::out); // | std::ofstream::app);

	
	

	 for(int i=0; i< state.internal_k_list_.size(); i++){
		
		state.internal_k_list_[i][0]=M_PI/20.1+.51*i;
		state.internal_k_list_[i][1]==M_PI/10.1+.51*i;
			
	}
	//210; j++){
		
	for(int y=0; y<grid; y++){	
		
		
		
	for(int x=0; x< grid; x++){
		
		
		// g.zero_state(state,orde);
		
		for(int i=0; i< state.internal_k_list_.size(); i++){
		// std::cout<<"Setting k for "<<i<<std::endl;
		
		if(i==0){
		state.internal_k_list_[i][0]=0.0;
		state.internal_k_list_[i][1]=0.0;}
		// else{
			// state.internal_k_list_[i][0]=g.random_real(2*M_PI);
			// state.internal_k_list_[i][1]=g.random_real(2*M_PI);
		// }
	
	}	
		
		for(int yshift=0; yshift<y ; yshift++){
			g.shift_state(state,0,ky_shift);
		}
		for(int xshift=0; xshift<x ; xshift++){
			g.shift_state(state,0,kx_shift);
		}
		
		// std::cout<<"State Dimension is "<< state.dim_<<std::endl;
		
		
		// if(j==0){
		// g.shift_state(state,0,k_shift);
		// }
		// if(j==1){
		// g.shift_state(state,0,minus_k_shift);	
		// }
	
    // g.shift_state(state,0,k_shift);
	// g.shift_state(state,1,k_shift);
	// g.shift_state(state,2,k_shift);
	// g.shift_state(state,3,k_shift);
	// g.shift_state_pi(state,4);

    AmiCalc::ami_vars_list vars_list;
	g.ami.construct_ami_vars_list(GG_AMI_MATRIX[4][0][0].R0_,GG_AMI_MATRIX[4][0][0].prefactor_, state, extern_list, vars_list);

	std::vector<double> Re_results, Im_results;
	// actual evaluation line - answer is in results 
	//graph.ami.evaluate_solutions(results, AMI_MATRIX[ord][num], vars_list);
	g.ami.evaluate_solutions(Re_results,Im_results, GG_AMI_MATRIX[4][0][0], vars_list);

	
	/* std::complex<double> calc_result;
	for(int run=0; run<200000; run++){
	
	calc_result=g.ami.evaluate(test_amiparms,ggm[4][5].ss_vec[j].R_, ggm[4][5].ss_vec[j].P_, ggm[4][5].ss_vec[j].S_,  avars);
	
	} */
	
	// std::cout<<"Evaluated to "<< calc_result<<std::endl;
for(int i=0; i<Re_results.size(); i++){
	
file<<x<<" "<< y<<" "<< state.internal_k_list_[0][0]<<" "<< state.internal_k_list_[0][1] <<" "<< Re_results[i]<<" "<< Im_results[i]<<std::endl;	
}
	
	}
	}

 // g.gos(graph_matrix, ggm);

//std::cout<<"graph matrix size is "<< graph_matrix.size()<<std::endl;

// g.gos(graph_matrix, ggm);


file.close();

}


AmiCalc::g_prod_t r0=g.ami.construct_multipole_example();
AmiCalc::ami_vars avars=g.ami.construct_4ord_ext_multipole_example();//construct_ext_example_Y();//construct_ext_multipole_example();

// Timing testing 
// clock_t total;
// int top=10000;

// total=clock();
 // for(int i=0; i< top; i++)
// {
AmiCalc::R_t R_array;
AmiCalc::P_t P_array;
AmiCalc::S_t S_array;

AmiCalc::ami_parms test_amiparms(4, .00001);


// g.ami.minimal_construct(test_amiparms, r0, R_array, P_array, S_array);
g.ami.construct(test_amiparms, r0, R_array, P_array, S_array);
// g.ami.print_final(4, R_array, P_array, S_array);
std::complex<double> calc_result=g.ami.evaluate(test_amiparms,R_array, P_array, S_array,  avars);

// std::cout<<"Evaluated to "<< calc_result<<std::endl;
// g.ami.print_final(3, R_array, P_array, S_array);



// }

 // total=clock()-total;

 

// float tim=(float)total/(float)CLOCKS_PER_SEC/(float)top;
// std::cout<<"It took "<<total/(float)CLOCKS_PER_SEC<<" and "<< tim <<" to do the whole loop per attempt clocks="<< (float)CLOCKS_PER_SEC<<std::endl;
// std::cout<<"This could be done "<< 1.0/tim*60.0 <<" times per minute :)"<<std::endl;


// g.ami.print_final(2, R_array,P_array,S_array);

// std::cout<<"There are "<< R_array[4].size()<< " terms"<<std::endl;



// for(int j=0;j< graph_matrix.size(); j++){
// for(int i=0; i<graph_matrix[j].size();i++){
	// std::cout<<i<<std::endl;
// g.number_vertices(graph_matrix[j][i]);
// g.print_all_edge_info(graph_matrix[j][i]);


// }
// }

// these I thought were not both in the set of isomorphic graphs...


////////####################

/* 
std::vector< std::vector< AmiGraph::graph_t>> graph_sorted;
// resize the order container 
graph_sorted.resize(10);

for (int k=0; k<= max; k++){
for(int i=0; i< graph_matrix[k].size(); i++){

if(i % 10000 ==0){
std::cout<<"On graph "<< i <<" sorted is now "<< graph_sorted[k].size()<<std::endl;
}
bool add=true;
for(int j=0; j<graph_sorted[k].size(); j++){

if(g.is_isomorphic(graph_matrix[k][i], graph_sorted[k][j])){ add=false; continue;}	
	
}

if(add){ graph_sorted[k].push_back(graph_matrix[k][i]);}

}

std::cout<<"Remaining graphs were "<<k<<" "<< graph_sorted[k].size()<<std::endl;
}
 */


////////####################

// for( int k=0; k<max;k++){
// std::cout<<"Remaining graphs were "<<k<<" "<< graph_sorted[k].size()<<std::endl;

// }
// for(int j=0;j< graph_sorted.size(); j++){
// for(int i=0; i<graph_sorted[j].size();i++){
	// std::cout<<i<<std::endl;
// g.number_vertices(graph_sorted[j][i]);
// g.print_all_edge_info(graph_sorted[j][i]);


// }
// }

// std::cout<<"Principle linne is length  "<<g.get_length_principle_line(graph_matrix[3][18])<<std::endl;
// std::cout<<g.is_isomorphic(graph_matrix[3][18],graph_matrix[3][19])<<std::endl;


// Timing testing 
// clock_t total;
// int max=100000;

// total=clock();
 // for(int i=0; i< max; i++){
	
	 // g.is_isomorphic(graph_matrix[3][18],graph_matrix[3][19]);
 // }
// total=clock()-total;
// std::cout<<"It took "<< (float)total/(float)CLOCKS_PER_SEC/(float)max<<" to check isomorphic per attempt"<<std::endl;


/* 
AmiGraph::barcode_t B1, B2;
AmiGraph::tree_t T1, T2;
AmiGraph::vertex_vector_t V1,V2;
// std::cout<<"Constructing barcodes"<<std::endl;

 g.symmetrize_bosonic_lines(graph_matrix[3][12]);
 g.symmetrize_bosonic_lines(graph_matrix[3][13]);
g.construct_graph_barcode(graph_matrix[3][12], B1, T1);
g.construct_graph_barcode(graph_matrix[3][13], B2, T2);

g.label_via_tree(T1,graph_matrix[3][12]);
g.label_via_tree(T2,graph_matrix[3][13]);

g.compare_trees(T1, graph_matrix[3][12], T2, graph_matrix[3][13]);
  */
 
 
 
 

 


// std::cout<<std::endl;
// std::cout<<g.is_1P_reducible_graph(graph_matrix[3][5])<<std::endl;

	


// std::vector< std::vector< std::vector<int> >> ForB1, ForB2;
//g.symmetrize_bosonic_lines(graph_matrix[2][1]);
//g.symmetrize_bosonic_lines(graph_matrix[2][2]);
// std::cout<<std::endl;
// g.print_all_edge_info(graph_matrix[2][1]);
// std::cout<<std::endl;
// g.print_all_edge_info(graph_matrix[2][2]);

// AmiGraph::barcode_t B1, B2;
 // g.construct_graph_barcode(graph_matrix[2][1], B1);
 // g.construct_graph_barcode(graph_matrix[2][2], B2);

 // std::cout<<g.is_isomorphic(graph_matrix[2][1],graph_matrix[2][2])<<std::endl;


// std::cout<<"Result"<<std::endl;

// for(int i=0; i<graph_matrix[3].size(); i++){
// g.print_all_edge_info(graph_matrix[3][i]);
// std::cout<<std::endl;

// }


//g.print_all_edge_info();




/* 

bool iso;
// bool iso=g.is_isomorphic(g.current_graph, g.current_graph);
// std::cout<<"Finished"<<std::endl;
// std::cout<<"Isomorphism check returned "<< iso <<std::endl;

g.proposed_graph=g.current_graph;
g.AB_rand(g.current_graph);
AmiGraph::graph_t original=g.current_graph;

g.AT_rand(g.proposed_graph);
g.AT_rand(g.current_graph);
iso=g.is_isomorphic(g.current_graph, g.proposed_graph);

std::cout<<"Isomorphism check returned "<< iso <<std::endl;

double result_sum=0;
double count=0;
int top=1;//00000;

clock_t total;
total=clock();

for(int i=0; i< top; i++){

g.proposed_graph=original;
g.current_graph=original;

g.AT_rand(g.proposed_graph);
g.AT_rand(g.current_graph);
iso=g.is_isomorphic(g.current_graph, g.proposed_graph);


result_sum+=double(iso);
count=count+1.0;
 
}

 total=clock()-total;

 
std::cout<<"Results are "<<std::endl;
std::cout<<result_sum<<" "<< count <<" "<< result_sum/count<<std::endl;


float tim=(float)total/(float)CLOCKS_PER_SEC/(float)top;
std::cout<<"It took "<<total/(float)CLOCKS_PER_SEC<<" and "<< tim <<" to do the whole loop per attempt clocks="<< (float)CLOCKS_PER_SEC<<std::endl;
std::cout<<"This could be done "<< 1.0/tim*60.0 <<" times per minute :)"<<std::endl;
 */


/////////
/*
/////////

int MAX_ORDER=5;
AmiCalc::solution_set_matrix_t AMI_MATRIX;
AMI_MATRIX.resize(MAX_ORDER);
std::string top_dir;
top_dir="../../../S_P_R_F_up_to_5th_order/S_P_R_F_up_to_order_5/";
double ereg=1e-8;
g.ami.load_solutions(top_dir, AMI_MATRIX, MAX_ORDER, ereg);

std::cout<<AMI_MATRIX.size()<<std::endl;
std::cout<<AMI_MATRIX[2].size()<<std::endl;

//// Read External Variables////
std::string infile="testfile.txt";
std::cout<<"Reading external parameters from "<< infile<<std::endl;
AmiCalc::external_variable_list extern_list;
g.ami.read_external(infile, extern_list);

// Pick an order and num to evaluate - 0== first order, 1== second order etc
int ord=4;
int num=6;
ord--;
num--;


std::complex<double> result_sum(0,0);
double result_resum2, result_imsum2; 
result_resum2=0;
result_imsum2=0;
double count=0;
int top=1000000;//00000;

clock_t total;
total=clock();

for(int i=0; i< top; i++){

// Do something to the internal state - here we randomize the internal state of 'g'
g.randomize_state(g.current_state, AMI_MATRIX[ord][num].ami_parms_.N_INT_, AMI_MATRIX[ord][num].prefactor_);

// given that state - package together everything so that ami can evaluate it 
AmiCalc::ami_vars_list vars_list;
g.ami.construct_ami_vars_list(AMI_MATRIX[ord][num].R0_, g.current_state, extern_list, vars_list);
//
//AmiCalc::ami_vars test_ami_eval_vars=g.construct_ami_vars( AMI_MATRIX[ord][num].R0_, g.current_state, extern_list[0]);
// evaluate that solution for ord, num, for ALL extern_list, external variables
// define results
std::vector<std::complex<double>> results;
// actual evaluation line - answer is in results 
g.ami.evaluate_solutions(results, AMI_MATRIX[ord][num], vars_list);

// for(int i=0; i< vars_list.size(); i++){
// std::cout<<"Result from direct import is "<< results[i].real()<<" + "<<  results[i].imag()<<"i"<<std::endl;
// }





result_sum+=results[0];
result_resum2+=pow(results[0].real(),2.);
result_imsum2+=pow(results[0].imag(),2.);
count=count+1.0;
 

}

 total=clock()-total;

 
double re_err=sqrt(result_resum2/count - result_sum.real()/count*result_sum.real()/count)/sqrt(count);
double im_err=sqrt(result_imsum2/count - result_sum.imag()/count*result_sum.imag()/count)/sqrt(count);
std::cout<<"Results are "<<std::endl;
std::cout<<result_sum<<" "<< count <<" "<< result_sum/count<<" with errors "<< re_err<<" "<<im_err<<std::endl;


float tim=(float)total/(float)CLOCKS_PER_SEC/(float)top;
std::cout<<"It took "<<total/(float)CLOCKS_PER_SEC<<" and "<< tim <<" to do the whole loop per attempt clocks="<< (float)CLOCKS_PER_SEC<<std::endl;
std::cout<<"This could be done "<< 1.0/tim*60.0 <<" times per minute :)"<<std::endl;


/////////
*/
/////////


// AmiCalc::S_t test_S;
// AmiCalc::P_t test_P;
// AmiCalc::R_t test_R;
// AmiCalc::g_prod_t test_R0;
// double test_prefactor;
// double test_ereg;
// int test_graph_order;
// AmiCalc::external_variable_list test_ext;

// std::string file_loc;
// std::string S_file, P_epsfile, P_alphafile, R_epsfile, R_alphafile, R0_file, f_file;
// file_loc="../../../S_P_R_F_up_to_5th_order/S_P_R_F_up_to_order_5/3_order/"; 
// S_file=file_loc+"S_m_3"+"_txt_files/"+"S_m_3_num_1.txt";
// P_epsfile=file_loc+"P_mnta_m_3"+"_txt_files/"+"P_mnta_m_3_num_1.txt";
// P_alphafile=file_loc+"P_freq_m_3"+"_txt_files/"+"P_freq_m_3_num_1.txt";
// R_epsfile=file_loc+"R_mnta_m_3"+"_txt_files/"+"R_mnta_m_3_num_1.txt";
// R_alphafile=file_loc+"R_freq_m_3"+"_txt_files/"+"R_freq_m_3_num_1.txt";
// R0_file=file_loc+"alpha_m_3"+"_txt_files/"+"alpha_m_3_num_1.txt";
// f_file=file_loc+"f_m_3"+"_txt_files/"+"f_m_3_num_1.txt";

// //std::cout<<S_file<<std::endl;

// test_ereg=1e-5;
// test_graph_order=3;

// g.ami.read_text_S_solutions(S_file, test_S);
// std::cout<<"Test_S has size "<< test_S.size()<<std::endl;
// //g.ami.write_S_readable(test_S);


// g.ami.read_text_P_solutions(P_epsfile, P_alphafile, test_P);
// g.ami.write_P_readable(test_P);
// //g.ami.print_P(4,test_P);


// g.ami.read_text_R_solutions(R_epsfile, R_alphafile, test_R, test_graph_order);
// g.ami.write_R_readable(test_R);


// g.ami.read_text_R0(R0_file, test_R0);

// test_prefactor=g.ami.load_prefactor(f_file, test_graph_order);

// // at this point everything is constructed?

// AmiCalc::ami_parms test_amiparms(test_graph_order, test_ereg);
// //
// test_ext=extern_list; // me being lazy
// AmiCalc::solution_set solution(test_R0, test_S, test_P, test_R, test_amiparms,  test_prefactor);

// do calculation with solution struct 
// g.randomize_state(g.current_state, solution.ami_parms_.N_INT_, solution.prefactor_);
	
// AmiCalc::ami_vars test_ami_eval_vars=g.construct_ami_vars(solution.R0_, g.current_state, extern_list[0]);

// std::complex<double> calc_result=g.ami.evaluate(solution.ami_parms_, solution.R_, solution.P_, solution.S_,  test_ami_eval_vars);

// std::cout<<"Result from direct import is "<< calc_result.real()<<" + "<< calc_result.imag()<<"i"<<std::endl;



/*
// don't need to construct, since test_R, P and S are already created
//g.ami.construct(testamiparms, R0, R_array, P_array, S_array);

std::complex<double> result_sum(0,0);
double result_resum2, result_imsum2; 
result_resum2=0;
result_imsum2=0;
double count=0;
int top=10000000;//00000;

clock_t total;
total=clock();


// evaluate 

for(int i=0; i< top; i++){
 //This creates a random set of internal k values and puts them into current state
g.randomize_state(g.current_state, test_graph_order, test_prefactor);
	
AmiCalc::ami_vars test_ami_eval_vars=g.construct_ami_vars(test_R0, g.current_state, extern_list[0]);

std::complex<double> calc_result=g.ami.evaluate(test_amiparms, test_R, test_P, test_S,  test_ami_eval_vars);

//std::cout<<"Result from direct import is "<< calc_result.real()<<" + "<< calc_result.imag()<<"i"<<std::endl;


result_sum+=calc_result;
result_resum2+=pow(calc_result.real(),2.);
result_imsum2+=pow(calc_result.imag(),2.);
count=count+1.0;
 

}

 total=clock()-total;

 
double re_err=sqrt(result_resum2/count - result_sum.real()/count*result_sum.real()/count)/sqrt(count);
double im_err=sqrt(result_imsum2/count - result_sum.imag()/count*result_sum.imag()/count)/sqrt(count);
std::cout<<"Results are "<<std::endl;
std::cout<<result_sum<<" "<< count <<" "<< result_sum/count<<" with errors "<< re_err<<" "<<im_err<<std::endl;


float tim=(float)total/(float)CLOCKS_PER_SEC/(float)top;
std::cout<<"It took "<<total/(float)CLOCKS_PER_SEC<<" and "<< tim <<" to do the whole loop per attempt clocks="<< (float)CLOCKS_PER_SEC<<std::endl;
std::cout<<"This could be done "<< 1.0/tim*60.0 <<" times per minute :)"<<std::endl;


/*/



// This is a block comment for graph related construction and evaluation 
/*


// check that graph is connected
//boost::add_vertex(g.current_graph);
// bool connected=g.is_connected_graph(g.current_graph);
// std::cout<<"Graph is connected bool = "<< connected <<std::endl;

// zero is false - 1 is true 

//g.AB_rand(g.current_graph);
//g.AIL_rand(g.current_graph);


 bool result;
 // g.label_systematic_redone(g.current_graph);
 // g.label_systematic_redone(g.current_graph);
 // g.label_systematic_redone(g.current_graph);
 // g.label_systematic_redone(g.current_graph);
 g.repeated_labelling(g.current_graph, result);
 bool reducible=g.is_one_particle_reducible(g.current_graph);
 std::cout<<"Graph is reducible bool = "<< reducible <<std::endl;
 
  bool skel=g.is_skeleton(g.current_graph);
 std::cout<<"Graph is skeleton graph bool = "<< skel <<std::endl;
   bool hubb=g.is_hubbard(g.current_graph);
 std::cout<<"Graph is hubbard graph bool = "<< hubb <<std::endl;



// g.AB_rand(g.current_graph);
g.number_vertices(g.current_graph);
g.print_all_edge_info();

std::complex<double> result_sum(0,0);
double result_resum2, result_imsum2;
result_resum2=0;
result_imsum2=0;
double count=0;
int top=4;//00000;

clock_t total;
total=clock();


//g.random_int(0,6);	
	
	//total=clock()-total;

// for(int i=0; i< top; i++){
// modify
//if(i%2==0){g.AB_rand(g.current_graph);}else{g.DB_rand(g.current_graph);}


// label 

// g.label_systematic(g.current_graph);
// g.check_momentum_conservation(g.current_graph, result);
// if(result==false){
// g.repeated_labelling(g.current_graph, result);
// }

// Construct
AmiCalc::g_prod_t R0=g.graph_to_R0(g.current_graph);
// g.ami.print_g_prod_info(R0);

AmiCalc::R_t R_array;
AmiCalc::P_t P_array;
AmiCalc::S_t S_array;

// some logic around this perhaps
int graph_order=g.graph_order(g.current_graph);
double prefactor=g.get_prefactor(g.current_graph, graph_order);

double ereg=1e-5;
//std::cout<<"Looks like a "<< graph_order<<" order graph "<<std::endl;
AmiCalc::ami_parms amiparms(graph_order, ereg);

g.ami.construct(amiparms, R0, R_array, P_array, S_array);

///

// g.ami.print_final(graph_order, R_array,P_array,S_array);
// std::cout<<"Total number of terms was "<<R_array[graph_order].size()<< std::endl;

// evaluate 

for(int i=0; i< top; i++){
 g.randomize_state(g.current_state, graph_order, prefactor);
// g.print_state(g.current_state);

// g.construct_state(g.current_state);	
// g.print_state(g.current_state);
	
AmiCalc::ami_vars ami_eval_vars=g.construct_ami_vars(R0, g.current_state, extern_list[0]);

// g.print(ami_eval_vars);

std::complex<double> calc_result=g.ami.evaluate(amiparms, R_array, P_array, S_array,  ami_eval_vars);
//std::cout<< calc_result.real()<<" "<<calc_result.imag()<<std::endl;
result_sum+=calc_result;
result_resum2+=pow(calc_result.real(),2.);
result_imsum2+=pow(calc_result.imag(),2.);
count=count+1.0;
 

 // bool reducible=g.is_one_particle_reducible(g.current_graph);
// std::cout<<"Graph is reducible bool = "<< reducible <<std::endl;
 
  // bool skel=g.is_skeleton(g.current_graph);
// std::cout<<"Graph is skeleton graph bool = "<< skel <<std::endl;
   // bool hubb=g.is_hubbard(g.current_graph);
 // std::cout<<"Graph is hubbard graph bool = "<< hubb <<std::endl;

}

 total=clock()-total;

 
 double re_err=sqrt(result_resum2/count - result_sum.real()/count*result_sum.real()/count)/sqrt(count);
double im_err=sqrt(result_imsum2/count - result_sum.imag()/count*result_sum.imag()/count)/sqrt(count);
std::cout<<"Results are "<<std::endl;
std::cout<<result_sum<<" "<< count <<" "<< result_sum/count<<" with errors "<< re_err<<" "<<im_err<<std::endl;


float tim=(float)total/(float)CLOCKS_PER_SEC/(float)top;
std::cout<<"It took "<<total/(float)CLOCKS_PER_SEC<<" and "<< tim <<" to do the whole loop per attempt clocks="<< (float)CLOCKS_PER_SEC<<std::endl;
std::cout<<"This could be done "<< 1.0/tim*60.0 <<" times per minute :)"<<std::endl;
//std::cout<< calc_result.real()<<" "<<calc_result.imag()<<std::endl;

//*/ // end block comment 

//g.ami.construct(amiparms, R0, R_array, P_array, S_array);


//g.print_edgeson_vert(12);
//g.print_edgeson_vert(13);
//g.print_edgeson_vert(14);

// Timing testing 
// clock_t total;
// int max=100000;

// total=clock();
 // for(int i=0; i< max; i++){
	
	 // g.label_systematic(g.current_graph);
 // }
// total=clock()-total;
// std::cout<<"It took "<< (float)total/(float)CLOCKS_PER_SEC/(float)max<<" to do the labelling part per attempt"<<std::endl;
// no O3 857k per minute. with O3 27 million times per minute  

//std::cout<<"Final"<<std::endl;
//g.print_all_edge_info();

//label_systematic




// for( int i=0; i<2000000; i++){
// std::cout<<i<<std::endl;
// g.RIB_rand(g.current_graph);
// g.number_vertices(g.current_graph);
// g.print_all_vertex_info();
// g.RBI_rand(g.current_graph);

// }


// std::cout<<"Printing info"<<std::endl;
// g.number_vertices(g.current_graph);
// g.print_all_vertex_info();
// g.print_all_edge_info();



// std::cout<<"Printint info Before AT"<<std::endl;
 // g.number_vertices(g.current_graph);
// g.print_all_vertex_info();
// g.print_all_edge_info();

// g.AT_rand(g.current_graph);
// std::cout<<"Printint info After AT"<<std::endl;
 // g.number_vertices(g.current_graph);
// g.print_all_vertex_info();
// g.print_all_edge_info();

// g.DT_rand(g.current_graph);
// std::cout<<"Printint info after DT"<<std::endl;
 // g.number_vertices(g.current_graph);
// g.print_all_vertex_info();
// g.print_all_edge_info();
//g.DIL_rand(g.current_graph);





//for (int i=0; i<10; i++){ g.insert_int_rand(g.current_graph);}


// g.number_vertices(g.current_graph);
// g.print_all_vertex_info();
// g.print_all_edge_info();


//out.open("figs/final.dot"); 
//write_graphviz(out, g.current_graph); 



    return 0;
}


// this was checking in_out edge stuff
/*
boost::graph_traits<AmiGraph::graph_t>::vertex_iterator vi, v_end;
boost::graph_traits<AmiGraph::graph_t>::in_edge_iterator ei, ei_end;
boost::graph_traits<AmiGraph::graph_t>::out_edge_iterator eo, eo_end;
//set_all_vertex_freq(g);

//
std::cout<<"Current graph Contains the following vertices and info"<<std::endl;
std::cout<<num_vertices(g.current_graph)<<std::endl;
for (boost::tie(vi,v_end) = vertices(g.current_graph); vi != v_end; ++vi){
std::cout << "Vertex " << *vi << " of vertex type "<< g.current_graph[*vi].type_ << std::endl;

std::cout<<"In edges are "<<std::endl;
for (boost::tie(ei,ei_end) = in_edges(*vi, g.current_graph); ei != ei_end; ++ei){
				
				std::cout<<*ei<<std::endl;
				
			}
std::cout<<"oun edges are "<<std::endl;			
			for (boost::tie(eo,eo_end) = out_edges(*vi,g.current_graph); eo != eo_end; ++eo){
std::cout<<*eo<<std::endl;
				
			 }


   }



*/



// bubble debug 

/*

std::cout<<"Final vector sizes are "<< bubble_vertex_list[0].size()<<" "<<bubble_edge_list[0].size()<<" "<<legs_edge_list[0].size()<<std::endl;

std::cout<< "Found bubbles Nb = "<< bubble_vertex_list.size()<<std::endl;
for (int i=0; i<bubble_edge_list[0].size(); i++){
	std::cout<<"Here"<<std::endl;
	std::cout<<bubble_edge_list[0][i]<<std::endl;
	
}

for (int i=0; i<bubble_vertex_list[0].size(); i++){
	std::cout<<"Here"<<std::endl;
	std::cout<<bubble_vertex_list[0][i]<<std::endl;
	
}

for (int i=0; i<legs_edge_list[0].size(); i++){
	std::cout<<"Here"<<std::endl;
	std::cout<<legs_edge_list[0][i]<<std::endl;
	
}


for( int i=0; i<500000; i++){
	// std::cout<<i<<std::endl;
// std::cout<<"***Before add***"<<std::endl;
// g.number_vertices(g.current_graph);
// g.print_all_edge_info();
// g.print_all_vertex_info();
// std::cout<<std::endl;	
// bubble_vertex_list.clear();
// bubble_edge_list.clear();
// legs_edge_list.clear();
// g.bubble_finder(g.current_graph, bubble_vertex_list,  bubble_edge_list, legs_edge_list);
// std::cout<< "Found bubbles before AB Nb = "<< bubble_vertex_list.size()<<std::endl;
	
	
	// std::cout<<"Failure on i= "<<i<<std::endl;
// std::cout<<"add bubble"<<std::endl;
g.AB_rand(g.current_graph);

// bubble_vertex_list.clear();
// bubble_edge_list.clear();
// legs_edge_list.clear();
// g.bubble_finder(g.current_graph, bubble_vertex_list,  bubble_edge_list, legs_edge_list);
// std::cout<< "Found bubbles After AB Nb = "<< bubble_vertex_list.size()<<std::endl;

// std::cout<<"***Before DB***"<<std::endl;
// g.number_vertices(g.current_graph);
// g.print_all_edge_info();
// g.print_all_vertex_info();
// std::cout<<std::endl;	
// std::cout<<" Starting DB ";
g.DB_rand(g.current_graph);

// std::cout<< " ending DB"<< std::endl;
// std::cout<<"***After DB***"<<std::endl;
// g.number_vertices(g.current_graph);
// g.print_all_edge_info();
// g.print_all_vertex_info();
// bubble_vertex_list.clear();
// bubble_edge_list.clear();
// legs_edge_list.clear();
// g.bubble_finder(g.current_graph, bubble_vertex_list,  bubble_edge_list, legs_edge_list);
// std::cout<< "Found bubbles after DB Nb = "<< bubble_vertex_list.size()<<std::endl;
}


*/

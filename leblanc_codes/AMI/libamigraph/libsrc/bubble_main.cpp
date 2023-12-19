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
  int pairs=0;
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
	if(argc==7){
		max=atoi(argv[1]);
		write=atoi(argv[2]);
		read=atoi(argv[3]);
		filter=atoi(argv[4]);
		reduce=atoi(argv[5]);
		pairs=atoi(argv[6]);
	}

bool ppvertex=true;//false;//true;//false;  
  
//int seed=0; // but this can be set to anything
auto now =std::chrono::high_resolution_clock::now();
auto seed = std::chrono::duration_cast<std::chrono::microseconds>(now.time_since_epoch()).count();
// int seed=0;

std::cout<< seed <<" " <<std::endl;
// AmiGraph g(AmiBase::FORCE,seed);//g(AmiBase::Pi_phuu, seed);//g(AmiCalc::Pi,seed);//g(AmiCalc::Hartree,seed);
//AmiGraph g(AmiBase::Greens,seed);
// AmiGraph g(AmiBase::Pi_phuu,seed);
 AmiGraph g(AmiBase::Pi_ppud,seed);
// AmiGraph g(AmiBase::Pi_phuu,seed);
// AmiGraph g(AmiBase::Sigma,seed); // TODO: This never works I don't think - it misses some graphs. 

g.ami_parameters.int_type_=AmiBase::hubbard;//AmiBase::hubbard;//AmiBase::coulomb;

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
AmiGraph::gg_matrix_t ggm_uu, ggm_ud;
AmiGraph::gg_matrix_t ggm_1Pirr, ggm_1Pred;
AmiGraph::gg_matrix_t ggm_even, ggm_odd;

if(max!=0){
std::vector< std::vector< AmiGraph::graph_t>> graph_matrix;
if(read==0){
	std::cout<<"Generating Graphs"<<std::endl;
// g.generate_bubble_graphs(graph_matrix,max,mpi_rank);
g.generate_graphs(graph_matrix,max,mpi_rank);
g.print_ggm(ggm);


for(int i=0; i< graph_matrix.size(); i++){
	for(int m=0; m< graph_matrix[i].size(); m++){
	

std::cout<<"--- Order "<<i<<"--- Num "<<m<<"----"<<std::endl;	
g.number_vertices(graph_matrix[i][m]);		
	std::cout<<"-----------"<<std::endl;
g.print_all_edge_info(graph_matrix[i][m]);
std::cout<<"-----------"<<std::endl;
	}
}

}
if(read==1){
g.read_ggmp("last_ggm",preggm, max);
g.print_ggm(preggm);
// g.ggm_remove_pairs(preggm);
// g.print_ggm(preggm);
g.ggm_to_gm(preggm,graph_matrix, 0);
}


///



if(filter){
std::cout<<"Beginning filtering"<<std::endl;
// std::cout<<"Labelling"<<std::endl;



// g.label_graphs(graph_matrix,0,max);
// g.min_ext_label_counts(graph_matrix);
// g.reduce_gm_rf(graph_matrix,mpi_rank);
// g.reduce_gm_oneleg(graph_matrix, mpi_rank);

// Don't want this filter for pi bubbles. it cuts out RPA chains
// 1PBose filter removes tadpoles also 

// g.reduce_gm_1PBose(graph_matrix, mpi_rank,1);
// need to label for skeleton to work 
// don't label the g0 graph 
// g.label_graphs(graph_matrix,0,max);
if(g.ami_parameters.TYPE_!=AmiBase::FORCE){
g.label_graphs(graph_matrix,1,max);
}
// g.reduce_gm_skeleton(graph_matrix,mpi_rank);
g.reduce_gm_tp(graph_matrix, mpi_rank,0);
if(g.ami_parameters.TYPE_!=AmiBase::density ){

// g.reduce_gm_tp(graph_matrix, mpi_rank,0);
// g.reduce_gm_oneleg(graph_matrix,mpi_rank);
// g.reduce_gm_hubbard(graph_matrix, mpi_rank,0);






for(int i=0; i< graph_matrix.size(); i++){
	for(int m=0; m< graph_matrix[i].size(); m++){
	

std::cout<<"--- Order "<<i<<"--- Num "<<m<<"----"<<std::endl;	
g.number_vertices(graph_matrix[i][m]);		
	std::cout<<"-----------"<<std::endl;
g.print_all_edge_info(graph_matrix[i][m]);
std::cout<<"-----------"<<std::endl;
	}
}

// std::cout<<"Isomorphism debugging for force "<<std::endl;

// g.is_isomorphic(graph_matrix[4][4], graph_matrix[4][5]);

int this_min=1;
if(g.ami_parameters.TYPE_==AmiBase::FORCE){this_min=2;

std::cout<<"Before probes raph sizes are "<<std::endl;

for(int i=0; i<graph_matrix.size(); i++){
	
	std::cout<<i<<" "<< graph_matrix[i].size()<<std::endl;
}


std::vector< std::vector< AmiGraph::graph_t>> gm_force;
g.gm_attach_probe_lines(graph_matrix,gm_force);
graph_matrix=gm_force;

std::cout<<"After Graph sizes are now "<<std::endl;

for(int i=0; i<graph_matrix.size(); i++){
	
	std::cout<<i<<" "<< graph_matrix[i].size()<<std::endl;
}

}


for(int i=0; i< graph_matrix.size(); i++){
	for(int m=0; m< graph_matrix[i].size(); m++){
	

std::cout<<"--- Order "<<i<<"--- Num "<<m<<"----"<<std::endl;	
g.number_vertices(graph_matrix[i][m]);		
	std::cout<<"-----------"<<std::endl;
g.print_all_edge_info(graph_matrix[i][m]);
std::cout<<"-----------"<<std::endl;
	}
}

// exit(0);

// reattach things and then

// if(g.ami_parameters.TYPE_==AmiBase::FORCE){
// g.label_graphs(graph_matrix,1,max);
// }


//g.reduce_gm_1PBose(graph_matrix, mpi_rank,this_min);
// g.reduce_gm_1PFermi(graph_matrix, mpi_rank,this_min);




// if(g.ami_parameters.TYPE_!=AmiCalc::doubleocc && g.ami_parameters.TYPE_!=AmiCalc::Pi_phuu && g.ami_parameters.TYPE_!=AmiCalc::Pi_phud && g.ami_parameters.TYPE_!=AmiCalc::Pi_ppuu && g.ami_parameters.TYPE_!=AmiCalc::Pi_ppud ){
// g.reduce_gm_1PBose(graph_matrix, mpi_rank,1);
// }

// if(g.ami_parameters.TYPE_==AmiCalc::Pi_phud){
// g.reduce_gm_RPA_chains(graph_matrix, mpi_rank,1);
	
// }

}

// if(g.ami_parameters.TYPE_== AmiCalc::doubleocc){
// g.reduce_gm_tp(graph_matrix, mpi_rank,0);	
// }

if(ppvertex){
g.reduce_gm_ppvertex(graph_matrix, mpi_rank,0);


}
// g.reduce_gm_fock(graph_matrix, mpi_rank,1);
//
// g.reduce_gm_ladder(graph_matrix, mpi_rank,1);

/* if(g.ami_parameters.TYPE_==AmiCalc::density || g.ami_parameters.TYPE_==AmiCalc::doubleocc){

g.gm_close_loops(graph_matrix);	
	
} */

ggm.resize(graph_matrix.size());
}


// 99% sure we only need to number the vertices up to this point. maybe fix this?
// g.number_vertices(Each Graph);
// g.print_all_edge_info(graph_matrix[0][0]);
// g.label_graphs(graph_matrix,0,max);

int start=0;
if(g.ami_parameters.TYPE_==AmiBase::Greens){ // || g.ami_parameters.TYPE_==AmiCalc::density){
	start=1;
}
if(ppvertex){start=1;}



std::cout<<"Writing prelabelled to file ggm"<<std::endl;
g.gm_to_ggm(graph_matrix, ggm);
g.write_ggm("prelabel",ggm);

std::cout<<"Labelling"<<std::endl;
g.label_graphs_sys(graph_matrix,start,max);

std::cout<<"Converting to ggm"<<std::endl;
g.gm_to_ggm(graph_matrix, ggm);
g.write_ggm("postlabel",ggm);


// for(int ord=0; ord< graph_matrix.size(); ord++){
	
// for(int i=0; i< graph_matrix[ord].size(); i++){
// std::vector< AmiGraph::graph_t> ct_vec;	
// g.generate_bubble_ct(graph_matrix[ord][i], ct_vec);

// std::cout<<"On ord "<< ord<<" gnum "<<i<<std::endl;
// g.print_all_edge_info(graph_matrix[ord][i]);
// std::cout<<"Found counter-terms n="<<ct_vec.size()<<std::endl;	
	
// }
	
// }



// return 0;

// g.label_graphs(graph_matrix,0,max);
// g.print_all_edge_info(graph_matrix[0][0]);


// for(int i=0; i< graph_matrix.size(); i++){
	// for(int m=0; m< graph_matrix[i].size(); m++){
	
	// std::cout<<"Ord and num "<< i<<" "<<m<<std::endl<<std::endl;
	// std::cout<<"-------------"<<std::endl;
	// g.print_all_edge_info(graph_matrix[i][m]);
	// std::cout<<"-------------"<<std::endl<<std::endl;
	
	// }
// }

// std::cout<<"Closing loops in bubble main"<<std::endl;
// g.gm_close_loops(graph_matrix);

// for(int i=0; i< graph_matrix.size(); i++){
	// for(int m=0; m< graph_matrix[i].size(); m++){
	
	// std::cout<<"Ord and num "<< i<<" "<<m<<std::endl<<std::endl;
	// std::cout<<"-------------"<<std::endl;
	// g.print_all_edge_info(graph_matrix[i][m]);
	// std::cout<<"-------------"<<std::endl<<std::endl;
	
	// }
// }


// g.reduce_gm_nested(graph_matrix,mpi_rank,0);
// exit(0);




bool success;
if (write==1){
g.write_ggm("ggm_orig",ggm);}

if(read && !filter){
ggm=preggm;
}	



// std::cout<<"Compare sizes "<< graph_matrix[0].size()<<" "<< ggm[0].size()<<std::endl;

bool scramble=true;

if(reduce){
std::cout<<"Label sets"<<std::endl;	
g.construct_label_sets(ggm, mpi_rank,start, 3600, 1800);

// std::cout<<"Zero label size is "<< ggm[0][0].labels.size()<<std::endl;
// std::cout<<"test"<<std::endl;
// g.trim_label_sets(ggm,mpi_rank,200);
std::cout<<"starting reduce"<<std::endl;
if(g.graph_type!=AmiBase::FORCE){
g.check_label_sets(ggm, mpi_rank, start);
}

// g.reduce_ggm_git_odd_random(ggm,mpi_rank,2);
g.reduce_ggm_git(ggm,mpi_rank, 30);
g.write_ggm("ggm_reduced",ggm);



g.print_ggm(ggm);

std::cout<<"Final isomorphism check "<<std::endl;
g.print_ggm(ggm);
for(int thisord=0; thisord<ggm.size(); thisord++){
	std::cout<<"Checking isomorphisms on "<<thisord<<std::endl;
for(int g1=0; g1< ggm[thisord].size(); g1++){
	for(int g2=0; g2<ggm[thisord].size();g2++){
    for(int n=0; n< ggm[thisord][g1].graph_vec.size(); n++){
				for(int n2=0; n2<ggm[thisord][g2].graph_vec.size(); n2++){

std::cout<<"Checking isomorphism "<<g1<<" "<<n<<" "<<g2<<" "<<n2<<" = "<<g.is_isomorphic(ggm[thisord][g1].graph_vec[n], ggm[thisord][g2].graph_vec[n2])<<std::endl;
		
		}
		}
	}
}
}


g.reduce_ggm_isomorphic(ggm);



g.write_ggm("ggm_iso",ggm);


// g.print_all_edge_info(ggm[3][0].graph_vec[0]);

// bool test=g.is_ppvertex_graph(ggm[3][0].graph_vec[0]);

// std::cout<<"Returned "<<test<<std::endl;

/* 
bool test=g.spin_iso(ggm[3][0].graph_vec[0], ggm[3][0].graph_vec[2]);

std::cout<<"Spin iso test returned "<< test<<std::endl;

g.symmetrize_internal_bosonic_lines(ggm[3][0].graph_vec[0]);
g.symmetrize_internal_bosonic_lines(ggm[3][0].graph_vec[2]);

bool test3=g.spin_iso(ggm[3][0].graph_vec[0], ggm[3][0].graph_vec[2]);

std::cout<<"Bose symmetrized Spin iso test returned "<< test3<<std::endl;


 std::vector< std::vector<int>> e1, e2;

g.get_labelled_bose_lists(e1,ggm[3][0].graph_vec[0]);
g.get_labelled_bose_lists(e2, ggm[3][0].graph_vec[2]);


bool first=g.bose_lists_equal(e1,e2);

std::cout<<"First test returns "<<first<<std::endl;



AmiGraph::barcode_t B1, B2;
AmiGraph::tree_t T1, T2;
// std::cout<<"Constructing barcodes"<<std::endl;
g.construct_graph_barcode(ggm[3][0].graph_vec[0], B1, T1);
g.construct_graph_barcode(ggm[3][0].graph_vec[2], B2, T2);

// std::cout<<"Checking equal"<<std::endl;
// bool second=true;
bool second=g.barcodes_are_equal(B1, B2);

std::cout<<"Second test results "<< second<<std::endl;



g.print_all_edge_info(ggm[3][0].graph_vec[0]);
std::cout<<"--"<<std::endl;
g.print_all_edge_info(ggm[3][0].graph_vec[2]);
std::cout<<"--"<<std::endl;

AmiGraph::pp_barcode_t ppb1, ppb2;

std::cout<<"Generating bc1"<<std::endl;
g.pp_barcode(ggm[3][0].graph_vec[0], ppb1);

std::cout<<"Generating bc2"<<std::endl;
g.pp_barcode(ggm[3][0].graph_vec[2], ppb2);

bool pptest=g.compare_pp_barcodes(ppb1,ppb2);

std::cout<<"pp barcode test resulted in "<< pptest<<std::endl;



bool test2=g.is_isomorphic(ggm[2][2].graph_vec[0], ggm[2][2].graph_vec[1]);

std::cout<<"Is_isomorphic function returned "<<test2<<std::endl;
 */

}

if(g.graph_type==AmiBase::Pi_phuu){
g.ggm_eo_split(ggm, ggm_even, ggm_odd);

g.write_ggm("ggm_even", ggm_even);
g.write_ggm("ggm_odd", ggm_odd);

}

g.ggm_1P_split(ggm,ggm_1Pirr,ggm_1Pred);

g.write_ggm("ggm_irr", ggm_1Pirr);
g.write_ggm("ggm_red", ggm_1Pred);



if(g.graph_type==AmiBase::Pi_phuu || g.graph_type==AmiBase::Pi_ppuu){

g.ggm_split(ggm, ggm_uu, ggm_ud);	

g.write_ggm("ggm_uu_plus", ggm_uu);
g.write_ggm("ggm_ud_plus", ggm_ud);

std::cout<<"Full ggm was "<<std::endl;
g.print_ggm(ggm);

std::cout<<std::endl<<"Split was uu "<<std::endl;
g.print_ggm(ggm_uu);

std::cout<<std::endl<<"Split was ud "<<std::endl;
g.print_ggm(ggm_ud);

if(pairs){
	
	std::cout<<"Trying to pair uu"<<std::endl;
	
g.construct_label_sets(ggm_uu, mpi_rank, start, 3600, 1200);
g.construct_label_sets(ggm_ud, mpi_rank, start, 3600, 1200);
g.ggm_populate_git_pairs(ggm_uu,mpi_rank,start);

std::cout<<"Trying to pair ud"<<std::endl;
g.ggm_populate_git_pairs(ggm_ud,mpi_rank,start);	
	

g.write_ggmp("ggm_uu_plus_paired", ggm_uu,0);
g.write_ggmp("ggm_ud_plus_paired", ggm_ud,0);	
	
	
}


	
}

if(g.graph_type==AmiBase::doubleocc){

g.ggm_split(ggm, ggm_uu, ggm_ud);	

// g.write_ggm("ggm_uu_plus", ggm_uu);
g.write_ggm("ggm_doubleocc_ud", ggm_ud);

std::cout<<"Full ggm was "<<std::endl;
g.print_ggm(ggm);

std::cout<<std::endl<<"Split was uu "<<std::endl;
g.print_ggm(ggm_uu);

std::cout<<std::endl<<"Split was ud "<<std::endl;
g.print_ggm(ggm_ud);
	
	
// set ggm for doubleocc to the ggm_ud partial_sort
ggm=ggm_ud;	
	
}




// g.reduce_ggm_gos(ggm,mpi_rank, 80);
// g.reduce_ggm_gos(ggm,mpi_rank, 200);

// g.print_ggm(ggm);
// g.write_ggm("ggm_reduced",ggm);
// g.read_ggm("ggm",ggm, max);
// g.print_ggm(ggm);


// g.print_all_edge_info(ggm[2][0].graph_vec[0]);
// std::cout<<"Test systematic"<<std::endl;
// std::cout<<num_vertices(ggm[2][0].graph_vec[0])<<std::endl;
// g.systematic_vertex_label(ggm[2][0].graph_vec[0]);

//TODO: What is the status of the graph_vec graphs? Are they modified by reduce_ggm_git? or no?



// g.ggm_label(ggm,start); // don't think this labelling is necessary


// g.construct_label_sets(ggm, mpi_rank, 0, 100, 100); // this is needed if not reduced

if(pairs && !reduce){g.construct_label_sets(ggm, mpi_rank, start, 3600, 1200);}
if(pairs){

// std::cout<<"Debugging"<<std::endl;

// bool done;
// AmiGraph::graph_t g2_out;
// AmiGraph::git_perm_set_t pst;
// double prefactor=1;
// g.git_compare_v2(ggm[2][0].graph_vec[0], ggm[2][0].graph_vec[1],done, prefactor,pst, g2_out);

// std::cout<<"Done, prefactor "<< done<<" "<<prefactor<<std::endl;

g.ggm_populate_git_pairs(ggm,mpi_rank,start);
}

// after all possible labelling is done we can finally close the frickin graphs 
// TODO: It is easier to create the loop from the ggm? Will the prefactor be wrong? 
//

// 1) take open greens or Pi bubble graphs - label them and THEN convert to gg_AMI_MATRIX 
// 2) then modify the LABELLED graph to close the loops.  need to copy all properties including labels of the edges.
// std::cout<<"Pre closing"<<std::endl;
// for(int i=0; i< ggm.size(); i++){
	// for(int m=0; m< ggm[i].size(); m++){
	// for(int k=0; k< ggm[i][m].graph_vec.size(); k++){
	// std::cout<<"Ord and num "<< i<<" "<<m<<std::endl<<std::endl;
	// std::cout<<"-------------"<<std::endl;
	// g.print_all_edge_info(ggm[i][m].graph_vec[k]);
	// std::cout<<"-------------"<<std::endl<<std::endl;
	// }
	// }
// }




if(g.graph_type==AmiBase::density || g.graph_type==AmiBase::doubleocc){
g.close_ggm(ggm);
}

g.print_ggm(ggm);
// std::cout<<"Post Closing"<<std::endl;
// for(int i=0; i< ggm.size(); i++){
	// for(int m=0; m< ggm[i].size(); m++){
	// for(int k=0; k< ggm[i][m].graph_vec.size(); k++){
	// std::cout<<"Ord and num "<< i<<" "<<m<<std::endl<<std::endl;
	// std::cout<<"-------------"<<std::endl;
	// g.print_all_edge_info(ggm[i][m].graph_vec[k]);
	// std::cout<<"-------------"<<std::endl<<std::endl;
	// }
	// }
// }
/* 
g.number_vertices(ggm[2][2].graph_vec[0]);
g.number_vertices(ggm[2][2].graph_vec[1]);

 AmiGraph::barcode_t B1, B2;
AmiGraph::tree_t T1, T2;
// std::cout<<"Constructing barcodes"<<std::endl;
// g.construct_graph_barcode(ggm[2][2].graph_vec[0], B1, T1);
// g.construct_graph_barcode(ggm[2][2].graph_vec[1], B2, T2);

// std::cout<<"Checking equal"<<std::endl;
// bool second=true;
// bool second=g.barcodes_are_equal(B1, B2);
// std::cout<<"Barcodes equal? -> "<< second<<std::endl;
// g.print_barcodes(B1,B2);

g.symmetrize_internal_bosonic_lines(ggm[2][2].graph_vec[0]);
g.symmetrize_internal_bosonic_lines(ggm[2][2].graph_vec[1]);

g.print_all_edge_info(ggm[2][2].graph_vec[0]);
g.print_all_edge_info(ggm[2][2].graph_vec[1]);

g.reconstruct_graph_barcode(ggm[2][2].graph_vec[0], B1, T1);
g.reconstruct_graph_barcode(ggm[2][2].graph_vec[1], B2, T2);
g.print_barcodes(B1,B2); */

// g.compare_trees(T1,ggm[2][2].graph_vec[0] , T2,ggm[2][2].graph_vec[1]);

bool printmapping=false;//true;

if(printmapping){
	
	g.construct_label_sets(ggm, mpi_rank, start, 600, 200);

std::cout<<"Git mapping"<<std::endl;




std::ofstream file;
file.open("outfile.dat");

for(int ord =2; ord< ggm.size(); ord++){
	for(int gr=0; gr<ggm[ord].size(); gr++){
	std::vector<int> map;
	g.generate_git_map(ggm[ord][gr], map);
	
	file<< ord <<" "<<gr<<" ";
	for(int i=0; i< map.size(); i++){
		file<<map[i]<<" ";
	}
	file<<std::endl;
	
	
	}
}

file.close();

}



g.ggm_construct_ami_sol(ggm, 1e-8, mpi_rank);

g.ggm_check_pair_momenta(ggm, mpi_rank);


// for(int i=0; i< ggm.size(); i++){
	// for(int m=0; m<ggm[i].size();m++){
	
	// std::cout<<"ord "<< i<<" and num "<< m <<" have "<< ggm[i][m].graph_vec.size()<<" graphs and "<< ggm[i][m].gp_vec.size()<<" pairs"<<std::endl;	
		
	// }	
// }
// g.print_all_edge_info(ggm[2][0].graph_vec[0]);
// g.print_all_edge_info(ggm[2][2].graph_vec[0]);
// g.print_all_edge_info(ggm[2][1].graph_vec[0]);
// g.print_all_edge_info(ggm[2][1].graph_vec[1]);

g.write_ggmp("ggm_paired", ggm,0);
// g.read_ggmp("ggm_paired",ggm,max);

//////////////////


int check=0;

  for(int i=0; i< ggm[check].size(); i++){
for(int j=0; j< ggm[check][i].graph_vec.size(); j++){	
	std::cout<<std::endl;
	std::cout<<"Printing group "<< i<< " graph "<<j<<std::endl;
	std::cout<<std::endl;
g.print_all_edge_info(ggm[check][i].graph_vec[j]);	

std::cout<<"Loop count "<< g.count_fermi_loops(ggm[check][i].graph_vec[j])<<" and order "<< g.graph_order(ggm[check][i].graph_vec[j])<<std::endl;
std::cout<<"Get prefactor is "<< g.get_prefactor(ggm[check][i].graph_vec[j],g.graph_order(ggm[check][i].graph_vec[j]))<<std::endl;

g.ami.print_final(check, ggm[check][i].ss_vec[j].R_, ggm[check][i].ss_vec[j].P_, ggm[check][i].ss_vec[j].S_);
	
}

for(int j=0; j< ggm[check][i].gp_vec.size(); j++){	
	std::cout<<std::endl;
	std::cout<<"Printing group "<< i<< " pair "<<j<<std::endl;
	std::cout<<std::endl;
g.print_all_edge_info(ggm[check][i].gp_vec[j].g1_);	

std::cout<<"Loop count "<< g.count_fermi_loops(ggm[check][i].gp_vec[j].g1_)<<" and order "<< g.graph_order(ggm[check][i].gp_vec[j].g1_)<<std::endl;
std::cout<<"Get prefactor is "<< g.get_prefactor(ggm[check][i].gp_vec[j].g1_,g.graph_order(ggm[check][i].gp_vec[j].g1_))<<std::endl;

	std::cout<<std::endl;
g.print_all_edge_info(ggm[check][i].gp_vec[j].g2_);	

std::cout<<"Loop count "<< g.count_fermi_loops(ggm[check][i].gp_vec[j].g2_)<<" and order "<< g.graph_order(ggm[check][i].gp_vec[j].g2_)<<std::endl;
std::cout<<"Get prefactor is "<< g.get_prefactor(ggm[check][i].gp_vec[j].g2_,g.graph_order(ggm[check][i].gp_vec[j].g2_))<<std::endl;

// g.ami.print_final(check, ggm[check][i].ss_vec[j].R_, ggm[check][i].ss_vec[j].P_, ggm[check][i].ss_vec[j].S_);
	
}



}  







////////////////////
std::cout<<"Converting to ggamim"<<std::endl;

NewAmiCalc::gg_solution_set_matrix_t GG_AMI_MATRIX;
g.ggm_to_ggamim(ggm, GG_AMI_MATRIX,0); 

std::cout<<"Done"<<std::endl;



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


// std::cout<<"Loop count "<< g.count_fermi_loops(ggm[0][0].graph_vec[0])<<" and order "<< g.graph_order(ggm[0][0].graph_vec[0])<<std::endl;
// std::cout<<"Get prefactor is "<< g.get_prefactor(ggm[0][0].graph_vec[0],g.graph_order(ggm[0][0].graph_vec[0]))<<std::endl;

//// Read External Variables////
std::string infile="ext_vars.dat";
std::cout<<"Reading external parameters from "<< infile<<std::endl;
NewAmiCalc::external_variable_list extern_list;
g.ami.read_external(infile, extern_list);

int orde=4;
int dim=2;

    NewAmiCalc::internal_state state(orde,dim);
	g.zero_state(state,orde);
	g.randomize_state(state,orde);
	

	 for(int i=0; i< state.internal_k_list_.size(); i++){
		
		state.internal_k_list_[i][0]=M_PI/20.1+.51*i;
		state.internal_k_list_[i][1]==M_PI/10.1+.51*i;
			
	}
	
		
	
		
	

   NewAmiCalc::ami_vars_list vars_list;
	g.ami.construct_ami_vars_list(GG_AMI_MATRIX[4][0][0].R0_,GG_AMI_MATRIX[4][0][0].prefactor_, state, extern_list, vars_list);

	std::vector<double> Re_results, Im_results;
	// actual evaluation line - answer is in results 
	//graph.ami.evaluate_solutions(results, AMI_MATRIX[ord][num], vars_list);
	g.ami.evaluate_solutions(Re_results,Im_results, GG_AMI_MATRIX[4][0][0], vars_list);

	
	


}





    return 0;
}


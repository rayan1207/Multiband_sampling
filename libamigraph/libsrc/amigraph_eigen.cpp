#include "amigraph.hpp"
#include <Eigen/Dense>
#include <Eigen/QR>    

// Eigen::MatrixXd A = ... // fill in A
// Eigen::MatrixXd pinv = A.completeOrthogonalDecomposition().pseudoInverse();

// simple solution to AX=Y

bool AmiGraph::equal(Eigen::MatrixXd &A1, Eigen::MatrixXd &A2){

if(A1.rows() != A2.rows() || A1.cols()!= A2.cols()){ return false;}

for(int i=0; i< A1.rows(); i++){
for(int j=0; j< A1.cols(); j++){

if(A1(i,j)!= A2(i,j)){return false;}	
	
}

}	
	
return true;	
}

bool AmiGraph::neg_equal(Eigen::MatrixXd &A1, Eigen::MatrixXd &A2){

if(A1.rows() != A2.rows() || A1.cols()!= A2.cols()){ return false;}

for(int i=0; i< A1.rows(); i++){
for(int j=0; j< A1.cols(); j++){

if(A1(i,j)!= -A2(i,j)){return false;}	
	
}

}	
	
return true;	
}

bool AmiGraph::abs_equal(Eigen::MatrixXd &A1, Eigen::MatrixXd &A2){

if(A1.rows() != A2.rows() || A1.cols()!= A2.cols()){ return false;}

for(int i=0; i< A1.rows(); i++){
for(int j=0; j< A1.cols(); j++){

if(std::abs(A1(i,j))!= std::abs(A2(i,j))){return false;}	
	
}

}	
	
return true;	
}

void AmiGraph::construct_transform(Eigen::MatrixXd &A, Eigen::MatrixXd &X, Eigen::MatrixXd &Y){
	
	
	A=X.transpose().fullPivLu().solve(Y.transpose());
	A.transposeInPlace();
	

}

void AmiGraph::prune(Eigen::MatrixXd &inout, float threshold)
{
    inout = (threshold < inout.array().abs()).select(inout, 0.0f);
}

void AmiGraph::full_alpha_to_matrix(graph_t &g, Eigen::MatrixXd &A_out){
		
edge_vector_t internal;
find_internal_fermionic_edges(g, internal);
// std::cout<<internal.size()<<std::endl;


// remove_simple_alphas(g,internal);
// std::cout<<internal.size()<<std::endl;
// A_out.resize(internal.size()+1,g[internal[0]].g_struct_.alpha_.size()-1);
A_out.resize(g[internal[0]].g_struct_.alpha_.size(),internal.size());

// for(int i=0; i< internal.size()+1; i++){
	for(int i=0; i< internal.size(); i++){

for(int m=0; m< g[internal[i]].g_struct_.alpha_.size(); m++){

A_out(m,i)=g[internal[i]].g_struct_.alpha_[m];
 // if(i==internal.size()){A_out(i,m)=0;}
		
}

}
A_out.transposeInPlace();
}


//TODO: This function completely bungles the actual graph structure.  Can only be used right before construction/evaluation.... I think... testing in progress.
void AmiGraph::replace_graph_with_matrix(Eigen::MatrixXd &A_in, graph_t &g ){
		
// presumes that g already exists and is the appropriate dimentions

edge_vector_t internal;
find_internal_fermionic_edges(g, internal);


for(int i=0; i< internal.size(); i++){
	for(int j=0; j< A_in.cols(); j++){
	
	g[internal[i]].g_struct_.alpha_[j]=A_in(i,j);
	
	}
	// this resets the epsilons 
g[internal[i]].g_struct_.eps_.clear();
g[internal[i]].g_struct_.eps_.resize(internal.size(),0);
g[internal[i]].g_struct_.eps_[i]=1;
	
	}

// then fix epsilons based on alpha?
// reset epsilons

// this fixes the reset epsilons for multipoles.
fix_epsilons(g);

}



void AmiGraph::alpha_to_matrix(graph_t &g, Eigen::MatrixXd &A_out){
		
edge_vector_t internal;
find_internal_fermionic_edges(g, internal);
// std::cout<<internal.size()<<std::endl;


remove_simple_alphas(g,internal);
// std::cout<<internal.size()<<std::endl;
// A_out.resize(internal.size()+1,g[internal[0]].g_struct_.alpha_.size()-1);
A_out.resize(internal.size(),g[internal[0]].g_struct_.alpha_.size()-1);

// for(int i=0; i< internal.size()+1; i++){
	for(int i=0; i< internal.size(); i++){

for(int m=0; m< g[internal[i]].g_struct_.alpha_.size()-1; m++){

A_out(i,m)=g[internal[i]].g_struct_.alpha_[m];
 // if(i==internal.size()){A_out(i,m)=0;}
		
}

}
A_out.transposeInPlace();
}

void AmiGraph::remove_simple_alphas(graph_t &g, edge_vector_t &internal){

// print_all_edge_info(g);
// std::cout<<std::endl;
// for (int i=0; i< internal.size(); i++){
// print_edge_info(internal[i],g);	
	
// }

bool add=false;
edge_vector_t temp;
for(int i=0; i< internal.size(); i++){
	int count=0;
	int count2=0;
for(int m=0; m< g[internal[i]].g_struct_.alpha_.size(); m++){
	if(g[internal[i]].g_struct_.alpha_[m]==1){count++;}
		count2+=std::abs(g[internal[i]].g_struct_.alpha_[m]);	
}
// std::cout<<"Count is "<<count<<" for i "<<i <<std::endl;

if(count!=1){temp.push_back(internal[i]);}


}	

/* edge_vector_t temp;
for(int i=0; i< internal.size(); i++){
	int count=0;
for(int m=0; m< g[internal[i]].g_struct_.alpha_.size(); m++){
	count+=std::abs(g[internal[i]].g_struct_.alpha_[m]);	
}
std::cout<<"Count is "<<count<<" for i "<<i <<std::endl;

if(count!=1){temp.push_back(internal[i]);}
// if(count==1){
// internal.erase(internal.begin()+i);
// i--;	
// }

}	 */
	
internal=temp;
	
}





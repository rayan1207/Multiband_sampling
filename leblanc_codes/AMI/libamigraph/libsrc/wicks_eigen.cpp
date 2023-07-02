#include "amigraph.hpp"




void AmiGraph::extract_ABC(symDET_element &SE, Eigen::MatrixXi &A, Eigen::MatrixXi &B, Eigen::MatrixXi &C){

int n_indep=(SE.plist.size()-1)/2;

Eigen::MatrixXi M=Eigen::MatrixXi::Zero(n_indep, 2*n_indep+1);

std::cout<<"Empty M is "<<std::endl;
std::cout<<M<<std::endl;

int shift=1;
int count=0;

for (int row=0; row<M.rows(); row++){
count=0;	
for(int item=0; item<SE.plist.size(); item++){	
	
	if(SE.plist[item].first/2==row){
		
		// M(row,item)+=1;
	
	if(SE.plist[item].second==2*n_indep){
	M(row,M.cols()-shift)+=1;
	shift++;
	// continue;
	}	
	else {
	M(row,count)+=1;
	count++;
	}

	}	
 if(SE.plist[item].second/2==row){
		
	if(SE.plist[item].first==2*n_indep){
	M(row,M.cols()-shift)+=-1;
	shift++;
// continue;	
	}

	else{M(row,count)+=-1;
	count++;
	}
		
		// M(row,item)+=-1;
	}

}

}


std::cout<<M<<std::endl;

A=M.block(0,0,M.rows(),n_indep-1);

std::cout<<A<<std::endl;

M(0,0)=5;

std::cout<<A<<std::endl;
std::cout<<M<<std::endl;



}	
#include "amigraph.hpp"

void AmiGraph::extract_M(symDET_element &SE, std::vector< std::vector< int > > &M_OUT){
	
M_OUT.clear();


int n_indep=(SE.plist.size()-1)/2;
std::cout<<"Number of independent variables is "<< n_indep<<std::endl;

// determine size of M matrix 
M_OUT.resize(n_indep);
for(int i=0; i<M_OUT.size(); i++){

M_OUT[i].resize(SE.plist.size(),0);
	
}

for( int row=0; row< M_OUT.size(); row++){

for( int item=0; item< SE.plist.size(); item++){

if(SE.plist[item].first/2==row){
M_OUT[row][item]=1;	
}
else if(SE.plist[item].second/2==row){
M_OUT[row][item]=-1;		
}

}
}
	
	
}
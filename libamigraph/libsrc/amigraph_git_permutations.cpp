#include "amigraph.hpp"




void AmiGraph::get_git_rc_permutations( Eigen::VectorXd &r1, Eigen::VectorXd &c1, std::vector< std::vector<int>> &row_perm_set, std::vector< std::vector<int>> &col_perm_set){
	
row_perm_set.clear();
col_perm_set.clear();	
	
std::vector<int> temp;
bool add=true;

int max=r1.maxCoeff();

for(int i=2; i<max+1;i++){
	
for (int j=0; j< r1.rows(); j++){

if( r1[j]==i){ temp.push_back(j);}

}	
if(temp.size()>1){
row_perm_set.push_back(temp);	
}	
temp.clear();
}



// columns 
max=c1.maxCoeff();
for(int i=2; i<max+1;i++){
	
for (int j=0; j< c1.rows()-1; j++){

if( c1[j]==i){ temp.push_back(j);}

}	
if(temp.size()>1){
col_perm_set.push_back(temp);	
}	
temp.clear();
}

	
	
}


bool AmiGraph::git_sign_permutations_v2(Eigen::MatrixXd &A1, Eigen::MatrixXd &A2, Eigen::VectorXd &energy, Eigen::VectorXd &kshifts, git_perm_set_t &pst_out){

// type 0 is k_sign_change that has options from k=0, k<A2.cols()-1
// type 1 is swap col signs, same range of options.

// first, get all options

bool success=false;

git_perm_set_t pst;
git_perm_set_t pst_scrambled;
git_perm_t gpt;	

gpt.resize(2,-1);
// three kinds of permutations


// k-shift permutations
for(int i=0; i< A2.cols()-1; i++){
gpt[0]=0;
gpt[1]=i;
pst.push_back(gpt);
}

//col sign swap permutations 
for(int i=0; i< A2.cols()-1; i++){
gpt[0]=1;
gpt[1]=i;
pst.push_back(gpt);
}	
// for now we don't flip any row signs

// for(int i=0; i< A2.rows(); i++){

// if(A2(i,A2.cols()-1)!=0){	
// gpt[0]=2;
// gpt[1]=i;
// pst.push_back(gpt);
// }
// }

// TODO: Right here we should scramble the pst vector. This way the permutations sampled always change when you run the code. 
pst_scrambled=pst;
std::shuffle(pst_scrambled.begin(), pst_scrambled.end(), rand_gen);
pst=pst_scrambled;

// int max=pst.size(); // Ideally would use the entire length, but the permutations take too long. So truncate at 7
int max=std::min(pst.size(),(size_t)7);
std::vector<int> v;

for(int i=0; i< pst.size(); i++){
v.push_back(i);	
}


// print_permutation_set(pst);
	
Eigen::MatrixXd A2c, D;

Eigen::VectorXd ec, kc;

int count=0;

for(int i=1; i<max; i++){
	
	
	
if(success){break;}	
	// std::cout<<"Working on sign Permutation Length "<<i<<std::endl;
do{
	count++;
	
	A2c=A2;
	ec=energy;
	kc=kshifts;
	
D=A1-A2c;
if(D.isZero(0)){
	// std::cout<<"No change required"<<std::endl;
success=true;
break;
}	
	
	for(int j=0; j< i; j++){
		// std::cout<<"Listing j values for i with pst.size()="<<pst.size()<<std::endl;
		// std::cout<<j<<" "<< i<<" "<<v[j]<<std::endl;
		 // std::cout<<pst[v[j]][0]<<" "<<pst[v[j]][1]<<" "<<std::endl;
git_apply_sign_permutation(A2c,ec, kc, pst[v[j]]);

	}
	
	// std::cout<<"Energy and k shifts "<<std::endl<< ec<< std::endl<<std::endl<<kc.transpose()<<std::endl;
	
D=A1-A2c;

if(D.isZero(0)){
	
// duplicate the pst vector for pst_out
for(int j=0; j< i; j++){
pst_out.push_back(pst[v[j]]);
}	
	
	// std::cout<<"Diff is zero, ALL DONE!"<<std::endl; 
// std::cout<<"Total checked is count="<<count<<std::endl;
success=true; 
A2=A2c;
energy=ec;
kshifts=kc;
break;}

	// if(i>3){ // this is for safety. something wrong with this if pst size is close to i 
	 std::reverse(v.begin()+i, v.end());
	// }
}while(std::next_permutation(v.begin(),v.end()));
	
	

}	

if(!success){
	// std::cout<<"Exhausted all permutations. no sign match."<<std::endl;
	}

return success;


	
}	

void AmiGraph::git_get_perm_set(std::vector< std::vector<int>> &row_perm_set, std::vector< std::vector<int>> &col_perm_set, git_perm_set_t &pst){
	

git_perm_t temp_pair;

temp_pair.resize(3,-1);
//temp_rc_set.resize(2);

// add empty no-swap cases

for(int set=0; set< row_perm_set.size(); set++){

for(int i=0; i< row_perm_set[set].size(); i++){
for(int j=i+1; j< row_perm_set[set].size(); j++){

temp_pair[0]=0;
temp_pair[1]=row_perm_set[set][i];
temp_pair[2]=row_perm_set[set][j];


pst.push_back(temp_pair);
// std::cout<<"Pst is size "<<pst.size()<<" with last entry of size "<<pst[pst.size()-1].size()<<std::endl;
temp_pair.clear();
temp_pair.resize(3,-1);

	
}
}	
	
}


// add all single permutations of cols
for(int set=0; set< col_perm_set.size(); set++){

for(int i=0; i< col_perm_set[set].size(); i++){
for(int j=i+1; j< col_perm_set[set].size(); j++){

temp_pair[0]=1;
temp_pair[1]=col_perm_set[set][i];
temp_pair[2]=col_perm_set[set][j];

pst.push_back(temp_pair);
temp_pair.clear();
temp_pair.resize(3,-1);

	
}
}	
	
}

// from those, create all combined permutations 


	
	
}

bool AmiGraph::git_sign_permutations(Eigen::MatrixXd &A1, Eigen::MatrixXd &A2, Eigen::VectorXd &energy, Eigen::VectorXd &kshifts){

// type 0 is k_sign_change that has options from k=0, k<A2.cols()-1
// type 1 is swap col signs, same range of options.

// first, get all options

bool success=false;

git_perm_set_t pst;
git_perm_t gpt;	

gpt.resize(2,-1);
// three kinds of permutations

//
//col sign swap permutations 
for(int i=0; i< A2.cols()-1; i++){
gpt[0]=1;
gpt[1]=i;
pst.push_back(gpt);
}
// k-shift permutations
for(int i=0; i< A2.cols()-1; i++){
gpt[0]=0;
gpt[1]=i;
pst.push_back(gpt);
}	


for(int i=0; i< A2.rows(); i++){

if(A2(i,A2.cols()-1)!=0){	
gpt[0]=2;
gpt[1]=i;
pst.push_back(gpt);
}
}

// int max=pst.size(); // Ideally would use the entire length, but the permutations take too long. So truncate at 7
int max=std::min(pst.size(),(size_t)7);
std::vector<int> v;

for(int i=0; i< pst.size(); i++){
v.push_back(i);	
}


// print_permutation_set(pst);
	
Eigen::MatrixXd A2c, D;

Eigen::VectorXd ec, kc;

int count=0;

for(int i=1; i<max; i++){
	
	
	
if(success){break;}	
	// std::cout<<"Working on sign Permutation Length "<<i<<std::endl;
do{
	count++;
	
	A2c=A2;
	ec=energy;
	kc=kshifts;
	
D=A1-A2c;
if(D.isZero(0)){
	// std::cout<<"No change required"<<std::endl;
success=true;
break;
}	
	
	for(int j=0; j< i; j++){
		// std::cout<<"Listing j values for i with pst.size()="<<pst.size()<<std::endl;
		// std::cout<<j<<" "<< i<<" "<<v[j]<<std::endl;
		 // std::cout<<pst[v[j]][0]<<" "<<pst[v[j]][1]<<" "<<std::endl;
git_apply_sign_permutation(A2c,ec, kc, pst[v[j]]);

	}
	
	// std::cout<<"Energy and k shifts "<<std::endl<< ec<< std::endl<<std::endl<<kc.transpose()<<std::endl;
	
D=A1-A2c;

if(D.isZero(0)){
	// std::cout<<"Diff is zero, ALL DONE!"<<std::endl; 
// std::cout<<"Total checked is count="<<count<<std::endl;
success=true; 
A2=A2c;
energy=ec;
kshifts=kc;
break;}

	// if(i>3){ // this is for safety. something wrong with this if pst size is close to i 
	 std::reverse(v.begin()+i, v.end());
	// }
}while(std::next_permutation(v.begin(),v.end()));
	
	

}	

if(!success){
	// std::cout<<"Exhausted all permutations. no sign match."<<std::endl;
	}

return success;


	
}	




bool AmiGraph::git_compare_n_permutations(int n, Eigen::MatrixXd A1, Eigen::MatrixXd &A2, git_perm_set_t &pst){

// std::cout<<"Pst length is "<< pst.size()<<std::endl;

// print_permutation_set(pst);

bool success=false;
int count=0;

std::vector<int> v;

for(int i=0; i< pst.size(); i++){
v.push_back(i);	
}

	
Eigen::MatrixXd A2c;
for(int i=1; i<n; i++){
	
if(success){break;}		
	// std::cout<<"Swap combinations "<<i<<std::endl;
do{
	
	count++;
	
	A2c=A2;
	// std::cout<<"Start loop"<<std::endl;
	for(int j=0; j< i; j++){
		 // std::cout<<j<<" "<<v[j]<<std::endl;
		 // std::cout<<pst[v[j]][0]<<" "<<pst[v[j]][1]<<" "<<pst[v[j]][2]<<" "<<std::endl;
		git_apply_permutation(A2c,pst[v[j]]);

	}
	// std::cout<<"End loop"<<std::endl;
	// std::cout<<std::endl;
	
	
if(abs_equal(A1,A2c)){
	A2=A2c;
	// std::cout<<"Found equal permutation after count="<<count<<std::endl;
success=true;
break;}	
	
	
	 // if(i>3){
	 std::reverse(v.begin()+i, v.end());
	 // }
}while(std::next_permutation(v.begin(),v.end()));
	
	

}	
	


return success;
	
}

bool AmiGraph::git_compare_rc_permutations(Eigen::MatrixXd A1, Eigen::MatrixXd &A2, git_perm_set_t &pst){

// std::cout<<"Pst length is "<< pst.size()<<std::endl;

// print_permutation_set(pst);

bool success=false;
int count=0;

std::vector<int> v;

int n=pst.size()+1;

//////// This is really important. if the pst permutations are too big, just assume they will never match.
if(pst.size()>9){
return false;	
}
//////////

for(int i=0; i< pst.size(); i++){
v.push_back(i);	
}

	
Eigen::MatrixXd A2c;
for(int i=1; i<n; i++){
	
if(success){break;}		
	// std::cout<<"Swap combinations "<<i<<std::endl;
do{
	
	count++;
	
	A2c=A2;
	// std::cout<<"Start loop"<<std::endl;
	for(int j=0; j< i; j++){
		 // std::cout<<j<<" "<<v[j]<<std::endl;
		 // std::cout<<pst[v[j]][0]<<" "<<pst[v[j]][1]<<" "<<pst[v[j]][2]<<" "<<std::endl;
		git_apply_permutation(A2c,pst[v[j]]);

	}
	// std::cout<<"End loop"<<std::endl;
	// std::cout<<std::endl;
	
	
if(abs_equal(A1,A2c)){
	A2=A2c;
	// std::cout<<"Found equal permutation after count="<<count<<std::endl;
success=true;
break;}	
	
	
	 // if(i>3){
	 std::reverse(v.begin()+i, v.end());
	 // }
}while(std::next_permutation(v.begin(),v.end()));
	
	

}	
	


return success;
	
}

void AmiGraph::print_permutation_set( git_perm_set_t &pst){

std::cout<<"Printing pst set of length "<<pst.size()<<std::endl;
// std::cout<<pst[1].size()<<std::endl;	
for(int i=0; i< pst.size(); i++){
	// std::cout<<"pst has size "<< pst[i].size();
// std::cout<<" for i="<<i<<" ";
if(pst[i][0]==0){std::cout<<" Row (type 1): ";}
if(pst[i][0]==1){std::cout<<" Col (type 2): ";}

for(int j=1; j<pst[i].size(); j++){
std::cout<<pst[i][j]<<" ";
}
std::cout<<std::endl;

// std::cout<<"Next"<<std::endl;
}	
	
	
}

// bool AmiGraph::git_compare_single_permutation(Eigen::MatrixXd A1, Eigen::MatrixXd &A2, git_perm_set_t &pst){

// for(int i=0; i< pst.size(); i++){
// Eigen::MatrixXd A2c=A2;
// git_apply_permutation(A2c, pst[i]);

// if(abs_equal(A1,A2c)){
	// A2=A2c;
	// std::cout<<"Found equal permutation"<<std::endl;
// return true;}

// }	
	
// return false;	
// }

void AmiGraph::git_apply_permutation(Eigen::MatrixXd &A2, git_perm_t &gpt){
	
if(gpt[0]==0){
A2.row(gpt[1]).swap(A2.row(gpt[2]));
// std::cout<<"Applying row swap "<<gpt[1]<<" "<<gpt[2]<<std::endl;
}	

if(gpt[0]==1){
A2.col(gpt[1]).swap(A2.col(gpt[2]));
// std::cout<<"Applying col swap "<<gpt[1]<<" "<<gpt[2]<<std::endl;	
reorder_matrix(A2); // not sure if this would help or harm	// Necessary because we don't include the row swaps for the first 'n' rows, that can always be trivially done like this reordering.
}
	
	
}


void AmiGraph::git_apply_sign_permutation(Eigen::MatrixXd &A2, Eigen::VectorXd &energy, Eigen::VectorXd &kshifts, git_perm_t &gpt){
	
if(gpt[0]==0){
git_shift_kp(gpt[1], A2, energy, kshifts);
}	

if(gpt[0]==1){
git_swap_col_signs(gpt[1], A2);
}
if(gpt[0]==2){
git_swap_row_signs(gpt[1], A2, energy);
}

	
	
}




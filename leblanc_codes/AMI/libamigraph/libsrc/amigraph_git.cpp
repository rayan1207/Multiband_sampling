#include "amigraph.hpp"

void AmiGraph::git_apply_T1(Eigen::MatrixXd &A2, Eigen::RowVectorXd &perm, int first, int second){
// This function just swaps the two columns of A2 - first and second - not allowed to swap with the last column.

A2.col(first).swap(A2.col(second));
perm.col(first).swap(perm.col(second));

}

void AmiGraph::git_apply_T2(Eigen::MatrixXd &A2, Eigen::RowVectorXd &perm, int first){

A2.col(first)=-A2.col(first);
perm.col(first)=-perm.col(first);	
	
}


void AmiGraph::git_apply_T3(int index, Eigen::MatrixXd &A2, Eigen::VectorXd &energy, Eigen::VectorXd &kshifts){
	
kshifts(index)++;
// adjust all the alphas
for(int i=0; i< A2.rows(); i++){

if(A2(i,index)!=0){	
energy(i)=energy(i)*(-1); // flip the energy index once

for(int j=0; j< A2.cols(); j++){	

A2(i,j)=-A2(i,j); // flip each element of the i'th row
 

}
}	
}	

	
	
}

void AmiGraph::permute_state(NewAmiCalc::internal_state &state, NewAmiCalc::k_vector_t &k_shift, git_perm_set_t &pst){

// std::cout<<"Permuting state"<<std::endl;

// print_permutation_set(pst);

// 0 is k_shift
// 1 is col sign swap - flips sign of momentum
for (int i=0; i< pst.size(); i++){

if(pst[i][0]==1){
	// std::cout<<"Sign shift on col "<< pst[i][1]<<std::endl;
	
	for(int dim=0; dim< state.dim_; dim++){
		state.internal_k_list_[pst[i][1]][dim]=-state.internal_k_list_[pst[i][1]][dim];
	}
	
}

if(pst[i][0]==0){
	// std::cout<<"k_shift on col "<< pst[i][1]<<std::endl;
shift_state(state,pst[i][1],k_shift);	
}

}
	
	
}

void AmiGraph::repeated_git_compare( graph_t &g1, graph_t &g2, bool &result, double &prefactor,git_perm_set_t &pst, graph_t &g2_out){

std::vector<int> count1,count2;	

			
label_counts(g1,count1);
label_counts(g2,count2);

if(count1.back() == count2.back()){
			
			std::sort(count1.begin(), count1.end());
			std::sort(count2.begin(), count2.end());
			
			if (count1==count2){
bool junk=false;

for(int i=0; i<200; i++){
repeated_labelling(g1,junk);
repeated_labelling(g2,junk);
			
// git_compare(g1, g2,result, prefactor);
git_compare_v2(g1, g2,result, prefactor,pst, g2_out);
 if(result){break;}
}
}
// }



}

	
	
	
	
}

void AmiGraph::git_compare_v2( graph_t &g1, graph_t &g2, bool &result, double &prefactor,git_perm_set_t &pst, graph_t &g2_out){


g2_out=g2;
	
result=false;

// print_all_edge_info(g1);
// std::cout<<std::endl;

// print_all_edge_info(g2);
// std::cout<<std::endl;

Eigen::MatrixXd A1,A2;

full_alpha_to_matrix(g1,A1);
full_alpha_to_matrix(g2,A2);

// std::cout<<"Starting with A1"<<std::endl<<A1<<std::endl<<std::endl<<"And A2"<<std::endl<<A2<<std::endl;

// Immediately reject if ext columns don't have the same number of ext frequencies

bool first=false;
if(A1.col(A1.cols()-1).count() != A2.col(A2.cols()-1).count()){
first=false;
// std::cout<<"Ext frequency counts not the same"<<std::endl;
result=false;
return;	
}
/////

// Next we allow trivial re-ordering of rows (trivial permutations)
// in order to put the independent lines at the top for easy comparison
reorder_matrix(A1);
reorder_matrix(A2);

Eigen::VectorXd energy=Eigen::VectorXd::Ones(A1.rows());
Eigen::VectorXd kshifts=Eigen::VectorXd::Zero(A1.cols());

// First - rearrange to be abs_equal via only row and column permutations. then replace g2 with that result that will be then fed into ami.

Eigen::VectorXd r1,r2,c1,c2;

row_col_counts(A1,r1,c1);
row_col_counts(A2,r2,c2);

first= row_col_compare(r1,r2,c1,c2);

// std::cout<<"Row/col compare returned "<< first<<std::endl;

// if either the row or the column ID doesn't match then exit 
if(!first){result=false; 
// std::cout<<"Pole-IDs are not permutations different"<<std::endl;
return;}

// Now that we know they match. we can sort the rows (trivial) and columns (not trivial) 

Eigen::MatrixXd row_Perm=Eigen::MatrixXd::Identity(A1.rows(), A1.rows());
Eigen::MatrixXd col_Perm=Eigen::MatrixXd::Identity(A1.cols(), A1.cols());

git_sort_rows(A1, A2, r1, r2, row_Perm);
git_sort_cols(A1, A2, c1, c2, col_Perm);

// discard these permutations
row_Perm=Eigen::MatrixXd::Identity(A1.rows(), A1.rows());
col_Perm=Eigen::MatrixXd::Identity(A1.cols(), A1.cols());



reorder_matrix(A1);
reorder_matrix(A2);

// std::cout<<"Working with A1"<<std::endl<<A1<<std::endl<<std::endl<<"And A2"<<std::endl<<A2<<std::endl;

// now A2 has been permuted such that the row and column pol IDs are the same.  But they may not be abs_equal - r1 and c1 should look like r2 and c2 now - but don't re-evaluate just use r1 and c1

// std::cout<<"R1 c1 ";
// std::cout<< r1<<" | "<< c1 << std::endl;
// std::cout<<"R2 c2 ";
// std::cout<< r2 <<" | "<< c2 <<std::endl;
row_col_counts(A1,r1,c1);
row_col_counts(A2,r2,c2);
// std::cout<<"R1 c1 ";
// std::cout<< r1<<" | "<< c1 << std::endl;
// std::cout<<"R2 c2 ";
// std::cout<< r2 <<" | "<< c2 <<std::endl;
// if they are already abs_equal then skip this part. if not manipulate them
	if(!abs_equal(A1,A2)){
		std::vector< std::vector<int>> row_perm_set, col_perm_set;

		get_git_rc_permutations(r1,c1, row_perm_set, col_perm_set);

		git_perm_set_t pst;
		git_get_perm_set(row_perm_set,col_perm_set,pst);

		// print_permutation_set(pst);

		bool p_compare= git_compare_rc_permutations( A1,A2, pst);

		// exit if cannot find absolute equal 
		if(!p_compare){	
		// std::cout<<"Could not find ABS_EQUAL during permutations, exiting"<<std::endl;
		return;
		}
	}
		

// std::cout<<"Abs Equal A1"<<std::endl<<A1<<std::endl<<std::endl<<"And A2"<<std::endl<<A2<<std::endl;
// The result of changes to A2 must therefore be put back onto g2 from this stage IF they are found to be equal
Eigen::MatrixXd A2_for_g2=A2;


// then from new_g2 we will find the T2 and T3 transformations that map from one to the other. once those are found, store the permutation set - we will use that to manipulate a state to construct the relevant cancelling state

// std::cout<<"Found equal within prefactor "<<std::endl;

bool success=false;
Eigen::MatrixXd  A2c;
A2c=A2; 
// first fix the ext_signs by defining a prefactor// Q what impact does this have on k's?  It is equivalent to flipping all the k's but only in one line... maybe should not do this? Would require that epsilon(-k)=epsilon(k) lets try skipping this

git_perm_set_t pst_result;

success=git_sign_permutations_v2(A1,A2c,energy,kshifts, pst_result);

// std::cout<<"Final Equal A1"<<std::endl<<A1<<std::endl<<std::endl<<"And A2"<<std::endl<<A2c<<std::endl;

if(!success){return;}

result=success;

// std::cout<<"SUCCESS!!"<<std::endl; 
// std::cout<<"The resulting Energy and k shifts "<<std::endl<< energy<< std::endl<<std::endl<<kshifts.transpose()<<std::endl;


double x=0;

for(int i=0; i<A2.cols(); i++){
x+=kshifts(i)*A2.col(i).count();
}	
//TODO: Feb 4 2020 uncommented this. no idea if it's correct or not.
// x+= energy.prod();

prefactor=energy.prod();//std::pow(-1.0, x);
// prefactor=energy.prod();//std::pow(-1.0, x)*energy.prod();

// std::cout<<"Prefactor is "<<prefactor<<" with x "<<x<<std::endl; 

// std::cout<<"---------PERMUTATION SET-------------"<<std::endl;
// print_permutation_set(pst_result);
pst=pst_result;

if(success){
replace_graph_with_matrix(A2_for_g2,g2_out); 
}

}

// typedef std::vector<graph_group> gg_vec_t;
// typedef std::vector< gg_vec_t > gg_matrix_t;

// todo this just spits out a list of graphs to remove 
// todo there is an issue that if there are only two in a group and they are isomorphic then neither gets added 
// todo this is questionable if it is working or not 
void AmiGraph::reduce_ggm_isomorphic(gg_matrix_t &ggm){
	
gg_matrix_t new_ggm;
gg_vec_t this_ggv;
graph_group this_gg;

for(int ord=0; ord<ggm.size(); ord++){
	// for each order. pick a graph and cycle through all others. 
	
	for(int gg=0; gg<ggm[ord].size(); gg++){
		
		for(int n1=0; n1< ggm[ord][gg].graph_vec.size(); n1++){
			
			bool add=true;
			
			// if(n1!=0){
	
	for(int gg2=0; gg2<ggm[ord].size(); gg2++){
	
		
	
		for(int n2=0; n2<ggm[ord][gg2].graph_vec.size(); n2++){
	
			bool result=is_isomorphic(ggm[ord][gg].graph_vec[n1],ggm[ord][gg2].graph_vec[n2]);
			
			if(result){
				if(gg!=gg2 || n1!=n2){
					if(n1>n2){
					add=false;
					}// don't want to remove both, just the second time this is found... maybe?
				}
				
				
				
			}
	
		}
				}// end gg2
			// } // end initial add 
				
				
				if(add){
			
			this_gg.graph_vec.push_back(ggm[ord][gg].graph_vec[n1]);
			
		}
		
			
		
	}
	/// These are commented out - since empty order groups need to not get mixed up
	/// Maybe this first one needs to be there so that we skip empty groups...
	if(this_gg.graph_vec.size()!=0){
		this_ggv.push_back(this_gg);
		this_gg.graph_vec.clear();
	}
	
	}// end gg 
	
	/// Even if there are no graphs of a given order still need to add them. 
	// if(this_ggv.size()!=0){
		new_ggm.push_back(this_ggv);
		this_ggv.clear();
	// }
}


	
	ggm=new_ggm;
	
	
	
	
}

void AmiGraph::reduce_ggm_git( gg_matrix_t &ggm, int mpi_rank, int cutoff){
	
int initial;
int fin;	
std::cout<<"Starting git reduction"<<std::endl;
	
for (int i=0; i< ggm.size(); i++){	
int count=0;
std::cout<<"Starting Order="<<i<<std::endl;
do{

// if(count%10==0 && count>0){	
// std::cout<<"Order ="<<i<<" count="<<count<<std::endl;}

initial=ggm[i].size();
reduce_ggv_git(ggm[i]);	
fin=ggm[i].size();
if(fin==initial){count++;}else{count=0; 
if(mpi_rank==0){
std::cout<<"Order "<<i<<" reduced to "<< fin<<std::endl;}
}
}while(count<cutoff);
}


}

void AmiGraph::reduce_ggm_git_odd_random( gg_matrix_t &ggm, int mpi_rank, int cutoff){
	
int initial;
int fin;	
std::cout<<"Starting git reduction"<<std::endl;
	
for (int i=0; i< ggm.size(); i++){	
int count=0;
std::cout<<"Starting Order="<<i<<std::endl;
do{

// if(count%10==0 && count>0){	
// std::cout<<"Order ="<<i<<" count="<<count<<std::endl;}

initial=ggm[i].size();
reduce_ggv_git_odd_random(ggm[i]);	
fin=ggm[i].size();
if(fin==initial){count++;}else{count=0; 
if(mpi_rank==0){
std::cout<<"Order "<<i<<" reduced to "<< fin<<std::endl;}
}
}while(count<cutoff);
}


}

void AmiGraph::reduce_ggv_git_odd_random( gg_vec_t &ggv){

//gos_compare(graph_vector[j], graph_vector[k],result, prefactor);
bool result=false;
double prefactor;
for(int i=0; i< ggv.size(); i++){
	if(ggv[i].graph_vec.size()%2==0){continue;}
	// if(i%10==0){std::cout<<"Comparing graph i="<<i<<" to others"<<std::endl;}
for(int j=i+1; j< ggv.size(); j++){
	if(ggv[j].graph_vec.size()%2==0){continue;}
	bool result=false;

if(j%20==0){std::cout<<"Made it to i,j="<<i<<","<<j<<" with sizes "<< ggv[i].graph_vec.size()<<" "<< ggv[j].graph_vec.size()<<std::endl;}
	
	
// for each group - pick random graph from a group to compare 
int rand1, rand2;
// do{
// get a random graph	
rand1=random_int(0, ggv[i].graph_vec.size()-1);
rand2=random_int(0, ggv[j].graph_vec.size()-1);
// std::cout<<"Checking labels at random "<<rand1<<" vs "<< rand2<<std::endl;
// now get a random label from each graph : this is done inside 'single_git_compare'

// }while(rand1==rand2);
//single_gos_compare(ggv[i].graph_vec[rand1], ggv[j].graph_vec[rand2], result, prefactor);
// std::cout<<"Comparing groups "<<i<<" "<<j<<std::endl;
// std::cout<<"Rand graphs  "<< rand1<<" of "<< ggv[i].graph_vec.size()<<" and "<< rand2 <<" of "<< ggv[j].graph_vec.size()<<std::endl;

// create a random label for each 
graph_t l1=ggv[i].graph_vec[rand1];
graph_t l2=ggv[j].graph_vec[rand2];
bool junk=true;
repeated_labelling(l1,junk);
repeated_labelling(l2,junk);

git_perm_set_t pst;
graph_t g2_out;

repeated_git_compare(l1,l2,result, prefactor,pst, g2_out);
// single_git_compare(l1,l2, result, prefactor, pst, g2_out);

// TODO: unclear where the pst should live. needs to know who it is in reference to... 

if( result){
	
prefactor=prefactor;
// std::cout<<"Collapse with prefactor "<< prefactor<<std::endl;
collapse_ggv(i, j, ggv, prefactor);
j--;
}

}
}
	
}


void AmiGraph::reduce_ggv_git( gg_vec_t &ggv){

//gos_compare(graph_vector[j], graph_vector[k],result, prefactor);
bool result=false;
double prefactor;
for(int i=0; i< ggv.size(); i++){
	// if(i%10==0){std::cout<<"Comparing graph i="<<i<<" to others"<<std::endl;}
for(int j=i+1; j< ggv.size(); j++){
	bool result=false;

if(j%20==0){std::cout<<"Made it to i,j="<<i<<","<<j<<std::endl;}	
	
// for each group - pick random graph from a group to compare 
int rand1, rand2;
// do{
// get a random graph	
rand1=random_int(0, ggv[i].graph_vec.size()-1);
rand2=random_int(0, ggv[j].graph_vec.size()-1);
// std::cout<<"Checking labels at random "<<rand1<<" vs "<< rand2<<std::endl;
// now get a random label from each graph : this is done inside 'single_git_compare'

// }while(rand1==rand2);
//single_gos_compare(ggv[i].graph_vec[rand1], ggv[j].graph_vec[rand2], result, prefactor);
// std::cout<<"Comparing groups "<<i<<" "<<j<<std::endl;
// std::cout<<"Rand graphs  "<< rand1<<" of "<< ggv[i].graph_vec.size()<<" and "<< rand2 <<" of "<< ggv[j].graph_vec.size()<<std::endl;
git_perm_set_t pst;
graph_t g2_out;
single_git_compare(ggv[i].labels[rand1], ggv[j].labels[rand2], result, prefactor, pst, g2_out);

// TODO: unclear where the pst should live. needs to know who it is in reference to... 

if( result){
	
prefactor=prefactor;
// std::cout<<"Collapse with prefactor "<< prefactor<<std::endl;
collapse_ggv(i, j, ggv, prefactor);
j--;
}

}
}
	
}

//dec 4 - rewrite this to systematically compare all labels rather than a random one.
void AmiGraph::single_git_compare(labels_t &L1, labels_t &L2, bool &result, double &prefactor, git_perm_set_t &pst, graph_t &g2_out){
bool out;
result=false;
std::vector<int> count1, count2;

graph_t g1, g2;
int rand1=random_int(0, L1.size()-1);
int rand2=random_int(0, L2.size()-1);
g1=L1[rand1];

// TODO: Figure out what to do with this pst...
// git_perm_set_t pst;

// for(int j=0; j<L2.size();j++){

// g2=L2[j];
g2=L2[rand2];

// std::cout<<"Git compare "<< rand1<<" "<<rand2<<std::endl;

// if(rand1==3 && rand2==116){
// print_all_edge_info(g1);
// std::cout<<std::endl;	
// print_all_edge_info(g2);
// }

			
label_counts(g1,count1);
label_counts(g2,count2);

if(count1.back() == count2.back()){
			
			std::sort(count1.begin(), count1.end());
			std::sort(count2.begin(), count2.end());
			
			if (count1==count2){
			
// git_compare(g1, g2,result, prefactor);
git_compare_v2(g1, g2,result, prefactor,pst, g2_out);
// if(result){break;}
			
}
// }



}
// if( result){
// std::cout<<"Found equivalent graphs from labels "<<rand1<<" vs "<< rand2<<std::endl;
// }

// std::cout<<"Result was "<<result<<std::endl;	
	
}



void AmiGraph::git_compare( graph_t &g1, graph_t &g2, bool &result, double &prefactor){
	
result=false;

// print_all_edge_info(g1);
// std::cout<<std::endl;

// print_all_edge_info(g2);
// std::cout<<std::endl;

Eigen::MatrixXd A1,A2;
// alpha_to_matrix(g1,A1);
// alpha_to_matrix(g2,A2);
full_alpha_to_matrix(g1,A1);
full_alpha_to_matrix(g2,A2);

bool first=false;
if(A1.col(A1.cols()-1).count() != A2.col(A2.cols()-1).count()){
first=false;
std::cout<<"Ext frequency counts not the same"<<std::endl;
result=false;
return;	
}



// std::cout<<"Starting with A1"<<std::endl<<A1<<std::endl<<std::endl<<"And A2"<<std::endl<<A2<<std::endl;
// these pull the independent rows to the top in order
reorder_matrix(A1);
reorder_matrix(A2);

// std::cout<<"Before Reorder A1"<<std::endl<<A1<<std::endl<<std::endl<<"And A2"<<std::endl<<A2<<std::endl;
	
Eigen::VectorXd r1,r2,c1,c2;

Eigen::VectorXd energy=Eigen::VectorXd::Ones(A1.rows());
Eigen::VectorXd kshifts=Eigen::VectorXd::Zero(A1.cols());


// energy.resize(A1.rows());
// kshifts.resize(A1.cols());
// for(int i=0; i<

// for (int i=0; i< energy.rows(); i++){
// energy[i]=1;	
// }
	
row_col_counts(A1,r1,c1);
row_col_counts(A2,r2,c2);

first= row_col_compare(r1,r2,c1,c2);

if(!first){result=false; 
std::cout<<"Pole-ID different"<<std::endl;
return;}

// std::cout<<"Comparison has same row column character "<<std::endl;


// make permutation tracking identity matrices
Eigen::MatrixXd row_Perm=Eigen::MatrixXd::Identity(A1.rows(), A1.rows());
Eigen::MatrixXd col_Perm=Eigen::MatrixXd::Identity(A1.cols(), A1.cols());

git_sort_rows(A1, A2, r1, r2, row_Perm);
git_sort_cols(A1, A2, c1, c2, col_Perm);
reorder_matrix(A1);
reorder_matrix(A2);

std::cout<<"After Reorder A1"<<std::endl<<A1<<std::endl<<std::endl<<"And A2"<<std::endl<<A2<<std::endl;

// define a k-permuation ex. (0 1 2 3...) -> (1 3 2 0)
Eigen::RowVectorXd k_perm=Eigen::VectorXd::Zero(A1.cols());
Eigen::VectorXd row_perm=Eigen::VectorXd::Zero(A1.rows());
for(int i=0; i< k_perm.cols(); i++){
k_perm(i)=i;	
}
k_perm=k_perm*col_Perm;

// std::cout<<"K-perm is now "<<k_perm<<std::endl;

for(int i=0; i< row_perm.rows(); i++){
row_perm(i)=i;	
}
std::cout<<"row_perm is now "<<row_perm<<std::endl;
row_perm=row_Perm*row_perm;

std::cout<<"row_perm is now "<<row_perm<<std::endl;
std::cout<<"row_Perm is now "<<row_Perm<<std::endl;

// reorder_matrix(A1);
// reorder_matrix(A2);

// std::cout<<"After sorting we get A1"<<std::endl<<A1<<std::endl<<std::endl<<"And A2"<<std::endl<<A2<<std::endl;

bool success=false;
Eigen::MatrixXd D, A2c;

if(!abs_equal(A1,A2)){
// next need to check for equal permutations
// std::cout<<"Not absolutely equal so looking at permutations "<<std::endl;
std::vector< std::vector<int>> row_perm_set, col_perm_set;

get_git_rc_permutations(r1,c1, row_perm_set, col_perm_set);

git_perm_set_t pst;
git_get_perm_set(row_perm_set,col_perm_set,pst);

bool p_compare= git_compare_n_permutations(pst.size()+1, A1,A2, pst);

if(!p_compare){	
std::cout<<"Could not find ABS_EQUAL during permutations, exiting"<<std::endl;
return;
}
}

// std::cout<<"After permutations we get A1"<<std::endl<<A1<<std::endl<<std::endl<<"And A2"<<std::endl<<A2<<std::endl;
	



if(abs_equal(A1,A2)){ 
// std::cout<<"Found equal within prefactor "<<std::endl;


// result=true;
// prefactor=get_prefactor(A1,A2);
 // std::cout<<"Prefactor is "<< prefactor<<std::endl;
 // std::cout<<std::endl<<A1<<std::endl<<std::endl<<A2<<std::endl<< std::endl;
 
 // std::cout<<"Energy and k shifts "<<std::endl<< energy<< std::endl<<std::endl<<kshifts.transpose()<<std::endl;

// work with a copy of A2 that is A2c 

A2c=A2; 
D=A1-A2;


git_fix_ext_signs(A1, A2c,energy );
success=git_sign_permutations(A1,A2c,energy,kshifts);

/* 
for(int run=0; run<max; run++){

std::cout<<"Now have "<< prefactor<<std::endl;
 std::cout<<std::endl<<A1<<std::endl<<std::endl<<A2c<<std::endl<< std::endl;
D=A1-A2c;

if(D.isZero(0)){std::cout<<"Diff is zero, ALL DONE! That took "<< run<<" git sweeps."<<std::endl; success=true; break;}

std::cout<<"Diff map is "<<std::endl<<std::endl<<D<<std::endl; 
// ext signs
git_fix_ext_signs(A1, A2c,energy );
// next flip any trivial signs in rows by pulling out prefactors
// git_swap_easy_rows(A1, A2c,energy );
// swap a row or column as necessary
git_decide_swap(A1,A2c,energy,kshifts);

std::cout<<"The Energy and k shifts right now are "<<std::endl<< energy<< std::endl<<std::endl<<kshifts.transpose()<<std::endl;

} */


} 

// if successful, replace A2 with A2c 
if(success){
replace_graph_with_matrix(A2,g2); // replace graph with A2 BEFORE the transformations that give you A2C.
A2=A2c;// unclear that this has any role to play at this point.

result=success;

// get prefactor here from epsilon
std::cout<<"SUCCESS!!"<<std::endl; 
std::cout<<"The resulting Energy and k shifts "<<std::endl<< energy<< std::endl<<std::endl<<kshifts.transpose()<<std::endl;

std::cout<<"Replacing graph with new unsorted A2 matrix"<<std::endl;


}

// std::cout<<"Exiting compare with final"<<std::endl;
// std::cout<<std::endl<<A1<<std::endl<<std::endl<<A2<<std::endl<< std::endl;


// compute prefactor - pretty easy at this point.
double x=0;

for(int i=0; i<A2.cols(); i++){
x+=kshifts(i)*A2.col(i).count();
}	
//x+= energy.prod();

std::cout<<"Energy prod "<< energy.prod() <<std::endl;
std::cout<<energy<<std::endl;

prefactor=energy.prod();//std::pow(-1.0, x)*energy.prod();

std::cout<<"Prefactor is "<<prefactor<<" with x "<<x<<std::endl; 
 
return;
	
}




void AmiGraph::reorder_matrix(Eigen::MatrixXd &A1){

Eigen::MatrixXd copy=A1;

// std::cout<<"Reordering"<<std::endl;

int num=copy.cols()-1;

int count=0;
int row_count=0;

for(int count=0; count<num; count++){
for(int i=0; i< A1.rows(); i++){
// std::cout<<i<<" "<< A1.rows()<< std::endl;	
// std::cout<<A1.row(i).count()<<std::endl;

if(A1.row(i).count()==1 && A1.row(i)(count)==1){
copy.row(row_count)=A1.row(i);
row_count++;

// std::cout<<"Copied "<< A1.row(i)<<std::endl;
}	

}	
}
// for(int i=0	; i< num; i++){
	// Eigen::VectorXd one=Eigen::VectorXd::Zero(A1.cols());
	// one(i)=1;
	// copy.row(i)=one;

// }

for(int i=0; i< A1.rows(); i++){

if(A1.row(i).count() >1){
copy.row(row_count)=A1.row(i);
row_count++;	
// std::cout<<"Copied "<< A1.row(i)<<std::endl;

}
}

A1=copy;	
	
}


void AmiGraph::git_decide_swap(Eigen::MatrixXd &A1, Eigen::MatrixXd &A2, Eigen::VectorXd &energy, Eigen::VectorXd &kshifts){
Eigen::MatrixXd D;
D=A1-A2;

// std::cout<<"Deciding swap "<<std::endl;
// std::cout<<"Diff matrix is "<<std::endl<<D<<std::endl;

// Eigen::VectorXd r1c1;
// get row and column counts of D
// row_col_counts(D,r1,c1);

int rc=0;
int diff=-1;
int index=-1;

int diffc=-1;
int indexc=-1;
int indexr=-1;
int diffr=-1;


for(int i=0; i< D.cols(); i++){
if (D.col(i).count()>diffc){
diffc=D.col(i).count();
rc=1;
indexc=i;
} 	
// std::cout<<"Index and diff "<<index<<" "<<diff<<" "<<rc<<std::endl;	
	
}

for(int i=0; i< D.rows(); i++){
if (D.row(i).count()>diffr){
diffr=D.row(i).count();
rc=0;
indexr=i;
} 	
// std::cout<<"Index and diff and rc "<<index<<" "<<diff<<" "<<rc<<std::endl;	
}

// std::cout<<"Diff r and c are "<< diffr<<" "<< diffc<<std::endl;

if(diffr>diffc){
diff=diffr;
rc=0;
index=indexr;	
}
if(diffc>diffr){
diff=diffc;
rc=1;
index=indexc;	
}
if(diffc==diffr){

int s=random_int(0,1);
// std::cout<<"Rolled random "<<s<<std::endl;
switch(s){

case 0:{
	diff=diffr;
rc=0;
index=indexr;
}
break;
case 1:{
diff=diffc;
rc=1;
index=indexc;
}
}	
	
	
}



if(diff==0){
	// std::cout<<"Decided to do nothing "<<std::endl;
return;
}


if(rc==0){
	// std::cout<<"Decided to swap row "<<index<<std::endl;
// git_swap_easy_rows(Eigen::MatrixXd &A1, Eigen::MatrixXd &A2, Eigen::VectorXd &energy );
git_k_sign_change(A1,A2, energy, kshifts );	

}else{
	// std::cout<<"Decided to swap col "<<index<<std::endl;
git_swap_col_signs(	index, A2 );
}

	
	
	
	
}

void AmiGraph::git_swap_col_signs(int index, Eigen::MatrixXd &A2){

for(int i=0; i<A2.rows(); i++){

A2(i,index)=A2(i,index)*(-1);

}	
	
}

void AmiGraph::git_shift_kp(int index, Eigen::MatrixXd &A2, Eigen::VectorXd &energy, Eigen::VectorXd &kshifts){
	
kshifts(index)++;
// adjust all the alphas
for(int i=0; i< A2.rows(); i++){

if(A2(i,index)!=0){	
energy(i)=energy(i)*(-1); // flip the energy index once
// Still really unsure about the role of energy(i) signs here 

for(int j=0; j< A2.cols(); j++){	

A2(i,j)=-A2(i,j); // flip each element of the i'th row
 

}
}	
}	
	
	
}

void AmiGraph::git_k_sign_change(Eigen::MatrixXd &A1, Eigen::MatrixXd &A2, Eigen::VectorXd &energy, Eigen::VectorXd &kshifts ){
// std::cout<<"Entering sign change..."<<std::endl;

Eigen::MatrixXd D;
D=A1-A2;

// std::cout<<"Entering sign change..."<<std::endl;
// std::cout<<"Diff matrix is "<<std::endl<<D<<std::endl;


bool is_zero=D.isZero(0);
if(is_zero){ return;} // do nothing if D is all zero

// decide what is the row to swap - one with the most diff 
int index=-1;
int count=0;

for(int i=0; i< D.rows(); i++){

if( D.row(i).count() > count){
count=D.row(i).count();
index=i;
}
	
}		
if(index==-1){ 
std::cout<<"Something isn't working correctly"<<std::endl;
return;} // again do nothing if D is all zero - this should never trigger.

// decide what is the minimal swap for a row 
// std::cout<<"Deciding..."<<std::endl;
count=100; // this is sloppy
int col_index=-1;
// std::cout<<A2.cols()<<std::endl;
for(int j=0; j< A2.cols()-1; j++){
// std::cout<<"j is "<<std::endl;
	if(A2(index,j)!=0){
		// if(j==0){count=A2.col(j).count();}else{
			if(A2.col(j).count()<count){
				
				count=A2.col(j).count();
				col_index=j;
				
			}
		
		// }
	
	}

}


// if col_index=-1 exit

if(col_index==-1){
	// std::cout<<"Col index is "<<col_index<<std::endl;
	return;} // this should never trigger.

// so now we are changing k_col_index 
kshifts(col_index)++;
// adjust all the alphas
for(int i=0; i< A2.rows(); i++){

if(A2(i,col_index)!=0){	

for(int j=0; j< A2.cols(); j++){	

A2(i,j)=-A2(i,j);

}
}	
}

// std::cout<<"Exiting"<<std::endl;

	
	
	
	return;
}


void AmiGraph::git_swap_row_signs(int row_index, Eigen::MatrixXd &A2,Eigen::VectorXd &energy ){

for(int j=0; j< A2.cols(); j++){	
A2(row_index,j)=A2(row_index,j)*(-1);

}	

// std::cout<<"Before "<<energy(row_index)<<std::endl;
energy(row_index)=energy(row_index)*(-1);
// std::cout<<"After "<<energy(row_index)<<std::endl;	
	
return;	
}

void AmiGraph::git_swap_easy_rows(Eigen::MatrixXd &A1, Eigen::MatrixXd &A2, Eigen::VectorXd &energy ){

Eigen::MatrixXd D;
D=A1-A2;

// std::cout<<"Diff matrix is "<<std::endl<<D<<std::endl;

std::vector<int> row_list,col_list;
for(int i=0; i< D.rows(); i++){
for(int j=0; j< D.cols(); j++){

if(D(i,j)!=0){
row_list.push_back(i);
col_list.push_back(j);	
}

}	
}

for (int item=0; item< row_list.size(); item++){

int count=A2.row(row_list[item]).count();

// if only one entry in a row then flip signs 
if(count==1){
// std::cout<<"Easy flip for row "<< row_list[item]<<std::endl;
git_swap_row_signs(row_list[item], A2, energy);
	
}
	
	
}

	
	
return;	
}

void AmiGraph::git_fix_ext_signs(Eigen::MatrixXd &A1, Eigen::MatrixXd &A2,Eigen::VectorXd &energy ){

// std::cout<<"Fixing ext signs "<<std::endl;

Eigen::MatrixXd D;
D=A1-A2;

// std::cout<<"Diff matrix is "<<std::endl<<D<<std::endl;
	
std::vector<int> list;
// find any ext that need to swap sign and pull out prefactor
for (int i=0; i< A2.rows(); i++){

if( D.col(D.cols()-1)(i) !=0){

// std::cout<<i<<" "<<A2.col(A2.cols()-1)(i)<<std::endl;
list.push_back(i);	

git_swap_row_signs(i, A2, energy);
// std::cout<<"Swapped row signs "<<i<<std::endl;

}

}		
	
	
	
	return;
}


// void AmiGraph::create_diff_map(Eigen::MatrixXd &A1, Eigen::MatrixXd &A2, Eigen::MatrixXd &D){

// D=A1-A2;

// return;
// }	

void AmiGraph::git_sort_rows(Eigen::MatrixXd &A1, Eigen::MatrixXd &A2, Eigen::VectorXd &r1, Eigen::VectorXd &r2, Eigen::MatrixXd &perm){

for(int i=0; i< A1.rows(); i++){
// std::cout<<"sorting on row "<<i<<std::endl;
if(r1(i)!= r2(i)){

for(int m=i+1; m< A1.rows(); m++){

if(r1(i)==r2(m)){
r2.row(m).swap(r2.row(i));
perm.row(m).swap(perm.row(i));	
continue;
}

}// m	
}// if	
}// i

Eigen::MatrixXd A2c=perm*A2;

 
// check explicitly if rows of A2c are equal to A1 after permutation 


if(r1==r2){
	
A2=perm*A2;	
}else{std::cout<<"Rows not equal"<<std::endl;} // I think this probably doesn't work correctly. Don't think I've ever seen this printed...



}

void AmiGraph::git_sort_cols(Eigen::MatrixXd &A1, Eigen::MatrixXd &A2, Eigen::VectorXd &r1, Eigen::VectorXd &r2, Eigen::MatrixXd &perm){
// std::cout<<"Beginning column sort with c1 and c2"<<std::endl;
// std::cout<<r1<<std::endl<<std::endl<<r2<<std::endl;

// not allowed to swap last column?

for(int i=0; i< A1.cols()-1; i++){
// std::cout<<"sorting on col "<<i<<std::endl;
if(r1(i)!= r2(i)){

for(int m=i+1; m< A1.cols()-1; m++){
// std::cout<<"on m "<<m<<std::endl;

if(r1(i)==r2(m)){
// std::cout<<"1"<<std::endl;
// std::cout<<r1<<std::endl<<std::endl<<r2<<std::endl;
r2.row(m).swap(r2.row(i));
// std::cout<<"1"<<std::endl;
perm.col(m).swap(perm.col(i));
// std::cout<<"1"<<std::endl;	
continue;
}

}// m	
}// if	
}// i

if(r1==r2){
	
A2=A2*perm;	
}else{
	// std::cout<<"cols not equal"<<std::endl;
}



}	



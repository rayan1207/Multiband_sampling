#include "amigraph.hpp"



void AmiGraph::reduce_ggm_gos( gg_matrix_t &ggm, int mpi_rank, int cutoff){
	
int initial;
int fin;	
	
for (int i=0; i< ggm.size(); i++){	
int count=0;
do{
initial=ggm[i].size();
reduce_ggv_gos(ggm[i]);	
fin=ggm[i].size();
if(fin==initial){count++;}else{count=0; 
if(mpi_rank==0){
std::cout<<"Order "<<i<<" reduced to "<< fin<<std::endl;}
}
}while(count<cutoff);
}


}

void AmiGraph::ggm_AID_reduce(int attempts, gg_matrix_t &ggm, std::vector< std::vector<int>> &gg_lists, int mpi_rank){
	if(mpi_rank==0){
	std::cout<<"---------"<<std::endl<<"Starting AID ggm reduction using existing labels"<<std::endl<<"---------"<<std::endl;	}	
	
for(int ord=1; ord< ggm.size(); ord++){
	if(mpi_rank==0){
	std::cout<<"Working on Order "<< ord<<std::endl;}
	
	int start=get_start_index(gg_lists[ord]);
	
	// std::vector<int> list;

int current=0;
int index=0;
	
//gos_compare(graph_vector[j], graph_vector[k],result, prefactor);
bool result=false;
double prefactor;
for(int i=0; i< ggm[ord].size(); i++){
	
	// bool inot= !in_list(i,list);
if(i==current){current=current+gg_lists[ord][index]; index++;}	
// if(inot){

for(int j=i+1; j< ggm[ord].size(); j++){
	
if( i<current && j<start && j>current){j=start;}
// std::cout<<i<<" "<<j<<std::endl;
// if(!in_list(j,list)){	
	
	bool result=false;
// for each group - pick random graph to compare 
int rand1=random_int(0, ggm[ord][i].graph_vec.size()-1);
int rand2=random_int(0, ggm[ord][j].graph_vec.size()-1);

for(int m=0; m<attempts;m++){
	// std::cout<<m<<std::endl;
	// single_gos_compare_labels(ggv[i].labels[rand1], ggv[j].labels[rand2], result, prefactor);
single_gos_compare(ggm[ord][i].labels[rand1], ggm[ord][j].labels[rand2], result, prefactor);
if(result){break;}
}

if( result){
	// std::cout<<"Found equivalent graphs "<<i<<" "<<rand1<<" vs "<< j<<" "<<rand2<<std::endl;
prefactor=prefactor*ggm[ord][i].prefactor[rand1];
collapse_ggv(i, j,ggm[ord], prefactor);
j--;
start--;// this may have bug in it when j> start 
}

}

// }
// }

}	
}
	
	
}


void AmiGraph::gm_to_ggm_all( std::vector< std::vector< AmiGraph::graph_t>> &gm, gg_matrix_t &ggm, int mpi_rank){
ggm.clear();
ggm.resize(gm.size());

for(int ord=0; ord< gm.size(); ord++){
	graph_group gg;
for(int m=0; m< gm[ord].size(); m++){	


// gg_vec_t ggv;

gg.graph_vec.push_back(gm[ord][m]);
gg.prefactor.push_back(1);
// ggv.push_back(gg);

// ggv.clear();
// gg.graph_vec.clear();
// gg.prefactor.clear();	
}
if(gm[ord].size()>0){
ggm[ord].push_back(gg);	
}
}



}	


void AmiGraph::ggm_to_gm(  gg_matrix_t &ggm, std::vector< std::vector< AmiGraph::graph_t>> &gm, int min){

gm.clear();
gm.resize(ggm.size());


for(int ord=min; ord< gm.size(); ord++){
// int index=0;	
for(int m=0; m< ggm[ord].size(); m++){	
for(int graph=0; graph<ggm[ord][m].graph_vec.size(); graph++){

// ggv.push_back(gg);
gm[ord].push_back(ggm[ord][m].graph_vec[graph]);	
// ggv.clear();
// gg.graph_vec.clear();
// gg.prefactor.clear();
// index++;
}

for(int pair=0; pair<ggm[ord][m].gp_vec.size(); pair++){

gm[ord].push_back(ggm[ord][m].gp_vec[pair].g1_);
gm[ord].push_back(ggm[ord][m].gp_vec[pair].g2_);

}
}
}



}	

void AmiGraph::gm_to_ggm( std::vector< std::vector< AmiGraph::graph_t>> &gm, gg_matrix_t &ggm){
ggm.clear();
ggm.resize(gm.size());

for(int ord=0; ord< gm.size(); ord++){
	
for(int m=0; m< gm[ord].size(); m++){	

graph_group gg;
// gg_vec_t ggv;

gg.graph_vec.push_back(gm[ord][m]);
gg.prefactor.push_back(1);
// gg.prefactor.push_back(gm[ord][m].prefactor);
// ggv.push_back(gg);
ggm[ord].push_back(gg);	
// ggv.clear();
// gg.graph_vec.clear();
// gg.prefactor.clear();	
}
}



}	

void AmiGraph::collapse_ggv(int i, int j, gg_vec_t &ggv, double prefactor){

for( int k=0; k< ggv[j].graph_vec.size(); k++){	

ggv[i].graph_vec.push_back(ggv[j].graph_vec[k]);
ggv[i].prefactor.push_back(ggv[j].prefactor[k]*prefactor);
ggv[i].labels.push_back(ggv[j].labels[k]);	
	
}

ggv.erase(ggv.begin()+j);

}

void AmiGraph::gg_check(gg_matrix_t &ggm){
for(int ord=1; ord< ggm.size(); ord++){
std::vector<graph_t> gv;
linearize_ggv(ggm[ord], gv);

if(gv.size()>0){

for(int i=0; i< gv.size(); i++){
for(int j=i+1; j< gv.size(); j++){
std::cout<<ord<<" "<<i<<" "<<j <<std::endl;
if(full_iso(gv[i],gv[j])){
	std::cout<<"Seems to be duplicate isomorphic somehow"<<std::endl;
	std::cout<<std::endl;
	print_all_edge_info(gv[i]);
	std::cout<<std::endl;
	
	std::cout<<std::endl;
	print_all_edge_info(gv[j]);
	std::cout<<std::endl;
}

}
}
}
}


}

void AmiGraph::linearize_ggv(gg_vec_t &ggv_in, std::vector<graph_t> &gv_out){
	
	for(int i=0; i<ggv_in.size(); i++){
		for(int j=0; j< ggv_in[i].graph_vec.size(); j++){
			std::cout<<i<<" "<<j<<std::endl;
		gv_out.push_back(ggv_in[i].graph_vec[j]);
		}
	}
	
}



void AmiGraph::reduce_ggv_gos( gg_vec_t &ggv){

//gos_compare(graph_vector[j], graph_vector[k],result, prefactor);
bool result=false;
double prefactor;
for(int i=0; i< ggv.size(); i++){
for(int j=i+1; j< ggv.size(); j++){
	bool result=false;
// for each group - pick random graph to compare 
int rand1=random_int(0, ggv[i].graph_vec.size()-1);
int rand2=random_int(0, ggv[j].graph_vec.size()-1);

//single_gos_compare(ggv[i].graph_vec[rand1], ggv[j].graph_vec[rand2], result, prefactor);

single_gos_compare(ggv[i].labels[rand1], ggv[j].labels[rand2], result, prefactor);

if( result){
	// std::cout<<"Found equivalent graphs "<<i<<" "<<rand1<<" vs "<< j<<" "<<rand2<<std::endl;
prefactor=prefactor*ggv[i].prefactor[rand1];
collapse_ggv(i, j, ggv, prefactor);
j--;
}

}
}
	
}

void AmiGraph::single_gos_compare(labels_t &L1, labels_t &L2, bool &result, double &prefactor){
bool out;
result=false;
std::vector<int> count1, count2;

graph_t g1, g2;
int rand1=random_int(0, L1.size()-1);
int rand2=random_int(0, L2.size()-1);
g1=L1[rand1];
g2=L2[rand2];

			
label_counts(g1,count1);
label_counts(g2,count2);

if(count1.back() == count2.back()){
			
			std::sort(count1.begin(), count1.end());
			std::sort(count2.begin(), count2.end());
			
			if (count1==count2){
			
gos_compare(g1, g2,result, prefactor);
			
}
}
	
	
}

void AmiGraph::repeated_gos(int &max, bool &result, int &j, int &k , std::vector< AmiGraph::graph_t> &graph_vector,  graph_group &gg_temp, std::vector<int> &list){

std::vector<int> count1, count2;
// graph_group gg_temp;
// gg_vec_t gg_vec;
bool out;
double prefactor;
// bool result;

result=false;

int success=0;
int pass=0;

bool r1,r2;
r1=true;
r2=true;


for(int m=0; m< 6; m++){
	result=false;
	 // std::cout<<m<<std::endl;
	
			// if(m==0){
			// label_systematic(graph_vector[j]);	
			// check_momentum_conservation(graph_vector[j], r1);
			// label_systematic(graph_vector[k]);
			// check_momentum_conservation(graph_vector[k], r2);
// std::cout<<r1<<" - "<<r2<<std::endl;			
			// if(!r1 || !r2){continue;}	
			
			// }else{
			
			repeated_labelling(graph_vector[j],out);
			repeated_labelling(graph_vector[k],out);
			
			// }	

			label_counts(graph_vector[j],count1);
			label_counts(graph_vector[k],count2);			
			
			if(count1.back() == count2.back()){
			/////// sort counts
			pass++;
 
			std::sort(count1.begin(), count1.end());
			std::sort(count2.begin(), count2.end());
			///////

			if (count1==count2){
				success++;
				// std::cout<<"Found possible pair for graph "<<j<<" as "<< k<<std::endl;  
				
				
				gos_compare(graph_vector[j], graph_vector[k],result, prefactor);
				
			if(result){
				// std::cout<<"Found pair for graph "<<j<<" as "<< k<<std::endl; 

			if(gg_temp.graph_vec.size()==0){ gg_temp.graph_vec.push_back(graph_vector[j]);	
			gg_temp.graph_vec.push_back(graph_vector[k]);
			gg_temp.prefactor.push_back(1);gg_temp.prefactor.push_back(prefactor);
			list.push_back(k);

			// graph_vector.erase(graph_vector.begin()+k); k--; // remove in this order specifically
			// graph_vector.erase(graph_vector.begin()+j); j--; k--;

			}
			else{ gg_temp.graph_vec.push_back(graph_vector[k]); gg_temp.prefactor.push_back(prefactor);
			list.push_back(k);
			// graph_vector.erase(graph_vector.begin()+k); k--;

			}
			break;	
				
			}
			
if(success>max || pass > max*5){break;}			
				
			}	 
				
				
			}	
			
		}


}	

int AmiGraph::get_start_index(std::vector<int> list){
int result=0;
for(int i=0; i< list.size(); i++){
	
if(list[i]==1){return result;}
result=result+list[i];
}	
	
return -1;	
}


// TODO: add comparison of graphs between groups and merge groups as necessary

 void AmiGraph::gos(std::vector< std::vector< AmiGraph::graph_t>> &graph_matrix, gg_matrix_t &ggm, std::vector< std::vector<int>> &gg_lists){


// assign max labels // this doesn't seem to work for all cases
// std::cout<<"Finding max labels"<<std::endl;
// max_ext_label_counts(graph_matrix);
// std::cout<<"Finished labels"<<std::endl;
std::cout<<"-------------"<<std::endl;
std::cout<<"Beginning GOS graph reduction "<<std::endl;
std::cout<<"-------------"<<std::endl;


std::vector<int> list;
bool out;
int max=1;

// int current=0;

for(int ord=1; ord< graph_matrix.size(); ord++){

std::vector<int> count1, count2;
graph_group gg_temp;

std::vector<int> stash;

gg_vec_t gg_vec;
int start=get_start_index(gg_lists[ord]);
// std::cout<<"start is "<<start<<std::endl;
std::cout<<"On order "<<ord<<std::endl;


std::vector< AmiGraph::graph_t> graph_vector= graph_matrix[ord];
list.clear();


int index=0;
int current=0;

// now process the rest as before 
for(int j=0; j<graph_vector.size(); j++){
	
bool jnot= !in_list(j,list);

if(j==current){current=current+gg_lists[ord][index]; index++;
// std::cout<<"current is "<<current<<std::endl;
}

if(jnot){
	
for(int k=j+1; k<	graph_vector.size() ; k++){
	
	
	if(j<current && k < start && k>current){ k=start;}
	
	
	if(!in_list(k,list) && jnot){	
	// std::cout<<j<<" "<<k<<std::endl;

repeated_gos(max,out, j, k, graph_vector, gg_temp,  list);
		
		
		}
}// end k

// if(gg_temp.graph_vec.size()==0){
// gg_temp.graph_vec.push_back(graph_matrix[ord][j]);	
// gg_temp.prefactor.push_back(1);	
// }

if(gg_temp.graph_vec.size()==0){  // if this triggers it means one of the AID pairs didn't fit in the group. so stash it away and compare again against everything in gg_vec 
stash.push_back(j);
//gg_temp.graph_vec.push_back(graph_matrix[ord][j]);	
//gg_temp.prefactor.push_back(1);	
}


}
if(gg_temp.graph_vec.size()>0){
gg_vec.push_back(gg_temp);
gg_temp.graph_vec.clear();
gg_temp.prefactor.clear();
}


}// end j

// now try to deal with stashed entries

// TODO: this part is a total hackjob. even if it works it's a mess...
// std::cout<<"Stash size is "<< stash.size()<<std::endl;


gg_temp.graph_vec.clear();
gg_temp.prefactor.clear();
bool cont=true;


for(int j=0; j< stash.size(); j++){
	// std::cout<<"Working on stash element "<<j<<" of "<< stash.size()<<std::endl;
cont=true;
	for(int m=0; m< gg_vec.size(); m++){
		
		for(int r=0; r< gg_vec[m].graph_vec.size(); r++){
			// int r=0;
		std::vector< AmiGraph::graph_t> gv_temp;
		graph_t g1,g2;

		g1=graph_matrix[ord][stash[j]];
		g2=gg_vec[m].graph_vec[r];

		gv_temp.push_back(g1); 
		gv_temp.push_back(g2); 
		int a=0;
		int b=1;
		repeated_gos(max,out, b,a, gv_temp, gg_temp,  list);
		
				if(out){ 
				// std::cout<<"Stash item added to group "<<m<<std::endl;
				gg_vec[m].graph_vec.push_back(graph_matrix[ord][stash[j]]); gg_vec[m].prefactor.push_back(gg_temp.prefactor.back());
				cont=false;
				}
		if(!cont){break;}
			
		}
	// cont if

	if(!cont){break;}
	}

if(gg_temp.graph_vec.size()==0){  
gg_temp.graph_vec.push_back(graph_matrix[ord][j]);	
gg_temp.prefactor.push_back(1);	


gg_vec.push_back(gg_temp);

}

gg_temp.graph_vec.clear();
gg_temp.prefactor.clear();

} 

ggm[ord]=gg_vec;	
}	


}	


void AmiGraph::gos(std::vector< std::vector< AmiGraph::graph_t>> &graph_matrix, gg_matrix_t &ggm){


// assign max labels // this doesn't seem to work for all cases
// std::cout<<"Finding max labels"<<std::endl;
// max_ext_label_counts(graph_matrix);
// std::cout<<"Finished labels"<<std::endl;

std::vector<int> list;
bool out;
int max=1;

for(int ord=1; ord< graph_matrix.size(); ord++){

std::vector<int> count1, count2;
graph_group gg_temp;
gg_vec_t gg_vec;


// std::cout<<ord<<std::endl;


std::vector< AmiGraph::graph_t> graph_vector= graph_matrix[ord];
list.clear();
for(int j=0; j<graph_vector.size(); j++){
	bool jnot=!in_list(j,list);
	
if(jnot){	
for(int k=j+1; k<	graph_vector.size() ; k++){
	
	
	if(!in_list(k,list) && jnot){	
	// std::cout<<j<<" "<<k<<std::endl;

repeated_gos(max,out, j, k, graph_vector, gg_temp,  list);
		
		
		 /* for(int m=0; m< 50; m++){
			
			repeated_labelling(graph_vector[j],out);
			repeated_labelling(graph_vector[k],out);
			label_counts(graph_vector[j],count1);
			label_counts(graph_vector[k],count2);	
			
			if(count1.back() == count2.back()){
			/////// sort counts
			

			std::sort(count1.begin(), count1.end());
			std::sort(count2.begin(), count2.end());
			///////

			if (count1==count2){
				// std::cout<<"Found possible pair for graph "<<j<<" as "<< k<<std::endl;  
				double prefactor;
				bool result;
				
				gos_compare(graph_vector[j], graph_vector[k],result, prefactor);
				
			if(result){
				// std::cout<<"Found pair for graph "<<j<<" as "<< k<<std::endl; 

			if(gg_temp.graph_vec.size()==0){ gg_temp.graph_vec.push_back(graph_vector[j]);	
			gg_temp.graph_vec.push_back(graph_vector[k]);
			gg_temp.prefactor.push_back(1);gg_temp.prefactor.push_back(prefactor);
			list.push_back(k);

			// graph_vector.erase(graph_vector.begin()+k); k--; // remove in this order specifically
			// graph_vector.erase(graph_vector.begin()+j); j--; k--;

			}
			else{ gg_temp.graph_vec.push_back(graph_vector[k]); gg_temp.prefactor.push_back(prefactor);
			list.push_back(k);
			// graph_vector.erase(graph_vector.begin()+k); k--;

			}
			break;	
				
			}
				
				
			}	 
				
				
			}	
			
		} */ 

		}
}// end k

if(gg_temp.graph_vec.size()==0){
gg_temp.graph_vec.push_back(graph_matrix[ord][j]);	
gg_temp.prefactor.push_back(1);	
}

}

if(gg_temp.graph_vec.size()>0){
gg_vec.push_back(gg_temp);
gg_temp.graph_vec.clear();
gg_temp.prefactor.clear();
}



}// end j

ggm[ord]=gg_vec;	
}	
}

bool AmiGraph::in_list(int k, std::vector<int> list){
	
	for(int i=0; i< list.size(); i++){
		
	if(k==list[i]){return true;}	
		
	}
return false;	
}



void AmiGraph::max_ext_label_counts(std::vector< std::vector< AmiGraph::graph_t>> &graph_matrix){
	
std::vector<int> count1;
graph_t temp_g;



for(int ord=1; ord< graph_matrix.size(); ord++){
for(int i=0; i< graph_matrix[ord].size(); i++){
// std::cout<<ord<<" "<<i<<std::endl;
temp_g=graph_matrix[ord][i];
int max_count=0;
bool out;
// std::cout<<"Marker 1"<<std::endl;
int fail=0;


do{
	fail++;
	repeated_labelling(temp_g,out);
	// std::cout<<"Marker 2"<<std::endl;
	label_counts(temp_g,count1);
	// std::cout<<"Marker 3"<<std::endl;
	
	if(count1.back() > max_count){
		// std::cout<<"Updating max_count"<<count1.back()<<std::endl;
	graph_matrix[ord][i]=temp_g;	
	graph_matrix[ord][i][boost::graph_bundle].ext_counts=count1.back();
	max_count=count1.back();
	fail=0;
	}
}while(fail<10);
	
}
}	
	
	
	
}

//TODO: Seems to be a bug in this since swapping to signed char 

void AmiGraph::min_ext_label_counts(std::vector< std::vector< AmiGraph::graph_t>> &graph_matrix){
	
std::vector<int> count1;
graph_t temp_g;



for(int ord=1; ord< graph_matrix.size(); ord++){
for(int i=0; i< graph_matrix[ord].size(); i++){
// std::cout<<ord<<" "<<i<<std::endl;
temp_g=graph_matrix[ord][i];
int min_count=1000;
bool out;
// std::cout<<"Marker 1"<<std::endl;
int fail=0;


do{
	fail++;
	repeated_labelling(temp_g,out);
	// std::cout<<"Marker 2"<<std::endl;
	label_counts(temp_g,count1);
	// std::cout<<"Marker 3"<<std::endl;
	
	if(count1.back() < min_count){
		// std::cout<<"Updating min_count "<<count1.back() <<std::endl;
	graph_matrix[ord][i]=temp_g;	
	graph_matrix[ord][i][boost::graph_bundle].ext_counts=count1.back();
	min_count=count1.back();
	fail=0;
	}
}while(fail<30);
	
}
}	
	
	
	
}



void AmiGraph::gos_compare( graph_t &g1, graph_t &g2, bool &result, double &prefactor){

result=false;

Eigen::MatrixXd A1,A2;
alpha_to_matrix(g1,A1);
alpha_to_matrix(g2,A2);

// if(A1.cwiseAbs().isApprox(A2.cwiseAbs())){std::cout<<"Found equal within prefactor"<<std::endl;}else{std::cout<<"Not equal at the start"<<std::endl;}

	
Eigen::VectorXd r1,r2,c1,c2;
	
row_col_counts(A1,r1,c1);
row_col_counts(A2,r2,c2);


// first do a sort by both rows AND columns to see if this can be a candidate
bool first= row_col_compare(r1,r2,c1,c2);
if(!first){result=false; return;}
// next lets switch rows until ok
// start at i=0,if r1i != r2i, let  m=i+1, ri vs rm . if equal then swap.

// make permutation tracking identity matrices
Eigen::MatrixXd row_Perm=Eigen::MatrixXd::Identity(A1.rows(), A1.rows());
Eigen::MatrixXd col_Perm=Eigen::MatrixXd::Identity(A1.cols(), A1.cols());

gos_sort_rows(A1, A2, r1, r2, row_Perm);
gos_sort_cols(A1, A2, c1, c2, col_Perm);

//if(A1.array().abs()==A2.array().abs()){std::cout<<"Found equal within prefactor"<<std::endl;}
//if(A1.cwiseAbs().isApprox(A2.cwiseAbs())){std::cout<<"Found equal within prefactor "<<std::endl;
if(abs_equal(A1,A2)){ 
// std::cout<<"Found equal within prefactor "<<std::endl;
result=true;
prefactor=get_prefactor(A1,A2);
// std::cout<<"Prefactor is "<< prefactor<<std::endl;
// std::cout<<std::endl<<A1<<std::endl<<std::endl<<A2<<std::endl<<std::endl;
return;
}

// next need to check for equal permutations
std::vector< std::vector<int>> row_perm_set, col_perm_set;

get_rc_permutations(r1,c1, row_perm_set, col_perm_set);

permutation_set_t pst;
expand_permutations(row_perm_set,col_perm_set,pst);

double eval_pf;
bool p_compare= compare_permutation_set(A1,A2, pst, eval_pf);

if(p_compare){ 
// std::cout<<"Found equal within prefactor during permutations"<<std::endl;
result=true;
prefactor=eval_pf;
// std::cout<<"Prefactor is "<< prefactor<<std::endl;

// std::cout<<std::endl<<A1<<std::endl<<std::endl<<A2<<std::endl<<std::endl;
}
	
	
	
}


void AmiGraph::get_rc_permutations( Eigen::VectorXd &r1, Eigen::VectorXd &c1, std::vector< std::vector<int>> &row_perm_set, std::vector< std::vector<int>> &col_perm_set){
	
row_perm_set.clear();
col_perm_set.clear();	
	
std::vector<int> temp;
bool add=true;
for (int i=0; i< r1.rows(); i++){
	add=true;
for(int j=0; j< r1.rows(); j++){

for(int m=0; m<i; m++){ if(r1(i)==r1(m)){ add=false;} }

if(add==true && r1(i)==r1(j)){ temp.push_back(j);}


}
if(temp.size()>1){
row_perm_set.push_back(temp);}
temp.clear();	
	
}

// modified c1 -> c1-1 for git
// for for c1
for (int i=0; i< c1.rows()-1; i++){
	add=true;
for(int j=0; j< c1.rows()-1; j++){

for(int m=0; m<i; m++){ if(c1(i)==c1(m)){ add=false;} }

if(add==true && c1(i)==c1(j)){ temp.push_back(j);}


}
if(temp.size()>1){
col_perm_set.push_back(temp);}
temp.clear();	
	
}
	
	
	
	
}


void AmiGraph::expand_permutations(std::vector< std::vector<int>> &row_perm_set, std::vector< std::vector<int>> &col_perm_set, permutation_set_t &pst){
	
// first add all single permutations of rows 

pair_t temp_pair;
set_t temp_set;
rc_set_t temp_rc_set;

temp_pair.resize(2,-1);
//temp_rc_set.resize(2);

// add empty no-swap cases

for(int set=0; set< row_perm_set.size(); set++){

for(int i=0; i< row_perm_set[set].size(); i++){
for(int j=i+1; j< row_perm_set[set].size(); j++){

temp_pair[0]=row_perm_set[set][i];
temp_pair[1]=row_perm_set[set][j];

temp_set.push_back(temp_pair);
temp_rc_set.push_back(temp_set);
temp_set.clear();
temp_rc_set.push_back(temp_set); // this adds an empty column set 

pst.push_back(temp_rc_set);
temp_rc_set.clear();

	
}
}	
	
}


// add all single permutations of cols
for(int set=0; set< col_perm_set.size(); set++){

for(int i=0; i< col_perm_set[set].size(); i++){
for(int j=i+1; j< col_perm_set[set].size(); j++){

temp_pair[0]=col_perm_set[set][i];
temp_pair[1]=col_perm_set[set][j];

temp_set.clear();
temp_rc_set.push_back(temp_set); // this adds an empty row set 
temp_set.push_back(temp_pair);
temp_rc_set.push_back(temp_set);
 

pst.push_back(temp_rc_set);
temp_rc_set.clear();

	
}
}	
	
}

// from those, create all combined permutations 	
	
	
}

bool AmiGraph::compare_permutation_set(Eigen::MatrixXd &A1, Eigen::MatrixXd &A2, permutation_set_t &pst, double &pf){

bool result=false;
	
for(int i=0; i< pst.size(); i++){
// std::cout<<"On permutation "<<i<<std::endl;
result=compare_permutation(A1,A2, pst[i], pf);

if(result){
// apply_permutation(A2, pst[i]); // added nov 27 2019 - if successful then actually permute A2	
	return result;}

}	
	
	
return result;	
}


bool AmiGraph::compare_permutation(Eigen::MatrixXd A1, Eigen::MatrixXd &A2, rc_set_t &rcst, double &pf){

Eigen::MatrixXd A2c=A2;
apply_permutation(A2c, rcst);

if(abs_equal(A1,A2c)){
	A2=A2c;
pf=get_prefactor(A1,A2);	
	return true;}	
	
return false;	
}

// These permutations cannot be correct
void AmiGraph::apply_permutation(Eigen::MatrixXd &A2, rc_set_t &rcst){
	
for(int i=0; i<rcst[0].size(); i++){

A2.row(rcst[0][i][0]).swap(A2.row(rcst[0][i][1]));	

}

for(int i=0; i<rcst[1].size(); i++){

A2.col(rcst[1][i][0]).swap(A2.col(rcst[1][i][1]));	
	
}
	
	
}


// TODO: This prefactor is incorrect
double AmiGraph::get_prefactor(Eigen::MatrixXd &A1, Eigen::MatrixXd &A2){
// TODO: Ideally we would fix what follows, but for now, gos does nothing to prefactor is correct procedure	
	return 1.0;
	
	
double result=1;
double count=0;

if(equal(A1,A2)){return 1;}
if(neg_equal(A1,A2)){return -1;}


for(int i=0; i<A1.rows(); i++){

if(A1.row(i).isApprox(-1.0*A2.row(i))){count=count+1.0*(A1.row(i).count()+1);}

}

for(int i=0; i< A1.cols(); i++){
	
	if(A1.col(i).isApprox(-1.0*A2.col(i))){count=count+1*(A1.col(i).count());}
	
}


// for(int i=0; i<A1.rows(); i++){
// for(int j=0; j<A1.cols(); j++){

// if(A1(i,j)!=A2(i,j)){
	// count=count+1.0*(A1.row(i).count()+1);
// }


// }
// }

result=std::pow(-1.0,count);	
	
return result;	
	
}


void AmiGraph::gos_sort_rows(Eigen::MatrixXd &A1, Eigen::MatrixXd &A2, Eigen::VectorXd &r1, Eigen::VectorXd &r2, Eigen::MatrixXd &perm){

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

if(r1==r2){
	
A2=perm*A2;	
}else{std::cout<<"Rows not equal"<<std::endl;} // I think this probably doesn't work correctly. Don't think I've ever seen this printed...



}

void AmiGraph::gos_sort_cols(Eigen::MatrixXd &A1, Eigen::MatrixXd &A2, Eigen::VectorXd &r1, Eigen::VectorXd &r2, Eigen::MatrixXd &perm){
// std::cout<<"Beginning column sort with c1 and c2"<<std::endl;
// std::cout<<r1<<std::endl<<std::endl<<r2<<std::endl;

// not allowed to swap last column?

for(int i=0; i< A1.cols(); i++){
// std::cout<<"sorting on col "<<i<<std::endl;
if(r1(i)!= r2(i)){

for(int m=i+1; m< A1.cols(); m++){
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
}else{std::cout<<"cols not equal"<<std::endl;}



}	

bool AmiGraph::row_col_compare(Eigen::VectorXd r1, Eigen::VectorXd r2, Eigen::VectorXd c1, Eigen::VectorXd c2){
	
std::sort(r1.data(), r1.data()+r1.size());
std::sort(r2.data(), r2.data()+r2.size());
std::sort(c1.data(), c1.data()+c1.size());
std::sort(c2.data(), c2.data()+c2.size());


if(r1 != r2 || c1 != c2){ return false;}	
	
return true;	
}


void AmiGraph::row_col_counts(Eigen::MatrixXd &A_in, Eigen::VectorXd &row, Eigen::VectorXd &col){
	// row.clear();
	// col.clear();
	
	Eigen::VectorXd new_row;
	Eigen::VectorXd new_col;
	
new_row.resize(A_in.rows());
new_col.resize(A_in.cols());	
	
for(int i=0; i< A_in.rows(); i++){
new_row(i)=A_in.row(i).count();	
}

for(int i=0; i< A_in.cols(); i++){
new_col(i)= A_in.col(i).count();	
}	
	
row=new_row;
col=new_col;	
	
	
}


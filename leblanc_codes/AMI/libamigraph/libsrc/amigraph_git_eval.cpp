#include "amigraph.hpp"

void AmiGraph::ggm_populate_git_pairs(gg_matrix_t &ggm, int mpi_rank, int min){
if(mpi_rank==0){
std::cout<<"Forming GIT pairs"<<std::endl;
}

for (int ord=min; ord< ggm.size(); ord++){
	
	if(mpi_rank==0){
	std::cout<<"Working on Order "<< ord<<std::endl;
	std::cout<<ggm[ord].size()<<std::endl;
	}
for(int num=0; num< ggm[ord].size(); num++){
if(mpi_rank==0){
	std::cout<<"Working on group "<< num<<std::endl;}	

std::vector<int> map;



generate_git_map(ggm[ord][num],map);
std::cout<<"Finished generating map"<<std::endl;
if(map.size()==1){continue;}

// if(mpi_rank==0){
	// std::cout<<"Map is "<<std::endl;
	// for(int i=0; i<map.size(); ++i){
  // std::cout << map[i] << ' ';
  	// }
	// std::cout<<std::endl;
  // }

bool reduce= std::count(map.begin(), map.end(), 0)==0 && !std::equal(map.begin() + 1, map.end(), map.begin());


if ( reduce )
{
    //all equal
    // std::cout << "All elements are equal each other" << std::endl;
	//std::cout << "All elements are NOT equal each other" << std::endl;
	std::cout<<"Reducing ord num "<<ord<<" "<<num<<" with size "<<ggm[ord][num].graph_vec.size()<<std::endl;
	
	int counter=0;
	int start,end;
	do{
		start=ggm[ord][num].gp_vec.size();
	group_reduce_to_pairs(ggm[ord][num], map);
		end=ggm[ord][num].gp_vec.size();
		
		if(start==end){counter++;
		std::cout<<"Counter made it to "<<counter<<std::endl;}else{counter=0;}
	}while(counter<40); // TODO: this is just going through a couple of times. maybe something else should happen?
}


// group_reduce_to_pairs(ggm[ord][num]);

}
}	
	
}	

void AmiGraph::generate_git_map(graph_group &gg, std::vector<int> &map_out){
	
std::vector<int> map(gg.graph_vec.size(),0);
// seed the first to be the reference 'one'
// std::cout<<gg.graph_vec.size()<<" "<< map.size()<<std::endl;
// std::cout<<graph_order(gg.graph_vec[0])<<std::endl;
map[0]=get_prefactor(gg.graph_vec[0],graph_order(gg.graph_vec[0]));;

// if(gg.graph_vec.size()==1){return;}

std::cout<<"Starting map is "<<std::endl;
	for(int i=0; i<map.size(); ++i){
  std::cout << map[i] << ' ';
  	}
	std::cout<<std::endl;

bool done=false;

// std::count(map.begin(), map.end(), 0);

// get sign relative to first 	
int count=0;
do{
	count++;
	
for(int i=0; i < gg.graph_vec.size(); i++){	
for(int j=i+1; j< gg.graph_vec.size(); j++){
done=false;

std::cout<<i<<" "<<j<<std::endl;
if(map[i]!=0 && map[j]!=0){continue;}
if(map[i]==0 && map[j]==0){continue;} // if both zero not helpful

// pick a random label
int rand1,rand2;
rand1=random_int(0, gg.labels[i].size()-1);
rand2=random_int(0, gg.labels[j].size()-1);

std::cout<<"Rand 1 and 2 "<< rand1<<" "<<rand2<<std::endl;

double prefactor=1;

graph_t g1, g2;

g1=gg.labels[i][rand1];
g2=gg.labels[j][rand2];

graph_t g2_out;
git_perm_set_t pst;

int ord1=graph_order(g1);
int ord2=graph_order(g2);

double p1=get_prefactor(g1,ord1);
double p2=get_prefactor(g2,ord2);

git_compare_v2(g1, g2,done, prefactor,pst, g2_out);

std::cout<<"Done is "<<done<<" On i and j "<< i<<" "<< j<<" with p1 p2 and prefactor of "<< p1<<" "<<p2<<" "<<prefactor<<std::endl;


if(done){
prefactor=prefactor*p1*p2;

if(map[i]==0){
	// std::cout<<"Inserting value for i "<<map[j]*(int)(prefactor*p1*p2)<<std::endl;
	map[i]=map[j]*(int)(prefactor); prefactor=1;
	count=0;
	}
if(map[j]==0){
	// std::cout<<"Inserting value for j "<<map[i]*(int)(prefactor*p1*p2)<<std::endl;
	map[j]=map[i]*(int)(prefactor); prefactor=1;
	count=0;
	}
}


}


}

std::cout<<"Count is "<<count<<std::endl;
std::cout<<"Current map is "<<std::endl;
	for(int i=0; i<map.size(); ++i){
  std::cout << map[i] << ' ';
  	}
	std::cout<<std::endl;

}while(std::count(map.begin(), map.end(), 0) !=0 && count<30);

// std::cout<<"Zero count is "<< std::count(map.begin(), map.end(), 0)<<std::endl;

// if(std::count(map.begin(), map.end(), 1)==2){

// for(int i=0; i< gg.graph_vec.size(); i++){
// print_all_edge_info(gg.graph_vec[i]);
// std::cout<<"------"<<std::endl;
// }
// }	
	
	map_out=map;
	
}


void AmiGraph::group_reduce_to_pairs(graph_group &gg, std::vector<int> &map){

int pairs=0;
graph_t g2_out;
 git_perm_set_t pst;

for(int i=0; i<gg.graph_vec.size();i++){
	
bool done=false;	
	
for(int j=i+1; j< gg.graph_vec.size();j++){

if(map[i]*map[j]==1){continue;}


// now i and j are possible pairs. So cycle through their labels until you find one that worked.
graph_t g1, g2;

for(int l1=0; l1< gg.labels[i].size(); l1++){
if(done){break;}
for(int l2=0; l2< gg.labels[j].size(); l2++){
if(done){break;}

double prefactor=1;
// std::cout<<"L1 and L2 "<<l1<<" "<<l2<<std::endl;
g1=gg.labels[i][l1];
g2=gg.labels[j][l2];


git_compare_v2(g1, g2,done, prefactor,pst, g2_out);

}	
}


	
if(done){
	// std::cout<<"Found a pair"<<std::endl;
	pairs++;
g2=g2_out;
// put graphs into pair and push into gp vector
git_pair gp(g1,g2_out,pst);
gg.gp_vec.push_back(gp);

// remove all parts from the group 
gg_pair_remove(gg,i,j);
// need to remove map entries for the pairs
// std::cout<<"Map has size "<< map.size()<<std::endl;
map.erase(map.begin()+i);
map.erase(map.begin()+j-1);
// std::cout<<"Map has size "<< map.size()<<std::endl;

i--;
j--;
j--;	

break;		
	
}



}
}

}






/* git_perm_set_t pst;
double prefactor=1;

graph_t g2_out;

int pairs=0;
// first, take zeroth graph. generate signs relative to first graph


// for each i, check through
bool keep_do=true;
do{
int start=gg.graph_vec.size();	
for(int i=0; i < gg.graph_vec.size(); i++){	
std::cout<<"i  "<<i<<std::endl;



bool done=false;
for(int j=i+1; j< gg.graph_vec.size(); j++){	


// std::cout<<"i and j "<<i<<" "<<j<<std::endl;

// std::cout<<"Label sizes are "<<gg.labels[i].size()<<" "<<gg.labels[j].size()<<std::endl;

graph_t g1, g2;


// now for each label set in i and j, check for git compare - 
for(int l1=0; l1< gg.labels[i].size(); l1++){
if(done){break;}
for(int l2=0; l2< gg.labels[j].size(); l2++){
if(done){break;}

// std::cout<<"L1 and L2 "<<l1<<" "<<l2<<std::endl;
g1=gg.labels[i][l1];
g2=gg.labels[j][l2];

int ord1=graph_order(g1);
int ord2=graph_order(g2);

double p1=get_prefactor(g1,ord1);
double p2=get_prefactor(g2,ord2);

git_compare_v2(g1, g2,done, prefactor,pst, g2_out);

prefactor=prefactor*p1*p2;

}	
}

// if found pair AND they cancel. them make a cancelling pair and remove from the group.
if(done && prefactor==-1){
	std::cout<<"Found a pair"<<std::endl;
	pairs++;
g2=g2_out;
// put graphs into pair and push into gp vector
git_pair gp(g1,g2_out,pst);
gg.gp_vec.push_back(gp);

// remove all parts from the group 
gg_pair_remove(gg,i,j);
i--;
// j--;
// j--;	

break;	
}

		
}



	
}

int end=gg.graph_vec.size();
if(start==end){keep_do=false;}

}while(keep_do);

std::cout<<"Complete: Found "<<pairs<< " pairs!"<<std::endl; */

/* void AmiGraph::group_reduce_to_pairs(graph_group &gg){

git_perm_set_t pst;
double prefactor=1;

graph_t g2_out;

int pairs=0;


// for each i, check through
for(int i=0; i < gg.graph_vec.size(); i++){	
std::cout<<"i  "<<i<<std::endl;

for(int j=i+1; j< gg.graph_vec.size(); j++){	

bool done=false;
// std::cout<<"i and j "<<i<<" "<<j<<std::endl;

// std::cout<<"Label sizes are "<<gg.labels[i].size()<<" "<<gg.labels[j].size()<<std::endl;

graph_t g1, g2;


// now for each label set in i and j, check for git compare - 
for(int l1=0; l1< gg.labels[i].size(); l1++){
if(done){break;}
for(int l2=0; l2< gg.labels[j].size(); l2++){
if(done){break;}

// std::cout<<"L1 and L2 "<<l1<<" "<<l2<<std::endl;
g1=gg.labels[i][l1];
g2=gg.labels[j][l2];

int ord1=graph_order(g1);
int ord2=graph_order(g2);

double p1=get_prefactor(g1,ord1);
double p2=get_prefactor(g2,ord2);

git_compare_v2(g1, g2,done, prefactor,pst, g2_out);

prefactor=prefactor*p1*p2;

}	
}

// if found pair AND they cancel. them make a cancelling pair and remove from the group.
if(done && prefactor==-1){
	std::cout<<"Found a pair"<<std::endl;
	pairs++;
g2=g2_out;
// put graphs into pair and push into gp vector
git_pair gp(g1,g2_out,pst);
gg.gp_vec.push_back(gp);

// remove all parts from the group 
gg_pair_remove(gg,i,j);
i--;
// j--;
// j--;	

break;	
}

		
}	
}

std::cout<<"Complete: Found "<<pairs<< " pairs!"<<std::endl;

} */

// remove the i'th and j'th elements from the group 
void AmiGraph::gg_pair_remove(graph_group &gg, int i, int j){

if(gg.ss_vec.size()==gg.graph_vec.size()){
	// note that since i is always less than j, that once you remove i, j_. j-1 
gg.ss_vec.erase(gg.ss_vec.begin()+i);
gg.ss_vec.erase(gg.ss_vec.begin()+j-1);
}	

gg.graph_vec.erase(gg.graph_vec.begin()+i);
gg.graph_vec.erase(gg.graph_vec.begin()+j-1);

gg.labels.erase(gg.labels.begin()+i);
gg.labels.erase(gg.labels.begin()+j-1);
	
	
}

// void AmiGraph::git_group_labels(graph_group &gg, int i, int j, double prefactor, git_perm_set_t &pst){

	
	
	
	
// }



/* void AmiGraph::ggm_form_git_pairs(gg_matrix_t &ggm, int mpi_rank){
	
if(mpi_rank==0){
std::cout<<"Forming GIT pairs"<<std::endl;
}
	
for (int ord=1; ord< ggm.size(); ord++){
	
	if(mpi_rank==0){
	std::cout<<"Working on Order "<< ord<<std::endl;}

for(int num=0; num< ggm[ord].size(); num++){

for(int first=0; first< ggm[ord][num].graph_vec.size(); first++){
	
for(int second=first+1; second<ggm[ord][num].graph_vec.size(); second++){

// compare 1st and second, if finds negative, make a pair (remove from graph_vec) and continue on to the next (first) graph 

// non_destructive_git(ggm[ord][num].graph_vec[first],ggm[ord][num].graph_vec[second], l1, l2


}
}
}

}	
	
	
	
	
} */





/* 
void AmiGraph::ggm_assign_git_perms(gg_matrix_t &ggm, int mpi_rank){
if(mpi_rank==0){
std::cout<<"Assigning GIT permutations to groups"<<std::endl;
}
	
for (int ord=1; ord< ggm.size(); ord++){
	
	if(mpi_rank==0){
	std::cout<<"Working on Order "<< ord<<std::endl;}

for(int num=0; num< ggm[ord].size(); num++){

for(int m=0; m< ggm[ord][num].graph_vec.size(); m++){

ggm[ord][num].perm_vec.resize(ggm[ord][num].graph_vec.size());

if(m==0){continue;}
// systematic_label(
assign_git_perms(ggm[ord][num].graph_vec[0],ggm[ord][num].graph_vec[m], ggm[ord][num].perm_vec[m]);


}
}

}	

}


// TODO: this really only works if the graphs are loaded in, and have not previously been processed by git.
void AmiGraph::assign_git_perms(graph_t &g1, graph_t &g2, git_perm_set_t &pst){
bool result=false;
double prefactor;
graph_t g1c=g1;
graph_t g2c=g2;
int count=0;

do{
count++;
bool lab=true;
repeated_labelling(g2c, lab);

git_compare_v2( g1c, g2c, result, prefactor,pst);

}while(!result && count<10000);

if(!result){std::cout<<"Failed to find pst"<<std::endl;}

g2=g2c;

return;
} */


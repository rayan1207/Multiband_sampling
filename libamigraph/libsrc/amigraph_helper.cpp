#include "amigraph.hpp"


int AmiGraph::factorial(int n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

void AmiGraph::combinations_repetition_util(std::vector<int> &chosen, std::vector< std::vector<int>> &list, int index, int r, int start, int end){
	
if(index==r){

list.push_back(chosen);
return;
}	

for (int i=start; i<=end; i++){
	
	chosen[index]=i;
	combinations_repetition_util(chosen,list, index+1, r, i, end); 
	
}


	
	
}

void AmiGraph::combinations_repetition(int n, int r, std::vector< std::vector<int>> &list ){
	
std::vector<int> temp_chosen(r);

combinations_repetition_util( temp_chosen, list, 0, r, 0, n-1);  	
	
	
}

void AmiGraph::combinations(int max, int length, std::vector< std::vector<int>> &list){

list.clear();

int n, r;

n=max;
r=length;

   std::vector<bool> v(n);
   std::fill(v.begin(), v.begin() + r, true);
   
   // std::vector< std::vector<int>> list;

   do {
	   std::vector<int> current;
       for (int i = 0; i < n; ++i) {
           if (v[i]) {
               std::cout << (i + 1) << " ";
			   current.push_back(i);
           }
       }
       std::cout << "\n";
	   list.push_back(current);
   } while (std::prev_permutation(v.begin(), v.end()));	
	
	
	return;
}


bool AmiGraph::is_number(const std::string& s)
{
    return !s.empty() && std::find_if(s.begin(), 
        s.end(), [](unsigned char c) { return !std::isdigit(c); }) == s.end();
}

void AmiGraph::resize_U_matrix(int N, U_matrix_type &U){
U.clear();
	
U.resize(N);
for(int i=0; i< U.size(); i++){

U[i].resize(N);

for(int j=0; j< U[i].size(); j++){

U[i][j].resize(N);

for(int m=0; m<U[i][j].size(); m++){

U[i][j][m].resize(N,0.0);
}
}
}	
	
	
}

void AmiGraph::print_matrix(std::vector< std::vector <int>> &M){

for (int i=0; i< M.size(); i++){

std::cout<<"(";
for(int j=0; j< M[i].size(); j++){

std::cout<<M[i][j]<<",";
}

std::cout<<")"<<std::endl;

}	
	
	return;
}


// In AmiGraph
void AmiGraph::bubble_finder(graph_t &g, vertex_vector_list_t &bubble_vertex_list, edge_vector_list_t &bubble_edge_list, edge_vector_list_t &legs_edge_list){
	
edge_vector_t fermionic_edges;

find_fermionic_edges(g, fermionic_edges);

//std::cout<<"Found edges, Ne = " << fermionic_edges.size()<<std::endl; 

edge_vector_t loop_edges;
edge_vector_t leg_edges;
vertex_vector_t loop_vertices;

boost::graph_traits<graph_t>::in_edge_iterator ei, ei_end;
boost::graph_traits<graph_t>::out_edge_iterator eo, eo_end;


for (int i=0; i< fermionic_edges.size(); i++){
	for(int j=i; j< fermionic_edges.size(); j++){
		
		
	//	std::cout<< i<<" "<< j<<std::endl;
		//std::cout<<fermionic_edges[i]<<" "<< fermionic_edges[j]<<std::endl;
	//	std::cout<<"Probing edges "<< fermionic_edges[i] <<" " << fermionic_edges[j] <<std::endl;
	
// this only works if loop id's are rigorously maintained	TODO: Try to renable this later
	//if( g[fermionic_edges[i]].fermi_loop_id == g[fermionic_edges[j]].fermi_loop_id )
	//{
		
		if( source(fermionic_edges[i], g) == target(fermionic_edges[j],g) && target(fermionic_edges[i], g) == source(fermionic_edges[j],g)){
			
			// std::cout<<"logic triggered for "<< i<<" "<<j<<std::endl;
			// std::cout<<fermionic_edges[i]<<" "<< fermionic_edges[j]<<std::endl;
			// std::cout<<"for vertices "<< source(fermionic_edges[i], g)<<" "<< source(fermionic_edges[j], g)<<std::endl;
			
			// put edges and vertices (2 of each) in a small vector
			loop_edges.push_back(fermionic_edges[i]);
			loop_edges.push_back(fermionic_edges[j]);
			
			loop_vertices.push_back(source(fermionic_edges[i], g));
			loop_vertices.push_back(source(fermionic_edges[j], g));
			
			
			int k=0;
			// get bosonic adjacent edges 
			for (boost::tie(ei,ei_end) = in_edges(source(fermionic_edges[i], g), g); ei != ei_end; ++ei){
				k++;
				//std::cout<<k<<std::endl;
				if( g[*ei].g_struct_.stat_==AmiBase::Bose){
				//	std::cout<<"Found bosonic edge "<< *ei<<std::endl;
				leg_edges.push_back(*ei);
				}
			}
			
			for (boost::tie(eo,eo_end) = out_edges(source(fermionic_edges[i], g),g); eo != eo_end; ++eo){

				if( g[*eo].g_struct_.stat_==AmiBase::Bose){
				//	std::cout<<"Found bosonic edge "<< *eo<<std::endl;
				leg_edges.push_back(*eo);
				}
			}
			
			for (boost::tie(ei,ei_end) = in_edges(source(fermionic_edges[j], g),g); ei != ei_end; ++ei){

				if( g[*ei].g_struct_.stat_==AmiBase::Bose){
				//	std::cout<<"Found bosonic edge "<< *ei<<std::endl;
				leg_edges.push_back(*ei);
				}
			}
			
			for (boost::tie(eo,eo_end) = out_edges(source(fermionic_edges[j], g),g); eo != eo_end; ++eo){

				if( g[*eo].g_struct_.stat_==AmiBase::Bose){
			//		std::cout<<"Found bosonic edge "<< *eo<<std::endl;
				leg_edges.push_back(*eo);
				}
			}
			
			// std::cout<<"Final vector sizes are "<< loop_vertices.size()<<" "<<loop_edges.size()<<" "<<leg_edges.size()<<std::endl;
			
			// add those vectors to vector lists
			bubble_vertex_list.push_back(loop_vertices);
			bubble_edge_list.push_back(loop_edges);
			legs_edge_list.push_back(leg_edges);
			
			
			
			
			// immediately clear the inner loop vector lists before continuing
			loop_edges.clear();
			loop_vertices.clear();
			leg_edges.clear();
			
			
			
			
		}
		
		
		
		
	}

	
		
	//}
		
}
	
	
	
}

void AmiGraph::find_two_unique_edges(graph_t &g, AmiBase::stat_type requested_type, edge_t &one, edge_t &two){
	
	bool abort;
do{
abort=false;
one = random_edge(g, rand_gen); 
if ( g[one].g_struct_.stat_==AmiBase::Bose ){ abort=true;} // if it's a bosonic line, pick again.
}while(abort == true);

do{
abort=false;
two = random_edge(g, rand_gen); 
if ( g[two].g_struct_.stat_==AmiBase::Bose || one == two  ){ abort=true;} // if it's a bosonic line, pick again.
}while(abort == true);


	
}


// TODO : double check this works right.
void AmiGraph::find_two_edges(graph_t &g, AmiBase::stat_type requested_type, edge_t &one, edge_t &two){
	
	bool cont_do; // logic is 'continue_do_loop'
do{
cont_do=false;
one = random_edge(g, rand_gen); 
if ( g[one].g_struct_.stat_!= requested_type){ cont_do=true;} // if it's a bosonic line, pick again.
}while(cont_do == true);

do{
cont_do=false;
two = random_edge(g, rand_gen); 
if ( g[two].g_struct_.stat_!= requested_type ){ cont_do=true;} // if it's a bosonic line, pick again.
}while(cont_do == true);


	
}




void AmiGraph::find_fermionic_edges(graph_t &g, edge_vector_t &vector){

vector.clear();

//std::cout<< "Finding Fermionic edges" <<std::endl;
boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
// Looking through all edges to find source and targets : this could be improved?
for (boost::tie(ei,ei_end) = edges(g); ei != ei_end; ++ei){

if( g[*ei].g_struct_.stat_==AmiBase::Fermi){
vector.push_back(*ei);
}
//std::cout << "Edge " << *ei << " of stat_type "<< g[*ei].g_struct_.stat_ << std::endl;

//std::cout<<"edge "<<jg[*ei].edge_number_<<" loop Total is "<< total<<std::endl;
}




}

void AmiGraph::find_internal_vertices(graph_t &g, vertex_vector_t &vector){
	
vector.clear();

boost::graph_traits<graph_t>::vertex_iterator vi, vi_end;

for (boost::tie(vi,vi_end)=vertices(g); vi!= vi_end; ++vi){

if (degree(*vi,g) !=1){
	vector.push_back(*vi);
}

}	
	
	
}

void AmiGraph::find_internal_fermionic_edges(graph_t &g, edge_vector_t &vector){

vector.clear();

//std::cout<< "Finding Fermionic edges" <<std::endl;
boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
// Looking through all edges to find source and targets : this could be improved?
for (boost::tie(ei,ei_end) = edges(g); ei != ei_end; ++ei){

if( g[*ei].g_struct_.stat_==AmiBase::Fermi){
//if(source(*ei,g)==target(*ei,g)){vector.push_back(*ei);}
if(degree(source(*ei,g),g) != 1 && degree(target(*ei,g),g)!=1){	
vector.push_back(*ei);
}

}
//std::cout << "Edge " << *ei << " of stat_type "<< g[*ei].g_struct_.stat_ << std::endl;

//std::cout<<"edge "<<jg[*ei].edge_number_<<" loop Total is "<< total<<std::endl;
}




}







void AmiGraph::find_rand_internal_edge(graph_t &g, AmiBase::stat_type requested_type, edge_t &one){
	
	bool cont_do; // logic is 'continue_do_loop'
do{
cont_do=false;
one = random_edge(g, rand_gen); 
if ( g[one].g_struct_.stat_!= requested_type || degree(source(one,g),g)==1 || degree(target(one,g),g)==1){ cont_do=true;} // if it's not the requested stat_type , pick again.
}while(cont_do == true);



	
}

void AmiGraph::find_one_fermi_bose_edge(graph_t &g, edge_t &bose, edge_t &fermi){
	
boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
bool bbose,bfermi;
bbose=false;
bfermi=false;

for (boost::tie(ei,ei_end) = edges(g); ei != ei_end; ++ei){


if( g[*ei].g_struct_.stat_==AmiBase::Bose){
bose=*ei; //vector.push_back(*ei);
bbose=true;
}

if( g[*ei].g_struct_.stat_==AmiBase::Fermi){
fermi=*ei; //vector.push_back(*ei);
bfermi=true;
}

if( bbose && bfermi){return;}


}
	
	
	
}


void AmiGraph::mpi_print_ggm( gg_matrix_t &ggm, int mpi_rank){

std::cout<<"Graph group size is "<<std::endl;

for(int ord=0; ord< ggm.size(); ord++){

int check=0;
int gpcheck=0;

for(int m=0; m < ggm[ord].size();m++){
for(int i=0; i< ggm[ord][m].graph_vec.size(); i++){

// std::cout<<ord<<" "<<m<<" "<<i<<std::endl;
check++;
if(ord>2){
// std::cout<<"Graph ord "<<ord<<" group "<< m <<" num "<<i<<" has n_labels="<<ggm[ord][m].labels[i].size()<<std::endl;
}
}
for(int i=0; i< ggm[ord][m].gp_vec.size(); i++){

// std::cout<<ord<<" "<<m<<" "<<i<<std::endl;
gpcheck++;
if(ord>2){
// std::cout<<"Graph ord "<<ord<<" group "<< m <<" num "<<i<<" has n_labels="<<ggm[ord][m].labels[i].size()<<std::endl;
}
}


}	
	
std::cout<<"There were "<<check<<" graphs and "<< gpcheck <<" pairs at order "<< ord<< " in "<< ggm[ord].size()<<" groups on rank "<< mpi_rank <<std::endl;	
}	
	
	
	
}


void AmiGraph::print_ggm( gg_matrix_t &ggm){

std::cout<<"Graph group size is "<<std::endl;

for(int ord=0; ord< ggm.size(); ord++){

int check=0;
int gpcheck=0;

for(int m=0; m < ggm[ord].size();m++){
for(int i=0; i< ggm[ord][m].graph_vec.size(); i++){

// std::cout<<ord<<" "<<m<<" "<<i<<std::endl;
check++;
if(ord>2){
// std::cout<<"Graph ord "<<ord<<" group "<< m <<" num "<<i<<" has n_labels="<<ggm[ord][m].labels[i].size()<<std::endl;
}
}

for(int i=0; i< ggm[ord][m].gp_vec.size(); i++){

// std::cout<<ord<<" "<<m<<" "<<i<<std::endl;
gpcheck++;
if(ord>2){
// std::cout<<"Graph ord "<<ord<<" group "<< m <<" num "<<i<<" has n_labels="<<ggm[ord][m].labels[i].size()<<std::endl;
}
}


}	
	
std::cout<<"There were "<<check<<" graphs and "<< gpcheck <<" pairs at order "<< ord<< " in "<< ggm[ord].size()<<" groups." <<std::endl;	
// std::cout<<"There were "<<check<<" graphs of order "<< ord<< " in "<< ggm[ord].size()<<" groups."<<std::endl;	
}	
	
	
	
}


void AmiGraph::find_bose_fermi_edges(graph_t &g, edge_vector_t &bose, edge_vector_t &fermi){

bose.clear();
fermi.clear();

//std::cout<< "Finding bosonic edges" <<std::endl;
boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
// Looking through all edges to find source and targets : this could be improved?
for (boost::tie(ei,ei_end) = edges(g); ei != ei_end; ++ei){

if( g[*ei].g_struct_.stat_==AmiBase::Bose){
bose.push_back(*ei);
}

if( g[*ei].g_struct_.stat_==AmiBase::Fermi){
fermi.push_back(*ei);
}
//std::cout << "Edge " << *ei << " of stat_type "<< g[*ei].g_struct_.stat_ << std::endl;

//std::cout<<"edge "<<jg[*ei].edge_number_<<" loop Total is "<< total<<std::endl;
}


	
	
}


void AmiGraph::find_bosonic_edges(graph_t &g, edge_vector_t &vector){

vector.clear();

//std::cout<< "Finding bosonic edges" <<std::endl;
boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
// Looking through all edges to find source and targets : this could be improved?
for (boost::tie(ei,ei_end) = edges(g); ei != ei_end; ++ei){

if( g[*ei].g_struct_.stat_==AmiBase::Bose){
vector.push_back(*ei);
}
//std::cout << "Edge " << *ei << " of stat_type "<< g[*ei].g_struct_.stat_ << std::endl;

//std::cout<<"edge "<<jg[*ei].edge_number_<<" loop Total is "<< total<<std::endl;
}



}

void AmiGraph::find_non_tadpole_bosonic_edges(graph_t &g, edge_vector_t &vector){

vector.clear();

//std::cout<< "Finding bosonic edges" <<std::endl;
boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
// Looking through all edges to find source and targets : this could be improved?
for (boost::tie(ei,ei_end) = edges(g); ei != ei_end; ++ei){

if( g[*ei].g_struct_.stat_==AmiBase::Bose){

vertex_t v;
edge_t e=*ei;
v=source(e,g);	
bool one=edge(v,v,g).second;	
v=target(e,g);
bool two=edge(v,v,g).second;

if(!one && !two){
vector.push_back(*ei);
}


}
//std::cout << "Edge " << *ei << " of stat_type "<< g[*ei].g_struct_.stat_ << std::endl;

//std::cout<<"edge "<<jg[*ei].edge_number_<<" loop Total is "<< total<<std::endl;
}



}


void AmiGraph::find_non_entry_bosonic_edges(graph_t &g, edge_vector_t &vector){

vector.clear();

//std::cout<< "Finding bosonic edges" <<std::endl;
boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
// Looking through all edges to find source and targets : this could be improved?
for (boost::tie(ei,ei_end) = edges(g); ei != ei_end; ++ei){

if( (g[*ei].g_struct_.stat_==AmiBase::Bose) && degree(source(*ei,g),g)!=1){
vector.push_back(*ei);
}
//std::cout << "Edge " << *ei << " of stat_type "<< g[*ei].g_struct_.stat_ << std::endl;

//std::cout<<"edge "<<jg[*ei].edge_number_<<" loop Total is "<< total<<std::endl;
}




}

void AmiGraph::find_non_external_bosonic_edges(graph_t &g, edge_vector_t &vector){

vector.clear();

//std::cout<< "Finding bosonic edges" <<std::endl;
boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
// Looking through all edges to find source and targets : this could be improved?
for (boost::tie(ei,ei_end) = edges(g); ei != ei_end; ++ei){

if( (g[*ei].g_struct_.stat_==AmiBase::Bose) && degree(source(*ei,g),g)!=1 && degree(target(*ei,g),g)!=1){
vector.push_back(*ei);
}
//std::cout << "Edge " << *ei << " of stat_type "<< g[*ei].g_struct_.stat_ << std::endl;

//std::cout<<"edge "<<jg[*ei].edge_number_<<" loop Total is "<< total<<std::endl;
}




}

void AmiGraph::put_back_legs(graph_t &g, vertex_vector_t &in_vv, vertex_vector_t &out_vv){
int alpha_size=0;
int eps_size=0;
// std::cout<<"Legs"<<std::endl;
do{
	// std::cout<<"In do"<<std::endl;

edge_t this_edge=random_edge(g,rand_gen);
	
alpha_size=g[this_edge].g_struct_.alpha_.size();
eps_size=g[this_edge].g_struct_.eps_.size();

// std::cout<<alpha_size<<" "<<eps_size<<std::endl;
	
}while(alpha_size==0 || eps_size==0);	

// std::cout<<"Determined sizes "<< alpha_size<< " and "<<eps_size<<std::endl;	
for(int i=0; i< in_vv.size(); i++){
	
vertex_t new_vert = add_vertex(g);
edge_t this_edge=add_edge(new_vert,in_vv[i], edge_info(AmiBase::Bose),g).first;

g[this_edge].g_struct_.alpha_.resize(alpha_size,0);
g[this_edge].g_struct_.eps_.resize(eps_size,0);
g[this_edge].label=labelled;

}	


for(int i=0; i< out_vv.size(); i++){
	
vertex_t new_vert = add_vertex(g);
edge_t this_edge=add_edge(out_vv[i],new_vert, edge_info(AmiBase::Bose),g).first;

g[this_edge].g_struct_.alpha_.resize(alpha_size,0);
g[this_edge].g_struct_.eps_.resize(eps_size,0);
}	
	
number_vertices(g);	
}

void AmiGraph::find_path_between_vertices(graph_t &g, vertex_t &v1, vertex_t &v2, bool &success, edge_vector_t &final_ev, vertex_vector_t &final_vv){
 success=false;
 bool keep_going=true;
 
 
edge_vector_list_t next_evl; 
edge_vector_list_t evl; 
vertex_vector_list_t vvl; 
vertex_vector_list_t next_vvl; 


vertex_vector_t v1start;

v1start.push_back(v1);
vvl.push_back(v1start);
evl.resize(vvl.size());

int steps=0;
do{
  
  steps++;
next_vvl.clear();
next_evl.clear();
	
for(int vv=0; vv< vvl.size(); vv++){

vertex_t current=vvl[vv].back();
vertex_t next;
edge_t e;

vertex_vector_t newvv=vvl[vv];
edge_vector_t newev=evl[vv];

get_adjacent_vertex_stat(g, current, next,e, AmiBase::Fermi);

// if(next==v1){keep_going=false; success=false;}
// std::cout<<"Next is "<<g[next].index_<<std::endl;
bool add=true;
for(int k=0; k<newvv.size(); k++){
  if(newvv[k]==next){ add=false;}
  
}

if(add){
newvv.push_back(next);
newev.push_back(e);

next_vvl.push_back(newvv);
next_evl.push_back(newev);
}

newvv=vvl[vv];
newev=evl[vv];

get_adjacent_vertex_stat(g, current, next,e, AmiBase::Bose);

// if(next==v1){keep_going=false; success=false;}
// if(next!=newvv[0]){
newvv.push_back(next);
newev.push_back(e);

next_vvl.push_back(newvv);
next_evl.push_back(newev);
// }


}

// std::cout<<next_vvl.size()<<std::endl;

vvl=next_vvl;
evl=next_evl;  

for(int i=0; i<vvl.size(); i++){
 
if (vvl[i].back()==v2){
keep_going=false;
final_vv=vvl[i];
final_ev=evl[i];
success=true;
}  
  
}

// std::cout<<steps<<" "<<vvl.size()<<std::endl;

if (steps>100){ keep_going=false; success=false;}

  
}while(keep_going); 


// std::cout<<"Keep going is "<< keep_going<<" success is "<< success<<std::endl; 

// std::cout<<"Vertex list is "<<std::endl;

// for (int i=0; i< final_vv.size(); i++){
// std::cout<<g[final_vv[i]].index_<<std::endl;
// }  

// std::cout<<"Edge list is "<<std::endl;

// for (int i=0; i< final_ev.size(); i++){
// std::cout<<g[source(final_ev[i],g)].index_ <<"-"<<g[target(final_ev[i],g)].index_ <<std::endl;
// }

 
  
}


void AmiGraph::fix_force_labels(graph_t &g, vertex_vector_t &in_vv, vertex_vector_t &out_vv){
 
// std::cout<<"Trying to find a path on graph "<<std::endl;
// print_all_edge_info(g);

edge_vector_t ev;
vertex_vector_t vv;
bool success=false;

for(int i=0; i< in_vv.size(); i++){
  for(int j=0; j<out_vv.size(); j++){
    // std::cout<<i<<" "<<j<<std::endl;
find_path_between_vertices(g, in_vv[i], out_vv[j], success, ev, vv);
if(success){break;}
  }
  if(success){break;}
}

if(success){

for(int i=0; i< ev.size(); i++){
 
g[ev[i]].g_struct_.alpha_.back()=1; 
  
}
}else{ throw std::runtime_error("Failed to find force line");}
 
  
}

/* void AmiGraph::fix_force_labels(graph_t &g, vertex_vector_t &in_vv, vertex_vector_t &out_vv){
  
  std::cout<<"Fixing for labels on graph "<<std::endl;
  print_all_edge_info(g);
  
  std::cout<<"In vertices are "<< g[in_vv[0]].index_<<" "<<g[in_vv[1]].index_<<std::endl;
  std::cout<<"out vertices are "<< g[out_vv[0]].index_<<" "<<g[out_vv[1]].index_<<std::endl;
	
bool keep_going=true;
bool success=false;

edge_vector_t ev;

for (int m=0; m< in_vv.size(); m++){
  keep_going=true;
  success=false;
ev.clear();
vertex_t current=in_vv[m];
vertex_t next;
edge_t e;

std::cout<<"Starting on "<<g[in_vv[m]].index_<<std::endl;

do{
	
get_adjacent_vertex_stat(g, current, next,e, AmiBase::Fermi);
ev.push_back(e);
std::cout<<"Moved to "<<g[next].index_<<std::endl;
// g[e].g_struct_.alpha_.back()=1;

for(int i=0; i< out_vv.size(); i++){
	
	if(next==out_vv[i]){
    success=true;
		keep_going=false;
		break;
	}
  
  if(next==in_vv[m]){
    success=false;
		keep_going=false;
		break;
	}
	
}

current=next;	
	
	
}while(keep_going);

if(success){break;}

} // end m loop 	

// At this point it was succesful so we found a path and can set the external lines=1

if(success){

for(int i=0; i< ev.size(); i++){
 
g[ev[i]].g_struct_.alpha_.back()=1; 
  
}
}else{ throw std::runtime_error("Failed to find force line");}	
	
	
}
 */
// void AmiGraph::fix_force_labels(graph_t &g, vertex_vector_t &in_vv, vertex_vector_t &out_vv){
	
// bool keep_going=true;

// vertex_t current=in_vv[0];
// vertex_t next;
// edge_t e;

// do{
	
// get_adjacent_vertex_stat(g, current, next,e, AmiBase::Fermi);
// g[e].g_struct_.alpha_.back()=1;

// for(int i=0; i< out_vv.size(); i++){
	
	// if(next==out_vv[i]){
		// keep_going=false;
		// break;
	// }
	
// }

// current=next;	
	
	
// }while(keep_going);
	
	
	
	
// }

void AmiGraph::find_force_LR_vertices(graph_t &g, vertex_vector_t &in_vv, vertex_vector_t &out_vv){
	
in_vv.clear();
out_vv.clear();

boost::graph_traits<graph_t>::edge_iterator ei, ei_end;

for(boost::tie(ei,ei_end)=edges(g); ei!= ei_end; ++ei){

if(g[*ei].g_struct_.stat_==AmiBase::Bose){


if(degree(source(*ei,g),g)==1){
	
	in_vv.push_back(target(*ei,g));
	
}

if(degree(target(*ei,g),g)==1){
	
	out_vv.push_back(source(*ei,g));
	
}





}	
	
}

	
	return;
}

// this function is specifically for force diagrams
void AmiGraph::find_mutable_bosonic_edges(graph_t &g, edge_vector_t &vector){

vector.clear();

// step 1 find  in and out vertices 

vertex_vector_t in_vv, out_vv;
find_force_LR_vertices(g,in_vv, out_vv);
// step 2 exclude any bosonic line whose source or target is those .


//std::cout<< "Finding bosonic edges" <<std::endl;
boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
// Looking through all edges to find source and targets : this could be improved?
for (boost::tie(ei,ei_end) = edges(g); ei != ei_end; ++ei){

bool add=true;
if( (g[*ei].g_struct_.stat_==AmiBase::Bose) ){

for(int i=0; i< in_vv.size(); i++){
	if(source(*ei,g)==in_vv[i] || target(*ei,g)==in_vv[i]){
		add=false;
	}
		
}


for(int i=0; i< out_vv.size(); i++){
	if(source(*ei,g)==out_vv[i] || target(*ei,g)==out_vv[i]){
		add=false;
	}
		
}


if(add){ 
vector.push_back(*ei);
}


}


//std::cout << "Edge " << *ei << " of stat_type "<< g[*ei].g_struct_.stat_ << std::endl;

//std::cout<<"edge "<<jg[*ei].edge_number_<<" loop Total is "<< total<<std::endl;
}




}

void AmiGraph::reset_visited(graph_t &g){

boost::graph_traits<graph_t>::vertex_iterator vi, vi_end;

for(boost::tie(vi,vi_end)=vertices(g); vi!=vi_end; ++vi){

g[*vi].visited=0;

}	
	
	
}

double AmiGraph::get_prefactor(graph_t &g, int order){

// std::cout<<"In get prefactor with order "<<order<<std::endl;
// print_all_edge_info(g);

double output;
int fermi_loops=fermi_connected_components(g);//count_fermi_loops(g);
if(ami_parameters.TYPE_==AmiBase::Pi_ppud || ami_parameters.TYPE_==AmiBase::Pi_ppuu){
  fermi_loops--;
}
if(ami_parameters.TYPE_==AmiBase::Sigma || ami_parameters.TYPE_==AmiBase::Greens || ami_parameters.TYPE_==AmiBase::density || ami_parameters.TYPE_==AmiBase::DOS || ami_parameters.TYPE_==AmiBase::Bare){
  if(order!=0){
  fermi_loops--;
  }
}

// if(ami_parameters.TYPE_==AmiBase::Greens ||  && order==0){
  // fermi_loops++;
// }


int sigma_ct=g[boost::graph_bundle].sigma_ct_count;

// std::cout<<"Found fermi loops"<< fermi_loops<<std::endl;

output=pow(-1, fermi_loops+order+sigma_ct);	
	
// std::cout<<"returning "<<output<<std::endl;
	
return output;	
}

void AmiGraph::find_start_end_vertices( graph_t &g, vertex_t &start, vertex_t &end){

boost::graph_traits<graph_t>::edge_iterator ei, ei_end;	
	
for (boost::tie(ei,ei_end)=edges(g); ei!= ei_end; ++ei){

if ( g[*ei].g_struct_.stat_== AmiBase::Bose &&  out_degree(target(*ei,g),g)==0  ){ 
 end=source(*ei,g);
}


if ( g[*ei].g_struct_.stat_== AmiBase::Bose &&  in_degree(source(*ei,g),g)==0  ){ 
 start=target(*ei,g);
}



}	
	
	return;
	
	
}


int AmiGraph::find_end_index(graph_t &g, vertex_t &v){

if(graph_type!=AmiBase::Pi_ppud){throw std::runtime_error("Function only works for ppud graphs");}	
		
boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
	
int result=-1;

bool found=false;


for (boost::tie(ei,ei_end) = edges(g); ei != ei_end; ++ei){
	if ( g[*ei].g_struct_.stat_== AmiBase::Bose &&  out_degree(target(*ei,g),g)==0  ){ 
	
	result=g[source(*ei,g)].index_;
	// std::cout<<"Returning result "<<result<<std::endl;
	v=source(*ei,g);

found=true;	

	}
	
	
}







if(!found){ throw std::runtime_error("Didn't find starting index for loop counting - in function");}
	
return result;	
	
	
}




// TODO: This may not work for all cases - needs a test 
// For bubbles this needs to happen BEFORE legs are removed 
int AmiGraph::find_start_index(graph_t &g, vertex_t &v){
	
	
// std::cout<<"Finding start index"<<std::endl;	
	// print_all_edge_info(g);
	
	boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
	
	// boost::graph_traits<graph_t>::vertex_iterator vi, vi_end;
int result=-1;

bool found=false;

// if a bubble then any starting vertex should work 
if(graph_type==AmiBase::FORCE || graph_type==AmiBase::doubleocc){
v=random_vertex(g,rand_gen);
result=g[v].index_;
return result;
	
}
// if( graph_type==AmiBase::Pi_phud || graph_type==AmiBase::Pi_phuu || graph_type==AmiBase::doubleocc ){
	
// v=random_vertex(g,rand_gen);
// result=	g[v].index_;

// return result;	
	
// }

if( graph_type==AmiBase::Pi_phud || graph_type==AmiBase::Pi_phuu || graph_type==AmiBase::doubleocc || graph_type==AmiBase::Pi_ppud || graph_type==AmiBase::Pi_ppuu ){


for (boost::tie(ei,ei_end) = edges(g); ei != ei_end; ++ei){
	// std::cout<<"On "<<g[source(*ei,g)].index_<<" with degrees "<< in_degree(source(*ei,g),g)<<" "<<degree(source(*ei,g),g)<<std::endl;
	if ( g[*ei].g_struct_.stat_== AmiBase::Bose &&  in_degree(source(*ei,g),g)==0  ){ 
	
	result=g[target(*ei,g)].index_;
	// std::cout<<"Returning result "<<result<<std::endl;
	v=target(*ei,g);

found=true;	
return result;

	}
	
	
}
}





		


for (boost::tie(ei,ei_end) = edges(g); ei != ei_end; ++ei){
	// std::cout<<"On "<<g[source(*ei,g)].index_<<" with degrees "<< in_degree(source(*ei,g),g)<<" "<<degree(source(*ei,g),g)<<std::endl;
	if (in_degree(source(*ei,g),g)==0  ){ 
	
	result=g[source(*ei,g)].index_;
	// std::cout<<"Returning result "<<result<<std::endl;
	v=source(*ei,g);

found=true;	
break;
	}
	
	
}

if(!found){
	
for (boost::tie(ei,ei_end) = edges(g); ei != ei_end; ++ei){
	std::cout<<"On edge with bmap "<< g[*ei].band_element.first <<" "<<  g[*ei].band_element.second<< " and source "<< g[source(*ei,g)].index_<<" and target "<< g[target(*ei,g)].index_ <<std::endl;
	if (degree(source(*ei,g),g)==1  ){ 
	
	std::cout<<"In here"<<std::endl;
	result=g[target(*ei,g)].index_;
	v=target(*ei,g);

found=true;
break;	
	}
	
	
}		
	
	
	
}

	

if(!found){ throw std::runtime_error("Didn't find starting index for loop counting - in function");}
	
return result;	
	
	
}

void AmiGraph::print_fullstate(NewAmiCalc::internal_state state){
	
std::cout<<"State dim="<<state.dim_<<std::endl;
std::cout<<"State order="<<state.order_<<std::endl;

std::cout<<"K bounds are "<< state.mink_<<" - "<<state.maxk_<<std::endl;
std::cout<<"Dispersion is "<<state.disp_<<std::endl;
	
std::cout<<"Printing internal k list "<<std::endl;

for (int i=0; i< state.internal_k_list_.size(); i++){
std::cout<<"k"<<i<<":(";
for(int x=0; x< state.internal_k_list_[i].size(); x++){
std::cout<<state.internal_k_list_[i][x]<<" ";

}
std::cout<<"), "<<std::endl;
}	
	
	
}


int AmiGraph::count_fermi_loops(graph_t &g){
	
// if(num_edges(g)==1 && graph_type==AmiCalc::density){

// return 1;

// }	

// std::cout<<"Counting loops - reset first "<<std::endl;
// reset visited
reset_visited(g);

graph_t graph; // copy graph for safety

graph=g;
// std::cout<<"Make a copy"<<std::endl;

vertex_t start_v;

// || (graph_type== AmiBase::doubleocc)
 if(graph_type==AmiBase::Pi_phuu || graph_type==AmiBase::Pi_phud || graph_type== AmiBase::doubleocc || ami_parameters.TYPE_==AmiBase::Pi_ppud || ami_parameters.TYPE_==AmiBase::Pi_ppuu || graph_type==AmiBase::FORCE ){
	
// need to remove external legs and vertices 

vertex_vector_t extern_vect_list;
edge_vector_t extern_edge_list;
find_external_vertices(graph, extern_vect_list, extern_edge_list);



// std::cout<<"Deleting legs"<<std::endl;
int index;
if(graph_type==AmiBase::FORCE){
	delete_legs(graph,extern_vect_list, extern_edge_list);
	index=find_start_index(graph, start_v);
}else{
	index=find_start_index(graph, start_v);

delete_legs(graph,extern_vect_list, extern_edge_list);
	
}

 
// std::cout<<"remaining vertices "<<std::endl;
// print_all_edge_info(graph);


// std::cout<<"Using start index of "<< index<<std::endl;

// boost::graph_traits<graph_t>::vertex_iterator vi, vi_end;

// for(boost::tie(vi,vi_end)=vertices(graph); vi!=vi_end; ++vi){

// std::cout<<"Index is "<< graph[*vi].index_<<std::endl;

// }	
 

if(graph[start_v].index_!=index && graph_type !=AmiBase::doubleocc){
throw std::runtime_error("Didn't find the starting index ");
}	

}	else{
	
find_start_index(graph, start_v);	
}
	
	
// std::cout<<"Start index is "<< graph[start_v].index_<<std::endl;	
	
// }

//print_all_edge_info(graph);
//int counter=0;

boost::graph_traits<graph_t>::edge_iterator ei, ei_end;

edge_vector_t edge_list;
for (boost::tie(ei,ei_end) = edges(graph); ei != ei_end; ++ei){

if(graph[*ei].g_struct_.stat_==AmiBase::Bose){
	
edge_list.push_back(*ei);
}

// if(graph[*ei].g_struct_.stat_==AmiBase::Fermi && graph[*ei].g_struct_.species_== 1){
	
// edge_list.push_back(*ei);
// }

	
}

// std::cout<<"Found all the bosonic edges"<<std::endl;

for(int i=0; i<edge_list.size();i++){
	
	remove_edge(edge_list[i],graph);
}
// std::cout<<"Graph with no bosonic lines now looks like "<<std::endl;
// print_all_edge_info(graph);

number_vertices(graph);
vertex_vector_t vert_list;

// std::cout<<"Starting with 
vertex_t current=start_v;//=vert_list[0]; // this will not work if vertex zero is no longer the external vertex 
vertex_t next;
bool instop=false;
bool outstop=false;
int counter=0;
bool first=true;
do{
	 // std::cout<<"Collecting non visited"<<std::endl;
collect_unvisited_vertices(graph, vert_list);
 // std::cout<<"Complete "<< vert_list.size()<<std::endl;
 
if(vert_list.size()==0){break;}else{

if(	!first){
current=vert_list[0];
}
}
instop=false;
first=false;

do{
 // std::cout<<"Source is "<<graph[current].index_<<std::endl;	

graph[current].visited=1;
get_adjacent_vertex_stat(graph, current, next, AmiBase::Fermi); 
// std::cout<<"Source is "<< graph[current].index_<<" next is "<<graph[next].index_<<std::endl;
 
 if( graph[next].visited==1){instop=true; graph[next].visited=1;};
 current=next;

} while(!instop);
 
 counter++;
} while(!outstop); 
	
g[boost::graph_bundle].num_loops=counter-1;	
	
// std::cout<<"Counted fermi loops is	"<< g[boost::graph_bundle].num_loops<<std::endl;

if((graph_type== AmiBase::Pi_phuu)	 ||  (graph_type== AmiBase::Pi_phud) || (graph_type== AmiBase::doubleocc)  || ami_parameters.TYPE_==AmiBase::Pi_ppud || ami_parameters.TYPE_==AmiBase::Pi_ppuu ||  ami_parameters.TYPE_==AmiBase::FORCE){// || graph_type==AmiCalc::density){

g[boost::graph_bundle].num_loops=counter;
return 	counter;
}
else{
	return counter-1;
}
	
}

void AmiGraph::delete_legs(graph_t &g, vertex_vector_t &v, edge_vector_t &edges){
	
for(int i=0; i< v.size(); i++){
clear_vertex(v[i],g);
remove_vertex(v[i],g);

}	
	
	
	
}


// this function looks at a bubble graph and pulls the k and k' vectors entering and leaving the vertex. 
// TODO: Checked that this does not work for pp vertex diagrams.  so needs fixing. - this is probably not the right one anyways...
void AmiGraph::get_kkp(std::vector<AmiBase::alpha_t> &kkp, AmiGraph::graph_t &g){

kkp.clear();

AmiGraph::vertex_t start_v;


int index=find_start_index(g, start_v);
	
boost::graph_traits<graph_t>::edge_iterator ei, ei_end;

for(boost::tie(ei,ei_end)=edges(g); ei!=ei_end; ++ei){

if(g[*ei].g_struct_.stat_==AmiBase::Fermi){

if(source(*ei,g)==start_v || target(*ei,g)==start_v){

kkp.push_back(g[*ei].g_struct_.alpha_);

}	

}

}	
	
// std::cout<<"Found kkp size "<< kkp.size()<<std::endl;	
	
}

void AmiGraph::get_pp_kkp(std::vector<AmiBase::alpha_t> &kkp, AmiGraph::graph_t &g){


kkp.clear();

AmiGraph::vertex_t start_v, end_v;

 
int index=find_start_index(g, start_v);
int end_index=find_end_index(g,end_v);

// std::cout<<"Start and end indices are "<< index <<" "<< end_index<<std::endl;	
	
boost::graph_traits<graph_t>::edge_iterator ei, ei_end;

for(boost::tie(ei,ei_end)=edges(g); ei!=ei_end; ++ei){

if(g[*ei].g_struct_.stat_==AmiBase::Fermi){

if(source(*ei,g)==start_v){
// std::cout<<"Found for "<< g[source(*ei,g)].index_<<"-"<<g[target(*ei,g)].index_<<std::endl;
kkp.push_back(g[*ei].g_struct_.alpha_);
break;
}	

// if( target(*ei,g)==end_v){
// std::cout<<"Found for "<< g[source(*ei,g)].index_<<"-"<<g[target(*ei,g)].index_<<std::endl;
// kkp.push_back(g[*ei].g_struct_.alpha_);
	
// }



}

}	
	
// now we have the start edge 

vertex_t current=start_v;
vertex_t next;
bool cont=true;
do{

edge_t e;

get_adjacent_vertex_stat(g, current, next,e, AmiBase::Fermi);

if(next==end_v){
 kkp.push_back(g[e].g_struct_.alpha_); 
 cont=false;
}

current=next;

}while(cont);	
  
  
	
	
	
}

// this function looks at a bubble graph and pulls the k entering the left of the bubble and k' vector that leaves the right vertex 
// in the case of pp bubble this will not work correctly!
void AmiGraph::get_cond_kkp(std::vector<AmiBase::alpha_t> &kkp, AmiGraph::graph_t &g){

kkp.clear();

AmiGraph::vertex_t start_v, end_v;


int index=find_start_index(g, start_v);
int end_index;
	
boost::graph_traits<graph_t>::edge_iterator ei, ei_end;

for(boost::tie(ei,ei_end)=edges(g); ei!=ei_end; ++ei){
if ( g[*ei].g_struct_.stat_== AmiBase::Bose &&  out_degree(target(*ei,g),g)==0  ){ 
	
	end_index=g[source(*ei,g)].index_;
	// std::cout<<"Returning result "<<result<<std::endl;
	end_v=source(*ei,g);

	}
	
}	

// std::cout<<"Start and end indices are "<< index <<" "<< end_index<<std::endl;	
	


for(boost::tie(ei,ei_end)=edges(g); ei!=ei_end; ++ei){

if(g[*ei].g_struct_.stat_==AmiBase::Fermi){

if(source(*ei,g)==end_v){
// std::cout<<"Found for "<< g[source(*ei,g)].index_<<"-"<<g[target(*ei,g)].index_<<std::endl;
kkp.push_back(g[*ei].g_struct_.alpha_);

}	

if( target(*ei,g)==start_v){
// std::cout<<"Found for "<< g[source(*ei,g)].index_<<"-"<<g[target(*ei,g)].index_<<std::endl;
kkp.push_back(g[*ei].g_struct_.alpha_);
	
}



}

}	
	
	
	
}




void AmiGraph::collect_unvisited_vertices(graph_t &g, vertex_vector_t &vertex_list){
	
vertex_list.clear();	
boost::graph_traits<graph_t>::vertex_iterator vi, vi_end;

for(boost::tie(vi,vi_end)=vertices(g); vi!=vi_end; ++vi){
	
if(g[*vi].visited==0){
vertex_list.push_back(*vi);
}	
	
}
	
	
}



int AmiGraph::graph_order(graph_t &g){
edge_vector_t edges;

find_internal_edges_stat(g, edges, AmiBase::Bose);

// edge_vector_t fedges;
// find_internal_edges_stat(g,fedges,AmiBase::Fermi);

// int eff_bose=0;
// for(int i=0; i<fedges.size(); i++){

// if(g[fedges[i]].g_struct_.species_==1){eff_bose++;}
  
// }
// std::cout<<"In graph order function we have "<<edges.size()<<" "<<eff_bose<<std::endl;

return edges.size(); //+eff_bose;
	
}

void AmiGraph::find_internal_edges_stat(graph_t &g, edge_vector_t &vector, AmiBase::stat_type requested){
	
	
	
vector.clear();

//std::cout<< "Finding Fermionic edges" <<std::endl;
boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
// Looking through all edges to find source and targets : this could be improved?
for (boost::tie(ei,ei_end) = edges(g); ei != ei_end; ++ei){

if( g[*ei].g_struct_.stat_==requested){
//if(source(*ei,g)==target(*ei,g)){vector.push_back(*ei);}
if(degree(source(*ei,g),g) != 1 && degree(target(*ei,g),g)!=1){	
vector.push_back(*ei);
}

}
//std::cout << "Edge " << *ei << " of stat_type "<< g[*ei].g_struct_.stat_ << std::endl;

//std::cout<<"edge "<<jg[*ei].edge_number_<<" loop Total is "<< total<<std::endl;
}


	
	
}


void AmiGraph::get_bubble_type(graph_t &g,bubble_type &b_type, AmiGraph::vertex_vector_t bubble_verts, AmiGraph::edge_vector_t bubble_edges){

bool one_in;
bool two_in;	

if(	target(bubble_edges[0],g)==bubble_verts[0] || target(bubble_edges[0],g)==bubble_verts[1]){one_in=true;}else{one_in=false;}
if(	target(bubble_edges[1],g)==bubble_verts[0] || target(bubble_edges[1],g)==bubble_verts[1]){two_in=true;}else{two_in=false;} 
	
if(	one_in && two_in){b_type=both_in;}else
if( one_in != two_in){b_type=directed;}else
if( one_in == false && two_in==false){b_type=both_out;}
	
	
}

void AmiGraph::find_tadpoles(graph_t &g, vertex_vector_t &tp_vec, vertex_vector_t &tp_conn_vec, edge_vector_t &tp_bose_edges, edge_vector_t &edge_a, edge_vector_t &edge_b ){

tp_vec.clear();
tp_conn_vec.clear();
tp_bose_edges.clear();
edge_a.clear();
edge_b.clear();
	
boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
boost::graph_traits<graph_t>::edge_iterator ei2, ei_end2;
edge_t bose;

for (boost::tie(ei,ei_end) = edges(g); ei != ei_end; ++ei){

if(source(*ei,g)==target(*ei,g)){
tp_vec.push_back(source(*ei,g));
find_bosonic_edge_on_three_pointvert(g, source(*ei,g), bose);
tp_bose_edges.push_back(bose);

// std::cout<<"found bose edge from "<< g[source(bose,g)].index_<<"-"<<g[target(bose,g)].index_<<std::endl;

if(source(bose,g)==source(*ei,g)){ tp_conn_vec.push_back(target(bose,g));}
if (target(bose,g)==source(*ei,g)){ tp_conn_vec.push_back(source(bose,g));}

for (boost::tie(ei2,ei_end2) = edges(g); ei2 != ei_end2; ++ei2){	

if(g[*ei2].g_struct_.stat_==AmiBase::Fermi){
	
	

  if(target(*ei2,g)==tp_conn_vec.back()){ edge_a.push_back(*ei2);}// std::cout<<"pushback for a"<< std::endl;}
  if(source(*ei2,g)==tp_conn_vec.back()){ edge_b.push_back(*ei2);}//std::cout<<"pushback for b"<< std::endl;}
	
}
	
}	

}	



	
	
}
	

	
}

void AmiGraph::find_tadpoles(graph_t &g, vertex_vector_t &tp_vec, vertex_vector_t &tp_conn_vec, edge_vector_t &tp_bose_edges){

tp_vec.clear();
tp_conn_vec.clear();
tp_bose_edges.clear();
	
boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
edge_t bose;

for (boost::tie(ei,ei_end) = edges(g); ei != ei_end; ++ei){

if(source(*ei,g)==target(*ei,g)){
tp_vec.push_back(source(*ei,g));
find_bosonic_edge_on_three_pointvert(g, source(*ei,g), bose);
tp_bose_edges.push_back(bose);

if(source(bose,g)==source(*ei,g)){ tp_conn_vec.push_back(target(bose,g));}else
if (target(bose,g)==source(*ei,g)){ tp_conn_vec.push_back(source(bose,g));}

}	
	
	
}
	
	
}

// NOTE: Due to force diagrams, need to not associate an external bosonic line with the actual tadpole 
void AmiGraph::find_bosonic_edge_on_three_pointvert(graph_t &g, vertex_t v, edge_t &bose_edge){


boost::graph_traits<graph_t>::in_edge_iterator ei, ei_end;
boost::graph_traits<graph_t>::out_edge_iterator eo, eo_end;

for (boost::tie(ei,ei_end) = in_edges(v,g); ei != ei_end; ++ei){

if( g[*ei].g_struct_.stat_==AmiBase::Bose){
	
	
	

bose_edge=*ei;



}
}

for (boost::tie(eo,eo_end) = out_edges(v,g); eo != eo_end; ++eo){

if( g[*eo].g_struct_.stat_==AmiBase::Bose){
	
bose_edge=*eo;
	
}
}


}



int AmiGraph::random_int( int min, int max){

std::uniform_int_distribution<int> int_dist(min,max);

return int_dist(rand_gen);

}

// random generator function:
// int AmiGraph::myrandom (int i) { return random_int(0,i-1);}


double AmiGraph::random_real(double max){ // produce real number from 0 to max 


return rand_dist(rand_gen)*max;

}

double AmiGraph::random_real(double min, double max){ // produce real number from 0 to max 


return rand_dist(rand_gen)*(max-min)+min;

}

std::vector<double>  AmiGraph::sobol_random_real(int size, int order){

	std::vector<double> sample(size);
	quasi_random_gen_t gen(engine_vec[order], boost::uniform_01<double>());

    // At this point you can use std::generate, generate member f-n, etc.
    std::generate(sample.begin(), sample.end(), gen);
    // engine.generate(sample.begin(), sample.end());
return 	sample;
	
}

std::vector<double>  AmiGraph::sobol_random_real(){

	std::vector<double> sample(12);
	quasi_random_gen_t gen(engine, boost::uniform_01<double>());

    // At this point you can use std::generate, generate member f-n, etc.
    std::generate(sample.begin(), sample.end(), gen);
    // engine.generate(sample.begin(), sample.end());
return 	sample;
	
}

void AmiGraph::find_in_out_edges_oftype(graph_t &g, vertex_t &v, edge_t &in_edge_output,edge_t &out_edge_output, AmiBase::stat_type stat_value){
	
	
boost::graph_traits<graph_t>::in_edge_iterator ei, ei_end;
boost::graph_traits<graph_t>::out_edge_iterator eo, eo_end;


for (boost::tie(ei,ei_end) = in_edges(v, g); ei != ei_end; ++ei){
				//std::cout<<"working"<<std::endl;
				
					// std::cout<< "vertex is "<<v <<std::endl;
				  // std::cout<< *ei << " "<< g[*ei].g_struct_.stat_ <<std::endl;
				// std::cout<< "source is "<< source(*ei,g) <<std ::endl;
				// std::cout<< "Target is "<< target(*ei,g) <<std ::endl;
				
				// std::cout<< "vertex is "<<v <<std::endl;
				
				// std::cout<< "Boolean comparison "<< bool(int(v)==int(target(*ei,g))) <<std::endl;
				
				
				if(target(*ei,g)==v ){
					// std::cout<< "vertex is "<<v <<std::endl;
				  // std::cout<< *ei << " "<< g[*ei].g_struct_.stat_ <<std::endl;
				if( g[*ei].g_struct_.stat_==stat_value){
				in_edge_output=*ei;
				}
				}
				
				// if(int(source(*ei,g))==int(v) ){
				// if( g[*ei].g_struct_.stat_==stat_value){
				// out_edge_output=*ei;
				// }
				// }
				
				// if(source(g[*ei],g)==g[v]){
				// if( g[*ei].g_struct_.stat_==stat_value){
				// out_edge_output=*ei;
				// }
				// }
				
				
}
			
 for (boost::tie(eo,eo_end) = out_edges(v,g); eo != eo_end; ++eo){
	 
	 	// std::cout<< "vertex is "<<v <<std::endl;
				  // std::cout<< *eo << " "<< g[*eo].g_struct_.stat_ <<std::endl;
				// std::cout<< "source is "<< source(*eo,g) <<std ::endl;
				// std::cout<< "Target is "<< target(*eo,g) <<std ::endl;

		 if(source(*eo,g)==v ){
			 	// std::cout<< "vertex is "<<v <<std::endl;
				  // std::cout<< *eo << " "<< g[*eo].g_struct_.stat_ <<std::endl;
				if( g[*eo].g_struct_.stat_==stat_value){
				out_edge_output=*eo;
				}
				}
		
 }


	
	
	
}











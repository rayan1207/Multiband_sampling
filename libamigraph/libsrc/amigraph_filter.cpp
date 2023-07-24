#include "amigraph.hpp"

// needed filters

// is diagram connected? this one is a bit tricky - dfs?

void AmiGraph::graph_matrix_update(graph_t &g, std::vector< std::vector<graph_t>> &graph_matrix){

// a bunch of do nothing conditions

// std::cout<<"In update function"<<std::endl;

// If the graph is representing a hubbard interaction graph then apply the hubbard filter 
if(ami_parameters.int_type_==AmiBase::hubbard){
	
if( !is_hubbard(g)){
// std::cout<<"Rejecting due to non-Hubbard"<<std::endl;
// number_vertices(g);	
// print_all_edge_info(g);
	return;}
 
}
 

// with the exception of green's function and density expansions - filter out one particle reducible graphs
bool isg=(ami_parameters.TYPE_==AmiBase::Greens)|| (ami_parameters.TYPE_==AmiBase::density) || (ami_parameters.TYPE_==AmiBase::DOS)|| (ami_parameters.TYPE_==AmiBase::ENERGY);
if( !isg ){ 
if( is_1P_reducible_graph(g)){
	
std::cout<<"Rejecting due to reducibility"<<std::endl;
// number_vertices(g);	
// print_all_edge_info(g);
return;	
}
}

	
int ord=AmiGraph::graph_order(g);

bool add=true;


// TODO: is_isomorphic function is not generalized 

if(graph_matrix[ord].size() !=0){
	
for(int i=0; i< graph_matrix[ord].size(); i++){
// std::cout<<"Checking isomorphism against "<<i<<" of total "<<graph_matrix[ord].size()<<std::endl;
if( is_isomorphic(g, graph_matrix[ord][i])){ add=false;
// std::cout<<"Found isomorphic "<<std::endl;
break;
}
	
}

	
}


if(add){graph_matrix[ord].push_back(g);
// std::cout<<"added!"<<std::endl;
}

	
}


int AmiGraph::count_bubbles(graph_t &g){

vertex_vector_list_t bubble_vertex_list;
edge_vector_list_t bubble_edge_list;
edge_vector_list_t legs_edge_list;
bubble_finder(g, bubble_vertex_list,  bubble_edge_list, legs_edge_list);	

return bubble_vertex_list.size();	
	
	
	
}


void AmiGraph::reorder_AID(std::vector< std::vector<graph_t>> &graph_matrix, std::vector< std::vector<int>> &gg_lists, int mpi_rank){
	if(mpi_rank==0){
	std::cout<<"---------"<<std::endl<<"Reordering AID Pairs"<<std::endl<<"---------"<<std::endl;}	
gg_lists.clear();
std::vector< std::vector<graph_t>> temp_graph_matrix;
temp_graph_matrix.resize(graph_matrix.size());
gg_lists.resize(graph_matrix.size());

std::vector<int> list;
std::vector<int> stash;

int count=0;
int index=0;

int numpairs=0;

// just duplicate the first couple of orders
for(int i=0; i< 2; i++){
	
for(int m=0; m< graph_matrix[i].size(); m++){	
temp_graph_matrix[i].push_back(	graph_matrix[i][m]);
}	
}

// sort remaining orders
for(int i=2; i< graph_matrix.size(); i++){
	list.clear();
	stash.clear();
	if(mpi_rank==0){
	std::cout<<"On order "<<i<<" ";}
	numpairs=0;

for(int m=0; m< graph_matrix[i].size(); m++){	
	count=0;
	if(!in_list(m,list)){ 

   for(int k=m+1; k< graph_matrix[i].size(); k++){
	
	if(!in_list(k,list)){
	if(is_AID(graph_matrix[i][m], graph_matrix[i][k])){ numpairs++;
		
		if(count==0){
			temp_graph_matrix[i].push_back(graph_matrix[i][m]); list.push_back(m); count++;
			temp_graph_matrix[i].push_back(graph_matrix[i][k]); list.push_back(k); count++;
		}else{
		temp_graph_matrix[i].push_back(graph_matrix[i][k]); list.push_back(k); count++;
		}
		// std::cout<<"Found AID for "<<m<<" in graph "<<k<<std::endl;
	}
	
	
	}
	
	
   }
   if(count!=0){
   gg_lists[i].push_back(count);
   }
   if(count==0){stash.push_back(m);}
}	


}

// std::cout<<"stash has size "<<stash.size()<<std::endl;

for(int m=0;m<stash.size(); m++){
temp_graph_matrix[i].push_back(graph_matrix[i][stash[m]]);
gg_lists[i].push_back(1);
}
if(mpi_rank==0){
std::cout<<"Found numpairs = "<< numpairs<<std::endl;}
}
	
	
// set equal to eachother	
graph_matrix=temp_graph_matrix;	
	
}

bool AmiGraph::edges_equal(graph_t &g1, graph_t &g2){

if(num_edges(g1)!= num_edges(g2)){return false;}

std::vector< std::vector<int>> e1, e2;

get_labelled_edge_list(e1,g1);
get_labelled_edge_list(e2,g2);

for(int i=0; i< e1.size(); i++){

if(e1[i]!=e2[i]){ return false;} 	
	
}

	
return true;	
}


// this function will not work for pp type pi graphs - since the principle loop is not closed. therefore we just number the vertices in whatever order. Not sure why this was necessary really other than for debugging. 
void AmiGraph::systematic_vertex_label( graph_t &g){
	
// 1) Label ingoing and outgoing lines first.
// 2) follow fermionic loops 	

if(graph_type==AmiBase::Pi_ppud || graph_type == AmiBase::Pi_ppuu || graph_type== AmiBase::FORCE){

number_vertices(g);	
	
	
return;	
}
	

// reset visited
reset_visited(g);


number_vertices(g);
// print_all_edge_info(g);
vertex_vector_t fermi_vert_list;

vertex_vector_t vert_list;
collect_unvisited_vertices(g, vert_list);
//
//vertex_vector_t current_list;
//vertex_vector_t next_list;

// TODO: By default lets try picking any vertex - and the next few sections will fix this if it can 

vertex_t start;

int index=find_start_index(g,start);

vertex_t current=start; //vert_list[0]; // this will not work if vertex zero is no longer the external vertex 
vertex_t next;

// first check if bosonic lines are symmetrized
boost::graph_traits<graph_t>::edge_iterator ei, ei_end;

// bool sym=false;
int deg=1;

for (boost::tie(ei,ei_end)=edges(g); ei!=ei_end; ++ei){
	
if(	g[*ei].g_struct_.stat_==AmiBase::Bose){ 

// &&( edge(source(*ei,g),target(*ei,g),g).second && edge(target(*ei,g),source(*ei,g),g).second ) 
if(degree(source(*ei,g),g)==1 || degree(source(*ei,g),g)==2 ){
	deg=degree(source(*ei,g),g);
	}

}
	
	
}





// check if ingoing and outgoing lines are fermionic or bosonic 



bool isbose=false;
for (boost::tie(ei,ei_end)=edges(g); ei!=ei_end; ++ei){
	// std::cout<<g[source(*ei,g)].index_<<" "<< g[target(*ei,g)].index_<<std::endl;
if(degree(source(*ei,g),g)==deg && (source(*ei,g) != target(*ei,g) )){
	// std::cout<<"Found start"<<std::endl;
	if(g[*ei].g_struct_.stat_==AmiBase::Bose){ isbose=true;}
	
	current=source(*ei,g); 
	
	// std::cout<<"Found current"<<std::endl;
	break;}
	
}

// std::cout<<"IS Bose is "<< isbose<<std::endl;

// if(isbose){
	
	
	
// }

// TODO: everything above does not work if it is a density or double occupancy graph 
// TODO: it would be best to never use this for density graphs 

// if( graph_type==AmiBase::density || graph_type==AmiBase::doubleocc){

// current=source(random_edge(g,rand_gen),g);	// pick the source of a random edge 
	
// }





bool instop=false;
bool outstop=false;
int counter=0;
do{
	 // std::cout<<"Collecting non visited"<<std::endl;
collect_unvisited_vertices(g, vert_list);
 // std::cout<<"Complete "<< vert_list.size()<<std::endl;
if(vert_list.size()==0){break;}else{
get_next_fermi_loop(fermi_vert_list, current, g);
//current=vert_list[0];
}
instop=false;
do{
 // std::cout<<"Source is "<<g[current].index_<<std::endl;	

g[current].visited=1;
g[current].index_=counter;
counter++;
fermi_vert_list.push_back(current);
// TODO: This only works if the main line is fermionic ??
bool success=false;
get_adjacent_vertex_stat(g, current, next, AmiBase::Fermi,success);
if(!success){break;}
// std::cout<<"Success was "<< success<<std::endl; 
// std::cout<<"Source is "<< g[current].index_<<" next is "<<g[next].index_<<std::endl;
 
 if( g[next].visited==1){instop=true; g[next].visited=1;};
 current=next;

} while(!instop);
// std::cout<<"Exit instop"<<std::endl;
 
 //counter++;
} while(!outstop); 
// std::cout<<"Exit outstop"<<std::endl;

	
//g[boost::graph_bundle].num_loops=counter-1;	
	
// std::cout<<"Counted fermi loops is	"<< g[boost::graph_bundle].num_loops<<std::endl;
	

	
	
}

void AmiGraph::get_next_fermi_loop(vertex_vector_t &list, vertex_t &current, graph_t &g){

vertex_t next;	
for(int i=0; i< list.size(); i++){
bool success=false;
get_adjacent_vertex_stat(g, list[i], next, AmiBase::Bose, success);
if( success){

if( g[next].visited==0){

current=next;
break;
}	
	
}

}	
	
	
	
}





// label all visited=0
// start at 0 vertex
// get adjacent vertices
// label visited 
// for each of those, get adjacent, label visited - add to a list
// repeat
// end condition on each thread is hitting a visted adjacent vertex 

// list of adjacents to add - which becomes a list to find adjacents of 
// if any vertex not visited then not connected

// this is basically a breadth first search 



// 1p reducible in ANY channel? Was this intentional?
// update it removes only fermionic lines so it is fine 
bool AmiGraph::is_1P_reducible_graph(graph_t g){
	// std::cout<<"In 1p red"<<std::endl;

// edge_t in, out;
edge_vector_t in, out;
	
boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
// Step through all edges and collect the alpha and epsilon labels of the 'g' graph
for (boost::tie(ei,ei_end) = edges(g); ei != ei_end; ++ei){

//g[*ei].g_struct_.stat_==AmiBase::Bose){

// if(g[*ei].g_struct_.stat_==AmiBase::Fermi){


if ( degree( source(*ei,g),g)==1 && (source(*ei,g) != target(*ei,g) )){ 
in.push_back(*ei);

}

if( degree( target(*ei,g),g)==1 && (source(*ei,g) != target(*ei,g))){
	out.push_back(*ei);
}

// }

}

for(int i=0; i<in.size(); i++){
remove_edge(in[i],g);
remove_vertex(source(in[i],g),g);
}
for(int i=0; i<out.size(); i++){
remove_edge(out[i],g);
remove_vertex(target(out[i],g),g);
}


	// number_vertices(g);
	// print_all_edge_info(g);

	
edge_vector_t internal_fermi, internal_copy;

find_fermionic_edges(g, internal_fermi);

graph_t copy;


for (int i=0; i<internal_fermi.size(); i++){
	// std::cout<<i<<std::endl;
copy=g;
find_fermionic_edges(copy, internal_copy);
// std::cout<<"Test removing edge "<<print_edge(internal_copy[i],g)<<std::endl;	
remove_edge(internal_copy[i],copy);

// std::cout<<std::endl;
// print_all_edge_info(copy);
// std::cout<<std::endl;

// std::cout<<"Checking if connected"<<std::endl;

// NOTE: Either !is_connected_graph or connected components can work here. there are comparible in time complexity it seems 
if(!is_connected_graph(copy)){

// if(connected_components(copy)!=1){
	// std::cout<<std::endl;
	// number_vertices(copy);
	// print_all_edge_info(copy);
	return true;}	
	
}


// std::cout<<"Out 1p red"<<std::endl;
	
return false;	
	
}


// working here 
// checks if spin of first and last bubbles are different: true, or the same :false.
bool AmiGraph::is_hubW_graph(graph_t g){

std::cout<<"Getting start and end"<<std::endl;

vertex_t start_v, end_v;
find_start_end_vertices(g, start_v, end_v);

spin_type start_spin, end_spin;

boost::graph_traits<graph_t>::edge_iterator ei, ei_end;

for(boost::tie(ei,ei_end)=edges(g); ei!=ei_end; ++ei){
	
if(source(*ei,g)==start_v || target(*ei,g)==start_v){
	start_spin=g[*ei].spin;
}
if(target(*ei,g)==end_v || source(*ei,g)== end_v){
	end_spin=g[*ei].spin;
}
	
	
}


if(start_spin == end_spin){return true;}


return false;

	
}

bool AmiGraph::is_ppvertex_graph(graph_t g){

std::cout<<"Getting start and end"<<std::endl;

vertex_t start_v, end_v;

int start_index=find_start_index(g, start_v);
int end_index=find_end_index(g, end_v);	

int extra_components=0;

boost::graph_traits<graph_t>::edge_iterator ei, ei_end;

for(boost::tie(ei,ei_end)=edges(g); ei!=ei_end; ++ei){
	
if(source(*ei,g)==start_v && target(*ei,g)==end_v){
	extra_components++;
}
	
	
}


std::cout<<"Symmetrize"<<std::endl;
symmetrize_internal_bosonic_lines(g);

vertex_vector_t extern_vect_list;
edge_vector_t extern_edge_list;
find_external_vertices(g, extern_vect_list, extern_edge_list);

std::cout<<"Deleting legs"<<std::endl;

delete_legs(g,extern_vect_list, extern_edge_list);

clear_vertex(start_v,g);
clear_vertex(end_v,g);
remove_vertex(start_v,g);
remove_vertex(end_v,g);

std::cout<<"Checking connectedness"<<std::endl;

// bool result=is_connected_graph(g);

int components=connected_components(g)+extra_components;

std::cout<<"Got result "<< components<<std::endl;

if(components==1){return true;}

return false;

/* 	
vertex_vector_list_t pv;
edge_vector_list_t pe;

get_principle_pp_lines(	g, pv, pe);
	
for(int i=0; i<pv[0].size(); i++){
std::cout<<pv[0][i]<<" ";
}	
std::cout<<std::endl;

for(int i=0; i<pv[1].size(); i++){
std::cout<<pv[1][i]<<" ";
}	
std::cout<<std::endl;

 */
	
}

void AmiGraph::get_principle_pp_lines(graph_t &g, vertex_vector_list_t &vvl, edge_vector_list_t&evl){

vvl.clear();
evl.clear();
vvl.resize(2);
evl.resize(2);

vertex_t start_v, end_v;

int start_index=find_start_index(g, start_v);
int end_index=find_end_index(g, end_v);	

vertex_t next, current;
vertex_t next2, current2;


current=start_v;
current2=start_v;

bool continue1=true;
bool continue2=true;

// get first 
vertex_vector_t nextv;
vertex_vector_t currentv;
nextv.resize(2);
currentv.resize(2);

currentv[0]=start_v;
currentv[1]=start_v;

std::vector<int> cont(2,0);
std::vector<int> nextcont(2,0);

int count=0;

boost::graph_traits<graph_t>::edge_iterator ei, ei_end;

for(boost::tie(ei,ei_end)=edges(g); ei!=ei_end; ++ei){

if( g[*ei].g_struct_.stat_==AmiBase::Fermi){


if(source(*ei,g)==currentv[count] && target(*ei,g)!=end_v){
nextv[count]=target(*ei,g);
cont[count]=1;
count++;	
}

	
}
	
}

for(int i=0; i< nextv.size(); i++){
if(cont[count]==1){
vvl[i].push_back(nextv[i]);	
}
}

// reset cont 
cont[0]=0;
cont[1]=0;




do{
	
	std::cout<<"Cont values "<< cont[0]<<" "<< cont[1]<<std::endl;
	
	currentv=nextv;
	
	nextcont.clear();
	nextcont.resize(2,0);

for(int i=0; i< currentv.size(); i++){

if(cont[i]==0){
bool found=get_next(g, currentv[i], nextv[i], end_v);
std::cout<<"Found was "<<found<<std::endl;
if(found){
vvl[i].push_back(nextv[i]);

}else{
	nextcont[i]=1;
}

 
}
	
	
}

for (int i=0 ; i< nextcont.size(); i++){
if(nextcont[i]==1){
cont[i]=nextcont[i];
}

}
	

}while(cont[0]!=1 || cont[1] !=1);
	
	
}


bool AmiGraph::get_next(graph_t &g, vertex_t &current, vertex_t &next, vertex_t &end_v){

std::cout<<"Looking for next starting at index "<<g[current].index_<<std::endl;

boost::graph_traits<graph_t>::edge_iterator ei, ei_end;

for(boost::tie(ei,ei_end)=edges(g); ei!=ei_end; ++ei){

if( g[*ei].g_struct_.stat_==AmiBase::Fermi){


if(source(*ei,g)==current && target(*ei,g)!=end_v){
next=target(*ei,g);	
return true;
}

	
}
	
}


return false;	
} 


// TODO: this does not appear to be correctly removing RPA diagrams?
bool AmiGraph::is_1P_bose_reducible_graph(graph_t g){
	
// print_all_edge_info(g);

edge_vector_t in, out;
	
boost::graph_traits<graph_t>::edge_iterator ei, ei_end;

AmiBase::stat_type external;

if( (graph_type == AmiBase::Pi_phuu) || (graph_type == AmiBase::Pi_phud) || (graph_type== AmiBase::doubleocc) || ami_parameters.TYPE_==AmiBase::Pi_ppud || ami_parameters.TYPE_==AmiBase::Pi_ppuu|| ami_parameters.TYPE_==AmiBase::FORCE){
	
external=AmiBase::Bose;	
}
else{
	external=AmiBase::Fermi;
}


// REMOVE external lines if they exist
// this only allows for a single external line

if(ami_parameters.TYPE_!=AmiBase::density &&   graph_type!=AmiBase::Greens&& ami_parameters.TYPE_!=AmiBase::doubleocc && ami_parameters.TYPE_!=AmiBase::DOS){

for (boost::tie(ei,ei_end) = edges(g); ei != ei_end; ++ei){

//g[*ei].g_struct_.stat_==AmiBase::Bose){


if(g[*ei].g_struct_.stat_==external){


if ( degree( source(*ei,g),g)==1 && (source(*ei,g) != target(*ei,g) )){ 
in.push_back(*ei);

}

if( degree( target(*ei,g),g)==1 && (source(*ei,g) != target(*ei,g))){
	out.push_back(*ei);
}

}

}

for(int i=0; i< in.size(); i++){
remove_edge(in[i],g);
remove_vertex(source(in[i],g),g);
}

for(int i=0; i<out.size(); i++){
remove_edge(out[i],g);
remove_vertex(target(out[i],g),g);
}

}

// std::cout<<"In theory have removed external lines"<<std::endl;

	// number_vertices(g);
	// print_all_edge_info(g);

	
edge_vector_t internal_fermi, internal_copy, internal_bose;

find_fermionic_edges(g, internal_fermi);
// find_bosonic_edges(g,internal_bose);

graph_t copy;

find_non_tadpole_bosonic_edges(g,internal_bose);

for (int i=0; i<internal_bose.size(); i++){
	// std::cout<<i<<std::endl;
copy=g;
// find_bosonic_edges(copy, internal_copy);
find_non_tadpole_bosonic_edges(copy, internal_copy);

// std::cout<<"Test removing edge "<<print_edge(internal_copy[i],g)<<std::endl;	
remove_edge(internal_copy[i],copy);
// std::cout<<"Removed"<<std::endl;


// std::cout<<std::endl;
// print_all_edge_info(copy);
// std::cout<<std::endl;


if(!is_connected_graph(copy)){
	// std::cout<<std::endl;
	// number_vertices(copy);
	// print_all_edge_info(copy);
	return true;}	
	
}

	
return false;	
	
}

bool AmiGraph::disconnected_loop(graph_t g){

bool output=false;
	
vertex_vector_t vv;
edge_vector_t ev;
	
get_principle_loop(g, vv, ev);	

boost::graph_traits<graph_t>::edge_iterator ei, ei_end;

vertex_t v1,v2;

	for (boost::tie(ei,ei_end) = edges(g); ei != ei_end; ++ei){

		if(g[*ei].g_struct_.stat_==AmiBase::Bose){
			
		if(degree(source(*ei,g),g)==1){ v1=target(*ei,g);} 
		if(degree(target(*ei,g),g)==1){ v2=source(*ei,g);} 	
			
			
		}

	}


bool foundv1=false;
bool foundv2=false;

for(int i=0; i< vv.size(); i++){

	if(vv[i]==v1){foundv1=true;}
	if(vv[i]==v2){foundv2=true;}
	
}

output=!(foundv1&&foundv2);

	
	
return output;
	
}


bool AmiGraph::has_fock_insertion(graph_t g){
	
// print_all_edge_info(g);

edge_t in, out;
	
boost::graph_traits<graph_t>::edge_iterator ei, ei_end;

AmiBase::stat_type external;

if( (graph_type == AmiBase::Pi_phuu) || (graph_type == AmiBase::Pi_phud)|| (graph_type== AmiBase::doubleocc) || ami_parameters.TYPE_==AmiBase::Pi_ppud || ami_parameters.TYPE_==AmiBase::Pi_ppuu){
	
external=AmiBase::Bose;	
}
else{
	external=AmiBase::Fermi;
}


for (boost::tie(ei,ei_end) = edges(g); ei != ei_end; ++ei){

//g[*ei].g_struct_.stat_==AmiBase::Bose){


if(g[*ei].g_struct_.stat_==external){


if ( degree( source(*ei,g),g)==1 && (source(*ei,g) != target(*ei,g) )){ 
in=*ei;

}

if( degree( target(*ei,g),g)==1 && (source(*ei,g) != target(*ei,g))){
	out=*ei;
}

}

}

remove_edge(in,g);
remove_vertex(source(in,g),g);
remove_edge(out,g);
remove_vertex(target(out,g),g);

	// number_vertices(g);
	// print_all_edge_info(g);

	
	
// boost::graph_traits<graph_t>::edge_iterator ei, ei_end;

boost::graph_traits<graph_t>::edge_iterator ei2, ei2_end;

for (boost::tie(ei,ei_end) = edges(g); ei != ei_end; ++ei){
	
if( g[*ei].g_struct_.stat_==AmiBase::Bose){

vertex_t sour, tar;

sour=source(*ei,g);
tar=target(*ei,g);

// for each bosonic edge. check if there is any fermionic edge that goes from the target to source or source to target 
for (boost::tie(ei2,ei2_end) = edges(g); ei2 != ei2_end; ++ei2){

if( g[*ei2].g_struct_.stat_==AmiBase::Fermi){
	
	
	// return true if find such an edge 
if(source(*ei2,g)==sour && target(*ei2,g)==tar){

return true;
}	

if(source(*ei2,g)==tar && target(*ei2,g)==sour){
	return true;
}
	
	
}



}

}	
	
	
}	

	
return false;	
	
}



// TODO: this will falsely pick up hartree tadpole diagrams ? maybe not, because it would delete the whole tadpole 
bool AmiGraph::is_ladder_graph(graph_t &g){

	
edge_vector_t internal_fermi, internal_copy, internal_bose;
edge_vector_t external;

find_principle_bose_lines(g, external);

// std::cout<<"Found list of internal bose lines of size "<< external.size()<<std::endl;

bool result=false;

graph_t copy;
for(int i=0; i< external.size();i++){

copy=g;
find_principle_bose_lines(copy, internal_bose);
split_graph(internal_bose[i],copy);

if(is_connected_graph(copy)==false){
result=true;	
}
}

return result;
	
}

void AmiGraph::find_principle_bose_lines(graph_t &g, edge_vector_t &ev){
	
ev.clear();	

// find a vertex on the principle line 

vertex_t start;

// reset visited 
reset_visited(g);

vertex_vector_t prin_vv;
edge_vector_t prin_ev;

get_principle_loop(g,prin_vv,prin_ev);


boost::graph_traits<graph_t>::edge_iterator ei, ei_end;

for (boost::tie(ei,ei_end)=edges(g); ei!=ei_end; ++ei){
	
if( g[*ei].g_struct_.stat_==AmiBase::Bose){

vertex_t sour, tar;
edge_t edge;
edge=*ei;

sour=source(*ei,g);
tar=target(*ei,g);

bool add=false;
	
for(int i=0; i< prin_vv.size(); i++){

if(	sour==prin_vv[i] || tar==prin_vv[i]){
	add=true;
	break;
}

}

if(add){
ev.push_back(edge);	
}



	
}

	
}

	
	
	
}

// TODO this is coded only for particle hole loops 
void AmiGraph::get_principle_loop(graph_t &g, vertex_vector_t &vv, edge_vector_t &ev){
vv.clear();
ev.clear();

vertex_t current, next;

if(graph_type == AmiBase::Pi_phuu || graph_type == AmiBase::Pi_phud || (graph_type== AmiBase::doubleocc) || ami_parameters.TYPE_==AmiBase::Pi_ppud || ami_parameters.TYPE_==AmiBase::Pi_ppuu){
	
	
	// find external in vertex vertex 
boost::graph_traits<graph_t>::edge_iterator ei, ei_end;

for (boost::tie(ei,ei_end)=edges(g); ei!=ei_end; ++ei){
	
if(degree(source(*ei,g),g)==1){current= target(*ei,g); break;}
	
}

// vv.push_back(current);

bool instop=false;
do{
 // std::cout<<"Source is "<<graph[current].index_<<std::endl;	

g[current].visited=1;
//g[current].index_=counter;
//counter++;
vv.push_back(current);
get_adjacent_vertex_stat(g, current, next, AmiBase::Fermi); 
if(edge(current,next,g).second){ev.push_back(edge(current,next,g).first);}
// if(g[next].visited==0){ vv.push_back(next);}
if( g[next].visited==1){instop=true; g[next].visited=1;}
current=next;

} while(!instop);	

	
}

}

void AmiGraph::split_graph(edge_t &e, graph_t &g){

vertex_t v1,v2;

v1=source(e,g);
v2=target(e,g);

clear_vertex(v1,g);
clear_vertex(v2,g);

remove_vertex(v1,g);
remove_vertex(v2,g);


return;	
	
}

void AmiGraph::find_entry_vertex(graph_t &g, vertex_t &v){
	
boost::graph_traits<graph_t>::vertex_iterator vi, vi_end;

for(boost::tie(vi,vi_end)=vertices(g); vi!=vi_end; ++vi){

if(degree(*vi,g)==1 && in_degree(*vi,g)==0){
v=*vi;
}	
	
}

	
	
	
}


bool AmiGraph::is_connected_graph(graph_t &g){
	
reset_visited(g);

number_vertices(g);
vertex_vector_t vert_list;

collect_unvisited_vertices(g, vert_list);


// std::cout<<"Start of analysis graph is "<<std::endl;
// print_all_edge_info(g);


//current=vert_list[0]; // again only works if the zero is still the entry external leg 

vertex_vector_t current_list;
vertex_vector_t next_list;

vertex_t in=random_vertex(g,rand_gen);
// find_entry_vertex(g,in);
// std::cout<<"Entry vertex has index "<< g[in].index_<<std::endl;



current_list.push_back(in); // this is our starter list

// std::cout<<"Starting with vertex numbered "<< g[vert_list[0]].index_<<std::endl;
// find source of edge with degree equal to one 

// boost::graph_traits<graph_t>::edge_iterator ei, ei_end;

// for (boost::tie(ei,ei_end)=edges(g); ei!=ei_end; ++ei){
	
// if(degree(source(*ei,g),g)==1){current_list.push_back(source(*ei,g)); break;}
	
	
// }


// need a function that takes in a list of vertices - returns a list of adjacent vertices, and sets every current_list as visited 
bool exit_do=false;
do{
get_adjacent_lists_undirected(g, current_list, next_list);
// for(int i=0; i<next_list.size();i++){
	// std::cout<<"Next are vertices "<<i<<" "<<g[next_list[i]].index_<<std::endl;
	
// }
label_visited(g, current_list);

// set current_list to next _list if next_list.size()!=0
// std::cout<<"Next list size is "<<next_list.size()<<std::endl;
if(next_list.size()==0){
exit_do=true;	
}else{
current_list=next_list;
}	
	
}while(!exit_do);
	
collect_unvisited_vertices(g, vert_list);
// std::cout<<"Final vert_list is size "<<vert_list.size()<<std::endl;
// std::cout<<"Exited do"<<std::endl;

// for(int i=0; i< vert_list.size(); i++){
	// std::cout<<vert_list[i]<<std::endl;
	// std::cout<<g[vert_list[i]].index_<<std::endl;
// }
// std::cout<<"End of analysis graph is "<<std::endl;
// print_all_edge_info(g);

if(vert_list.size()==0){return true;}else{return false;}	
	
	
}	


void AmiGraph::get_adjacent_lists_undirected(graph_t &g, vertex_vector_t &current_list, vertex_vector_t &next_list){
	
next_list.clear();	

// std::cout<<"Current list has size "<< current_list.size()<<std::endl;
	
for (int i=0 ; i< current_list.size();i++){


boost::graph_traits<graph_t>::in_edge_iterator iei, iei_end;

for(boost::tie(iei,iei_end)=in_edges(current_list[i],g); iei!=iei_end; ++iei){
	edge_t edge=*iei;
 // std::cout<<"In edge is "<< print_edge(edge,g)<<std::endl;	
if(g[source(*iei,g)].visited==0){
next_list.push_back(source(*iei,g));
}	
	
}

boost::graph_traits<graph_t>::out_edge_iterator oei, oei_end;

for(boost::tie(oei,oei_end)=out_edges(current_list[i],g); oei!=oei_end; ++oei){
	edge_t edge=*oei;
	// std::cout<<"Out edge is "<< print_edge(edge,g)<<" with visited ("<<g[source(edge,g)].visited<<","<<g[target(edge,g)].visited<<")"<<std::endl;
	// std::cout<<"Out edge is "<< print_edge(edge,g)<<std::endl;
if(g[target(edge,g)].visited==0){
	// std::cout<<"Pushing into next list"<<std::endl;
next_list.push_back(target(edge,g));
}	
	
}


// boost::graph_traits<graph_t>::adjacency_iterator ai, ai_end;

// for(boost::tie(ai,ai_end)=adjacent_vertices(current_list[i],g); ai!=ai_end; ++ai){

// if (g[*ai].visited==0){
// next_list.push_back(*ai);
	
// }	
	
	
// }


}	

	
}

void AmiGraph::get_adjacent_lists(graph_t &g, vertex_vector_t &current_list, vertex_vector_t &next_list){
	
next_list.clear();	
	
for (int i=0 ; i< current_list.size();i++){


boost::graph_traits<graph_t>::adjacency_iterator ai, ai_end;

for(boost::tie(ai,ai_end)=adjacent_vertices(current_list[i],g); ai!=ai_end; ++ai){

if (g[*ai].visited==0){
next_list.push_back(*ai);
	
}	
	
	
}


}	

	
}

void AmiGraph::label_visited(graph_t &g, vertex_vector_t &current_list){
	
for (int i=0 ; i< current_list.size();i++){

g[current_list[i]].visited=1;	
}	
	
}


bool AmiGraph::is_oneleg(graph_t g){
  
 
vertex_vector_t extern_vect_list;
edge_vector_t extern_edge_list; 
  
find_external_vertices(g, extern_vect_list, extern_edge_list);
delete_legs(g,extern_vect_list, extern_edge_list);


boost::graph_traits<graph_t>::edge_iterator ei, ei_end;

edge_t thisedge;
for (boost::tie(ei,ei_end) = edges(g); ei != ei_end; ++ei){

thisedge=*ei;
if(is_labelled_bose_zero(thisedge,g)){return true;}


}


return false;
  
}





// is diagram one particle reducible - are any labelles equal to only the external label?

// find external edges - get its label
// cycle through all edges and pick out any that match labels
// TODO TEST THIS FUNCTION 

// TODO: this function will not think an unlabelled graph is reducible 
bool AmiGraph::is_one_particle_reducible(graph_t &g){

// find the external legs
vertex_vector_t extern_vect_list;
edge_vector_t extern_edge_list;
edge_t one;
find_external_vertices(g, extern_vect_list, extern_edge_list);	

boost::graph_traits<graph_t>::edge_iterator ei, ei_end;

bool is_ext=false;

for(boost::tie(ei,ei_end)=edges(g); ei!=ei_end; ++ei){

// check that it isn't an external leg already
is_ext=false;
for(int i=0; i<extern_edge_list.size();i++){
if(*ei==extern_edge_list[i]){is_ext=true; break;}	
	
}

// if edge has external label then it is reducible, so return true 
if(!is_ext){
one=*ei;
bool check=edge_labels_are_equal(one, extern_edge_list[0],g);
//std::cout<<"ON edge "<< print_edge(one,g)<<" with bool "<<check<<std::endl;
if(check){
return true;	
}
// if(edge_labels_are_equal(*ei, extern_edge_list[0],g)){
// return true;
// }	
	
}


}	
	
return false;	
	
}



// Filter - check if is skeleton diagram from labelled graph 
// for each fermionic edge - check if it's labels are the same as each other edge - except the external legs
// - get list of all fermionic edges that are not the external
// - from i to end j=i to end check if equal.  If any two are equal immediately break and report that it is NOT a skeleton diagram 

bool AmiGraph::is_skeleton(graph_t &g){
	
	

edge_vector_t internal_fermi;

find_internal_fermionic_edges(g,internal_fermi);	

for(int i=0; i< internal_fermi.size(); i++){
for(int j=i; j< internal_fermi.size(); j++){


if( i!=j){

bool check=edge_labels_are_equal(internal_fermi[i], internal_fermi[j],g);
if(check){
	return false;
}
	
	
}




}	
	
	
}

	

return true;	
	
	
}



bool AmiGraph::is_nested(graph_t &g){
	
graph_t myg=g;

// print_all_edge_info(myg);

edge_vector_t internal_fermi;
get_principle_line(myg,internal_fermi);	

// std::cout<<"Principle line is size "<< internal_fermi.size()<<std::endl;
if(internal_fermi.size()<4){return false;} // catch for short diagrams that can't possibly be nested
// THIS Won't work in general. 

internal_fermi.erase(internal_fermi.begin());
internal_fermi.pop_back();


for(int i=0; i< internal_fermi.size(); i++){
for(int j=i; j< internal_fermi.size(); j++){


if( i!=j){

bool check=edge_labels_are_equal(internal_fermi[i], internal_fermi[j],myg);
// std::cout<<"Edge check "<< i<<" "<< j<<" "<<check<<std::endl;
// print_edge_info(internal_fermi[i],myg);
// print_edge_info(internal_fermi[j],myg);
if(check){
	return true;
}
	
	
}




}	
	
	
}

	

return false;	
	
	
}




/// TODO: this does not work for mixed spin graphs. Needs updating.
bool AmiGraph::is_hubbard(graph_t &g){
	
	// std::cout<<"In is hubbard"<<std::endl;
	
edge_vector_t bosonic_lines;

find_bosonic_edges(g,bosonic_lines);

bool spin_issue_1=false;
bool spin_issue_2=false;

for(int i=0; i<bosonic_lines.size(); i++){
	
	vertex_t sour=source(bosonic_lines[i],g);
	vertex_t targ=target(bosonic_lines[i],g);

if(!(degree(sour,g)==1 || degree(targ,g)==1)){
	

if(spin_mismatch(sour,g) || spin_mismatch(targ,g)){
	return false;
}
	
	
if( return_spin_on_vertex(g, sour, spin_issue_1 ) == return_spin_on_vertex(g,targ, spin_issue_2 )){
	return false;
}
}
// this returns false if the bubble was attached in such a way that spin was not maintained
// if(spin_issue_1 || spin_issue_2){
	// return false;
// }


}	
	
// std::cout<<"Exit is hubbard"<<std::endl;	
	
return true;	
	
}

bool AmiGraph::spin_mismatch(vertex_t &v, graph_t &g){

bool result=false;
	
boost::graph_traits<graph_t>::in_edge_iterator iei, iei_end;
boost::graph_traits<graph_t>::out_edge_iterator oei, oei_end;	

std::vector<spin_type> spin_vector;
	
for(boost::tie(iei,iei_end)=in_edges(v,g); iei!=iei_end; ++iei){
	
if(g[*iei].g_struct_.stat_==AmiBase::Fermi){	
spin_vector.push_back(g[*iei].spin);	
}
}

for(boost::tie(oei,oei_end)=out_edges(v,g); oei!=oei_end; ++oei){
	
if(g[*oei].g_struct_.stat_==AmiBase::Fermi){	
spin_vector.push_back(g[*oei].spin);	
}
}

// TODO: Should really test this function. If not all the spins on this line are equal then throw an error
//for

// FOR NOW THIS IS HARDCODED
if( spin_vector[0]!= spin_vector[1]){ 

result=true;
// number_vertices(g);	
// print_all_edge_info(g);
// throw std::runtime_error("Something has gone wrong, two fermi lines have different spins"); }
}

return result;
	
}

AmiGraph::spin_type AmiGraph::return_spin_on_vertex(graph_t &g, vertex_t &v, bool &fail){
	
boost::graph_traits<graph_t>::in_edge_iterator iei, iei_end;
boost::graph_traits<graph_t>::out_edge_iterator oei, oei_end;	

std::vector<spin_type> spin_vector;
	
for(boost::tie(iei,iei_end)=in_edges(v,g); iei!=iei_end; ++iei){
	
if(g[*iei].g_struct_.stat_==AmiBase::Fermi){	
spin_vector.push_back(g[*iei].spin);	
}
}

for(boost::tie(oei,oei_end)=out_edges(v,g); oei!=oei_end; ++oei){
	
if(g[*oei].g_struct_.stat_==AmiBase::Fermi){	
spin_vector.push_back(g[*oei].spin);	
}
}

// TODO: Should really test this function. If not all the spins on this line are equal then throw an error
//for

// FOR NOW THIS IS HARDCODED
if( spin_vector[0]!= spin_vector[1]){ 

fail=true;
// number_vertices(g);	
// print_all_edge_info(g);
// throw std::runtime_error("Something has gone wrong, two fermi lines have different spins"); }
}


return spin_vector[0];	
	
	
}









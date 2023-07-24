#include "amigraph.hpp"


void AmiGraph::generate_sigma_ct( graph_t &g_in, std::vector< graph_t> &ct_vec, int maxdots){
// max is the maximum number of dots to add to one graph 

graph_t gc=g_in;

edge_vector_t fermi_edges;
// find_fermionic_edges(g_in, fermi_edges);

if(graph_type==AmiBase::density|| graph_type==AmiBase::DOS || graph_type==AmiBase::Greens){
  find_fermionic_edges(g_in, fermi_edges);
}else{
find_internal_fermionic_edges(g_in, fermi_edges);
}
// we now create all combinations of 
int n=fermi_edges.size();
for( int i=1; i<=maxdots; i++){

int r=i;
std::vector< std::vector<int>> list;

combinations_repetition(n,r,list);

// std::cout<<"Comb list has size "<<list.size()<<std::endl;
// for(int i=0; i< list.size(); i++){
	
	// for(int j=0; j<list[i].size(); j++){
	// std::cout<<list[i][j];	
	// }
	// std::cout<<std::endl;
	
	
// }
// now for each list item we will make a new graph that is then converted to a ct graph 

for (int comb=0; comb<list.size(); comb++){
	
// convert comb list to an index and length pair 
std::vector<std::pair<int,int>> chainlist;// vector of size n that contains how MANY dots to add to each current line for this particular combination 
	// std::cout<<"On Comb "<< comb<<std::endl;
	// for(int j=0; j<list[comb].size(); j++){
	// std::cout<<list[comb][j];	
	// }
	// std::cout<<std::endl;
	
// making chain lists 
			for(int i=0; i<n;i++){
				
			std::pair<int,int> temp(i,0);
				
				for (int j=0; j< list[comb].size(); j++){
				
			if(list[comb][j]==i){
			temp.second++;
			}	
					
				}

			if(temp.second>0){
			chainlist.push_back(temp);
			}	
				
			}


			// for(int i=0; i<chainlist.size(); i++){

			// std::cout<<"i="<<i<<" "<<chainlist[i].first<<" "<< chainlist[i].second<<std::endl;	
				
			// }
//////

// next for each chain 
graph_t gtemp=gc;
for(int i=0; i<chainlist.size(); i++){


// std::cout<<"i="<<i<<" "<<chainlist[i].first<<" "<< chainlist[i].second<<std::endl;	
insert_chain(gc,gtemp, fermi_edges[chainlist[i].first],chainlist[i].second); 

// print_all_edge_info(gc);


}	

ct_vec.push_back(gtemp);
	
	
}





}

	
}


void AmiGraph::swap_alphas(AmiGraph::graph_t &g, int first, int second){
	
boost::graph_traits<AmiGraph::graph_t>::edge_iterator ei, ei_end;

for (boost::tie(ei,ei_end)=edges(g); ei!=ei_end; ++ei){
	
int temp=g[*ei].g_struct_.alpha_[second];
g[*ei].g_struct_.alpha_[second]=g[*ei].g_struct_.alpha_[first];	
g[*ei].g_struct_.alpha_[first]=temp;
	
	
}	
	
	
	
}


void AmiGraph::zero_external(AmiGraph::graph_t &g, int index){
	
boost::graph_traits<AmiGraph::graph_t>::edge_iterator ei, ei_end;

for (boost::tie(ei,ei_end)=edges(g); ei!=ei_end; ++ei){
	

g[*ei].g_struct_.alpha_[index]=0;	
	
	
}
	
fix_epsilons(g);	
	
}

// TODO this is for hubbard only since bosonic lines won't get their Q's modified 
void AmiGraph::zero_external_Q(AmiGraph::graph_t &g){
	
// walk through and find all lines where the alphas are identical when ignoring the external, index (assume it is the last one)	
// if they are equal then set their epsilons to all be the same 	
	
	
boost::graph_traits<AmiGraph::graph_t>::edge_iterator ei, ei_end;
boost::graph_traits<AmiGraph::graph_t>::edge_iterator ei2, ei_end2;

for (boost::tie(ei,ei_end)=edges(g); ei!=ei_end; ++ei){
	for (boost::tie(ei2,ei_end2)=edges(g); ei2!=ei_end2; ++ei2){
	
	if((g[*ei].g_struct_.stat_!= AmiBase::Fermi) || (g[*ei2].g_struct_.stat_!= AmiBase::Fermi)){
		continue;
	}
	
	bool same=true;
		
		for(int i=0; i<g[*ei].g_struct_.alpha_.size()-1; i++){

			if(g[*ei].g_struct_.alpha_[i]!=g[*ei2].g_struct_.alpha_[i]){
				same=false;
				break;
			}	
		}
	
	if(same){
	
	g[*ei2].g_struct_.eps_=g[*ei].g_struct_.eps_;
		
	}
	
	
	}
}
	
// fix_epsilons(g);	
	
}




// TODO this is for hubbard only since bosonic lines won't get their Q's modified 
void AmiGraph::PIPI_external_Q(AmiGraph::graph_t &g){
	
// walk through and find all lines where the alphas are identical when ignoring the external, index (assume it is the last one)	
// if they are equal then set the epsilon of the one with non-zero alpha to be the negative of the other one.
 // std::cout<<"Checking pipi "<<std::endl;	
	
boost::graph_traits<AmiGraph::graph_t>::edge_iterator ei, ei_end;
boost::graph_traits<AmiGraph::graph_t>::edge_iterator ei2, ei_end2;

for (boost::tie(ei,ei_end)=edges(g); ei!=ei_end; ++ei){
	for (boost::tie(ei2,ei_end2)=edges(g); ei2!=ei_end2; ++ei2){
	
	if((g[*ei].g_struct_.stat_!= AmiBase::Fermi) || (g[*ei2].g_struct_.stat_!= AmiBase::Fermi)){
		continue;
	}
	
	bool same=true;
  bool negative=true;
	
  // if alpha is both non-zero and also 
	if(g[*ei].g_struct_.alpha_.back()==g[*ei2].g_struct_.alpha_.back() || g[*ei].g_struct_.alpha_.back()==-g[*ei2].g_struct_.alpha_.back()){
		continue;
	}
		
		for(int i=0; i<g[*ei].g_struct_.alpha_.size()-1; i++){

			if(g[*ei].g_struct_.alpha_[i]!=g[*ei2].g_struct_.alpha_[i]){
				same=false;
				break;
			}	
		}
    
    for(int i=0; i<g[*ei].g_struct_.alpha_.size()-1; i++){

			if(g[*ei].g_struct_.alpha_[i]!=-g[*ei2].g_struct_.alpha_[i]){
				negative=false;
				break;
			}	
		}
    
	
	if(same || negative){
    
    std::cout<<"Found a conflict"<<std::endl;
	
	if(g[*ei].g_struct_.alpha_.back()!=0){
		for(int m=0; m<g[*ei].g_struct_.eps_.size(); m++){
		g[*ei].g_struct_.eps_[m]=-g[*ei2].g_struct_.eps_[m];
    
		}
    g[*ei].g_struct_.alpha_.back()=0;
	}
	
	if(g[*ei2].g_struct_.alpha_.back()!=0){
		for(int m=0; m<g[*ei].g_struct_.eps_.size(); m++){
		g[*ei2].g_struct_.eps_[m]=-g[*ei].g_struct_.eps_[m];
    
		}
    g[*ei2].g_struct_.alpha_.back()=0;
	}
	
		
	}
	
	
	}
} 
	
// now go through and if any alpha is not zero, flip its epsilon and zero the alpha 
//set the alpha external indexes to 0. this will let us do pi,pi AND Omega=0	
// boost::graph_traits<AmiGraph::graph_t>::edge_iterator ei, ei_end;

for (boost::tie(ei,ei_end)=edges(g); ei!=ei_end; ++ei){
	
if(g[*ei].g_struct_.alpha_.back()!=0){
  for(int m=0; m<g[*ei].g_struct_.eps_.size(); m++){
		g[*ei].g_struct_.eps_[m]=-g[*ei].g_struct_.eps_[m];
    
		}
  
g[*ei].g_struct_.alpha_.back()=0;	

}	
	
}
		
	
}



void AmiGraph::save_ct_eff_orders(gg_matrix_t &ggm){

std::ofstream file;
file.open("eff_orders.dat");
file<<"# order num ext_index eff_ord ct_count sigma_ct_count "<<std::endl;	
	

for(int ord=0; ord< ggm.size(); ord++){

	for(int group=0; group<ggm[ord].size(); group++){
		// if(group%5==0){std::cout<<"On group "<< group<<std::endl;}
	for(int graph=0; graph< ggm[ord][group].graph_vec.size(); graph++){	
	
	int order=ord;//graph_order(ggm[ord][group].graph_vec[graph]);
	int bubble_ct=ggm[ord][group].graph_vec[graph][boost::graph_bundle].ct_count;
	int sigma_ct=ggm[ord][group].graph_vec[graph][boost::graph_bundle].sigma_ct_count;
	
	int sum=order+bubble_ct+sigma_ct+1;
	
	file<<ord<<" "<< group<<" "<<graph<<" "<< order+1<<" "<<bubble_ct<<" "<<sigma_ct<<" "<<sum<<std::endl;
	
	
	
	}
	}
}
	
	
	
}

// on first call presumes g_in==ctg and then modifies only the ctg 
void AmiGraph::insert_chain(graph_t &g_in, graph_t &ctg,  edge_t &e, int length){




// std::cout<<"adding chain to edge from "<<g_in[source(e,g_in)].index_<<" - to - "<<g_in[target(e,g_in)].index_<<std::endl;

// first find the edge and vertices in gc with the correct index values

vertex_t vA,vB;
edge_t cedge;

boost::graph_traits<graph_t>::edge_iterator ei, ei_end;

for(boost::tie(ei,ei_end)=edges(ctg); ei!=ei_end; ++ei){

if(	(ctg[source(*ei,ctg)].index_== g_in[source(e,g_in)].index_) && (ctg[target(*ei,ctg)].index_== g_in[target(e,g_in)].index_) && g_in[*ei].g_struct_.stat_==AmiBase::Fermi){
vA=	source(*ei,ctg);
vB= target(*ei,ctg);
	cedge=*ei;
}
		
}


// add vertices and give them new index_ values 
int num_vert=num_vertices(ctg);
vertex_vector_t vv(length);
for(int i=0; i< length; i++){
	vv[i]=add_vertex(ctg);
	ctg[vv[i]].index_=num_vert+i;
}


// connect inner vertices 
// std::cout<<"Connecting inner"<<std::endl;
for(int i=1; i< length; i++){

add_edge(vv[i-1],vv[i],edge_info(g_in[e].g_struct_.eps_, g_in[e].g_struct_.alpha_,  g_in[e].g_struct_.stat_,g_in[e].fermi_loop_id, g_in[e].spin), ctg);

}	

// connect chain to source and target 
// std::cout<<"Connecting outer"<<std::endl;

add_edge(vA,vv[0], edge_info(g_in[e].g_struct_.eps_, g_in[e].g_struct_.alpha_,  g_in[e].g_struct_.stat_,g_in[e].fermi_loop_id, g_in[e].spin)  ,ctg);

add_edge(vv[vv.size()-1],vB, edge_info(g_in[e].g_struct_.eps_, g_in[e].g_struct_.alpha_,  g_in[e].g_struct_.stat_,g_in[e].fermi_loop_id, g_in[e].spin)  ,ctg);

// remove original edge 
// std::cout<<"Removing original"<<std::endl;

remove_edge(cedge,ctg);
// std::cout<<"Edge removed?"<<std::endl;

ctg[boost::graph_bundle].sigma_ct_count+=length;

// std::cout<<"Exiting insert chain function"<<std::endl;
}	



void AmiGraph::generate_bubble_ct( graph_t &g_in, std::vector< graph_t> &ct_vec){
	
ct_vec.clear();	
	
if( graph_order(g_in)==0){
std::cerr<<"Can't generate counter-term for zeroth order graph"<<std::endl;
return;	
	
}
// first check that the graph is labelled - if not, return 
// either - size=0 for alpha or ALL alphas are zero
if(not_labelled(g_in)){

std::cerr<<"Can't generate counter-term from unlabelled graph. Label it first."<<std::endl;
return;	
	
}

// count bare bubbles
vertex_vector_list_t bv_list; // bubble vertex list 
edge_vector_list_t be_list; // bubble edge list
edge_vector_list_t leg_e_list;// bosonic legs of the bubble 

// problem - the vertex and edge descriptors become invalid if we change anything 
bubble_finder(g_in, bv_list, be_list, leg_e_list);

// std::cout<<"Found bubbles n="<< be_list.size()<<std::endl;

//// create vector of bubble structs 
std::vector< bubble > bubble_vector;

for( int i=0; i< bv_list.size(); i++){
bubble bub;
descriptor_to_index(g_in, bv_list[i], be_list[i], leg_e_list[i], bub);

bubble_vector.push_back(bub);	
	
}

if(bubble_vector.size()!= be_list.size()){
	throw std::runtime_error("Didn't convert bubbles correctly");
}
/////

// we now create all combinations of 
int max=bubble_vector.size();
for( int i=0; i<bubble_vector.size(); i++){

int length=i+1;

std::vector< std::vector<int>> list;

combinations(max,length,list);
//for each combination m - remove each bubble_vector from the graph 


for(int m=0; m< list.size(); m++){

//	

// make a copy of the graph
graph_t gc=g_in;
graph_t g_out;

bg_to_ct(gc,g_out, bubble_vector, list[m]);

// for(int n=0; n< list[m].size(); n++){

// std::cout<<"About to call bgtoct"<<std::endl;
// print_all_edge_info(gc);

// bg_to_ct(gc, g_out, bubble_vector[list[m][n]]);

// gc=g_out;

// }

// ct_vec.push_back(gc);
// Set size of CT
g_out[boost::graph_bundle].ct_count=list[m].size();
// std::cout<<"Set graph bundle to ct_count "<< g_out[boost::graph_bundle].ct_count<<std::endl;
	
ct_vec.push_back(g_out);	
	
}

	
}



	
}


void AmiGraph::bg_to_ct(graph_t gin, graph_t &g_out, std::vector< bubble > bubble_vector, std::vector<int> list){

// pre-process bubble_vector for overlaps

for(int n=0; n<list.size()-1; n++){
	for(int j=n+1; j<list.size(); j++){

if(	bubble_vector[list[n]].b==bubble_vector[list[j]].A){
bubble_vector[list[j]].A=bubble_vector[list[n]].a;	
	
}
	

if(	bubble_vector[list[n]].b==bubble_vector[list[j]].B){
bubble_vector[list[j]].B=bubble_vector[list[n]].a;	
	
}	
	
	}
	
}


	
for(int n=0; n< list.size(); n++){

// std::cout<<"About to call bgtoct"<<std::endl;
// print_all_edge_info(gin);

// std::cout<<"Trying to remove bubble "<<std::endl;
// std::cout<<bubble_vector[list[n]].A<<" "<<bubble_vector[list[n]].a<<" "<<bubble_vector[list[n]].b<<" "<<bubble_vector[list[n]].B<<std::endl;

bg_to_ct(gin, g_out, bubble_vector[list[n]]);

gin=g_out;

}


	
	
	
}


// three options.  both in. both out. or directed. if directed then the alphas for the legs will be the same 

void AmiGraph::descriptor_to_index(graph_t &g, vertex_vector_t &bv, edge_vector_t &be, edge_vector_t &le , bubble &bub){
// Find A, a, b and B first


if(le.size()!=2 || be.size()!=2){
	throw std::runtime_error("This is not a bubble");
}

if (g[target(le[0],g)].index_ == g[bv[0]].index_ || g[target(le[0],g)].index_ == g[bv[1]].index_){
	
bub.a=	g[target(le[0],g)].index_;
bub.A=g[source(le[0],g)].index_;
	
}

if (g[source(le[0],g)].index_ == g[bv[0]].index_ || g[source(le[0],g)].index_ == g[bv[1]].index_){
	
bub.a=	g[source(le[0],g)].index_;
bub.A=g[target(le[0],g)].index_;
	
}

//

if (g[target(le[1],g)].index_ == g[bv[0]].index_ || g[target(le[1],g)].index_ == g[bv[1]].index_){
	
bub.b=	g[target(le[1],g)].index_;
bub.B=g[source(le[1],g)].index_;
	
}

if (g[source(le[1],g)].index_ == g[bv[0]].index_ || g[source(le[1],g)].index_ == g[bv[1]].index_){
	
bub.b=	g[source(le[1],g)].index_;
bub.B=g[target(le[1],g)].index_;
	
}




// for(int i=0; i< le.size(); i++){
	
// for(int j=0; j< bv.size(); j++){	

// if( target(le[i],g)==bv[j]){
	
// bub.A=g[source(le[i],g)].index_;
// bub.a=g[target(le[i],g)].index_;	
	
// }

// if(source(le[i],g)==bv[j]){
	
// bub.b=g[source(le[i],g)].index_;
// bub.B=g[target(le[i],g)].index_;
	
// }
	
	
// }	
// }

std::pair<int,int> l1(bub.A,bub.a);
std::pair<int,int> l2(bub.b,bub.B);
std::vector< std::pair<int,int>> leg_edges;
leg_edges.push_back(l1);
leg_edges.push_back(l2);

bub.leg_edges=leg_edges;


// get bubble edge directions right 

// std::pair<int,int> b1;
// std::pair<int,int> b2;
std::vector< std::pair<int,int>> bub_edges;

for(int i=0; i<be.size(); i++){

std::pair<int,int> b1;
b1.first=g[source(be[i],g)].index_;
b1.second=g[target(be[i],g)].index_;

bub_edges.push_back(b1);
	
	
}



bub.bubble_edges=bub_edges;

	
	
}

// bubble graph to counter-term
void AmiGraph::bg_to_ct(graph_t &g, graph_t &ctg, bubble &bub ){

// set counter-term to be a copy 
ctg=g;	

// print_all_edge_info(ctg);

edge_t l1,l2;
edge_vector_t bev;
vertex_t A,a,b,B;

boost::graph_traits<graph_t>::edge_iterator ei, ei_end;

for(boost::tie(ei,ei_end)=edges(ctg); ei!=ei_end; ++ei){

// find legs first 

if( ctg[*ei].g_struct_.stat_==AmiBase::Bose){
	
	// std::cout<<"On bosonic edge from "<<ctg[source(*ei, ctg)].index_<<" to "<<ctg[target(*ei, ctg)].index_<<std::endl;
	
if( ctg[source(*ei, ctg)].index_== bub.A &&	ctg[target(*ei, ctg)].index_== bub.a){

l1=*ei;	
A=source(*ei, ctg);
a=target(*ei, ctg);

// std::cout<<"Set A1's"<<std::endl;
	
}

if( ctg[source(*ei, ctg)].index_== bub.a &&	ctg[target(*ei, ctg)].index_== bub.A){

l1=*ei;	
a=source(*ei, ctg);
A=target(*ei, ctg);
	// std::cout<<"Set A2's"<<std::endl;
}

if( ctg[source(*ei, ctg)].index_== bub.b &&	ctg[target(*ei, ctg)].index_== bub.B){

l2=*ei;	
b=source(*ei, ctg);
B=target(*ei, ctg);
	// std::cout<<"Set B1's"<<std::endl;
}

if( ctg[source(*ei, ctg)].index_== bub.B &&	ctg[target(*ei, ctg)].index_== bub.b){

l2=*ei;	
B=source(*ei, ctg);
b=target(*ei, ctg);
	// std::cout<<"Set B2's"<<std::endl;
}

	
}

/// done with bosonic. just need the b1 and b2 edges 

if (ctg[*ei].g_struct_.stat_==AmiBase::Fermi){
	
bool condition=	 (ctg[source(*ei,ctg)].index_ ==bub.a && ctg[target(*ei,g)].index_== bub.b) ||
					(ctg[source(*ei,ctg)].index_ ==bub.b && ctg[target(*ei,ctg)].index_== bub.a);
if (condition){

edge_t one_edge=*ei;

bev.push_back(one_edge);

}	
	
	
}


	
	
	
	
}



// std::cout<<"bev vector of size "<<bev.size()<<std::endl;

// print_all_edge_info(ctg);

///// now I have the appropriate descriptors for the ctg graph copy 


// Next I need to know what is the independent variable attached to the bubble.  Presumably that variable does not exist elsewhere.  

// Presume that one of the bubble edges must have an independent label - so find its index 

int indep_index=-1;
bool success=false;
std::vector<int> eps_ind;
for(int n=0; n< bev.size(); n++){
int count=0;
int temp=-1;
for (int i=0; i< ctg[bev[n]].g_struct_.alpha_.size(); i++){

// std::cout<<"On alpha "<<ctg[bev[n]].g_struct_.alpha_[i]<<std::endl;

if(	ctg[bev[n]].g_struct_.alpha_[i]!=0){

count++;	
temp=i;	
}

	
}

if(count==1){
	
indep_index=temp;
success=true;	
	
}

for(int i=0; i<ctg[bev[n]].g_struct_.eps_.size(); i++){

if(ctg[bev[n]].g_struct_.eps_[i]==1){
eps_ind.push_back(i);	
}


}

}

if(!success){throw std::runtime_error("Failed to identify independent variable in counter-term diagram");}

// now iterate through and remove that indep_index from every alpha
// boost::graph_traits<graph_t>::edge_iterator ei, ei_end;

// std::cout<<"Found indep index of "<< indep_index<<std::endl;
// print_all_edge_info(ctg);

// sort the epsilon indices
std::sort(eps_ind.begin(), eps_ind.end());

for(boost::tie(ei,ei_end)=edges(ctg); ei!=ei_end; ++ei){

ctg[*ei].g_struct_.alpha_.erase( ctg[*ei].g_struct_.alpha_.begin()+indep_index  );	
	
for(int i=0; i< eps_ind.size(); i++){

ctg[*ei].g_struct_.eps_.erase(ctg[*ei].g_struct_.eps_.begin()+ eps_ind[i]-i);

}	
	
}

// std::cout<<"Removed indep_index from all edges "<<std::endl;

// all done that part!

// print_all_edge_info(ctg);






// First need to make new edges and give them the info of the bosonic edge
// lets choose form 
// A -->-- a -->--B 
// this means leave A to a alone.  But then make a new edge from a to B with the same info as is in A to a. 
//edge_info(AmiCalc::epsilon_t epsilon, AmiCalc::alpha_t alpha, AmiCalc::stat_type type, int loop_id, spin_type spin_val)

// std::cout<<"Adding edge to "<< g[a].index_<<" - "<< g[B].index_<<std::endl;
edge_t newedge=add_edge(a,B,edge_info(g[l1].g_struct_.eps_, g[l1].g_struct_.alpha_,  g[l1].g_struct_.stat_, g[l1].fermi_loop_id, g[l1].spin), ctg).first;

// std::cout<<"Added edge"<<std::endl;
g[newedge].label=labelled; // mark the edge as labelled. 

// std::cout<<"Added an edge "<<std::endl;

// now we will remove the edges 

// remove_edge(l1,ctg); // don't remove l1, since it was fine. 
remove_edge(l2,ctg);
for(int i=0; i< bev.size(); i++){

remove_edge(bev[i],ctg);	
	
}

// need to delete the extra vertex .... 
clear_vertex(b,ctg);
remove_vertex(b,ctg);




	
	
	
	
} 

bool AmiGraph::not_labelled(graph_t &g){
// std::cout<<"In not labelled function"<<std::endl;
boost::graph_traits<graph_t>::edge_iterator ei, ei_end;

bool allzero=true;

for( boost::tie(ei,ei_end)=edges(g); ei!=ei_end; ++ei){
	
	// std::cout<<"on edge "<< g[source(*ei,g)].index_<<"-"<< g[target(*ei,g)].index_<<std::endl;
	
if( g[*ei].g_struct_.alpha_.size()==0){
	// std::cout<<"Alpha is size zero"<<std::endl;
	return true;
}

if( g[*ei].g_struct_.eps_.size()==0){
	// std::cout<<"Eps is size zero"<<std::endl;
	return true;
}

if(allzero){

for(int i=0; i< g[*ei].g_struct_.alpha_.size(); i++){

if(g[*ei].g_struct_.alpha_[i]!=0){ allzero=false;}

}	

// std::cout<<"eps size is "<< g[*ei].g_struct_.eps_.size()<<std::endl;
for(int i=0; i< g[*ei].g_struct_.eps_.size(); i++){

// std::cout<<"eps is "<<g[*ei].g_struct_.eps_[i]<<std::endl;
if(g[*ei].g_struct_.eps_[i]!=0){ allzero=false;}

}
	
	
}

	
	
	
}

if(allzero){return true;}
	
	
return false;	
	
}
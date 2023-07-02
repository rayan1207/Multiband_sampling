#include "amigraph.hpp"




void AmiGraph::symM_make_graphs(symM &SM){
	
for(int ord=0; ord< SM.size(); ord++){
std::cout<<"On ord "<<ord<<std::endl;

construct_graphs(SM[ord]);


// for(int n=0; n< SM[ord].size(); n++){
// std::cout<<"On n="<<n<<std::endl;
// for(int p=0; p<SM[ord][n].plist.size(); p++){
	
// std::cout<<"("<<SM[ord][n].plist[p].first<<","<<SM[ord][n].plist[p].second<<")"<<std::endl;
	
	
// }

// this seems redundent but is subtle pointer for struct issue in c++ I think 
// symDET_element se=SM[ord][n];
// convert_symelement_to_graph(&SM[ord][n]);	
	
// }
}	
	

}

void AmiGraph::symM_bandmaps(symM &SM){
	
for(int ord=0; ord< SM.size(); ord++){

for(int g=0; g< SM[ord].size(); g++){
	
generate_bandmap(SM[ord][g].graph, SM[ord][g].bmap);	

}
}	
		
}

// this should probably use R0_ that is already populated
void AmiGraph::generate_bandmap(graph_t &g, band_map &bm){

bm.clear();

// bm.resize(g[random_edge(g, rand_gen)].g_struct_.eps_.size());

boost::graph_traits<graph_t>::edge_iterator ei, ei_end;

// THis only gets the internal lines 
for( boost::tie(ei,ei_end)=edges(g); ei!=ei_end; ++ei){

if( g[*ei].g_struct_.stat_==AmiBase::Fermi){
	
for(int i=0; i< g[*ei].g_struct_.eps_.size(); i++){

if(g[*ei].g_struct_.eps_[i]==1){
	
	bm.push_back(g[*ei].band_element);
	
// std::cout<<"Band element is "<< g[*ei].band_element.first<<" "<<g[*ei].band_element.second<<std::endl;	
	
// std::cout<<"Set bm to "<<bm[i].first <<" "<< bm[i].second<<std::endl;
}

}	
	
	
	
}



}	


// iterate through again and toss external lines to the end of the bmap 
for( boost::tie(ei,ei_end)=edges(g); ei!=ei_end; ++ei){

if( g[*ei].g_struct_.stat_==AmiBase::Fermi){

bool allzero=true;

	
for(int i=0; i< g[*ei].g_struct_.eps_.size(); i++){
if(g[*ei].g_struct_.eps_[i]==1){
allzero=false;
}
}	
	
if(allzero){

bm.push_back(g[*ei].band_element);
}	
	
	
}
}
	
	
}

void AmiGraph::construct_rtss_ami_sol(symDET_element &SE){
	
if(SE.use_rtss==1){

ami.amibase.construct(SE.rtss.ami_parms_, SE.rtss.R0_, SE.rtss.R_, SE.rtss.P_, SE.rtss.S_);
	
}	
	
}


 
void AmiGraph::symDET_element_construct_ami_sol(symDET_element &SE, double ereg){

graph_to_R0(SE.graph, SE.ss.R0_);
extract_bose_alphas(SE.graph, SE.ss.bose_alphas_);
extract_Ulist(SE.graph, SE.Ulist);
int loops=count_fermi_loops(SE.graph);
int order = graph_order(SE.graph);


	// SM[ord][g].ss.prefactor_= SM[ord][g].sign*std::pow(-1.0/2.0,order)*std::pow(2,loops)/(double)factorial(order);
SE.ss.prefactor_=SE.sign*std::pow(-1.0/2.0,order)*std::pow(2,loops)/(double)factorial(order);
	
AmiBase::ami_parms amiparms(order, ereg, ami_parameters.TYPE_);
	
SE.ss.ami_parms_=amiparms;	
ami.amibase.construct(SE.ss.ami_parms_, SE.ss.R0_, SE.ss.R_, SE.ss.P_, SE.ss.S_);
// all we need to do is reconstruct	
	
	
	
} 


void AmiGraph::symM_construct_ami_sol(symM &SM, double ereg){

std::cout<<"SM is size "<<SM.size()<<std::endl;
for(int ord=0; ord<SM.size(); ord++){
	std::cout<<"On order "<<ord<<std::endl;
	std::cout<<"With size "<< SM[ord].size()<<std::endl;
	for(int g=0; g<SM[ord].size(); g++){
	// std::cout<<"On graph "<<g<<std::endl;
	// print_all_edge_info(SM[ord][g].graph);
	graph_to_R0(SM[ord][g].graph, SM[ord][g].ss.R0_);
	
	// std::cout<<"Extract bose "<<std::endl;
	extract_bose_alphas(SM[ord][g].graph, SM[ord][g].ss.bose_alphas_);
	
	extract_Ulist(SM[ord][g].graph, SM[ord][g].Ulist);
	
	// TODO: check if this loop counting is accurate or not 
	int loops=count_fermi_loops(SM[ord][g].graph);
	
	
	int order=graph_order(SM[ord][g].graph);
	// unsure if missing minus sign 
	SM[ord][g].ss.prefactor_= SM[ord][g].sign*std::pow(-1.0/2.0,order)*std::pow(2,loops)/(double)factorial(order);
	// std::cout<<"Found loops="<<loops<<std::endl;

// SE.ss.prefactor_=SE.sign*std::pow(-1,order)*std::pow(2,loops);
	
	
	// int order=graph_order(SM[ord][g].graph);
	// SM[ord][g].ss.prefactor_=SM[ord][g].sign*std::pow(-1,order);
	
	AmiBase::ami_parms amiparms(order, ereg, ami_parameters.TYPE_);
	
	SM[ord][g].ss.ami_parms_=amiparms;
	
	
	
	ami.amibase.construct(SM[ord][g].ss.ami_parms_, SM[ord][g].ss.R0_, SM[ord][g].ss.R_, SM[ord][g].ss.P_, SM[ord][g].ss.S_);
	
	
	
	}
}
}

void AmiGraph::construct_minimal_avbv_sets(U_matrix_type &UM, std::vector<std::vector<int>> &avlist,std::vector<std::vector<int>> &bvlist, int order, int band_num, int ext_index){
	
	
std::vector< std::vector<int>> avbvlist;	
	// avbvlist.clear();

// Generate avbv list except for external 
int maxval=band_num;
int maxlength=4*order;
do{
	
std::vector<std::vector<int>> outlist;

avbv_add_col( avbvlist, outlist, maxval, maxlength);

avbvlist=outlist;	
	
}while(avbvlist[0].size()<maxlength);

avlist.clear();
bvlist.clear();

avlist.resize(avbvlist.size());
bvlist.resize(avbvlist.size());

for(int i=0; i< avbvlist.size(); i++){
	
avlist[i].resize(2*order);
bvlist[i].resize(2*order);

for(int j=0; j< 2*order; j++){

avlist[i][j]=avbvlist[i][2*j];
bvlist[i][j]=avbvlist[i][2*j+1]; 	
	
}

avlist[i].push_back(ext_index);
bvlist[i].push_back(ext_index);

double result= eval_Ulist(UM, avlist[i], bvlist[i]);

if (result==0){
	
	avlist.erase(avlist.begin()+i);
	bvlist.erase(bvlist.begin()+i);
	avbvlist.erase(avbvlist.begin()+i);
	i--;	
	
}


	
	
}

	
	
	
}



// These functions do not work 
/* 
void AmiGraph::construct_offdiag_delta_sets(std::vector< std::vector<int>> &dlist, int order, int band_num){

dlist.clear();

int maxval=band_num;
int maxlength=2*order;	

do{
std::vector<std::vector<int>> outlist;	
avbv_add_col(dlist, outlist, maxval, maxlength);

dlist=outlist;

} while( dlist[0].size()<maxlength);

// for(int i=0; i< dlist.size(); i++){

// dlist[i].push_back(in_band);
// dlist[i].push_back(out_band);

// }	
	
	
}


void AmiGraph::get_avbv_from_offdiag_delta( symDET_element &SD,  std::vector<int> &delta, std::vector<int> &ai, std::vector<int> &bi, int in_band, int out_band){

// this is hard coded for sigma graphs? or G graphs?
int ext_index= (SD.plist.size()-1)/2;

ai=delta;
bi.resize(ai.size(),0);

ai.push_back(out_band);
bi.push_back(in_band);

for( int i=0; i< SD.plist.size(); i++){
	
if (SD.plist[i].first==ext_index){
ai[SD.plist[i].second]=bi[SD.plist[i].first];
}	else

if(SD.plist[i].second==ext_index){
	bi[SD.plist[i].first]=ai[SD.plist[i].second];
}
else {
	bi[SD.plist[i].first]=ai[SD.plist[i].second];	
}
	
}


	
} */




void AmiGraph::construct_delta_sets(std::vector< std::vector<int>> &dlist, int order, int band_num, int ext_index){

dlist.clear();

int maxval=band_num;
int maxlength=2*order;	

do{
std::vector<std::vector<int>> outlist;	
avbv_add_col(dlist, outlist, maxval, maxlength);

dlist=outlist;

} while( dlist[0].size()<maxlength);

for(int i=0; i< dlist.size(); i++){

dlist[i].push_back(ext_index);

}	
	
	
}

void AmiGraph::get_avbv_from_delta( symDET_element &SD,  std::vector<int> &delta, std::vector<int> &ai, std::vector<int> &bi){

ai=delta;
bi.resize(ai.size(),0);

for( int i=0; i< SD.plist.size(); i++){
	
bi[SD.plist[i].first]=ai[SD.plist[i].second];	
	
	
}


	
}



void AmiGraph::general_construct_avbv_sets(std::vector<std::vector<int>> &avlist,std::vector<std::vector<int>> &bvlist, int order, int band_num, int ext_index){
	
	
std::vector< std::vector<int>> avbvlist;	
	// avbvlist.clear();


// Generate avbv list except for external 
int maxval=band_num;
int maxlength=4*order;
do{
	
std::vector<std::vector<int>> outlist;

avbv_add_col( avbvlist, outlist, maxval, maxlength);

avbvlist=outlist;	
	
}while(avbvlist[0].size()<maxlength);

avlist.clear();
bvlist.clear();

avlist.resize(avbvlist.size());
bvlist.resize(avbvlist.size());

for(int i=0; i< avbvlist.size(); i++){
	
avlist[i].resize(2*order);
bvlist[i].resize(2*order);

for(int j=0; j< 2*order; j++){

avlist[i][j]=avbvlist[i][2*j];
bvlist[i][j]=avbvlist[i][2*j+1]; 	
	
}

avlist[i].push_back(ext_index);
bvlist[i].push_back(ext_index);
	
	
}

	
	
	
}


void AmiGraph::general_construct_avbv_sets(std::vector<std::vector<int>> &avlist,std::vector<std::vector<int>> &bvlist, int order, int band_num, int in_index, int out_index){
	
	
std::vector< std::vector<int>> avbvlist;	
	// avbvlist.clear();


// Generate avbv list except for external 
int maxval=band_num;
int maxlength=4*order;
do{
	
std::vector<std::vector<int>> outlist;

avbv_add_col( avbvlist, outlist, maxval, maxlength);

avbvlist=outlist;	
	
}while(avbvlist[0].size()<maxlength);

avlist.clear();
bvlist.clear();

avlist.resize(avbvlist.size());
bvlist.resize(avbvlist.size());

for(int i=0; i< avbvlist.size(); i++){
	
avlist[i].resize(2*order);
bvlist[i].resize(2*order);

for(int j=0; j< 2*order; j++){

avlist[i][j]=avbvlist[i][2*j];
bvlist[i][j]=avbvlist[i][2*j+1]; 	
	
}

avlist[i].push_back(in_index);
bvlist[i].push_back(out_index);
	
	
}

	
	
	
}



void AmiGraph::avbv_add_col( std::vector<std::vector<int> > &inlist, std::vector<std::vector<int> > &outlist, int maxval, int maxlength){

outlist.clear();

if(inlist.size()==0){ 

std::vector<int> entry;
for(int i=0; i<maxval; i++){
	entry.push_back(i);
	outlist.push_back(entry);
	entry.clear();
}



return;}

if(inlist[0].size()>maxlength){return;}

std::vector<int> vals;

for( int i=0; i< maxval; i++){
	
	vals.push_back(i);
	
}


for( int i=0; i< inlist.size(); i++){

for( int j=0; j< vals.size(); j++){

std::vector<int> entry;

for(int m=0; m< inlist[i].size(); m++){
	
entry.push_back(inlist[i][m]);
	
}
entry.push_back(vals[j]);

outlist.push_back(entry);
entry.clear();
	
	
	
}
	
}



/* 

std::cout<<"Resulting outlist is now "<<std::endl;

for(int i=0; i< outlist.size(); i++){
	
	for(int j=0; j< outlist[i].size(); j++){
		
	std::cout<<outlist[i][j]<<" ";	
		
	}
	std::cout<<std::endl;
	
} */
	
	
	
	
	
}


/* void AmiGraph::construct_avbv_sets(std::vector< std::vector<int>> &avlist, std::vector< std::vector<int>> &bvlist, int order, int num_bands, int ext_index){

avlist.clear();
bvlist.clear();
// for (int ord=0; ord<order; ord++){
std::vector<int> v(4*order);

for(int ones=0; ones< v.size()+1; ones++){
	// std::cout<<"Size is "<< avlist.size()<<std::endl;

for(int i=0; i< v.size(); i++){
	
v[i]=0;	
	
}

for(int j=ones; j<v.size(); j++){
v[j]=1;	
}

do {
	
	// std::cout<<"On permutation"<<std::endl;
	// for(int i=0; i<v.size(); i++){
		// std::cout<<v[i]<<" ";
		
	// }
	// std::cout<<std::endl;
	
	std::vector<int> av,bv;
	av.resize(2*order);
	bv.resize(2*order);
	
	for(int aind=0; aind<v.size()/2; aind++){
	av[aind]=v[aind*2];	
	bv[aind]=v[aind*2+1];	
	}
	
	

av.push_back(ext_index);
bv.push_back(ext_index);

// std::cout<<"Pushing into avlist ";
// for(int i=0; i< av.size(); i++){
// std::cout<<av[i]<<" ";
// }	
// std::cout<<std::endl;

avlist.push_back(av);
bvlist.push_back(bv);
av.clear();
bv.clear();		
	
}while(std::next_permutation(v.begin(),v.end()));







}
// }	
	
	
	
} */

void AmiGraph::random_avbv(std::vector<int> &av, std::vector<int> &bv, int order, int num_bands, int ext_index){
	
av.resize(order*2);
bv.resize(order*2);
randomize(av, 0, num_bands-1);
randomize(bv, 0, num_bands-1);

av.push_back(ext_index);
bv.push_back(ext_index);	
	
	
	
}

void AmiGraph::randomize(std::vector<int> &v, int min, int max){
	
for(int i=0; i< v.size(); i++){	
	
v[i]=random_int(min,max);	
}
}

bool AmiGraph::is_diagonal(band_map &bm, std::vector<int> &av, std::vector<int> &bv){
	

for(int i=0; i< bm.size(); i++){

// std::cout<<"Checking equality between "<<std::endl;
// std::cout<<bv[bm[i].first]<<"?=?"<<av[bm[i].second]<<std::endl;

if (av[ bm[i].second] != bv[ bm[i].first]){
	return false;
}

}	

return true;	
}

void AmiGraph::assign_solution(symDET_element &SD){
	
// set rtss to ss // TODO: this might be expensive. should consider alternative approaches
SD.use_rtss=0;
SD.rtss=SD.ss;	

// std::cout<<"Working on rtss with R0"<<std::endl;

// ami.print_g_prod_info(SD.rtss.R0_);
// std::cout<<"R0_ Species list"<<std::endl;

// for(int m=0; m< SD.rtss.R0_.size(); m++){

// std::cout<<SD.rtss.R0_[m].species_<<std::endl;	
// std::cout<<SD.rtss.R0_[m].species_<<std::endl;	
	
// }
	
// first check if species are the same in R0_ of ss.  if no, then use_rtss=0 and return. if yes - let 	

for(int i=0; i< SD.ss.R0_.size(); i++){
for(int j=i+1; j< SD.ss.R0_.size(); j++){

// std::cout<<"eps For i: ";
// for(int m=0; m< SD.rtss.R0_[i].eps_.size(); m++){

// std::cout<<	SD.rtss.R0_[i].eps_[m]<<" ";
	
// }
// std::cout<<std::endl;

// std::cout<<"eps For j: ";
// for(int m=0; m< SD.rtss.R0_[j].eps_.size(); m++){

// std::cout<<	SD.rtss.R0_[j].eps_[m]<<" ";
	
// }
// std::cout<<std::endl;


if( SD.ss.R0_[j].species_==SD.ss.R0_[i].species_){
// std::cout<<"Species equal on i="<<i<<"    j="<<j<<std::endl;

SD.use_rtss=1;	

SD.rtss.R0_[j].eps_=SD.rtss.R0_[i].eps_;
	
}

}	
	
}





// for(int i=0; i< SD.ss.R0_.size(); i++){

// for(int eps=0; eps< SD.ss.R0_[i].eps_.size(); eps++){

// for(int eps2=eps+1; eps2< SD.ss.R0_[i].eps_.size(); eps2++){

// if( SD.ss.R0_[i].eps_

// }	
	
// }

// }	
	
	
}



// this presumes we already determined it to be diagonal 
/// dup_map would nominally be a vector of length max_species and have entries dup_map[i]=i. but if there is a duplicate energy entry then they get merged, so only one species label exists, but the duplicate remains in av/bv to be summed over. 
void AmiGraph::assign_species(symDET_element &SD, std::vector<int> &av, std::vector<int> &bv, std::vector<int> &dup_map){

// std::cout<<"Assigning species at runtime "<<std::endl;
// std::cout<< SD.bmap.size()<<" "<< SD.ss.R0_.size()<<std::endl;



for(int i=0; i< SD.ss.R0_.size(); i++){

for(int eps=0; eps< SD.ss.R0_[i].eps_.size(); eps++){

if (SD.ss.R0_[i].eps_[eps]==1){
	//TODO: Somehow this is very subtle in how contractions are assigned - the bv list is ordered while av list is not 
SD.ss.R0_[i].species_=dup_map[bv[SD.bmap[i].first]]; //bv[SD.bmap[eps].first];
// std::cout<<"Species assigned "<<SD.ss.R0_[i].species_<<" for bmap ("<<SD.bmap[i].first<<","<<SD.bmap[i].second<<")"<<std::endl; 
}

}




}	
	
	
	
}

void AmiGraph::read_U_matrix(std::string filename, U_matrix_type &UM, int maxval){
std::ifstream infile_stream;

infile_stream.open(filename);	

if(infile_stream.fail()) // checks to see if file opended 
    { 
	std::cout<<filename<<std::endl;
      throw std::runtime_error("Could not open U file");
    } 		
	
std::string line;

while (std::getline(infile_stream,line)){

std::stringstream ss(line);

int a,b,c,d;
double value;

bool read= bool( ss>> a);

if(read){
ss>>b>>c>>d>> value;

// std::cout<<"Read line and found "<<a<<" "<<b<<" "<<c<<" "<<d<<" "<<value<<std::endl;

if( a<=maxval &&b<=maxval && c<=maxval && d<= maxval){
a--;
b--;
c--;
d--;

// std::cout<<UM.size()<<std::endl;

UM[a][b][c][d]=value;
}


}

}	

infile_stream.close();
	
	
}

void AmiGraph::reduce_avbv_sets(U_matrix_type &UM, std::vector< std::vector< int>> &avlist,std::vector< std::vector< int>> &bvlist){
	
std::vector< std::vector< int>> newa,newb;
	
for (int i=0; i< avlist.size(); i++){

double result= eval_Ulist(UM, avlist[i], bvlist[i]);


if(result!=0){
	
	newa.push_back(avlist[i]);
	newb.push_back(bvlist[i]);
	
}


// if(result==0){
	
	// avlist.erase(avlist.begin()+i);
	// bvlist.erase(bvlist.begin()+i);
	// i--;	
	
// }


}	

avlist=newa;
bvlist=newb;
	
	
}
	
	
double AmiGraph::eval_Ulist(U_matrix_type &Uvals, std::vector<int> &av, std::vector<int> &bv){

double result=1.0;

int length=av.size()-1;
length=length/2;

// std::cout<<"===="<<std::endl;
for(int item=0; item< length; item++){
// int a,b,c,d;
// a=

//from the U matrix - form the antisymmetrize Uabcd=Vabcd - Vadcb
// std::cout<<"Mult U_"<<av[0+item*2]<<bv[0+item*2]<<av[1+item*2]<<bv[1+item*2]<<"="<<Uvals[av[0+item*2]][bv[0+item*2]][av[1+item*2]][bv[1+item*2]]<<std::endl;

double uval=Uvals[av[0+item*2]][bv[0+item*2]][av[1+item*2]][bv[1+item*2]];///4.0;

// double uval=(Uvals[av[0+item*2]][bv[0+item*2]][av[1+item*2]][bv[1+item*2]] - 
				// Uvals[av[0+item*2]][bv[1+item*2]][av[1+item*2]][bv[0+item*2]])/4.0;
// std::cout<<"Uval is "<<uval<<std::endl;	
if(uval==0){return 0.0;}			
result=result*uval;	
	
	
}

return result;

}	


void AmiGraph::extract_Ulist(graph_t &g, Ulist_type &U){
	
// std::vector<int> Uitem(4);

boost::graph_traits<graph_t>::edge_iterator ei, ei_end;

		for( boost::tie(ei,ei_end)=edges(g); ei!=ei_end; ++ei){

			if( g[*ei].g_struct_.stat_== AmiBase::Bose){
			U.push_back(g[*ei].fourindex);
				
			}
		
	
	
}

}


void AmiGraph::symM_filter_mpi(symM &SM, int mpi_rank, int comm_size){

	
if(mpi_rank==0){std::cout<<"Creating linearized index for filtering"<<std::endl;}	
	
// create linearized index with values:  ord g keep==1
std::vector<  std::vector <int>  > linearized_index, bc_linearized_index;
// linearized_index.resize(SM.size());
// bc_linearized_index.resize(SM.size());

std::vector<int> throw_away(3);
for(int ord=0; ord< SM.size(); ord++){
	for(int num=0; num<SM[ord].size(); num++){

// std::cout<<"Appending to index "<< ord<<" "<< num<<std::endl;

	throw_away[0]=ord;
	throw_away[1]=num;
	throw_away[2]=0;  // 0 will mean to keep it. anything else to discard it 
	
linearized_index.push_back(throw_away);
bc_linearized_index.push_back(throw_away);	
	}
}

std::cout<<"Linearized index size "<<linearized_index.size()<<std::endl;



// for(int i=0; i< linearized_index.size(); i++){
	
// std::cout<<i<<" "<<linearized_index[i][0]<<" "<<linearized_index[i][1]<<" "<<linearized_index[i][2]<<std::endl;	
	
	
// }






// filter each item and store result in linearized_index 
// if(mpi_rank==0){std::cout<<"Filtering each item "<<std::endl;}
for(int i=0; i< linearized_index.size(); i++){
	// if(mpi_rank==0){std::cout<<i<<" "; std::cout<<linearized_index[i][0]<<" ";std::cout<<linearized_index[i][1]<<" ";std::cout<<linearized_index[i][2]<<std::endl;}

if(i%comm_size!=mpi_rank){ continue;}

if(i>1000 && i % 500 ==0){
std::cout<<"On i=:"<< i<< " of "<< linearized_index.size()<<std::endl;
}


if(!is_connected_graph(SM[linearized_index[i][0]][linearized_index[i][1]].graph)){

	linearized_index[i][2]=1;
	continue;
	
}	



if(is_1P_reducible_graph(SM[linearized_index[i][0]][linearized_index[i][1]].graph)){
	linearized_index[i][2]=1;
	continue;
	
}
	
	
}

// now each process has their own linearized index 
// if(mpi_rank==0){std::cout<<"MPI Call"<<std::endl;}
for(int i=0; i< linearized_index.size(); i++){	
	
MPI_Allreduce(&linearized_index[i][2], &bc_linearized_index[i][2], 1, MPI_INT, MPI_SUM,  MPI_COMM_WORLD);
	
}


symM outSM;
symDET term;

int lastord=0;

outSM.resize(SM.size());

// if(mpi_rank==0){std::cout<<"Collecting results "<<std::endl;}
for (int i=0; i< bc_linearized_index.size(); i++){
int ord=bc_linearized_index[i][0];
int g=bc_linearized_index[i][1];

// std::cout<<ord<<" "<<g<<" "<<lastord<<std::endl;

if(ord!=lastord){
	if (term.size()!=0){
		// std::cout<<"Pushing term into outSM element "<<lastord<<std::endl;
	outSM[lastord]=term;
	}
term.clear();	
}
	
if(	bc_linearized_index[i][2]==0){

term.push_back(SM[ord][g]);	
	
}
	
lastord=ord;	
}
outSM[lastord]=term;
// if(mpi_rank==0){std::cout<<"end  "<<std::endl;}

SM.clear();
SM.shrink_to_fit();
SM=outSM;
	 
	
}


void AmiGraph::symM_filter(symM &SM){



for(int ord=0; ord< SM.size(); ord++){
// int count=0;
std::cout<<"On order "<<ord<<std::endl;
for(int g=0; g< SM[ord].size(); g++){
// std::cout<<count<<std::endl;

// count++;	
// std::cout<<"g"<<g<<std::endl;	
if(g>1000 && g % 500 ==0){
std::cout<<"On g=:"<< g<< " of "<< SM[ord].size()<<std::endl;
}	
	
// std::cout<<"Considering graph "<<std::endl;
// print_all_edge_info(SM[ord][g].graph);	

// graph_t copy;
// copy=SM[ord][g].graph;
// symmetrize_bosonic_lines(copy);



// TODO the connected components filter seems to fail for bubble graphs - consider it dangerous for now 
if(!is_connected_graph(SM[ord][g].graph)){

	
	
// if(connected_components(SM[ord][g].graph)!=1){
// std::cout<<"Apparently it is not connected"<<std::endl;

	SM[ord].erase(SM[ord].begin() + g);
	g--;
	continue;
	
}

	
// this won't work on disconnected diagrams - so much filter those out first
// if(has_fock_insertion(SM[ord][g].graph)){
	// SM[ord].erase(SM[ord].begin() + g);
	// g--;
	// continue;
	
// }

// if(count_tadpoles(SM[ord][g].graph)!=0){
	
	// SM[ord].erase(SM[ord].begin() + g);
	// g--;
	// continue;
	
// }


if(is_1P_reducible_graph(SM[ord][g].graph)){
	SM[ord].erase(SM[ord].begin() + g);
	g--;
	continue;
	
	
}

// TODO: unsure if this should be here:


// if(is_1P_bose_reducible_graph(SM[ord][g].graph)){
	// SM[ord].erase(SM[ord].begin() + g);
	// g--;
	// continue;
	
	
// }



}
}	
	
	
	
}




void AmiGraph::symM_label(symM &SM){
	
std::cout<<"Labelling"<<std::endl;
for(int ord=0; ord< SM.size(); ord++){
std::cout<<"On order"<<std::endl;
for(int g=0; g< SM[ord].size(); g++){

if(g%200==0){
std::cout<<"On graph "<< g<<" out of "<< SM[ord].size()<<std::endl;
}
bool success=false;
sys_label(SM[ord][g].graph, success);

// std::cout<<"Exited sys label function"<<std::endl;
if(!success){
	std::cerr<<"Warning: Failed to label systematically graph ord"<<ord<<" num "<<g<<std::endl;
std::cerr<<"Falling back to repeated half-random labelling - does this graph have a tadpole?"<<std::endl;
// systematic_vertex_label(SM[ord][g].graph);	
repeated_labelling(SM[ord][g].graph, success);

if(!success){throw std::runtime_error("Failed to label graph - exiting"); }
if(success){std::cerr<<"Warning Resolved: Random labelling succeeded for graph ord "<<ord<<" num "<<g<<std::endl;}
}	

}
}	
	
	
	
	
}




void AmiGraph::construct_graphs(symDET &SD){

for( int n=0; n<SD.size(); n++){

// std::cout<<"Converting plist to graph for n "<<n<<std::endl;


SD[n].graph=return_graph(SD[n]);

number_vertices(SD[n].graph);
	
}

	
	
	
	
}

void AmiGraph::print_plist(std::vector< std::pair<int,int> > &plist){
 
for(int i=0; i<plist.size(); i++){
std::cout<<"("<<plist[i].first<<","<<plist[i].second<<") ";

}  
std::cout<<std::endl;
  
}

std::vector< std::pair<int,int> > AmiGraph::return_contraction(graph_t g){
 
// number_vertices(g);
std::vector< std::pair<int,int> > output;

boost::graph_traits<graph_t>::edge_iterator ei, ei_end;

// THis only gets the internal lines 
int current=0;
for( boost::tie(ei,ei_end)=edges(g); ei!=ei_end; ++ei){

if(g[*ei].g_struct_.stat_==AmiBase::Bose){
  
g[source(*ei,g)].index_=current;
current++;
g[target(*ei,g)].index_=current;
current++;  
  
}

}

vertex_vector_t vv;
edge_vector_t ev;

find_external_vertices(g,vv,ev);

if(vv.size()!=2){throw std::runtime_error("Expected to find 2 external vertices");}

for(int i=0; i< vv.size(); i++){

g[vv[i]].index_=current;
  
}

// at this point the vertices have in principle been relabelled and so just collecting the pairs should work 

for( boost::tie(ei,ei_end)=edges(g); ei!=ei_end; ++ei){

if(g[*ei].g_struct_.stat_==AmiBase::Fermi){

std::pair<int,int> this_pair(g[source(*ei,g)].index_,g[target(*ei,g)].index_); 

output.push_back(this_pair);
  
}}
 
  
 return output; 
}

AmiGraph::graph_t AmiGraph::return_graph(symDET_element se){

graph_t gout;
	

int num_vert=se.plist.size()+1;
int in=se.plist.size()-1;
int out=se.plist.size();

// std::cout<<"Num vertices, in, out "<<num_vert<<" "<<in<<" "<<out<<std::endl;

vertex_vector_t vertices(num_vert);

// for(int i=0; i<se.plist.size(); i++){
	
	// std::cout<<"Checking PLIST "<<i<<" ("<< se.plist[i].first<<","<<se.plist[i].second<<")"<<std::endl;
// }



// std::cout<<"Converting to graph - Creating vertices "<<std::endl;
// make disconnected vertices
for (int i=0; i< num_vert; i++){

vertices[i]=	add_vertex(gout);	
	
}

// std::cout<<"Adding Fermionic lines "<<se.plist.size()<<std::endl;
// add all of the fermionic edges unless it is connecting the external to itself 
for(int i=0; i<se.plist.size(); i++){
	
	// std::cout<<"On plist pair "<<i<<" ("<< se.plist[i].first<<","<<se.plist[i].second<<")"<<std::endl;
int sour,tar;

sour=se.plist[i].first;
tar=se.plist[i].second;

if(se.plist[i].first==se.plist.size()-1){
sour=in;
}
if(se.plist[i].second==se.plist.size()-1){
tar=out;	
}

// if this would connect the external lines instead connect in to out 
if(se.plist[i].first==in && se.plist[i].second==in){
sour=in;
tar=out;	
}

// std::pair<int,int> thispair(sour,tar);
// edge_t test=
add_edge(vertices[sour],vertices[tar], edge_info(AmiBase::Fermi, se.plist[i]), gout);	
	
// std::cout<<"Testing band_element "<< gout[test].band_element.first<<" "<<	gout[test].band_element.second<<std::endl;
	
}

// check if 

// add bosonic lines 
// std::cout<<"Adding bosonic lines"<<std::endl;
for(int i=0; i<se.plist.size()-1; i=i+2){
	
std::vector<int> four;

four.push_back(se.plist[i].second);	
four.push_back(se.plist[i].first);	
four.push_back(se.plist[i+1].second);	
four.push_back(se.plist[i+1].first);	

add_edge(vertices[se.plist[i].first], vertices[se.plist[i+1].first],edge_info(AmiBase::Bose,four), gout);	
	
}
		
	
	
return gout;	
	
}

void AmiGraph::convert_symelement_to_graph( symDET_element &se){

int num_vert=se.plist.size()+1;
int in=se.plist.size()-1;
int out=se.plist.size();

// std::cout<<"Num vertices, in, out "<<num_vert<<" "<<in<<" "<<out<<std::endl;

vertex_vector_t vertices(num_vert);

// for(int i=0; i<se.plist.size(); i++){
	
	// std::cout<<"Checking PLIST "<<i<<" ("<< se.plist[i].first<<","<<se.plist[i].second<<")"<<std::endl;
// }



// std::cout<<"Converting to graph - Creating vertices "<<std::endl;
// make disconnected vertices
for (int i=0; i< num_vert+1; i++){

vertices[i]=	add_vertex(se.graph);	
	
}

// std::cout<<"Adding Fermionic lines "<<se.plist.size()<<std::endl;
// add all of the fermionic edges unless it is connecting the external to itself 
for(int i=0; i<se.plist.size(); i++){
	
	// std::cout<<"On plist pair "<<i<<" ("<< se.plist[i].first<<","<<se.plist[i].second<<")"<<std::endl;
int sour,tar;

sour=se.plist[i].first;
tar=se.plist[i].second;

if(se.plist[i].first==se.plist.size()-1){
sour=in;
}	
if(se.plist[i].second==se.plist.size()-1){
tar=out;	
}

// if this would connect the external lines, skip 
if(se.plist[i].first==in && se.plist[i].second==in){
continue;	
}

// std::pair<int,int> thispair(sour,tar);
	
add_edge(vertices[sour],vertices[tar], edge_info(AmiBase::Fermi, se.plist[i]), se.graph);	
	
}

// check if 

// add bosonic lines 
// std::cout<<"Adding bosonic lines"<<std::endl;
for(int i=0; i<se.plist.size()-1; i=i+2){
	
std::vector<int> four;

four.push_back(se.plist[i].second);	
four.push_back(se.plist[i].first);	
four.push_back(se.plist[i+1].second);	
four.push_back(se.plist[i+1].first);	

add_edge(vertices[se.plist[i].first], vertices[se.plist[i+1].first],edge_info(AmiBase::Bose,four), se.graph);	
	
}
	
	
	
}



/* 

// here we will remove all 1P reducible and disconnected graphs - and keep tadpoles
void AmiGraph::symM_filter_sigma(symM &SM){
	
for(int ord=0; ord< SM.size(); ord++){




}	
	
	
} */


void AmiGraph::construct_symM_full(symM &SM, int min, int max){
	
SM.resize(max+1);	
	
for (int i=min ; i<max+1; i++){

int length=i; //2*i+1;

symDET SD;
make_symDET(length,SD);

SM[i]=SD;


}




// even entries contain connected graphs 
// std::vector< std::vector< int>> collection;

// get_coll(collection, 5);

for(int i=1; i<max; i=i+2){
symDET SD_OUT;
int shift=i;
SD_multiply(SM[i], SM[max-i], SD_OUT, shift);

std::cout<<"Multiplication gave n="<<SD_OUT.size()<<std::endl;


for(int i=0; i< SD_OUT.size(); i++){
	std::cout<<i<<std::endl;
std::cout<<"---"<<std::endl;	
for(int j=0; j< SD_OUT[i].plist.size(); j++){
std::cout<<"("<<SD_OUT[i].plist[j].first<<","<<SD_OUT[i].plist[j].second<<")"<<std::endl;
}	
std::cout<<"---"<<std::endl;	
	
}

}
 

// for(int i=0; i< collection.size(); i++){
	// std::cout<<"{ ";
	// for(int j=0; j<collection[i].size(); j++){
	// std::cout<<collection[i][j]<<" ";
	// }
	// std::cout<<" }"<<std::endl;
// }




	
}

void AmiGraph::SD_multiply(symDET &SD1, symDET &SD2, symDET &SD_out, int shift){

// SD_out=SD1;

for(int one=0; one< SD1.size(); one++){
	for (int two=0; two<SD2.size();two++){
	
	symDET_element sde;
	
	sde.plist=SD1[one].plist;
	
	std::pair<int,int> p(0,0);
	for(int i=0; i< SD2[two].plist.size(); i++){
		
		p.first=SD2[two].plist[i].first+shift;
		p.second=SD2[two].plist[i].second+shift;
		sde.plist.push_back(p);
		
	}
	
	// sde.plist.insert(sde.plist.end(), SD2[two].plist.begin(),SD2[two].plist.end());
	sde.sign=SD1[one].sign*SD2[two].sign;
	
	SD_out.push_back(sde);
	
	}
}
	
}


void AmiGraph::get_coll(std::vector< std::vector<int>> &cv, int val){
cv.clear();

int outer=val;
if(outer %2==0){exit;} 
// get the set of integers that sum to outer 
// ex for 3, need 2 and 1
// ex for 5 need 2, 2, 1 - 3, 2 - 4, 1 

// need set of translations. 2, 2 to 5 could be 1,2 and 3,4 or 3,4 and 1,2 (ignore this for now)


std::vector<int> set;

for(int odd=1; odd< outer; odd=odd+2){
set.push_back(odd);
int sum=odd;
int start=val-1;
do{
if(sum+start<=outer	){
sum+=start;
set.push_back(start);	
}else{start=start-2;}
	

}while(sum!=outer);

cv.push_back(set);
set.clear();


}



}

void AmiGraph::construct_sigma_symM(symM &SM, int min, int max){
	
SM.resize(max+1);	
	
for (int i=min ; i<max+1; i++){

int length=2*i+1;

symDET SD;
make_symDET(length,SD);

SM[i]=SD;


}	
} 



void AmiGraph::construct_sigma_HF_symM(symM &SM, int min, int max){
	
SM.resize(max+1);	
	
for (int i=min ; i<max+1; i++){

std::cout<<"On ord i="<<i<<std::endl;

int length=2*i+1;

symDET SD;
make_symDET_sigma_HF(length,SD);

SM[i]=SD;


}	
}



void AmiGraph::make_symDET_sigma_HF(int L, symDET &SD){
SD.clear();

std::vector<int> v;

for( int i=0; i< L; i++){
	
v.push_back(i);	
}


int pos=0;
int neg=0;
int count=0;


do{
	count++;
	
symDET_element se;

std::pair<int,int> entry;
int x=0;

bool add;

add=true;

for(int i=0; i< v.size(); i++){	



// std::cout<<v[i]<<" ";	
// this removes tadpoles
if(i==v[i]){add=false;}

// this removes fock insertions
if((i)/2==(v[i])/2){add=false;

// std::cout<<"("<<i<<","<<v[i]<<")"<<std::endl;
// std::cout<<(i)/2<<" "<<(v[i])/2<<std::endl;
}

entry.first=i;
entry.second=v[i];

se.plist.push_back(entry);

for(int j=i+1; j< v.size(); j++){
	
if(v[i]>v[j]){x++;}	
	
}

}

se.sign=std::pow(-1,x);
// std::cout<<"Sign is "<< se.sign<<" vs lexigographical "<< std::pow(-1,count) <<std::endl;

// std::cout<<std::endl;
// std::cout<<"Add is "<<add<<std::endl;

graph_t g;
g=return_graph(se);

// std::cout<<"Checking connected "<<std::endl;
if(!is_connected_graph(g)){
	add=false;
// continue;	
}	
// std::cout<<"Add is currently "<<add<<std::endl;

// std::cout<<"Checking 1P R  "<<std::endl;
if(add){
if(is_1P_reducible_graph(g)){
	add=false;
	// continue;
	
}
}

// std::cout<<"Add is currently "<<add<<std::endl;



if(add){
// se.graph=g;
// number_vertices(se.graph);	
	
SD.push_back(se);
}
se.plist.clear();

if(se.sign==1){pos++;}
if(se.sign==-1){neg++;}



	
	}while(std::next_permutation(v.begin(),v.end()));



std::cout<<"Pos and neg were "<< pos<<" "<<neg<<std::endl;


}	



void AmiGraph::construct_HF_symM(symM &SM, int min, int max){
	
SM.resize(max+1);	
	
for (int i=min ; i<max+1; i++){

int length=2*i+1;

symDET SD;
make_symDET_HF(length,SD);

SM[i]=SD;


}	
}

void AmiGraph::make_symDET_HF(int L, symDET &SD){
SD.clear();

std::vector<int> v;

for( int i=0; i< L; i++){
	
v.push_back(i);	
}


int pos=0;
int neg=0;
int count=0;


do{
	count++;
	
symDET_element se;

std::pair<int,int> entry;
int x=0;

bool add;

add=true;

for(int i=0; i< v.size(); i++){	



// std::cout<<v[i]<<" ";	
// this removes tadpoles
if(i==v[i]){add=false;}

// this removes fock insertions
if((i)/2==(v[i])/2){add=false;

// std::cout<<"("<<i<<","<<v[i]<<")"<<std::endl;
// std::cout<<(i)/2<<" "<<(v[i])/2<<std::endl;
}

entry.first=i;
entry.second=v[i];

se.plist.push_back(entry);

for(int j=i+1; j< v.size(); j++){
	
if(v[i]>v[j]){x++;}	
	
}

}

se.sign=std::pow(-1,x);
// std::cout<<"Sign is "<< se.sign<<" vs lexigographical "<< std::pow(-1,count) <<std::endl;

// std::cout<<std::endl;
// std::cout<<"Add is "<<add<<std::endl;
if(add){
SD.push_back(se);
}
se.plist.clear();

if(se.sign==1){pos++;}
if(se.sign==-1){neg++;}



	
	}while(std::next_permutation(v.begin(),v.end()));



std::cout<<"Pos and neg were "<< pos<<" "<<neg<<std::endl;


}	


void AmiGraph::make_symDET(int L, symDET &SD){
SD.clear();

std::vector<int> v;

for( int i=0; i< L; i++){
	
v.push_back(i);	
}


int pos=0;
int neg=0;


do{
	
symDET_element se;

std::pair<int,int> entry;
int x=0;

for(int i=0; i< v.size(); i++){	



// std::cout<<v[i]<<" ";	

entry.first=i;
entry.second=v[i];

se.plist.push_back(entry);

for(int j=i+1; j< v.size(); j++){
	
if(v[i]>v[j]){x++;}	
	
}

}

se.sign=std::pow(-1,x);

// std::cout<<std::endl;

SD.push_back(se);
se.plist.clear();

if(se.sign==1){pos++;}
if(se.sign==-1){neg++;}



	
	}while(std::next_permutation(v.begin(),v.end()));



std::cout<<"Pos and neg were "<< pos<<" "<<neg<<std::endl;


}	

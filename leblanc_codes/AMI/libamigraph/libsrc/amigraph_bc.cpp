
#include "amigraph.hpp"


// typedef std::vector< solution_set_matrix_t > gg_solution_set_matrix_t;

// typedef std::vector< std::vector<solution_set> > solution_set_matrix_t;

/* {
solution_set(){}
	
g_prod_t R0_;
S_t S_;
P_t P_;
R_t R_;
ami_parms ami_parms_;
double prefactor_;
} */



/* 

void AmiGraph::broadcast_ggamim( NewAmiCalc::gg_solution_set_matrix_t & gg,  NewAmiCalc::gg_solution_set_matrix_t & gg_bc){
	
// int mpir;

// MPI_Comm_rank(MPI_COMM_WORLD, &mpir);	
	

// &matrix_one, MATRIX_SIZE * MATRIX_SIZE	
// MPI_Bcast(data, num_elements, MPI_INT, 0, MPI_COMM_WORLD);	

// std::vector<int> test;
// int nbit=0;
// if(mpi_rank==0){
// test.push_back(2);
// test.push_back(5);

// nbit=test.size();
// MPI_Bcast(&nbit,1,MPI_INT,0, MPI_COMM_WORLD);
// MPI_Bcast(test.data(), nbit,MPI_INT,0, MPI_COMM_WORLD);
// }
// else{
// MPI_Bcast(&nbit,1,MPI_INT,0, MPI_COMM_WORLD);
// test.resize(nbit);
// MPI_Bcast(test.data(),nbit,MPI_INT,0,MPI_COMM_WORLD);	
	
// }
	
	
// std::cout<<"This is rank "<< mpi_rank<<" with size "<< test.size()<<" and entry "<< test[1]<<std::endl;	


// std::vector<int> test;
// int nbit=0;
// if(mpi_rank==0){
	// test.push_back(2);
// test.push_back(4);
// }

// broadcast_stdvec_i(test,0);
	// std::cout<<"This is rank "<< mpi_rank<<" with size "<< test.size()<<" and entry "<< test[1]<<std::endl;



	
	
	
	
} */

// typedef std::vector<g_prod_t> Ri_t;
// typedef std::vector<pole_array_t> Pi_t;
// typedef std::vector<sign_t> Si_t;

// typedef std::vector<g_prod_array_t> R_t;  // TODO: is there a reason this is not a vector of Ri_t ???
// typedef std::vector<Pi_t> P_t;
// typedef std::vector<Si_t> S_t;


void AmiGraph::broadcast_ggm(AmiGraph::gg_matrix_t &ggm, int from){
	
int len=0;

if(mpi_rank==from){
	len=ggm.size();
	MPI_Bcast(&len, 1, MPI_INT, from, MPI_COMM_WORLD);
}else{
	MPI_Bcast(&len, 1, MPI_INT, from, MPI_COMM_WORLD);
	ggm.resize(len);
}
	
for (int i=0; i< len; i++){
	std::cout<<"On "<< mpi_rank<<" order "<< i<<std::endl;
broadcast_ggv(ggm[i],from);
}			
	
}

void AmiGraph::broadcast_ggv(AmiGraph::gg_vec_t &ggv, int from){
	
int len=0;

if(mpi_rank==from){
	len=ggv.size();
	MPI_Bcast(&len, 1, MPI_INT, from, MPI_COMM_WORLD);
}else{
	
	MPI_Bcast(&len, 1, MPI_INT, from, MPI_COMM_WORLD);
	ggv.resize(len);
}
	
for (int i=0; i< len; i++){
	// std::cout<<"On "<< mpi_rank<<" group "<< i<<std::endl;
broadcast_graph_group_pairs(ggv[i],from);
}			
	
}

//this only broadcasts the solutions and pst parts of the gp_vec 
void AmiGraph::broadcast_graph_group_pairs(AmiGraph::graph_group &gg, int from){
	
	// std::cout<<"BC Group pairs "<<std::endl;
// broadcast_graph_vec(gg.graph_vec,from);
broadcast_gp_vec(gg.gp_vec,from);
// std::cout<<"BC prefactor"<<std::endl;
broadcast_stdvec_d(gg.prefactor,from);
// broadcast_labels(gg.labels,from); // typically labels is empty when broadcast 
// std::cout<<"BC group sols"<<std::endl;
broadcast_solution_set_vec(gg.ss_vec,from);	
	
}


void AmiGraph::broadcast_solution_set(NewAmiCalc::solution_set &s, int from){
broadcast_g_prod_t(s.R0_,from);
broadcast_S_t(s.S_,from);
broadcast_R_t(s.R_,from);
broadcast_P_t(s.P_,from);
broadcast_ami_parms(s.ami_parms_,from);
broadcast_bose_alphas(s.bose_alphas_,from);

broadcast_g_prod_t(s.Unique,from);
broadcast_Pi(s.Unique_poles,from);
broadcast_R_ref_t(s.Rref,from);
broadcast_R_ref_t(s.Eval_list,from);


MPI_Bcast(&s.prefactor_,1, MPI_DOUBLE,from, MPI_COMM_WORLD);
MPI_Bcast(&s.loops_,1, MPI_INT,from, MPI_COMM_WORLD);	

MPI_Bcast(&s.ct_count_,1, MPI_DOUBLE,from, MPI_COMM_WORLD);
MPI_Bcast(&s.sigma_ct_count_,1, MPI_INT,from, MPI_COMM_WORLD);	
MPI_Bcast(&s.order_shift,1, MPI_INT,from, MPI_COMM_WORLD);
	
}


void AmiGraph::broadcast_R_ref_t(AmiBase::R_ref_t &Rref, int from){

int len=0;

if(mpi_rank==from){
	len=Rref.size();
	MPI_Bcast(&len, 1, MPI_INT, from, MPI_COMM_WORLD);
}else{
	MPI_Bcast(&len, 1, MPI_INT, from, MPI_COMM_WORLD);
	Rref.resize(len);
}

for (int i=0; i< len; i++){
	
broadcast_ref_v_t(Rref[i],from);	
	
}	
}


void AmiGraph::broadcast_ref_v_t(AmiBase::ref_v_t &ref_v, int from){

int len=0;

if(mpi_rank==from){
	len=ref_v.size();
	MPI_Bcast(&len, 1, MPI_INT, from, MPI_COMM_WORLD);
}else{
	MPI_Bcast(&len, 1, MPI_INT, from, MPI_COMM_WORLD);
	ref_v.resize(len);
}

for (int i=0; i< len; i++){
	
broadcast_pair_int_int(ref_v[i],from);	
	
}	
}

void AmiGraph::broadcast_pair_int_int(std::pair<int,int> &ref_t, int from){


MPI_Bcast(&ref_t.first, 1, MPI_INT, from, MPI_COMM_WORLD);
MPI_Bcast(&ref_t.second, 1, MPI_INT, from, MPI_COMM_WORLD);	


}


void AmiGraph::broadcast_gp_vec(std::vector<AmiGraph::git_pair> &gp_vec, int from){
	
	
int len=0;
// std::cout<<"On "<<mpi_rank<<" length "<< gp_vec.size()<<std::endl;

if(mpi_rank==from){
	len=gp_vec.size();
	MPI_Bcast(&len, 1, MPI_INT, from, MPI_COMM_WORLD);
	
}else{
	MPI_Bcast(&len, 1, MPI_INT, from, MPI_COMM_WORLD);
	if(len>0){
	gp_vec.resize(len);
	}
}
	
for (int i=0; i< len; i++){
	// std::cout<<"BC Git pairs"<<std::endl;
broadcast_gp(gp_vec[i],from);
}		
		
}
/* 
git_pair(){}		
	
graph_t g1_;
graph_t g2_;
git_perm_set_t pst_;
NewAmiCalc::solution_set s1_, s2_;	
}; */

void AmiGraph::broadcast_gp(AmiGraph::git_pair &gp, int from){
broadcast_solution_set(gp.s1_,from);
broadcast_solution_set(gp.s2_,from);
broadcast_git_perm_set(gp.pst_,from);	
}


// typedef std::vector<int> git_perm_t;
// typedef std::vector<git_perm_t> git_perm_set_t;
void AmiGraph::broadcast_git_perm_set(AmiGraph::git_perm_set_t &pst, int from){
	
int len=0;

if(mpi_rank==from){
	len=pst.size();
	MPI_Bcast(&len, 1, MPI_INT, from, MPI_COMM_WORLD);
}else{
	MPI_Bcast(&len, 1, MPI_INT, from, MPI_COMM_WORLD);
	pst.resize(len);
}
	
for (int i=0; i< len; i++){
// broadcast_git_perm(pst[i],from);
broadcast_stdvec_i(pst[i],from);
}		
			
	
}

// void AmiGraph::broadcast_git_perm(git_perm_t &perm, int from){
	
// }


/* struct graph_group {
  
  std::vector<graph_t> graph_vec;
  std::vector<git_pair> gp_vec;
  
  std::vector<double> prefactor;
  std::vector<labels_t> labels;
  std::vector<NewAmiCalc::solution_set> ss_vec; // solution set vector 
  
  // git permutation set 
  // git_perm_list_t perm_vec;
     
}; */
/* 
void AmiGraph::broadcast_graph_vec(std::vector<AmiGraph::graph_t> &gv, int from){
	
int len=0;

if(mpi_rank==from){
	len=gv.size();
	MPI_Bcast(&len, 1, MPI_INT, from, MPI_COMM_WORLD);
}else{
	MPI_Bcast(&len, 1, MPI_INT, from, MPI_COMM_WORLD);
	gv.resize(len);
}
	
for (int i=0; i< len; i++){
broadcast_graph(gv[i],from);
}			
	
}
 */
// void AmiGraph::broadcast_graph(AmiGraph::graph_t &g, int from){
// MPI_Bcast(&g,1, MPI_GRAPH, from, MPI_COMM_WORLD);
// } 










void AmiGraph::broadcast_S_t(AmiBase::S_t &s, int from){
int len=0;

if(mpi_rank==from){

len=s.size();
MPI_Bcast(&len, 1, MPI_INT, from, MPI_COMM_WORLD);
}else{
MPI_Bcast(&len, 1, MPI_INT, from, MPI_COMM_WORLD);
s.resize(len);

}	

for(int i=0; i< len; i++){
	
broadcast_Si(s[i],from);	
	
}
	
	
	
}


void AmiGraph::broadcast_P_t(AmiBase::P_t &p, int from){
int len=0;

if(mpi_rank==from){

len=p.size();
MPI_Bcast(&len, 1, MPI_INT, from, MPI_COMM_WORLD);
}else{
MPI_Bcast(&len, 1, MPI_INT, from, MPI_COMM_WORLD);
p.resize(len);

}	

for(int i=0; i< len; i++){
	
broadcast_Pi(p[i],from);	
	
}
	
	
	
}



void AmiGraph::broadcast_R_t(AmiBase::R_t &r, int from){
int len=0;

if(mpi_rank==from){

len=r.size();
MPI_Bcast(&len, 1, MPI_INT, from, MPI_COMM_WORLD);
}else{
MPI_Bcast(&len, 1, MPI_INT, from, MPI_COMM_WORLD);
r.resize(len);

}	

for(int i=0; i< len; i++){
	
broadcast_Ri(r[i],from);	
	
}
	
	
	
}



void AmiGraph::broadcast_Si(AmiBase::Si_t &si, int from){

int len=0;

if(mpi_rank==from){
	len=si.size();
	MPI_Bcast(&len, 1, MPI_INT, from, MPI_COMM_WORLD);
}else{
	MPI_Bcast(&len, 1, MPI_INT, from, MPI_COMM_WORLD);
	si.resize(len);
}

for (int i=0; i< len; i++){
	
broadcast_stdvec_d(si[i],from);	
	
}	
}

void AmiGraph::broadcast_pole_array(AmiBase::pole_array_t &p, int from){

int len=0;
if(mpi_rank==from){
	len=p.size();
	MPI_Bcast(&len, 1, MPI_INT, from, MPI_COMM_WORLD);
}else{
	MPI_Bcast(&len, 1, MPI_INT, from, MPI_COMM_WORLD);
	p.resize(len);
}
for (int i=0; i< len; i++){
	
broadcast_pole_struct(p[i],from);	
	
}		
	
	
}

// void AmiGraph::broadcast_solution_set(NewAmiCalc::solution_set &s, int from){
// broadcast_g_prod_t(s.R0_,from);
// broadcast_S_t(s.S_,from);
// broadcast_R_t(s.R_,from);
// broadcast_P_t(s.P_,from);
// broadcast_ami_parms(s.ami_parms_,from);
// broadcast_bose_alphas(s.bose_alphas_,from);

// MPI_Bcast(&s.prefactor_,1, MPI_DOUBLE,from, MPI_COMM_WORLD);
// MPI_Bcast(&s.loops_,1, MPI_INT,from, MPI_COMM_WORLD);	

// MPI_Bcast(&s.ct_count_,1, MPI_DOUBLE,from, MPI_COMM_WORLD);
// MPI_Bcast(&s.sigma_ct_count_,1, MPI_INT,from, MPI_COMM_WORLD);	
	
// }

void AmiGraph::broadcast_kkpm(kkp_matrix_t &kkpm, int from ){
int len=0;	

if(mpi_rank==from){
len=kkpm.size();
MPI_Bcast(&len,1,MPI_INT,from,MPI_COMM_WORLD);
}else{
MPI_Bcast(&len,1,MPI_INT,from,MPI_COMM_WORLD);
kkpm.resize(len);
}	

for (int i=0; i<len; i++){
broadcast_kkpv(kkpm[i],0);
	
}
	
	
}

void AmiGraph::broadcast_kkpv(kkp_vec_t &kkpv, int from ){
int len=0;	

if(mpi_rank==from){
len=kkpv.size();
MPI_Bcast(&len,1,MPI_INT,from,MPI_COMM_WORLD);
}else{
MPI_Bcast(&len,1,MPI_INT,from,MPI_COMM_WORLD);
kkpv.resize(len);
}	

for (int i=0; i<len; i++){
broadcast_alpha_vec_list(kkpv[i],0);
	
}
	
	
}

void AmiGraph::broadcast_alpha_vec_list(alpha_vec_list_t &avl, int from ){
int len=0;	

if(mpi_rank==from){
len=avl.size();
MPI_Bcast(&len,1,MPI_INT,from,MPI_COMM_WORLD);
}else{
MPI_Bcast(&len,1,MPI_INT,from,MPI_COMM_WORLD);
avl.resize(len);
}	

for (int i=0; i<len; i++){
broadcast_alpha_vec(avl[i],0);
	
}
	
	
}

void AmiGraph::broadcast_alpha_vec(alpha_vec_t &av, int from ){
int len=0;	

if(mpi_rank==from){
len=av.size();
MPI_Bcast(&len,1,MPI_INT,from,MPI_COMM_WORLD);
}else{
MPI_Bcast(&len,1,MPI_INT,from,MPI_COMM_WORLD);
av.resize(len);
}	

for (int i=0; i<len; i++){
broadcast_stdvec_i(av[i],0);
	
}
	
	
}



void AmiGraph::broadcast_gg_solution_set_matrix_t(NewAmiCalc::gg_solution_set_matrix_t &s, int from){
int len=0;

if(mpi_rank==from){
	len=s.size();
	MPI_Bcast(&len, 1, MPI_INT, from, MPI_COMM_WORLD);
}else{
	MPI_Bcast(&len, 1, MPI_INT, from, MPI_COMM_WORLD);
	s.resize(len);
}
	
for (int i=0; i< len; i++){
broadcast_solution_set_matrix_t(s[i],from);
}		
}

void AmiGraph::broadcast_solution_set_matrix_t(NewAmiCalc::solution_set_matrix_t &s, int from){
int len=0;

if(mpi_rank==from){
	len=s.size();
	MPI_Bcast(&len, 1, MPI_INT, from, MPI_COMM_WORLD);
}else{
	MPI_Bcast(&len, 1, MPI_INT, from, MPI_COMM_WORLD);
	s.resize(len);
}
	
for (int i=0; i< len; i++){
broadcast_solution_set_vec(s[i],from);
}		
}

void AmiGraph::broadcast_solution_set_vec(NewAmiCalc::solution_set_vec_t &s, int from){
int len=0;

if(mpi_rank==from){
	len=s.size();
	MPI_Bcast(&len, 1, MPI_INT, from, MPI_COMM_WORLD);
}else{
	MPI_Bcast(&len, 1, MPI_INT, from, MPI_COMM_WORLD);
	s.resize(len);
}
	
for (int i=0; i< len; i++){
broadcast_solution_set(s[i],from);
}	
	
	
}


void AmiGraph::broadcast_Pi(AmiBase::Pi_t &p, int from){

int len=0;

if(mpi_rank==from){
	len=p.size();
	MPI_Bcast(&len, 1, MPI_INT, from, MPI_COMM_WORLD);
}else{
	MPI_Bcast(&len, 1, MPI_INT, from, MPI_COMM_WORLD);
	p.resize(len);
}

for (int i=0; i< len; i++){
	
broadcast_pole_array(p[i],from);	
	
}	
}

void AmiGraph::broadcast_Ri(AmiBase::Ri_t &ri, int from){

int len=0;

if(mpi_rank==from){
	len=ri.size();
	MPI_Bcast(&len, 1, MPI_INT, from, MPI_COMM_WORLD);
}else{
	MPI_Bcast(&len, 1, MPI_INT, from, MPI_COMM_WORLD);
	ri.resize(len);
}

for (int i=0; i< len; i++){
	
broadcast_g_prod_t(ri[i],from);	
	
}	
}



void AmiGraph::broadcast_state(NewAmiCalc::internal_state &s, int from){
MPI_Bcast(&s.dim_, 1, MPI_INT, from, MPI_COMM_WORLD);
MPI_Bcast(&s.order_, 1, MPI_INT, from, MPI_COMM_WORLD);	

broadcast_k_list(s.internal_k_list_,  from);	
broadcast_hopping_list(s.t_list_,from);
broadcast_hopping_list(s.tp_list_,from);
	
	
}

void AmiGraph::broadcast_ext_vars(NewAmiCalc::ext_vars &e, int from){
	
broadcast_k_list(e.external_k_list_, from);
MPI_Bcast(&e.KDIM_, 1, MPI_INT, from, MPI_COMM_WORLD);	
broadcast_stdvec_cd(e.external_freq_, from);
MPI_Bcast(&e.BETA_,1, MPI_DOUBLE, from, MPI_COMM_WORLD);
MPI_Bcast(&e.MU_,1, MPI_DOUBLE_COMPLEX, from, MPI_COMM_WORLD);	
	
}

/* 



k_vector_t external_k_vector_;	
int KDIM_;
frequency_t external_freq_;

double BETA_;
std::complex<double> MU_;
}; */

// void AmiGraph::broadcast_len(int &len, int from){

// if(mpi_rank==from){

// MPI_Bcast(&len,1,MPI_INT,from, MPI_COMM_WORLD);
// }else{
// MPI_Bcast(&len,1,MPI_INT,from, MPI_COMM_WORLD);
	
	
	
	
// }

// typedef std::vector< double> k_vector_t;
// typedef std::vector< k_vector_t> k_vect_list_t;

void AmiGraph::broadcast_hopping_list(NewAmiCalc::hopping_list_t &t, int from){
	
int nbit=0;
    if(mpi_rank==from){
        nbit=t.size();
        MPI_Bcast(&t, 1, MPI_INT, from, MPI_COMM_WORLD);
        MPI_Bcast(t.data(), nbit, MPI_INT, from, MPI_COMM_WORLD);
    }else{

		MPI_Bcast(&t, 1, MPI_INT, from, MPI_COMM_WORLD);
		t.resize(nbit);
        MPI_Bcast(t.data(), nbit, MPI_INT, from, MPI_COMM_WORLD);
        

    }

	
	
}

void AmiGraph::broadcast_k_list(NewAmiCalc::k_vect_list_t &k,  int from){
int len=0;
    if(mpi_rank==from){
        len=k.size();
        MPI_Bcast(&len, 1, MPI_INT, from, MPI_COMM_WORLD);
        // MPI_Bcast(dub.data(), nbit, MPI_DOUBLE, from, MPI_COMM_WORLD);
    }else{

		MPI_Bcast(&len, 1, MPI_INT, from, MPI_COMM_WORLD);
		k.resize(len);
        // MPI_Bcast(dub.data(), nbit, MPI_DOUBLE, from, MPI_COMM_WORLD);
    

    }
	
	
	for	(int i=0; i<len; i++){
	
	broadcast_stdvec_d(k[i],from);
    
    }
}

/* struct internal_state{

internal_state(int k_length, int dim){
internal_k_list_.assign(k_length, std::vector< double>(dim,0.0));	
dim_=dim;
order_=k_length;
}

internal_state(){}

void initialize(int k_length, int dim){
internal_k_list_.assign(k_length, std::vector< double>(dim,0.0));	
dim_=dim;
order_=k_length;
	
}

k_vect_list_t internal_k_list_;
int order_;
int dim_;
//double prefactor_;
// R_t R_array_;
// P_t P_array_;
// S_t S_array_;

}; */





// struct ami_parms{
// ami_parms(int N_INT,  double E_REG){
// N_INT_=N_INT;
// E_REG_=E_REG;
// TYPE_=static_cast<NewAmiCalc::graph_type>(0); /// by default sigma
// }


void AmiGraph::broadcast_ami_parms(AmiBase::ami_parms &p, int from){
MPI_Bcast(&p.N_INT_, 1, MPI_INT, from, MPI_COMM_WORLD);	
MPI_Bcast(&p.N_EXT_, 1, MPI_INT, from, MPI_COMM_WORLD);	
MPI_Bcast(&p.E_REG_, 1, MPI_DOUBLE, from, MPI_COMM_WORLD);	
MPI_Bcast(&p.TYPE_, 1, MPI_INT, from, MPI_COMM_WORLD);	
}

/* 
typedef std::vector<double> sign_t;
typedef std::vector<g_struct> g_prod_t;
typedef std::vector<pole_struct> pole_array_t;

typedef std::vector<sign_t> sign_array_t;
typedef std::vector<g_prod_t> g_prod_array_t;
 */

void AmiGraph::broadcast_g_prod_t(AmiBase::g_prod_t &gprod, int from){
	
	

int nbit=0;
    if(mpi_rank==from){
        nbit=gprod.size();
        MPI_Bcast(&nbit, 1, MPI_INT, from, MPI_COMM_WORLD);
        //MPI_Bcast(dub.data(), nbit, MPI_DOUBLE, from, MPI_COMM_WORLD);
    }else{

		MPI_Bcast(&nbit, 1, MPI_INT, from, MPI_COMM_WORLD);
		gprod.resize(nbit);
       // MPI_Bcast(dub.data(), nbit, MPI_DOUBLE, from, MPI_COMM_WORLD);
    }

for(int i=0; i< gprod.size(); i++){
broadcast_g_struct(gprod[i], from);
}	
		
	
}



/* 
struct pole_struct{


epsilon_t eps_;
alpha_t alpha_;
int index_;
int multiplicity_=1;
std::vector<int> which_g_;

}; */



void AmiGraph::broadcast_pole_struct( AmiBase::pole_struct &p, int from){

broadcast_stdvec_i(p.eps_,from);
broadcast_stdvec_i(p.alpha_,from);	
broadcast_stdvec_i(p.which_g_,from);	
MPI_Bcast(&p.index_, 1, MPI_INT, from, MPI_COMM_WORLD);	
MPI_Bcast(&p.multiplicity_, 1, MPI_INT, from, MPI_COMM_WORLD);	
MPI_Bcast(&p.der_, 1, MPI_INT, from, MPI_COMM_WORLD);	

}


/* 
struct g_struct{
g_struct(epsilon_t eps, alpha_t alpha, stat_type stat){
eps_=eps;
alpha_=alpha;
stat_=stat;

} */
// epsilon_t and alpha_t are std::vector int, stat is int 

void AmiGraph::broadcast_bose_alphas( std::vector<AmiBase::alpha_t> &bose, int from){
	
	int len=0;
    if(mpi_rank==from){
        len=bose.size();
        MPI_Bcast(&len, 1, MPI_INT, from, MPI_COMM_WORLD);
        // MPI_Bcast(dub.data(), nbit, MPI_DOUBLE, from, MPI_COMM_WORLD);
    }else{

		MPI_Bcast(&len, 1, MPI_INT, from, MPI_COMM_WORLD);
		bose.resize(len);
        // MPI_Bcast(dub.data(), nbit, MPI_DOUBLE, from, MPI_COMM_WORLD);
    

    }
	
	
	for	(int i=0; i<len; i++){
	
	broadcast_stdvec_i(bose[i],from);
    
    }
	
	
}

void AmiGraph::broadcast_g_struct( AmiBase::g_struct &g, int from){

broadcast_stdvec_i(g.eps_,from);
broadcast_stdvec_i(g.alpha_,from);	
// broadcast_stdvec_i(g.eps_indices_,from);
MPI_Bcast(&g.stat_, 1, MPI_INT, from, MPI_COMM_WORLD);
MPI_Bcast(&g.species_, 1, MPI_INT, from, MPI_COMM_WORLD);		
	
}

void AmiGraph::broadcast_stdvec_d(std::vector<double> &dub, int from){


// MPI_Bcast(v.data(), v.size(), MPI_INT, 0, MPI_COMM_WORLD);


int nbit=0;
    if(mpi_rank==from){
        nbit=dub.size();
        MPI_Bcast(&nbit, 1, MPI_INT, from, MPI_COMM_WORLD);
        MPI_Bcast(dub.data(), nbit, MPI_DOUBLE, from, MPI_COMM_WORLD);
    }else{

		MPI_Bcast(&nbit, 1, MPI_INT, from, MPI_COMM_WORLD);
		dub.resize(nbit);
        MPI_Bcast(dub.data(), nbit, MPI_DOUBLE, from, MPI_COMM_WORLD);
        

    }

	
	
}

void AmiGraph::broadcast_stdvec_i(std::vector<int> &dub, int from){

int nbit=0;
    if(mpi_rank==from){
        nbit=dub.size();
        MPI_Bcast(&nbit, 1, MPI_INT, from, MPI_COMM_WORLD);
        MPI_Bcast(dub.data(), nbit, MPI_INT, from, MPI_COMM_WORLD);
    }else{

		MPI_Bcast(&nbit, 1, MPI_INT, from, MPI_COMM_WORLD);
		dub.resize(nbit);
        MPI_Bcast(dub.data(), nbit, MPI_INT, from, MPI_COMM_WORLD);
        

    }
	
	
}

void AmiGraph::broadcast_stdvec_cd(std::vector<std::complex<double>> &cdub, int from){
	
	int nbit=0;
    if(mpi_rank==from){
        nbit=cdub.size();
        MPI_Bcast(&nbit, 1, MPI_INT, from, MPI_COMM_WORLD);
        MPI_Bcast(cdub.data(), nbit, MPI_DOUBLE_COMPLEX, from, MPI_COMM_WORLD);
    }else{

		MPI_Bcast(&nbit, 1, MPI_INT, from, MPI_COMM_WORLD);
		cdub.resize(nbit);
        MPI_Bcast(cdub.data(), nbit, MPI_DOUBLE_COMPLEX, from, MPI_COMM_WORLD);
        

    }
	
	
	
}






/* 

 //now reading the file
    int nbitemread=0;
    float* buffer;
    if(rank==0){
        ifstream file ("example.txt",  ios::in |ios::binary);
        file.read ((char*)&nbitemread, sizeof(int));
        buffer=new float[nbitemread];
        file.read ((char*)buffer,nbitemread* sizeof(float));
        file.close();
        //communication
        MPI_Bcast(&nbitemread, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(buffer, nbitemread, MPI_FLOAT, 0, MPI_COMM_WORLD);
    }else{

        MPI_Bcast(&nbitemread, 1, MPI_INT, 0, MPI_COMM_WORLD);
        //nbitemread is meaningfull now
        buffer=new float[nbitemread];
        MPI_Bcast(buffer, nbitemread, MPI_FLOAT, 0, MPI_COMM_WORLD);

    } */






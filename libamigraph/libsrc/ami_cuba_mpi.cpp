//=======================================================================
// Copyright 2018 JPF LeBlanc
//=======================================================================
#include "amigraph.hpp"
#include <ctime>
#include <Eigen/Dense>
#include "cuba.h"
#include "mpi.h"

static int ami_integrand(const int *ndim, const double xx[],
  const int *ncomp, double ff[], void *userdata);

AmiGraph graph;
int global_integral=0;
// in amical simple_residue a lot of time is spent allocating memory? worth fixing?


// NOte, 3rd order translation is ( 4,2,1,3,5,6)

int main(int argc, char *argv[])
{
  
int mpi_rank;

int comm_size;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);	
  
  
  int max=2;
  int min=0;
  int maxeval=10000;
  int intseed=0;
  
  
  if(argc==1) 
        printf("\nNo Extra Command Line Argument Passed Other Than Program Name"); 
    if(argc==2) 
    { 
        max=atoi(argv[1]);
    }  
	if(argc==3) 
    { 
        min=atoi(argv[1]);
		max=atoi(argv[2]);
    } 
	if(argc==4) 
    { 
        min=atoi(argv[1]);
		max=atoi(argv[2]);
		maxeval=atoi(argv[3]);
				
    } 
	if(argc==5) 
    { 
        min=atoi(argv[1]);
		max=atoi(argv[2]);
		maxeval=atoi(argv[3]);
		intseed=atoi(argv[4]);
		
    } 
	if(argc==6) 
    { 
        min=atoi(argv[1]);
		max=atoi(argv[2]);
		maxeval=atoi(argv[3]);
		intseed=atoi(argv[4]);
		global_integral=atoi(argv[5]);
    } 


std::cout<<min<<" "<<max<<" "<<maxeval<<std::endl;
  
  
//int seed=0; // but this can be set to anything
auto now =std::chrono::high_resolution_clock::now();
auto seed = std::chrono::duration_cast<std::chrono::microseconds>(now.time_since_epoch()).count();
//int seed=0;

std::cout<< seed <<" " <<std::endl;
AmiGraph g(AmiBase::Sigma, seed);//g(AmiCalc::Pi,seed);//g(AmiCalc::Hartree,seed);



AmiGraph::vertex_vector_list_t bubble_vertex_list;
AmiGraph::edge_vector_list_t bubble_edge_list;
AmiGraph::edge_vector_list_t legs_edge_list;
g.bubble_finder(g.current_graph, bubble_vertex_list,  bubble_edge_list, legs_edge_list);

AmiGraph::edge_t one, two;




	
		std::cout<<"Reading external parameters from ext_vars.dat"<<std::endl;
		NewAmiCalc::external_variable_list extern_list;

		std::string infile("ext_vars.dat");
		g.ami.read_external(infile, extern_list);	
		std::cout<<"external parameters read on rank "<<mpi_rank<<" are"<<std::endl;
for(int i=0; i<extern_list.size();i++){
		std::cout<<extern_list[i].BETA_<<" "<<extern_list[i].MU_<<" "<< extern_list[i].H_<<" "<<extern_list[i].KDIM_<<" "<<extern_list[i].external_k_list_[0][0]<<" "<<extern_list[i].external_k_list_[0][1] <<" "<<extern_list[i].external_freq_[0]<<std::endl;
		}	
	
std::vector< std::vector< AmiGraph::graph_t>> graph_matrix;
// g.generate_bubble_graphs(graph_matrix,max,mpi_rank);

// bool success;


// g.label_graphs(graph_matrix,1,max);
// // g.min_ext_label_counts(graph_matrix);
// g.reduce_gm_rf(graph_matrix,mpi_rank);
// g.reduce_gm_tp(graph_matrix, mpi_rank, 1);
// //g.reduce_gm_oneleg(graph_matrix, mpi_rank);
// g.reduce_gm_1PBose(graph_matrix, mpi_rank,1);


// g.reduce_gm_ladder(graph_matrix, mpi_rank,1);

AmiGraph::gg_matrix_t ggm;
// ggm.resize(graph_matrix.size());
g.read_ggmp("../../../sigma_graphs/ggm_hf/",ggm, max);
std::cout<<"Completed read"<<std::endl;


// g.label_graphs_sys(graph_matrix,0,max);
// g.gm_to_ggm(graph_matrix, ggm);
// g.mpi_print_ggm(ggm,mpi_rank);



g.mpi_print_ggm(ggm, mpi_rank);

g.ggm_label(ggm, 0);


for(int m=0; m< ggm.size();m++){
for(int i=0; i< ggm[m].size();i++){
for(int j=0; j< ggm[m][i].graph_vec.size(); j++){
	std::cout<<"Printing graph "<<m<<" "<<i<<" "<<j<<std::endl;
g.print_all_edge_info(ggm[m][i].graph_vec[j]);
}
}}


g.ggm_construct_ami_sol(ggm, 0, mpi_rank);

NewAmiCalc::solution_set_matrix_t AMI_MATRIX;
g.ggm_to_amim(ggm, AMI_MATRIX);

std::cout<<"Completed graph content on rank "<< mpi_rank<<" of "<< comm_size<<std::endl;

std::vector< std::vector< std::vector< std::vector<double> > > > result_matrix, g_result_matrix;
// resize vectors
double last, g_last;

result_matrix.resize(max+1);
g_result_matrix.resize(max+1);
for(int ord=0; ord<result_matrix.size(); ord++){
	// std::cout<<"AMI_MATRIX of ord and size "<<ord<<" "<<AMI_MATRIX[ord].size() <<std::endl;
result_matrix[ord].resize(AMI_MATRIX[ord].size());
g_result_matrix[ord].resize(AMI_MATRIX[ord].size());
	for(int i=0; i< result_matrix[ord].size(); i++){
	result_matrix[ord][i].resize(extern_list.size());
	g_result_matrix[ord][i].resize(extern_list.size());
	for (int j=0; j< result_matrix[ord][i].size(); j++){
	result_matrix[ord][i][j].resize(4,0);	
	g_result_matrix[ord][i][j].resize(4,0);	
	}
		
	}
	
}
// std::cout<<"result_matrix[0].size is "<<result_matrix[0].size()<<std::endl;

std::cout<<"Finished setting up result matrix on rank "<< mpi_rank<<" of "<< comm_size<<std::endl;



// create a linearized index vector
std::vector<int> index(4);
std::vector< std::vector<int>> index_vec;
int count=0;

for (int ord=min; ord< max+1; ord++){
int ami_size=AMI_MATRIX[ord].size();
int extern_size=extern_list.size();
	for (int i=0; i< ami_size; i++){
	// #pragma omp parallel for 
		for (int j=0; j<extern_size; j++){
			// std::cout<<" ord is "<<ord<<std::endl;
		index[0]=ord;
		index[1]=i;
		index[2]=j;
		index[3]=count;
		count++;
		index_vec.push_back(index);
			
			
		}
	}
}

std::cout<<"Starting loop on rank "<< mpi_rank<<std::endl;


for(int i=0; i<index_vec.size(); i++){
			
if (i%comm_size != mpi_rank) continue;			

		//printf("ord = %d, i = %d, j= %d, threadId = %d \n",ord, i, j, omp_get_thread_num());
int ord=index_vec[i][0];
int num=index_vec[i][1];
int extvar=index_vec[i][2];


	std::cout<<"This is rank "<< mpi_rank <<" working on ord="<<ord<<", num="<< num <<", and extvar="<<extvar<<std::endl;
		// std::cout<<ord<<" "<< i <<" "<< j <<std::endl;

		int nsamples=maxeval;
		result_matrix[ord][num][extvar]=g.integrator(AMI_MATRIX[ord][num], extern_list[extvar], nsamples, intseed);
		
		// last=result_matrix[ord][num][extvar][2];
		
		std::cout<<result_matrix[ord][num][extvar][0]<<" "<<result_matrix[ord][num][extvar][1] <<" "<<result_matrix[ord][num][extvar][2]<<" "<<result_matrix[ord][num][extvar][3]<<std::endl;
		
		// result_list[j]=g.integrator(AMI_MATRIX[ord][i], extern_list[j], nsamples);
		// result_list.push_back(g.integrator(AMI_MATRIX[ord][i], extern_list[j], nsamples));

		// std::cout<<"Result was "<<std::endl;
		// std::cout<<"
		// std::cout<<"Re: "<< ord_list[i][j][0]<<" +/- "<<ord_list[i][j][1] <<" and Im: "<<ord_list[i][j][2]<<" +/- "<<ord_list[i][j][3]<<std::endl;

		
} 


for (int ord=min; ord< max+1; ord++){
	for(int i=0; i< g_result_matrix[ord].size(); i++){
		for(int j=0; j< g_result_matrix[ord][i].size(); j++){
			// MPI_Reduce(&last, &g_last, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			for(int m=0; m<4; m++){
		MPI_Reduce(&result_matrix[ord][i][j][m], &g_result_matrix[ord][i][j][m], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			}
		}
	}
}


MPI_Finalize();


if(mpi_rank==0){
	std::cout<<"This is master rank "<<mpi_rank<<" finalizing result."<<std::endl;
std::ofstream file;
file.open("outfile.dat");
file<<"# order num ext_index ext_freq Im_ext_freq Beta reMu imMu kx ky Re err_Re Im err_Im"<<std::endl;

for (int ord=min; ord< max+1; ord++){
	std::cout<<"Writing order "<<ord<<std::endl;
	for(int i=0; i< g_result_matrix[ord].size(); i++){
		// std::cout<<"Writing graph "<<i<<std::endl;
		
		for(int j=0; j< g_result_matrix[ord][i].size(); j++){
			// std::cout<<"Writing extern "<<j<<std::endl;

file<< ord<<" "<<i<<" "<<j<<" ";
file<< extern_list[j].external_freq_[0].real()<<" "<<extern_list[j].external_freq_[0].imag()<<" ";
file<<extern_list[j].BETA_<<" "<< extern_list[j].MU_.real()<<" "<<extern_list[j].MU_.imag()<<" ";
file<<extern_list[j].external_k_list_[0][0]<<" "<<extern_list[j].external_k_list_[0][1]<<" ";
file<< g_result_matrix[ord][i][j][0]<<" "<<g_result_matrix[ord][i][j][1] <<" "<<g_result_matrix[ord][i][j][2]<<" "<<g_result_matrix[ord][i][j][3]<<std::endl;

		}

	}
}
std::cout<<"Closing file "<<std::endl;

file.close();
std::cout<<"File Closed on Master "<<std::endl;
} 
 // exit(0);

return 0;
}


#define NCOMP 2
#define NVEC 1
#define EPSREL 1e-2
#define EPSABS 1e-8
#define VERBOSE 0
#define LAST 4
#define SEED 0
// #define MINEVAL 2000
#define NSTART 100
#define NINCREASE 50
#define NBATCH 100
#define GRIDNO 0

#define STATEFILE NULL
#define SPIN NULL

#define NNEW 500
#define NMIN 2
#define FLATNESS 0.25

#define KEY1 20//47
#define KEY2 4
#define KEY3 1
#define MAXPASS 5
#define BORDER 0.0
#define MAXCHISQ 0.2
#define MINDEVIATION 0.04
#define NGIVEN 0

#define NEXTRA 0

#define KEY 0

std::vector< double > AmiGraph::integrator(NewAmiCalc::solution_set &sol, NewAmiCalc::ext_vars &ext, int nsamples, int seed){
	
std::cout<<"In integrator with "<<ext.external_k_list_[0][0]	<<" "<< ext.external_k_list_[0][1]<<std::endl;
	
	std::vector<double> output;
// std::complex<double> output;	
int NDIM=ext.KDIM_*sol.ami_parms_.N_INT_;
std::cout<<"NDIM is "<< NDIM<<std::endl;
// int NCOMP=2;
//#define USERDATA NULL
// int NVEC=1;
// double EPSREL=1e-4;//5e-5;
// double EPSABS=1e-7;//5e-5;
// int VERBOSE=0;
// int LAST=0;//4
// int SEED=0;//(unsigned int)seedvar;//0
int MINEVAL=1000;//nsamples;
int MAXEVAL=nsamples*NDIM;//200000;

// int NSTART=1000;
// int NINCREASE=500;
// int NBATCH=1000;
// int GRIDNO=-1;
// char* STATEFILE=NULL;

// int NNEW=1000;
// double FLATNESS=1.0/2.0;//25.

// int KEY1=3000;
// int KEY2=1800;
// int KEY3=1;
// int MAXPASS=10;//5
// double BORDER=1e-6;
// double MAXCHISQ=4;
// double MINDEVIATION=.52500;
// int NGIVEN=0;
int LDXGIVEN=NDIM;
// int NEXTRA= 0;

// int KEY=0;  

// #define STATEFILE NULL
// #define SPIN NULL

int comp, nregions, neval, fail;
double integral[NCOMP], error[NCOMP], prob[NCOMP];

NewAmiCalc::evaluation_set eval(sol,ext);

//////////////////////////
// double complex saved;

if(global_integral==0){
  printf("\n------------------- Divonne test -------------------\n");

  // Divonne(NDIM, NCOMP, ami_integrand, (void*)&eval, NVEC,
    // EPSREL, EPSABS, VERBOSE, seed,
    // MINEVAL, MAXEVAL, KEY1, KEY2, KEY3, MAXPASS,
    // BORDER, MAXCHISQ, MINDEVIATION,
    // NGIVEN, LDXGIVEN, NULL, NEXTRA, NULL,
    // STATEFILE, SPIN,
    // &nregions, &neval, &fail, integral, error, prob);

 Divonne(NDIM, NCOMP, ami_integrand, (void*)&eval,
    EPSREL, EPSABS, VERBOSE, seed,
    MINEVAL, MAXEVAL, KEY1, KEY2, KEY3, MAXPASS,
    BORDER, MAXCHISQ, MINDEVIATION,
    NGIVEN, LDXGIVEN, NULL, NEXTRA, NULL,
    &nregions, &neval, &fail, integral, error, prob);
		
	

  printf("DIVONNE RESULT:\tnregions %d\tneval %d\tfail %d\n",  nregions, neval, fail);
  for( comp = 0; comp < NCOMP; ++comp ){
    printf("DIVONNE RESULT:\t%.6e +- %.6e\tp = %.3f\n",      (double)integral[comp], (double)error[comp], (double)prob[comp]);
	  
	  output.push_back((double)integral[comp]);
	  output.push_back((double)error[comp]);
  }
}
  
if(global_integral==1){  
   printf("-------------------- Vegas test --------------------\n");

  // Vegas(NDIM, NCOMP, ami_integrand, (void*)&eval, NVEC,
    // EPSREL, EPSABS, VERBOSE, SEED,
    // MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
    // GRIDNO, STATEFILE, SPIN,
    // &neval, &fail, integral, error, prob);

 Vegas(NDIM, NCOMP, ami_integrand, (void*)&eval, 
    EPSREL, EPSABS, VERBOSE, SEED,
    MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
    GRIDNO, STATEFILE,
    &neval, &fail, integral, error, prob);

  printf("VEGAS RESULT:\tneval %d\tfail %d\n",
    neval, fail);
  for( comp = 0; comp < NCOMP; ++comp ){
    printf("VEGAS RESULT:\t%.6e +- %.6e\tp = %.3f\n",
      (double)integral[comp], (double)error[comp], (double)prob[comp]);
	output.push_back((double)integral[comp]);
	  output.push_back((double)error[comp]);
  }
}

 
 if(global_integral==2){
 printf("\n-------------------- Cuhre test --------------------\n");

  // Cuhre(NDIM, NCOMP,  ami_integrand, (void*)&eval, NVEC,
    // EPSREL, EPSABS, VERBOSE | LAST,
    // MINEVAL, MAXEVAL, KEY,
    // STATEFILE, SPIN,
    // &nregions, &neval, &fail, integral, error, prob);

Cuhre(NDIM, NCOMP,  ami_integrand, (void*)&eval, 
    EPSREL, EPSABS, VERBOSE | LAST,
    MINEVAL, MAXEVAL, KEY,
    &nregions, &neval, &fail, integral, error, prob);



  printf("CUHRE RESULT:\tnregions %d\tneval %d\tfail %d\n",
    nregions, neval, fail);
  for( comp = 0; comp < NCOMP; ++comp ){
    printf("CUHRE RESULT:\t%.6e +- %.6e\tp = %.3f\n",
      (double)integral[comp], (double)error[comp], (double)prob[comp]);


	output.push_back((double)integral[comp]);
	  output.push_back((double)error[comp]);

  }
}

// output=integral[0]+I*integral[1];
//output=std::complex<double>{integral[0],integral[1]};

// printf("Returning output %e + %e I \n", creal(output),cimag(output));

// output.push_back(integral[0]);
// output.push_back(error[0]



return output ;///M_PI;
	
	
}

static int ami_integrand(const int *ndim, const double xx[],
  const int *ncomp, double ff[], void *userdata){
	  
NewAmiCalc ami;	  

double kmin,kmax;
double ktemp;
kmin=0.0;
kmax=2*M_PI;

NewAmiCalc::evaluation_set *eval;	 
NewAmiCalc::solution_set sol;
NewAmiCalc::ext_vars ext ;

eval=(NewAmiCalc::evaluation_set *) userdata;

sol=eval->sol_;
ext=eval->ext_vars_;

NewAmiCalc::k_vector_t k;
NewAmiCalc::k_vect_list_t k_list;

int kdim=ext.KDIM_;

int gsize=sol.R0_.size();

// kmax=ext.maxk_;
// kmin=ext.mink_;
// std::cout<<"Kdim is "<<kdim<<std::endl;


// first, translate the array of xx[i] into kx,ky sets 
// tehn, construct a state 

// new function - takes the array of xx[] and returns the state

int j=0;
do{
for(int i=j; i<kdim+j; i++){
	ktemp=kmin+(kmax-kmin)*xx[i];
	// std::cout<<"ktemp and xx are "<<ktemp<<" "<<xx[i]<<" "<<i<<std::endl;
	k.push_back(ktemp);
	
}
k_list.push_back(k);
k.clear();
j=j+kdim;
}while(j< *ndim);

double Jk=kmax-kmin;

// std::cout<<k_list.size()<<std::endl;
// std::cout<<k_list[0].size()<<std::endl;

// AmiCalc::ami_vars_list vars_list;
// graph.ami.construct_ami_vars_list(AMI_MATRIX[ord][num].R0_, graph.current_state, extern_list, vars_list);

NewAmiCalc::internal_state state(k_list.size(), k_list[0].size());
// reset all hoppings to 1
state.t_list_.clear();
state.t_list_.resize(gsize, 1);	
state.tp_list_.clear();
state.tp_list_.resize(gsize,0);
// std::cout<< "hopping is size "<<state.t_list_.size()<<std::endl;


state.internal_k_list_=k_list;
// state.prefactor_=sol.prefactor_; // TODO: depricated

// graph.print_fullstate(state);

AmiBase::ami_vars vars=ami.construct_ami_vars(sol.R0_,sol.prefactor_, state, ext);



std::complex<double> calc_result=ami.amibase.evaluate(sol.ami_parms_, sol.R_, sol.P_, sol.S_,  vars);

// std::cout<<"Integrand gave "<<calc_result<<std::endl;

if(abs(calc_result.real())>100000 || abs(calc_result.imag())>100000){ calc_result=(0,0);
}

// double factor=std::pow(Jk, *ndim);
// double factor=10000;
if(std::isnan(calc_result.real()) || std::isnan(calc_result.imag())){
	calc_result=(0,0);
}

ff[0]=calc_result.real();//*factor;
ff[1]=calc_result.imag();//*factor;
	  
	  
	return 0;  
  }


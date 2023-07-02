//=======================================================================
// Copyright 2018 JPF LeBlanc
//=======================================================================
#include "amigraph.hpp"
#include <ctime>
#include <Eigen/Dense>
#include "cuba.h"

static int ami_integrand(const int *ndim, const double xx[],
  const int *ncomp, double ff[], void *userdata);

// in amical simple_residue a lot of time is spent allocating memory? worth fixing?


// NOte, 3rd order translation is ( 4,2,1,3,5,6)

int main(int argc, char *argv[])
{
  int max=0;
  int min=2;
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
  
  
//int seed=0; // but this can be set to anything
auto now =std::chrono::high_resolution_clock::now();
auto seed = std::chrono::duration_cast<std::chrono::microseconds>(now.time_since_epoch()).count();
//int seed=0;

std::cout<< seed <<" " <<std::endl;
AmiGraph g(AmiCalc::Bare, seed);//g(AmiCalc::Pi,seed);//g(AmiCalc::Hartree,seed);



AmiGraph::vertex_vector_list_t bubble_vertex_list;
AmiGraph::edge_vector_list_t bubble_edge_list;
AmiGraph::edge_vector_list_t legs_edge_list;
g.bubble_finder(g.current_graph, bubble_vertex_list,  bubble_edge_list, legs_edge_list);

AmiGraph::edge_t one, two;

int mpi_rank=0;

if(max!=0){
	
	
		std::cout<<"Reading external parameters from ext_vars.dat"<<std::endl;
		AmiCalc::external_variable_list extern_list;

		std::string infile("ext_vars.dat");
		g.ami.read_external(infile, extern_list);	
	
	
std::vector< std::vector< AmiGraph::graph_t>> graph_matrix;
g.generate_graphs(graph_matrix,max,mpi_rank);

bool success;

g.label_graphs(graph_matrix,max);
g.min_ext_label_counts(graph_matrix);
g.reduce_gm_rf(graph_matrix,mpi_rank);
g.reduce_gm_tp(graph_matrix, mpi_rank);
g.reduce_gm_oneleg(graph_matrix, mpi_rank);


AmiGraph::gg_matrix_t ggm;
ggm.resize(graph_matrix.size());


g.gm_to_ggm(graph_matrix, ggm);

g.ggm_construct_ami_sol(ggm, 1e-8, mpi_rank);

AmiCalc::solution_set_matrix_t AMI_MATRIX;
g.ggm_to_amim(ggm, AMI_MATRIX);


// std::complex<double> mu(0,0);
// double beta=5;
// int dim=2;
// std::complex<double> freq(0,M_PI/beta);

// AmiCalc::ext_vars ext_vars(dim,beta,mu);
// ext_vars.external_k_vector_[0]=M_PI;
// ext_vars.external_k_vector_[1]=M_PI/3;	
// ext_vars.external_freq_[0]=freq;

// g.integrator(AMI_MATRIX[4][9], extern_list[1], 200000);

std::vector< std::vector< std::vector< std::vector<double> > > > result_matrix;
std::vector< std::vector<  std::vector<double>  > > ord_list;
std::vector< std::vector<double> > result_list;
// resize vectors
result_list.resize(extern_list.size());
for (int i=0; i< result_list.size(); i++){
result_list[i].resize(4);	
}
//

result_matrix.resize(max+1);
for(int ord=1; ord<result_matrix.size(); ord++){
result_matrix[ord].resize(AMI_MATRIX[ord].size());
	for(int i=0; i< result_matrix[ord].size(); i++){
	result_matrix[ord][i].resize(extern_list.size());
	for (int j=0; j< result_matrix[ord][i].size(); j++){
	result_matrix[ord][i][j].resize(4);	
	}
		
	}
	
}





for (int ord=min; ord< max+1; ord++){


	ord_list.resize(AMI_MATRIX[ord].size());
	for(int i=0; i< ord_list.size(); i++){
	ord_list[i].resize(extern_list.size());
	for (int j=0; j< ord_list[i].size(); j++){
	ord_list[i][j].resize(4);	
	}
		
	}
int ami_size=AMI_MATRIX[ord].size();
int extern_size=extern_list.size();
//collapse(2) shared(ord_list,AMI_MATRIX, extern_list)
#pragma omp parallel for default(none) shared(ord_list, AMI_MATRIX, extern_list, g, ord, ami_size, extern_size)  collapse(2)
	for (int i=0; i< ami_size; i++){
	// #pragma omp parallel for 
		for (int j=0; j<extern_size; j++){

		printf("ord = %d, i = %d, j= %d, threadId = %d \n",ord, i, j, omp_get_thread_num());

		// std::cout<<ord<<" "<< i <<" "<< j <<std::endl;

		int nsamples=10000*ord;
		ord_list[i][j]=g.integrator(AMI_MATRIX[ord][i], extern_list[j], nsamples);
		// result_list[j]=g.integrator(AMI_MATRIX[ord][i], extern_list[j], nsamples);
		// result_list.push_back(g.integrator(AMI_MATRIX[ord][i], extern_list[j], nsamples));

		// std::cout<<"Result was "<<std::endl;
		// std::cout<<"Re: "<< ord_list[i][j][0]<<" +/- "<<ord_list[i][j][1] <<" and Im: "<<ord_list[i][j][2]<<" +/- "<<ord_list[i][j][3]<<std::endl;

		

		}
	// ord_list[i][j]=result_list[j];
	// result_list.clear();
	}
	
for (int i=0; i< AMI_MATRIX[ord].size(); i++){
	// #pragma omp parallel for 
		for (int j=0; j< extern_list.size(); j++){
			
	result_matrix[ord][i][j]=ord_list[i][j];
			
}}
	
	
	
// result_matrix[ord]=ord_list;
ord_list.clear();

}



std::ofstream file;
file.open("outfile.dat");
file<<"# order num ext_index ext_freq Im_ext_freq Beta Mu kx ky Re err_Re Im err_Im"<<std::endl;

for (int ord=min; ord< max+1; ord++){
	std::cout<<"Writing order "<<ord<<std::endl;
	for(int i=0; i< result_matrix[ord].size(); i++){
		// std::cout<<"Writing graph "<<i<<std::endl;
		
		for(int j=0; j< result_matrix[ord][i].size(); j++){
			// std::cout<<"Writing extern "<<j<<std::endl;

file<< ord<<" "<<i<<" "<<j<<" ";
file<< extern_list[j].external_freq_[0].real()<<" "<<extern_list[j].external_freq_[0].imag()<<" ";
file<<extern_list[j].BETA_<<" "<< extern_list[j].MU_.real()<<" ";
file<<extern_list[j].external_k_vector_[0]<<" "<<extern_list[j].external_k_vector_[1]<<" ";
file<< result_matrix[ord][i][j][0]<<" "<<result_matrix[ord][i][j][1] <<"  "<<result_matrix[ord][i][j][2]<<"  "<<result_matrix[ord][i][j][3]<<std::endl;




		}

	}
}

file.close();

}






    return 0;
}


#define NCOMP 2
#define NVEC 1
#define EPSREL 1e-4
#define EPSABS 1e-8
#define VERBOSE 0
#define LAST 0
#define SEED 0
#define MINEVAL 0
#define NSTART 1000
#define NINCREASE 500
#define NBATCH 1000
#define GRIDNO -1

#define STATEFILE NULL
#define SPIN NULL

#define NNEW 1000
#define NMIN 2
#define FLATNESS 0.5

#define KEY1 3000
#define KEY2 1800
#define KEY3 1
#define MAXPASS 10
#define BORDER 1e-6
#define MAXCHISQ 4.
#define MINDEVIATION .5
#define NGIVEN 0

#define NEXTRA 0

#define KEY 0

std::vector< double > AmiGraph::integrator(AmiCalc::solution_set &sol, AmiCalc::ext_vars &ext, int nsamples){
	
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
// int MINEVAL=000;
int MAXEVAL=nsamples;//200000;

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

AmiCalc::evaluation_set eval(sol,ext);

//////////////////////////
// double complex saved;

  printf("\n------------------- Divonne test -------------------\n");

  Divonne(NDIM, NCOMP, ami_integrand, (void*)&eval, NVEC,
    EPSREL, EPSABS, VERBOSE, SEED,
    MINEVAL, MAXEVAL, KEY1, KEY2, KEY3, MAXPASS,
    BORDER, MAXCHISQ, MINDEVIATION,
    NGIVEN, LDXGIVEN, NULL, NEXTRA, NULL,
    STATEFILE, SPIN,
    &nregions, &neval, &fail, integral, error, prob);

  // printf("DIVONNE RESULT:\tnregions %d\tneval %d\tfail %d\n",  nregions, neval, fail);
  for( comp = 0; comp < NCOMP; ++comp ){
    // printf("DIVONNE RESULT:\t%.8f +- %.8f\tp = %.3f\n",      (double)integral[comp], (double)error[comp], (double)prob[comp]);
	  
	  output.push_back((double)integral[comp]);
	  output.push_back((double)error[comp]);
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
	  
AmiCalc ami;	  

double kmin,kmax;
double ktemp;
kmin=0.0;
kmax=2*M_PI;

AmiCalc::evaluation_set *eval;	 
AmiCalc::solution_set sol;
AmiCalc::ext_vars ext ;

eval=(AmiCalc::evaluation_set *) userdata;

sol=eval->sol_;
ext=eval->ext_vars_;

AmiCalc::k_vector_t k;
AmiCalc::k_vect_list_t k_list;

int kdim=ext.KDIM_;
// std::cout<<"Kdim is "<<kdim<<std::endl;


// first, translate the array of xx[i] into kx,ky sets 
// tehn, construct a state 

// new function - takes the array of xx[] and returns the state

int j=0;
do{
for(int i=j; i<kdim+j; i++){
	ktemp=kmin+(kmax-kmin)*xx[i];
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

AmiCalc::internal_state state(k_list.size(), k_list[0].size());
state.internal_k_list_=k_list;
state.prefactor_=sol.prefactor_;

AmiCalc::ami_vars vars=ami.construct_ami_vars(sol.R0_, state, ext);

std::complex<double> calc_result=ami.evaluate(sol.ami_parms_, sol.R_, sol.P_, sol.S_,  vars);

// std::cout<<"Integrand gave "<<calc_result<<std::endl;

if(abs(calc_result.real())>1 || abs(calc_result.imag())>1){ calc_result=(0,0);
}

// double factor=std::pow(Jk, *ndim);
// double factor=10000;
ff[0]=calc_result.real();//*factor;
ff[1]=calc_result.imag();//*factor;
	  
	  
	return 0;  
  }


//=======================================================================
// Copyright 2018 JPF LeBlanc
//=======================================================================
#include "amigraph.hpp"


int main(int argc, char *argv[])
{

int min=1;
int max=1;
double beta=1;	
	
auto now =std::chrono::high_resolution_clock::now();
auto seed = std::chrono::duration_cast<std::chrono::microseconds>(now.time_since_epoch()).count();	
	
if(argc==3) 
    { 
        min=atoi(argv[1]);
		max=atoi(argv[2]);
		
    }
if(argc==4){
min=atoi(argv[1]);
		max=atoi(argv[2]);
		beta=atof(argv[3]);

}  

int bands=2;//4;
int ext_band=0;
int sc_it=1;

if(argc==6){
min=atoi(argv[1]);
		max=atoi(argv[2]);
		beta=atof(argv[3]);
		bands=atoi(argv[4]);
		sc_it=atoi(argv[5]);
		// ext_band=atoi(argv[5]);


}	
	

// seed=37832;	
AmiGraph g(AmiBase::Sigma, seed);	

g.ami.not_molecule=0;
int isnan_count=0;
	
AmiGraph::symM SM;


auto first = std::chrono::steady_clock::now();	
// g.construct_sigma_symM(SM, min,max);

// g.construct_symM_full(SM,min,max);

g.construct_HF_symM(SM, min,max);



auto second  = std::chrono::steady_clock::now();

std::chrono::duration<double> elapsed = second-first;
    std::cout << "Time to Run symDET: " << elapsed.count() << "s\n"<<std::endl;






// for (int ord=0; ord< SM.size(); ord++){
// for(int i=0; i< SM[ord].size(); i++){

// std::cout<<"Possible graph with sign "<<SM[ord][i].sign<<" i="<<i<<std::endl;	
// for(int j=0; j< SM[ord][i].plist.size(); j++){
// std::cout<<"("<<SM[ord][i].plist[j].first<<","<<SM[ord][i].plist[j].second<<")"<<std::endl;
// }

// }}
// return 0;

// std::vector< std::vector< int>> M;

// std::cout<<"Extracted M matrix is"<<std::endl;
// g.extract_M(SM[ord][i],M);
// g.print_matrix(M); 

// std::cout<<"Using eigen"<<std::endl;

// EigenM_PI/beta*(2*n+1);MatrixXi A,B,C;

// g.extract_ABC(SM[ord][i], A,B,C);

// std::cout<<"Attempt to triangularize"<<std::endl;

// std::vector< std::vector< int>> T;

// triangularize_matrix(M,T);

// g.print_matrix(T);




// }	
// }
// return 0;




auto start = std::chrono::steady_clock::now();
g.symM_make_graphs(SM);	

auto filt=std::chrono::steady_clock::now();


std::cout<<"Prefiltering the last order is size "<< SM[max].size()<<std::endl;
g.symM_filter(SM);
std::cout<<"Postfiltering the remaining graphs are "<<std::endl;
for(int i=0; i< SM.size(); i++){
	std::cout<<"Ord "<<i<<": "<< SM[i].size()<<std::endl;
}

std::cout<<"Labelling graphs "<<std::endl;

auto lab=std::chrono::steady_clock::now();
 
g.symM_label(SM);
 auto end = std::chrono::steady_clock::now();

std::chrono::duration<double> make_time = filt-start;
std::chrono::duration<double> filt_time = lab-filt;
std::chrono::duration<double> lab_time = end-lab;
    std::cout << "Time to make graphs: " << make_time.count() << "s\n"<<std::endl;
	std::cout << "Time to filter graphs: " << filt_time.count() << "s\n"<<std::endl;
	std::cout << "Time to label graphs: " << lab_time.count() << "s\n"<<std::endl;
// g.count_fermi_loops(SM[2][15].graph);

// for(int i=0; i< SM[max].size(); i++){
// std::cout<<"Printing graph  "<<i<<" with loops="<<g.count_fermi_loops(SM[max][i].graph)<<std::endl;
// g.print_all_edge_info(SM[max][i].graph);	
// }

// return 0;

// int mpi_rank=0;

// std::cout<<"Constructing AMI sols"<<std::endl;

 // auto start = std::chrono::steady_clock::now();
  start = std::chrono::steady_clock::now();
   
double ereg=0;//1e-8;
g.symM_construct_ami_sol(SM, ereg);
g.symM_bandmaps(SM);

// TODO: test bandmaps

// std::cout<<"-------- term "<< beta<<"----"<<"--------"<<std::endl;
// g.print_all_edge_info(SM[2][beta].graph);
// std::cout<<std::endl;

// for(int i=0; i< SM[2][beta].bmap.size(); i++){
	
	// std::cout<<SM[2][beta].bmap[i].first<<","<<SM[2][beta].bmap[i].second<<std::endl;
// }

 end = std::chrono::steady_clock::now();
    std::chrono::duration<double>  elapsed_seconds = end-start;
    std::cout << "Time to construct AMI solutions: " << elapsed_seconds.count() << "s\n"<<std::endl;

// return 0;	



int dim=2;
	
// setup H2 
std::cout<<"Reading hii"<<std::endl;
std::string hii_filename("h2sto_e0.txt");	
// std::string hii_filename("H2_cc_pvtz_e0.txt");
// std::string hii_filename("H2_cc_pvdz_e0.txt");
g.ami.read_hii(hii_filename, bands); 

for(int i=0; i< g.ami.global_hii.size(); i++){

std::cout<<i<<" "<<	g.ami.global_hii[i]<<std::endl;
	
}

if(g.ami.global_hii.size()!= bands){ throw std::runtime_error("Not enough bands in hii file");}

std::complex<double> mu=(g.ami.global_hii[0]+g.ami.global_hii[1])/2.0;	

// TODO: This may be a subtle but very important point - that my epsilon function returns the negative of a result, not the energy. 
for(int i=0; i<g.ami.global_hii.size(); i++){

g.ami.global_hii[i]=(g.ami.global_hii[i]-mu);	// minus sign already a part of eval_epsilon in libami - line 632
	
}

std::vector<int> dup_map(g.ami.global_hii.size()); 
			// assign duplicate map 
			for(int i=0; i< dup_map.size(); i++){
				dup_map[i]=i;
			}
			
			for(int i=0; i< dup_map.size(); i++){
				for(int j=i+1; j< dup_map.size(); j++){
					
					if(g.ami.global_hii[i]==g.ami.global_hii[j]){
						dup_map[j]=i;
					}
					
				}
			}

// duplicate initial hii 
std::vector<std::complex<double>> initial_hii;
initial_hii=g.ami.global_hii;




std::cout<<"Reading U"<<std::endl;
std::string UM_filename("h2sto_U.txt");
// std::string UM_filename("H2_cc_pvtz_U.txt");
// std::string UM_filename("H2_cc_pvdz_U.txt");


AmiGraph::U_matrix_type UM;
g.resize_U_matrix(bands,UM);
g.read_U_matrix(UM_filename, UM, bands);
	
std::cout<<"Done"<<std::endl;

for(int a1=0; a1< UM.size(); a1++){
	for(int a2=0; a2< UM[a1].size(); a2++){
	for(int a3=0; a3< UM[a1][a2].size(); a3++){
	for(int a4=0; a4< UM[a1][a2][a3].size(); a4++){
	
	std::cout<< a1<<" "<< a2<<" "<< a3<<" "<< a4<<" "<< UM[a1][a2][a3][a4]<<" "<< UM[a1][a2][a3][a4]-UM[a1][a4][a3][a2]<<std::endl;
	}}}
}

/// try to evaluate something:

//int orde=2;

// AmiCalc::internal_state state(orde,dim);
// g.zero_state(state,orde);

// AmiCalc::ami_vars_list vars_list;
// g.ami.construct_ami_vars_list(GG_AMI_MATRIX[4][0][0].R0_,GG_AMI_MATRIX[4][0][0].prefactor_, state, extern_list, vars_list);
	
	
// g.print_all_edge_info(SM[2][0].graph);	

// g.ami.print_final(5,SM[2][0].ss.R_,SM[2][0].ss.P_,SM[2][0].ss.S_);

// AmiCalc::ext_vars ext;

// std::cout<<"Randomizing avbv"<<std::endl;
// std::vector<int> av,bv;
// g.random_avbv(av,bv, 2, 2, 0); //random_avbv(std::vector<int> &av, std::vector<int> &bv, int order, int num_bands, int ext_index)

// std::vector< std::vector<int>> avlist,bvlist;
// g.construct_avbv_sets(avlist, bvlist,orde, bands,ext_band);
// std::cout<<"avlist and bv list sizes are "<<avlist.size()<<" "<<bvlist.size()<<std::endl;

// TODO: probably don't need to reduce ahead of time. Instead just loop through and check if U==0 first and if it does then skip AMI evaluation parts 
// g.reduce_avbv_sets(UM,avlist,bvlist);// reduce based on actual U values provided only if the UM list contains zeros. 

// std::cout<<"avlist and bv list sizes are "<<avlist.size()<<" "<<bvlist.size()<<std::endl;

// std::cout<<"avone "<<avlist[0][0]<<" "<<avlist[0][1]<<" ";


// for (int comb=0; comb<avlist.size(); comb++){
// std::cout<<comb <<"av : ";
// for (int part=0; part< avlist[comb].size(); part++){
	
	// std::cout<<avlist[comb][part]<<" ";
	
// }
// std::cout<<std::endl;

// std::cout<<comb <<"bv : ";
// for (int part=0; part< avlist[comb].size(); part++){
	
	// std::cout<<bvlist[comb][part]<<" ";
	
// }
// std::cout<<std::endl;

// std::cout<<"U x U is : "<< g.eval_Ulist(UM, avlist[comb], bvlist[comb])<<std::endl;

// }

// int term=beta;
/* 
int can_count=0;
int all_count=0;

for(int comb=0; comb<avlist.size(); comb++){

av=avlist[comb];
bv=bvlist[comb];

for(int term=0; term<SM[2].size(); term++){

std::cout<<"-------- term "<< term<<"----"<<comb<<"--------"<<std::endl;
g.print_all_edge_info(SM[2][term].graph);
std::cout<<std::endl;

std::cout<<"av:";
for (int i=0; i< av.size(); i++){
	std::cout<<av[i]<<" ";
}
std::cout<<std::endl;

std::cout<<"bv:";
for (int i=0; i< bv.size(); i++){
	std::cout<<bv[i]<<" ";
}
std::cout<<std::endl;

// std::cout<<"av:"<< av[0]<<" "<< av[1]<<" "<<av[2]<<" "<<av[3]<<std::endl;
// std::cout<<"av:"<< bv[0]<<" "<< bv[1]<<" "<<bv[2]<<" "<<bv[3]<<std::endl;

// for(int i=0; i<SM[2][term].bmap.size(); i++){
// }


bool can_evaluate=g.is_diagonal(SM[2][term].bmap, av, bv);	
	if(can_evaluate){can_count++;
	std::cout<<"Can evaluate"<<std::endl;
	}else{
		std::cout<<"No evaluate"<<std::endl;
	}
	
	all_count++;
	
}

	
	
}


std::cout<<"All count "<<all_count<<std::endl;
std::cout<<"Can count "<<can_count<<std::endl;
  */

// return 0;
int count=0;
int totalcount=0;

int nmin=0;
int nmax=1024;//3000;//30;
std::ofstream file;
std::string fname="wicks_sc_sto6g_beta"+std::to_string(beta)+"_maxit"+std::to_string(sc_it)+".dat";
// "sto6g_james_"+std::to_string(max)+"_beta"+std::to_string(beta)+"_ext"+std::to_string(ext_band)+"_bc"+std::to_string(bands)+".dat";
file.open(fname,  std::ofstream::out | std::ofstream::trunc); //app);

std::ofstream termfile;
termfile.open("termfile.dat", std::ofstream::out | std::ofstream::trunc);


// Set up self_consistent sigma 
std::vector< std::vector< std::complex<double> >> sigma;
sigma.resize(bands);
std::vector< std::vector< std::complex<double> >> next_sigma;
// next_sigma.resize(bands);
for (int i=0; i< sigma.size(); i++){
	sigma[i].resize(nmax, std::complex<double>((0,0)));
	// next_sigma[i].resize(nmax, std::complex<double>((0,0)));
}

// set up freq and gamma values
std::vector<double> freq(nmax);
std::vector<double> gamma(nmax);

for(int i=0; i< freq.size(); i++){

freq[i]=-4+8.0/nmax*i;
gamma[i]=1e-2;	
	
}


for (int iteration=0; iteration< sc_it; iteration++){
	
// set up next_sigma and zero it 
next_sigma.clear();	

next_sigma.resize(bands);
for (int i=0; i< sigma.size(); i++){
	next_sigma[i].resize(nmax, std::complex<double>((0,0)));
}
	
	
	for (int ext_band=0; ext_band< bands; ext_band++){



for (int ord=min; ord<=max; ord++){

NewAmiCalc::internal_state state(ord,dim);
g.zero_state(state,ord);

NewAmiCalc::ext_vars ext;
ext.BETA_=beta;


auto t1 = std::chrono::steady_clock::now();

// std::cout<<"Randomizing avbv"<<std::endl;
std::vector<int> av,bv;
std::vector< std::vector<int>> avlist,bvlist, dlist;
// g.construct_minimal_avbv_sets(UM,avlist, bvlist,ord, bands,ext_band);
// g.general_construct_avbv_sets(avlist, bvlist,ord, bands,ext_band);
std::cout<<"av and bv-list sizes are "<<avlist.size()<<" "<<bvlist.size()<<std::endl;

g.construct_delta_sets(dlist, ord, bands, ext_band);
std::cout<<"delta-list size is "<<dlist.size()<<std::endl;
// return 0;

auto t2= std::chrono::steady_clock::now();

std::chrono::duration<double> avbv_time = t2-t1;
std::cout << "Time to construct delta sets: " << avbv_time.count() << "s\n"<<std::endl;

// auto t3= std::chrono::steady_clock::now();
// g.reduce_avbv_sets(UM,avlist,bvlist);
// std::cout<<"Reduced av and bv-list sizes are "<<avlist.size()<<" "<<bvlist.size()<<std::endl;
// auto t4= std::chrono::steady_clock::now();
// avbv_time = t4-t3;
// std::cout << "Time to reduce avbv sets: " << avbv_time.count() << "s\n"<<std::endl;

/* 
for(int i=0; i< avlist.size(); i++){
	
	std::cout<<"av:";
	for(int j=0; j< avlist[i].size(); j++){
		
		std::cout<<avlist[i][j]<<" ";
		
	}
	std::cout<<std::endl;
	
	std::cout<<"bv:";
	for(int j=0; j< bvlist[i].size(); j++){
		
		std::cout<<bvlist[i][j]<<" ";
		
	}
	std::cout<<std::endl;
	
	
} */

// return 0;


for(int n=nmin; n<nmax; n++){
	
// set global_hii for this particular n, for each band 

for(int i=0; i< g.ami.global_hii.size(); i++){
g.ami.global_hii[i]=initial_hii[i]-sigma[i][n];	
}
	

	// double freq=-10+20.0/nmax*n; //0;//-10.0+6.0/nmax*n;
    // double gamma=0;//1e-2;//M_PI/beta*(2*n+1);//1e-4;//M_PI/beta*(2*n+1);//1e-3;//M_PI/beta*(2*n+1);//1e-3;	
	
	std::cout<<"on n="<<n<<std::endl;
	
	std::complex<double> total_sum(0,0);

for(int term=0; term<SM[ord].size(); term++){

std::complex<double> term_sum(0,0);

// if(n==0){
// std::cout<<"On term "<<term <<std::endl;
// g.print_all_edge_info(SM[ord][term].graph);
// }

// for(int comb=0; comb<avlist.size(); comb++){
for(int comb=0; comb<dlist.size(); comb++){


totalcount++;
// av=avlist[comb];
// bv=bvlist[comb];

g.get_avbv_from_delta(SM[ord][term], dlist[comb], av,bv);

// if(comb!=116 && comb!=229){ continue;}

// std::cout<<"Evaluating term "<<term<<" for avbv: ";
// for(int e=0; e<av.size(); e++){

// std::cout<<av[e]<<" "<<bv[e]<<" ";	
	
	
// }
// std::cout<<std::endl;


// std::cout<<std::endl;

// std::cout<<"Checking evaluation criteria"<<std::endl;
// bool can_evaluate=g.is_diagonal(SM[ord][term].bmap, av, bv);

// std::cout<<"Returned "<< can_evaluate<<std::endl;

bool can_evaluate=true;

if(can_evaluate){
	
double umult=g.eval_Ulist(UM, av, bv);

if(umult==0){ continue;}

	
// std::cout<<"Can evaluate - av bv are "<<std::endl;

// for(int i=0; i< av.size(); i++){
	
// std::cout<<av[i]<< " ";	
	
// }
// std::cout<<std::endl;
// for(int i=0; i< bv.size(); i++){
	
// std::cout<<bv[i]<< " ";	
	
// }

g.assign_species(SM[ord][term],av,bv, dup_map);

// if(n==0){
// std::cout<<"On term "<<term <<std::endl;
// g.print_all_edge_info(SM[ord][term].graph);
// }

g.assign_solution(SM[ord][term]);
g.construct_rtss_ami_sol(SM[ord][term]);

// std::cout<<"On graph "<<term<<std::endl;

// g.print_all_edge_info(SM[ord][term].graph);

// std::cout<<"Standard solution set is "<<std::endl;
// g.ami.print_final(5,SM[ord][term].ss.R_,SM[ord][term].ss.P_,SM[ord][term].ss.S_);
// std::cout<<"------- rtss solution is ----------"<<std::endl;
// g.ami.print_final(5,SM[ord][term].rtss.R_,SM[ord][term].rtss.P_,SM[ord][term].rtss.S_);

// for(int i=0; i< SM[ord][term].ss.R0_.size(); i++){
	
	// int spec=SM[ord][term].ss.R0_[i].species_;
	// std::cout<<"Species assigned is "<< spec<<" with energy hii="<<g.ami.global_hii[spec]<<std::endl;
	
// }


// for(int i=0; i< av.size(); i++){
	
	// std::cout<<av[i]<<" ";
// }
// std::cout<<std::endl;

// for(int i=0; i< bv.size(); i++){
	
	// std::cout<<bv[i]<<" ";
// }
// std::cout<<std::endl;

ext.BETA_=beta;//1.0;//50.0;
// double freq=-1.0+1.0/nmax*n; //  M_PI/ext.BETA_*(2*n+1)
ext.external_freq_[0]=std::complex<double>(freq[n],gamma[n]);// M_PI/ext.BETA_*(2*n+1));			//std::complex<double>(0.0,freq);//M_PI/ext.BETA_);
// ext.MU_=(g.ami.global_hii[0]+g.ami.global_hii[1])/2.0;
	
std::complex<double> calc_result;
	if(!SM[ord][term].use_rtss){	
// if(true){	
// std::cout<<"Evaluating ss"<<std::endl;
// SM[ord][term].ss.prefactor_=1.0;
AmiBase::ami_vars amivar=g.ami.construct_ami_vars(SM[ord][term].ss.R0_, SM[ord][term].ss.prefactor_, state, ext);	
calc_result=g.ami.amibase.evaluate(SM[ord][term].ss.ami_parms_,SM[ord][term].ss.R_, SM[ord][term].ss.P_, SM[ord][term].ss.S_,  amivar);	
}else{

// SM[ord][term].rtss.prefactor_	
	
	// std::cout<<"Evaluating rtss with prefactor "<<SM[ord][term].rtss.prefactor_<<std::endl;
AmiBase::ami_vars amivar=g.ami.construct_ami_vars(SM[ord][term].rtss.R0_, SM[ord][term].rtss.prefactor_, state, ext);

// std::cout<<"Energy for evaluation is: ";
// for(int i=0; i< amivar.energy_.size(); i++){
// std::cout<<	amivar.energy_[i]<<" on species "<< SM[ord][term].rtss.R0_[i].species_<< ", ";
	
// }
// std::cout<<std::endl;

// g.ami.print_final(5,SM[ord][term].rtss.R_,SM[ord][term].rtss.P_,SM[ord][term].rtss.S_);
	

calc_result=g.ami.amibase.evaluate(SM[ord][term].rtss.ami_parms_,SM[ord][term].rtss.R_, SM[ord][term].rtss.P_, SM[ord][term].rtss.S_,  amivar);	

// g.ami.print_final(5,SM[ord][term].rtss.R_,SM[ord][term].rtss.P_,SM[ord][term].rtss.S_);
// g.print_all_edge_info(SM[ord][term].graph);
	
	
}


	
// std::cout<<"Result of ("<<SM[ord][term].use_rtss<<") for term "<<term<<" = "<< calc_result<<" with comb "<<comb<<" with Us "<<umult<<" ";	
// std::cout<<calc_result*umult<<std::endl;

termfile<<ord<<" "<<term<<" "<<comb<<" "<<n<<" "<<M_PI/ext.BETA_*(2*n+1)<<" "<<(calc_result*umult).real()<<" "<< (calc_result*umult).imag()<<std::endl;

if(!std::isnan(calc_result.real())){	
next_sigma[ext_band][n]+=calc_result*umult;
total_sum+=calc_result*umult;//g.eval_Ulist(UM, av, bv);
term_sum+=calc_result*umult;
}else{
	isnan_count++;
	
// std::cout<<"Found nan at ord and term "<<ord<<" "<< term <<" for combination "<<comb<<std::endl;	
// g.ami.print_final(5,SM[ord][term].rtss.R_,SM[ord][term].rtss.P_,SM[ord][term].rtss.S_);
// g.print_all_edge_info(SM[ord][term].graph);


// for(int i=0; i< av.size(); i++){
	
	// std::cout<<av[i]<<" ";
// }
// std::cout<<std::endl;

// for(int i=0; i< bv.size(); i++){
	
	// std::cout<<bv[i]<<" ";
// }
// std::cout<<std::endl;

	// return 0;
	
	
}
count++;

}else{
	
	
// std::cout<<"Can't evaluate, not diagonal"<<std::endl;
}



} // end comb loop

// std::cout<<"Beta is "<< beta <<" "<< ext.BETA_<<std::endl;
// termfile<<ord<<" "<<term<<" "<<n<<" "<<M_PI/ext.BETA_*(2*n+1)<<" "<< term_sum.real()<<" "<< term_sum.imag()<<std::endl;
 

} // end term loop 	

// std::cout<<"Total sum was "<< total_sum<<std::endl;
// std::cout<<"Total Count was "<<totalcount<<std::endl<<"Count was"<<count<<std::endl;

if(isnan_count!=0){
	std::cout<<"Some NAN were returned, count="<<isnan_count<<std::endl;
}

// std::complex<double> eta(0,0.05);
// std::complex<double> G=1.0/(freq-g.ami.global_hii[ext_band]-total_sum  );

// file<<ord<<" "<<n<<" "<<freq<<" "<< gamma<<" "<< total_sum.real()<<" "<< total_sum.imag()<<" "<<G.real()<<" "<<std::abs(G.imag())<<std::endl;

}// end matsubara loop (n loop)

}// end ord loop 

}// end ext_band loop



sigma=next_sigma;





}// end iteration loop


// write sigma to a file 
for(int i=0; i< sigma.size(); i++){
	for(int j=0; j< sigma[i].size(); j++){
	
	std::complex<double> G=1.0/(freq[j]-initial_hii[i]-sigma[i][j] );
		
	file<< i<< " "<< j <<" "<< freq[j]<<" "<< gamma[j]<<" "<< sigma[i][j].real() <<" "<<sigma[i][j].imag()<<" "<<G.real()<<" "<<std::abs(G.imag())<<std::endl;
	
	}
}


file.close();
termfile.close();







// g.ggm_construct_ami_sol(ggm, 1e-8, mpi_rank);


// std::cout<<"Printing Graphs"<<std::endl;
// for (int ord=0; ord< SM.size(); ord++){
	// std::cout<<"On ord "<<ord<<std::endl;
// for(int i=0; i< SM[ord].size(); i++){
// std::cout<<"On i "<<i<<std::endl;
// g.number_vertices(SM[ord][i].graph);
// g.print_all_edge_info(SM[ord][i].graph);	
	
// }
// }


	
// AmiGraph::symDET SD;

// g.make_symDET(5, SD);

// for(int i=0; i< SD.size(); i++){

// std::cout<<"Possible graph with sign "<<SD[i].sign<<std::endl;	
// for(int j=0; j< SD[i].plist.size(); j++){
// std::cout<<"("<<SD[i].plist[j].first<<","<<SD[i].plist[j].second<<")"<<std::endl;


// }	
	
	
// }

	

    return 0;
}
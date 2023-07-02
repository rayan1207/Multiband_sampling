#include "ami_spec.hpp"



void AmiSpec::find_closest_points_in_vector(double &closest_lt,double &closest_gt,double point, std::vector<double> vec){
  std::vector<double>  comparison_vector;
  comparison_vector=vec;
  for (int n=0;n<vec.size();n++){
    comparison_vector[n]=std::abs(comparison_vector[n]-point);
    //std::cout<<"comparison vector  "<<comparison_vector[n]<<"   "<<vec[n]<<std::endl;
  }
  //for(int n=0;n<comparison_vector.size();n++){std::cout << "comparison: " << comparison_vector[n]<<std::endl;}
  auto minimum_dif_iter = std::min_element(comparison_vector.begin(),comparison_vector.end());
  double minimum_dif=*minimum_dif_iter;//gives minimum difference between point and contents of vec
  //std::cout<< "\nminimum_dif:  " << minimum_dif << std::endl;
  int minimum_dif_index=0;


  for (int n=0;n<vec.size();n++){
    if (comparison_vector[n]==minimum_dif) minimum_dif_index=n;
    }
	double minimum2=comparison_vector[0];//what if index=0 is closest point?
	int  minimum2_dif_index=0;
	if(minimum_dif_index==0) {
		minimum2=comparison_vector[1];
		minimum2_dif_index=1;
	}
	for (int n=0;n<vec.size();n++){
    if(comparison_vector[n]<minimum2 && n!=minimum_dif_index) {minimum2=comparison_vector[n]; minimum2_dif_index=n;}

  }

  if (vec[minimum_dif_index]>vec[minimum2_dif_index]) {closest_lt=vec[minimum2_dif_index];closest_gt=vec[minimum_dif_index];}
  else {closest_gt=vec[minimum2_dif_index];closest_lt=vec[minimum_dif_index];}

}

// to be called inside the read function . 
void AmiSpec::initialize_linterp(){
	

// std::vector< std::vector<double>::iterator > grid_iter_list;
  grid_iter_list.push_back(AMI_spec_kx_vector_simplified.begin());
  grid_iter_list.push_back(AMI_spec_ky_vector_simplified.begin());
  grid_iter_list.push_back(AMI_spec_freq_vector_simplified.begin());

const int length[] = {(int)AMI_spec_kx_vector_simplified.size(), (int)AMI_spec_ky_vector_simplified.size(),(int)AMI_spec_freq_vector_simplified.size()};

  // array<int,3> grid_sizes;
  grid_sizes[0] = length[0];
  grid_sizes[1] = length[1];
  grid_sizes[2] = length[2];
  
num_elements = grid_sizes[0] * grid_sizes[1]*grid_sizes[2];	
	

// interp_ML_REAL()=InterpMultilinear<3, double>(grid_iter_list.begin(), grid_sizes.begin(), AMI_spec_se_Re_vector.data(), AMI_spec_se_Re_vector.data() + num_elements);
// interp_ML_IMAG()=InterpMultilinear<3, double>(grid_iter_list.begin(), grid_sizes.begin(), AMI_spec_se_Im_vector.data(), AMI_spec_se_Im_vector.data() + num_elements);	
	
	
	
}

void AmiSpec::initialize_W_linterp(){
	

// std::vector< std::vector<double>::iterator > grid_iter_list;
  W_grid_iter_list.push_back(AMI_W_kx_vector_simplified.begin());
  W_grid_iter_list.push_back(AMI_W_ky_vector_simplified.begin());
  W_grid_iter_list.push_back(AMI_W_freq_vector_simplified.begin());

const int length[] = {(int)AMI_W_kx_vector_simplified.size(), (int)AMI_W_ky_vector_simplified.size(),(int)AMI_W_freq_vector_simplified.size()};

  // array<int,3> grid_sizes;
  W_grid_sizes[0] = length[0];
  W_grid_sizes[1] = length[1];
  W_grid_sizes[2] = length[2];
  
W_num_elements = W_grid_sizes[0] * W_grid_sizes[1]*W_grid_sizes[2];	
	
	
}



std::complex<double> AmiSpec::linterp_W(NewAmiCalc::k_vector_t &k, std::complex<double> &X)
{
	
if(X.real()<AMI_W_freq_vector_simplified[0]) return std::complex<double> (0.0,0);
if(X.real()>AMI_W_freq_vector_simplified.back()) return std::complex<double> (0.0,0);


// std::cout<<"Input k is "<< k[0]<<" "<<k[1]<<std::endl;

NewAmiCalc::k_vector_t k_copy=remap_k(k);
// std::cout<<"Remapped k is "<< k_copy[0]<<" "<<k_copy[1]<<std::endl;

array<double,3> args = {k_copy[0], k_copy[1],X.real()};

// set up interpolator ( it would be ideal if this was done only once but I think it might be ok... )

InterpMultilinear<3, double> interp_ML_REAL(W_grid_iter_list.begin(), W_grid_sizes.begin(), AMI_W_se_Re_vector.data(), AMI_W_se_Re_vector.data() + W_num_elements);
InterpMultilinear<3, double> interp_ML_IMAG(W_grid_iter_list.begin(), W_grid_sizes.begin(), AMI_W_se_Im_vector.data(), AMI_W_se_Im_vector.data() + W_num_elements);


return std::complex<double>(interp_ML_REAL.interp(args.begin()),interp_ML_IMAG.interp(args.begin()));

	
}


std::complex<double> AmiSpec::linterp_sigma(NewAmiCalc::k_vector_t &k, std::complex<double> &X)
{
	
if(X.real()<AMI_spec_freq_vector_simplified[0]) return std::complex<double> (0.0,-global_gamma);
if(X.real()>AMI_spec_freq_vector_simplified.back()) return std::complex<double> (0.0,-global_gamma);

NewAmiCalc::k_vector_t k_copy=remap_k(k);

array<double,3> args = {k_copy[0], k_copy[1],X.real()};

// set up interpolator ( it would be ideal if this was done only once but I think it might be ok... )

InterpMultilinear<3, double> interp_ML_REAL(grid_iter_list.begin(), grid_sizes.begin(), AMI_spec_se_Re_vector.data(), AMI_spec_se_Re_vector.data() + num_elements);
InterpMultilinear<3, double> interp_ML_IMAG(grid_iter_list.begin(), grid_sizes.begin(), AMI_spec_se_Im_vector.data(), AMI_spec_se_Im_vector.data() + num_elements);


return std::complex<double>(interp_ML_REAL.interp(args.begin()),interp_ML_IMAG.interp(args.begin()));

	
}


// TODO this does not actually work correctly 
NewAmiCalc::k_vector_t AmiSpec::remap_k(NewAmiCalc::k_vector_t &k_in){
	
NewAmiCalc::k_vector_t k_copy=k_in;
double min=0.0;
double max=2.0*M_PI;
double range=max-min;
for(int for_counter=0;for_counter<k_copy.size();for_counter++){
	if(k_copy[for_counter]>max){
		double mod= std::floor(k_copy[for_counter]/range)*range;      //(std::fmod(abs(k_copy[for_counter]),range))*range;
		k_copy[for_counter]=k_copy[for_counter]-mod;
	}
	if(k_copy[for_counter]<min){
		double mod= (std::abs(std::floor(k_copy[for_counter]/range)))*range;   //(std::fmod(abs(k_copy[for_counter]),range)+1.0)*range;
		k_copy[for_counter]=k_copy[for_counter]+mod;
	}

}	
	
return k_copy;	
	
}


std::complex<double> AmiSpec::get_sigma(NewAmiCalc::k_vector_t &k, std::complex<double> &X)
{
//std::cout<<"entered get sigma"<<std::endl;
  //std::cout<<AMI_spec_freq_vector_simplified[0]<<std::endl;
  //std::cout<<AMI_spec_freq_vector_simplified.back()<<std::endl;

  //change this so that it defaults to the self energy for the 1st/last frequency point

if(X.real()<AMI_spec_freq_vector_simplified[0]) return std::complex<double> (0.0,0.0);
if(X.real()>AMI_spec_freq_vector_simplified.back()) return std::complex<double> (0.0,0.0);
 double k_min=AMI_spec_kx_vector_simplified[0];
 double k_max=AMI_spec_kx_vector_simplified.back();
 
 
NewAmiCalc::k_vector_t k_copy=k;
double min=0.0;
double max=2.0*M_PI;
double range=max-min;
for(int for_counter=0;for_counter<k_copy.size();for_counter++){

	if(k_copy[for_counter]>k_max||k_copy[for_counter]<k_min){
		double modulo_shifted=k_copy[for_counter]-k_min;
		k_copy[for_counter]=+k_min+modulo_shifted-floor(modulo_shifted/(2*M_PI))*2*M_PI;//fmod(k_copy[for_counter])

	}

}

// std::cout<<"k:  "<<k[0]<<"   "<<k[1]<<std::endl;
// std::cout<<"k copy:  "<<k_copy[0]<<"   "<<k_copy[1]<<std::endl;
// std::cout<<"X: "<< X<<std::endl;
 

//std::cout<<"completed boundary checks"<<std::endl;
//std::cout<<"where is the seg fault 1"<<std::endl;
  //std::vector<double> point_wanted_modulo;
  //if(X>M_PI) point_wanted_modulo.push_back(fmod())
  //identify 8 closest points
  double freq_lt=AMI_spec_freq_vector_simplified[0];
  double freq_gt=AMI_spec_freq_vector_simplified[0];
  double kx_lt=AMI_spec_kx_vector_simplified[0];
  double kx_gt=AMI_spec_kx_vector_simplified[0];
  double ky_lt=AMI_spec_ky_vector_simplified[0];
  double ky_gt=AMI_spec_ky_vector_simplified[0];

    //if(freq_vector_simplified[n] > freq_lt && abs(X - freq_vector_simplified[n]) < abs(X - freq_lt)) freq_lt = freq_vector_simplified[n];
    //if(freq_vector_simplified[n] < freq_gt && abs(X - freq_vector_simplified[n]) < abs(X - freq_gt)) freq_gt = freq_vector_simplified[n];

//std::cout<<"where is the seg fault 1.25"<<std::endl;
//std::cout<<X.real()<<"  "<<k[0]<<"  "<<k[0]<<std::endl;
  find_closest_points_in_vector(freq_lt,freq_gt, X.real(), AMI_spec_freq_vector_simplified);
//std::cout<<"found closest freq points"<<std::endl;
  find_closest_points_in_vector(kx_lt,kx_gt, k_copy[0], AMI_spec_kx_vector_simplified);
//std::cout<<"found closest kx points"<<std::endl;
  find_closest_points_in_vector(ky_lt,ky_gt, k_copy[1], AMI_spec_ky_vector_simplified);
 //// std::cout<<"found closest ky points"<<std::endl;
  //std::cout<<"where is the seg fault 1.5"<<std::endl;
  //find_closest_points_in_vector(double &closest_lt,double &closest_gt,double point, std::vector<double> vec)
//std::cout<<"found closest points"<<std::endl;
//std::cout<<"closest points in vector:  "<< X.real() << "   " << freq_lt <<"   "<<freq_gt<<"   "<<k_copy[0]<<"   "<<kx_lt <<"   "<<kx_gt<<"   "<<k_copy[1]<<"   "<<ky_lt<<"   "<<ky_gt<<std::endl;
  //std::cout<<"spec vectors: "<<std::endl;
  //for (int n=0;n<AMI_spec_kx_vector.size();n++){
//std::cout<< AMI_spec_freq_vector[n] << "   " << AMI_spec_kx_vector[n] << "   " <<  AMI_spec_ky_vector[n] <<std::endl;
  //}
  
  int lll_corner;
  int llg_corner;
  int lgl_corner;
  int lgg_corner;
  int gll_corner;
  int glg_corner;
  int ggl_corner;
  int ggg_corner;

assign_corners_indices( freq_lt, freq_gt, kx_lt, kx_gt, ky_lt, ky_gt, AMI_spec_freq_vector,  AMI_spec_kx_vector, AMI_spec_ky_vector, lll_corner, llg_corner, lgl_corner, lgg_corner, gll_corner, glg_corner, ggl_corner, ggg_corner);


  
//std::cout<<"set locations of corners"<<std::endl;
  //std::cout<<"where is the seg fault 2"<<std::endl;
  /*std::cout << "lll: "<<AMI_spec_freq_vector[lll_corner]<<"  "<<AMI_spec_kx_vector[lll_corner]<<"  "<<AMI_spec_ky_vector[lll_corner] << "  Real: " << AMI_spec_se_Re_vector[lll_corner] << "  Imag: " << AMI_spec_se_Im_vector[lll_corner] <<std::endl;
  std::cout << "ggg: "<<AMI_spec_freq_vector[ggg_corner]<<"  "<<AMI_spec_kx_vector[ggg_corner]<<"  "<<AMI_spec_ky_vector[ggg_corner] << "  Real: " << AMI_spec_se_Re_vector[ggg_corner] << "  Imag: " << AMI_spec_se_Im_vector[ggg_corner] <<std::endl;
	std::cout << "llg: "<<AMI_spec_freq_vector[llg_corner]<<"  "<<AMI_spec_kx_vector[llg_corner]<<"  "<<AMI_spec_ky_vector[llg_corner] << "  Real: " << AMI_spec_se_Re_vector[llg_corner] << "  Imag: " << AMI_spec_se_Im_vector[llg_corner] <<std::endl;
	std::cout << "ggl: "<<AMI_spec_freq_vector[ggl_corner]<<"  "<<AMI_spec_kx_vector[ggl_corner]<<"  "<<AMI_spec_ky_vector[ggl_corner] << "  Real: " << AMI_spec_se_Re_vector[ggl_corner] << "  Imag: " << AMI_spec_se_Im_vector[ggl_corner] <<std::endl;
	std::cout << "lgl: "<<AMI_spec_freq_vector[lgl_corner]<<"  "<<AMI_spec_kx_vector[lgl_corner]<<"  "<<AMI_spec_ky_vector[lgl_corner] << "  Real: " << AMI_spec_se_Re_vector[lgl_corner] << "  Imag: " << AMI_spec_se_Im_vector[lgl_corner] <<std::endl;
  std::cout << "gll: "<<AMI_spec_freq_vector[gll_corner]<<"  "<<AMI_spec_kx_vector[gll_corner]<<"  "<<AMI_spec_ky_vector[gll_corner] << "  Real: " << AMI_spec_se_Re_vector[gll_corner] << "  Imag: " << AMI_spec_se_Im_vector[gll_corner] <<std::endl;
	std::cout << "glg: "<<AMI_spec_freq_vector[glg_corner]<<"  "<<AMI_spec_kx_vector[glg_corner]<<"  "<<AMI_spec_ky_vector[glg_corner] << "  Real: " << AMI_spec_se_Re_vector[glg_corner] << "  Imag: " << AMI_spec_se_Im_vector[glg_corner] <<std::endl;
	std::cout << "lgg: "<<AMI_spec_freq_vector[lgg_corner]<<"  "<<AMI_spec_kx_vector[lgg_corner]<<"  "<<AMI_spec_ky_vector[lgg_corner] << "  Real: " << AMI_spec_se_Re_vector[lgg_corner] << "  Imag: " << AMI_spec_se_Im_vector[lgg_corner] <<std::endl;*/

  //linearinterpolation along freq axis
  //std::cout<<"vector sizes: " << AMI_spec_se_Im_vector.size()<< "   "<<AMI_spec_freq_vector.size()<<"   "<< AMI_spec_se_Re_vector.size()<<std::endl;
  //std::cout<<"corners:  "<<lll_corner<<"   "<<llg_corner<<"   "<<lgl_corner<<"   "<<lgg_corner<<"   "<<gll_corner<<"   "<<glg_corner<<"   "<<ggl_corner<<"   "<<ggg_corner<<std::endl;
  
 /* double se_Im_ll_corner = AMI_spec_se_Im_vector[lll_corner] + (X.real() - AMI_spec_freq_vector[lll_corner])*(AMI_spec_se_Im_vector[llg_corner] - AMI_spec_se_Im_vector[lll_corner])/(AMI_spec_freq_vector[llg_corner] - AMI_spec_freq_vector[lll_corner]);
  double se_Im_lg_corner = AMI_spec_se_Im_vector[lgl_corner] + (X.real() - AMI_spec_freq_vector[lgl_corner])*(AMI_spec_se_Im_vector[lgg_corner] - AMI_spec_se_Im_vector[lgl_corner])/(AMI_spec_freq_vector[lgg_corner] - AMI_spec_freq_vector[lgl_corner]);
  double se_Im_gl_corner = AMI_spec_se_Im_vector[gll_corner] + (X.real() - AMI_spec_freq_vector[gll_corner])*(AMI_spec_se_Im_vector[glg_corner] - AMI_spec_se_Im_vector[gll_corner])/(AMI_spec_freq_vector[glg_corner] - AMI_spec_freq_vector[gll_corner]);
  double se_Im_gg_corner = AMI_spec_se_Im_vector[ggl_corner] + (X.real() - AMI_spec_freq_vector[ggl_corner])*(AMI_spec_se_Im_vector[ggg_corner] - AMI_spec_se_Im_vector[ggl_corner])/(AMI_spec_freq_vector[ggg_corner] - AMI_spec_freq_vector[ggl_corner]);

	double se_Re_ll_corner = AMI_spec_se_Re_vector[lll_corner] + (X.real() - AMI_spec_freq_vector[lll_corner])*(AMI_spec_se_Re_vector[llg_corner] - AMI_spec_se_Re_vector[lll_corner])/(AMI_spec_freq_vector[llg_corner] - AMI_spec_freq_vector[lll_corner]);
	double se_Re_lg_corner = AMI_spec_se_Re_vector[lgl_corner] + (X.real() - AMI_spec_freq_vector[lgl_corner])*(AMI_spec_se_Re_vector[lgg_corner] - AMI_spec_se_Re_vector[lgl_corner])/(AMI_spec_freq_vector[lgg_corner] - AMI_spec_freq_vector[lgl_corner]);
	double se_Re_gl_corner = AMI_spec_se_Re_vector[gll_corner] + (X.real() - AMI_spec_freq_vector[gll_corner])*(AMI_spec_se_Re_vector[glg_corner] - AMI_spec_se_Re_vector[gll_corner])/(AMI_spec_freq_vector[glg_corner] - AMI_spec_freq_vector[gll_corner]);
	double se_Re_gg_corner = AMI_spec_se_Re_vector[ggl_corner] + (X.real() - AMI_spec_freq_vector[ggl_corner])*(AMI_spec_se_Re_vector[ggg_corner] - AMI_spec_se_Re_vector[ggl_corner])/(AMI_spec_freq_vector[ggg_corner] - AMI_spec_freq_vector[ggl_corner]);
  */
double se_Im_ll_corner;
double se_Im_lg_corner ;
double se_Im_gl_corner ;
double se_Im_gg_corner ;

double se_Re_ll_corner;
double se_Re_lg_corner;
double se_Re_gl_corner;
double se_Re_gg_corner;

  assign_corners_freq_interp(se_Im_ll_corner,se_Im_lg_corner,se_Im_gl_corner, se_Im_gg_corner, se_Re_ll_corner, se_Re_lg_corner, se_Re_gl_corner, se_Re_gg_corner, AMI_spec_se_Im_vector, AMI_spec_se_Re_vector, lll_corner, llg_corner, lgl_corner, lgg_corner, gll_corner, glg_corner, ggl_corner, ggg_corner,X);
  //std::cout<<"performed linear interpolation"<<std::endl;
  //std::cout<<"\nse_Re_ll_corner:  "<<se_Re_ll_corner<<std::endl;
  //std::cout<<"se_Re_gg_corner:  "<<se_Re_gg_corner<<std::endl;

  //bilinear interpolation in momentum plane
  double bilinear_prefactor=1/((kx_gt-kx_lt)*(ky_gt-ky_lt));
  //std::cout<<"where is the seg fault 3"<<std::endl;
  //std::cout<<"prefactor:  "<<bilinear_prefactor<<std::endl;
  //std::cout<<"se_Re_ll_corner*(kx_gt-k[0])*(ky_gt-k[1]):  "<<se_Re_ll_corner*(kx_gt-k[0])*(ky_gt-k[1])<<std::endl;
  //std::cout << "se_Re_gg_corner*(k[0]-kx_lt)*(k[1]-ky_lt):  " << se_Re_gg_corner*(k[0]-kx_lt)*(k[1]-ky_lt)<<std::endl;
  double bilinear_main_Re=se_Re_ll_corner*(kx_gt-k_copy[0])*(ky_gt-k_copy[1]) + se_Re_lg_corner*(kx_gt-k_copy[0])*(k_copy[1]-ky_lt) + se_Re_gl_corner*(k_copy[0]-kx_lt)*(ky_gt-k_copy[1]) + se_Re_gg_corner*(k_copy[0]-kx_lt)*(k_copy[1]-ky_lt);
  double  se_Re_interpolated=bilinear_prefactor*bilinear_main_Re;

	double bilinear_main_Im=se_Im_ll_corner*(kx_gt-k_copy[0])*(ky_gt-k_copy[1]) + se_Im_lg_corner*(kx_gt-k_copy[0])*(k_copy[1]-ky_lt) + se_Im_gl_corner*(k_copy[0]-kx_lt)*(ky_gt-k_copy[1]) + se_Im_gg_corner*(k_copy[0]-kx_lt)*(k_copy[1]-ky_lt);
  double  se_Im_interpolated=bilinear_prefactor*bilinear_main_Im;
  //std::cout<<"return"<<std::endl;
//std::cout<<"where is the seg fault 4"<<std::endl;

// std::cout<<"Returning "<< std::complex<double> (se_Re_interpolated,se_Im_interpolated) <<std::endl;

  return std::complex<double> (se_Re_interpolated,se_Im_interpolated);

}

std::vector<double> AmiSpec::return_simple_grid_vector(std::vector<double> &in_vector){
  std::vector<double> out_vector;
  out_vector.push_back(in_vector[0]);
  //std::cout<<"simplification in vector size: "<<in_vector.size()<<std::endl;
  for (int n=0;n<in_vector.size();n++){
    //std::cout<<"in vector "<<n<<"  "<<in_vector[n]<<std::endl;
    if (in_vector[n] != out_vector.back()){
      if (std::count(out_vector.begin(),out_vector.end(),in_vector[n] )) return out_vector;
      else {
        //std::cout<<"push back!"<<std::endl;
        out_vector.push_back(in_vector[n]);
      }
    }
  }
  return out_vector;
  //std::cout<<"problem! Return simple grid vector never returned a vector"<<std::endl;
}

//std::complex<double> AmiSpec::get_sigma(NewAmiCalc::k_vector_t &k, std::complex<double> &X){}

void AmiSpec::read_self_energy(std::string file_name){
  

  //expects first line of file to be column headings and throws away. Without heading the interpoolation will be messed up.
  //Also expects labelling to be low to high, behaviour other way around is unknown
  std::ifstream MyReadFile(file_name);
  std::string myText;
  std::stringstream ss;
  //std::ifstream MyReadFile(argv[1]);
  //std::cout << argv[1] << std::endl;
  //std::vector<float> frequency;
  //std::vector<float> gf;
  //int count = 0;
  //std::vector<double> se_Im_vector;
  //std::vector<double> se_Im_err_vector;
  //std::vector<double> se_Re_vector;
  //std::vector<double> se_Re_err_vector;
  //std::vector<double> freq_vector;
  //std::vector<double> ky_vector;
  //std::vector<double> kx_vector;


  std::vector<std::string> output_vector;
  std::string a_string;
  getline (MyReadFile, myText);
  while (getline (MyReadFile, myText)) {
    //std::cout <<"myText  " <<myText << std::endl;
    ss.str( myText);

   output_vector.clear();
   //   std::cout << "contents of ss: " << ss.str() << std::endl;
   //std::cout << "a string OUTSIDE: " << a_string << std::endl;
     a_string="";
   while(ss >> a_string){
     //std::cout << "a string: " << a_string << std::endl;
     output_vector.push_back(a_string);
   }
   //std::cout<<"got to here "<< output_vector.back()<<std::endl;
   AMI_spec_kx_vector.push_back(std::stod(output_vector[0]));
   AMI_spec_ky_vector.push_back(std::stod(output_vector[1]));
   AMI_spec_freq_vector.push_back(std::stod(output_vector[2]));
   AMI_spec_se_Re_vector.push_back(std::stod(output_vector[3]));
  AMI_spec_se_Im_vector.push_back(std::stod(output_vector[4]));
   AMI_spec_se_Re_err_vector.push_back(std::stod(output_vector[5]));
   AMI_spec_se_Im_err_vector.push_back(std::stod(output_vector[6]));
   ss.clear();
  }
  
  AMI_spec_ky_vector_simplified = return_simple_grid_vector(AMI_spec_ky_vector);
  AMI_spec_kx_vector_simplified = return_simple_grid_vector(AMI_spec_kx_vector);
 AMI_spec_freq_vector_simplified = return_simple_grid_vector(AMI_spec_freq_vector);

std::cout<<"Kx, ky and w vector sizes are "<<std::endl;
std::cout<<AMI_spec_kx_vector_simplified.size()<<" "<< AMI_spec_ky_vector_simplified.size()<<" "<<AMI_spec_freq_vector_simplified.size()<<std::endl;
// exit(0);

std::cout<<"Initializing linterp"<<std::endl;
initialize_linterp();



 //std::cout<<"simplified vectors"<<std::endl;
 //std::cout<<AMI_spec_freq_vector_simplified<<std::endl;
 //std::cout<<AMI_spec_kx_vector_simplified<<std::endl;
 //std::cout<<AMI_spec_ky_vector_simplified<<std::endl;
//for (std::vector<double>::const_iterator i = AMI_spec_ky_vector_simplified.begin(); i != AMI_spec_ky_vector_simplified.end(); ++i)
//    std::cout << *i << ' ';
//std::cout<<"\n";std::cout<<"\n";std::cout<<"\n";
/*for (std::vector<double>::const_iterator i = AMI_spec_kx_vector.begin(); i != AMI_spec_kx_vector.end(); ++i)    std::cout << *i << ' ';
std::cout<<"\n";std::cout<<"\n";std::cout<<"\n";

for (std::vector<double>::const_iterator i = AMI_spec_kx_vector_simplified.begin(); i != AMI_spec_kx_vector_simplified.end(); ++i)
    std::cout << *i << ' ';
std::cout<<"\n";std::cout<<"\n";std::cout<<"\n";
*/
//for (std::vector<double>::const_iterator i = AMI_spec_freq_vector_simplified.begin(); i != AMI_spec_freq_vector_simplified.end(); ++i)    std::cout << *i << ' ';
 //std::cout<<"\n";std::cout<<"\n";std::cout<<"\n";
 
 //std::cout<<"\n";
   //tk::spline spline_I(freq_vector,se_Re_vector);
   //std::cout<< spline_I(std::stod(argv[2]))<<std::endl;
 //std::vector<double> point_wanted={0.0,0.314159,0.0};//0.314159
 //std::cout <<   get_sigma(point_wanted, kx_vector, ky_vector, freq_vector, ky_vector_simplified,  kx_vector_simplified, freq_vector_simplified,se_Im_vector, se_Re_vector) << std::endl;
  MyReadFile.close();

}


void AmiSpec::read_W(std::string file_name){
  

  //expects first line of file to be column headings and throws away. Without heading the interpoolation will be messed up.
  //Also expects labelling to be low to high, behaviour other way around is unknown
  std::ifstream MyReadFile(file_name);
  std::string myText;
  std::stringstream ss;

  std::vector<std::string> output_vector;
  std::string a_string;
  getline (MyReadFile, myText);
  while (getline (MyReadFile, myText)) {
    //std::cout <<"myText  " <<myText << std::endl;
    ss.str( myText);

   output_vector.clear();
   //   std::cout << "contents of ss: " << ss.str() << std::endl;
   //std::cout << "a string OUTSIDE: " << a_string << std::endl;
     a_string="";
   while(ss >> a_string){
     //std::cout << "a string: " << a_string << std::endl;
     output_vector.push_back(a_string);
   }
   //std::cout<<"got to here "<< output_vector.back()<<std::endl;
   AMI_W_kx_vector.push_back(std::stod(output_vector[0]));
   AMI_W_ky_vector.push_back(std::stod(output_vector[1]));
   AMI_W_freq_vector.push_back(std::stod(output_vector[2]));
   AMI_W_se_Re_vector.push_back(std::stod(output_vector[3]));
  AMI_W_se_Im_vector.push_back(std::stod(output_vector[4]));
   AMI_W_se_Re_err_vector.push_back(std::stod(output_vector[5]));
   AMI_W_se_Im_err_vector.push_back(std::stod(output_vector[6]));
   ss.clear();
  }
  
  AMI_W_ky_vector_simplified = return_simple_grid_vector(AMI_W_ky_vector);
  AMI_W_kx_vector_simplified = return_simple_grid_vector(AMI_W_kx_vector);
 AMI_W_freq_vector_simplified = return_simple_grid_vector(AMI_W_freq_vector);

std::cout<<"Kx, ky and w vector sizes are "<<std::endl;
std::cout<<AMI_W_kx_vector_simplified.size()<<" "<< AMI_W_ky_vector_simplified.size()<<" "<<AMI_W_freq_vector_simplified.size()<<std::endl;
// exit(0);

std::cout<<"Initializing linterp"<<std::endl;
initialize_W_linterp();



  MyReadFile.close();

}



void AmiSpec::assign_corners_indices(double freq_lt,double freq_gt,double kx_lt,double kx_gt,double ky_lt,double ky_gt,std::vector<double>& AMI_spec_freq_vector, std::vector<double>& AMI_spec_kx_vector,std::vector<double>& AMI_spec_ky_vector,int& lll_corner,int& llg_corner,int& lgl_corner,int& lgg_corner,int& gll_corner,int& glg_corner,int& ggl_corner,int& ggg_corner){
int assignment_counter=0;
for (int n=0;n<AMI_spec_kx_vector.size();n++){
    if (AMI_spec_freq_vector[n]==freq_lt && AMI_spec_kx_vector[n]==kx_lt && AMI_spec_ky_vector[n]==ky_lt) {lll_corner=n; assignment_counter++;}
    if (AMI_spec_freq_vector[n]==freq_gt && AMI_spec_kx_vector[n]==kx_lt && AMI_spec_ky_vector[n]==ky_lt) {llg_corner=n; assignment_counter++;}
    if (AMI_spec_freq_vector[n]==freq_lt && AMI_spec_kx_vector[n]==kx_lt && AMI_spec_ky_vector[n]==ky_gt) {lgl_corner=n; assignment_counter++;}
    if (AMI_spec_freq_vector[n]==freq_gt && AMI_spec_kx_vector[n]==kx_lt && AMI_spec_ky_vector[n]==ky_gt) {lgg_corner=n; assignment_counter++;}
    if (AMI_spec_freq_vector[n]==freq_lt && AMI_spec_kx_vector[n]==kx_gt && AMI_spec_ky_vector[n]==ky_lt) {gll_corner=n; assignment_counter++;}
    if (AMI_spec_freq_vector[n]==freq_gt && AMI_spec_kx_vector[n]==kx_gt && AMI_spec_ky_vector[n]==ky_lt) {glg_corner=n; assignment_counter++;}
    if (AMI_spec_freq_vector[n]==freq_lt && AMI_spec_kx_vector[n]==kx_gt && AMI_spec_ky_vector[n]==ky_gt) {ggl_corner=n; assignment_counter++;}
    if (AMI_spec_freq_vector[n]==freq_gt && AMI_spec_kx_vector[n]==kx_gt && AMI_spec_ky_vector[n]==ky_gt) {ggg_corner=n; assignment_counter++;}
    if (assignment_counter==8) break;
  }


}

void AmiSpec::assign_corners_freq_interp(double& se_Im_ll_corner,double& se_Im_lg_corner,double& se_Im_gl_corner,double& se_Im_gg_corner,double& se_Re_ll_corner,double& se_Re_lg_corner,double& se_Re_gl_corner,double& se_Re_gg_corner,std::vector<double>& AMI_spec_se_Im_vector,std::vector<double>& AMI_spec_se_Re_vector,int lll_corner,int llg_corner,int lgl_corner,int lgg_corner,int gll_corner,int glg_corner,int ggl_corner,int ggg_corner,std::complex<double> &X){
  se_Im_ll_corner = AMI_spec_se_Im_vector[lll_corner] + (X.real() - AMI_spec_freq_vector[lll_corner])*(AMI_spec_se_Im_vector[llg_corner] - AMI_spec_se_Im_vector[lll_corner])/(AMI_spec_freq_vector[llg_corner] - AMI_spec_freq_vector[lll_corner]);
  se_Im_lg_corner = AMI_spec_se_Im_vector[lgl_corner] + (X.real() - AMI_spec_freq_vector[lgl_corner])*(AMI_spec_se_Im_vector[lgg_corner] - AMI_spec_se_Im_vector[lgl_corner])/(AMI_spec_freq_vector[lgg_corner] - AMI_spec_freq_vector[lgl_corner]);
  se_Im_gl_corner = AMI_spec_se_Im_vector[gll_corner] + (X.real() - AMI_spec_freq_vector[gll_corner])*(AMI_spec_se_Im_vector[glg_corner] - AMI_spec_se_Im_vector[gll_corner])/(AMI_spec_freq_vector[glg_corner] - AMI_spec_freq_vector[gll_corner]);
  se_Im_gg_corner = AMI_spec_se_Im_vector[ggl_corner] + (X.real() - AMI_spec_freq_vector[ggl_corner])*(AMI_spec_se_Im_vector[ggg_corner] - AMI_spec_se_Im_vector[ggl_corner])/(AMI_spec_freq_vector[ggg_corner] - AMI_spec_freq_vector[ggl_corner]);

	se_Re_ll_corner = AMI_spec_se_Re_vector[lll_corner] + (X.real() - AMI_spec_freq_vector[lll_corner])*(AMI_spec_se_Re_vector[llg_corner] - AMI_spec_se_Re_vector[lll_corner])/(AMI_spec_freq_vector[llg_corner] - AMI_spec_freq_vector[lll_corner]);
	se_Re_lg_corner = AMI_spec_se_Re_vector[lgl_corner] + (X.real() - AMI_spec_freq_vector[lgl_corner])*(AMI_spec_se_Re_vector[lgg_corner] - AMI_spec_se_Re_vector[lgl_corner])/(AMI_spec_freq_vector[lgg_corner] - AMI_spec_freq_vector[lgl_corner]);
	se_Re_gl_corner = AMI_spec_se_Re_vector[gll_corner] + (X.real() - AMI_spec_freq_vector[gll_corner])*(AMI_spec_se_Re_vector[glg_corner] - AMI_spec_se_Re_vector[gll_corner])/(AMI_spec_freq_vector[glg_corner] - AMI_spec_freq_vector[gll_corner]);
	se_Re_gg_corner = AMI_spec_se_Re_vector[ggl_corner] + (X.real() - AMI_spec_freq_vector[ggl_corner])*(AMI_spec_se_Re_vector[ggg_corner] - AMI_spec_se_Re_vector[ggl_corner])/(AMI_spec_freq_vector[ggg_corner] - AMI_spec_freq_vector[ggl_corner]);
  

}

/*void AmiSpec::assign_corners_freq_interp(double& se_Im_ll_corner,double& se_Im_lg_corner,double& se_Im_gl_corner,double& se_Im_gg_corner,double& se_Re_ll_corner,double& se_Re_lg_corner,double& se_Re_gl_corner,double& se_Re_gg_corner,std::vector<double>& AMI_spec_se_Im_vector,std::vector<double>& AMI_spec_se_Re_vector,double freq_lt,double freq_gt,double kx_lt,double kx_gt,double ky_lt,std::complex<double> &X){
  se_Im_ll_corner = AMI_spec_se_Im_vector[lll_corner] + (X.real() - AMI_spec_freq_vector[lll_corner])*(AMI_spec_se_Im_vector[llg_corner] - AMI_spec_se_Im_vector[lll_corner])/(AMI_spec_freq_vector[llg_corner] - AMI_spec_freq_vector[lll_corner]);
  se_Im_lg_corner = AMI_spec_se_Im_vector[lgl_corner] + (X.real() - AMI_spec_freq_vector[lgl_corner])*(AMI_spec_se_Im_vector[lgg_corner] - AMI_spec_se_Im_vector[lgl_corner])/(AMI_spec_freq_vector[lgg_corner] - AMI_spec_freq_vector[lgl_corner]);
  se_Im_gl_corner = AMI_spec_se_Im_vector[gll_corner] + (X.real() - AMI_spec_freq_vector[gll_corner])*(AMI_spec_se_Im_vector[glg_corner] - AMI_spec_se_Im_vector[gll_corner])/(AMI_spec_freq_vector[glg_corner] - AMI_spec_freq_vector[gll_corner]);
  se_Im_gg_corner = AMI_spec_se_Im_vector[ggl_corner] + (X.real() - AMI_spec_freq_vector[ggl_corner])*(AMI_spec_se_Im_vector[ggg_corner] - AMI_spec_se_Im_vector[ggl_corner])/(AMI_spec_freq_vector[ggg_corner] - AMI_spec_freq_vector[ggl_corner]);

	se_Re_ll_corner = AMI_spec_se_Re_vector[lll_corner] + (X.real() - AMI_spec_freq_vector[lll_corner])*(AMI_spec_se_Re_vector[llg_corner] - AMI_spec_se_Re_vector[lll_corner])/(AMI_spec_freq_vector[llg_corner] - AMI_spec_freq_vector[lll_corner]);
	se_Re_lg_corner = AMI_spec_se_Re_vector[lgl_corner] + (X.real() - AMI_spec_freq_vector[lgl_corner])*(AMI_spec_se_Re_vector[lgg_corner] - AMI_spec_se_Re_vector[lgl_corner])/(AMI_spec_freq_vector[lgg_corner] - AMI_spec_freq_vector[lgl_corner]);
	se_Re_gl_corner = AMI_spec_se_Re_vector[gll_corner] + (X.real() - AMI_spec_freq_vector[gll_corner])*(AMI_spec_se_Re_vector[glg_corner] - AMI_spec_se_Re_vector[gll_corner])/(AMI_spec_freq_vector[glg_corner] - AMI_spec_freq_vector[gll_corner]);
	se_Re_gg_corner = AMI_spec_se_Re_vector[ggl_corner] + (X.real() - AMI_spec_freq_vector[ggl_corner])*(AMI_spec_se_Re_vector[ggg_corner] - AMI_spec_se_Re_vector[ggl_corner])/(AMI_spec_freq_vector[ggg_corner] - AMI_spec_freq_vector[ggl_corner]);
  

}*/
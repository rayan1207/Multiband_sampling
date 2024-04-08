#include "mini_ami.hpp"
#include "algorithm"


void mband::print_assigned_species(std::vector<std::vector<std::vector<int>>> interaction_species){
	     for (int i = 0; i < interaction_species.size(); i++) {
        std::cout << "Interactions are  " << i << ":" << std::endl;
        // Iterate through each row
        for (int j = 0; j < interaction_species[i].size(); j++) {
            // Iterate through each column
            for (int k = 0; k < interaction_species[i][j].size(); k++) {
                std::cout << interaction_species[i][j][k] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
	
}
void mband::print_interactions(AmiGraph::graph_t &graph,AmiGraph::edge_vector_t b_vector, std::vector<AmiGraph::edge_vector_t> f_vector){
	std::cout  << " Printing interaction edges of each bosonic line " << std::endl;
	std::cout << " (1) incoming source, (2) outgoing source , (3)  incoming target , (4) outgoing target" << "\n ";
    //std::assert(b_vector.size()==f_vector.size());
	 for (int i=0;i <b_vector.size();i++){
		 std::cout << "Edge = ("  << graph[source(b_vector[i],graph)].index_<<"," << graph[target(b_vector[i],graph)].index_ <<") \n "  ;
		 
		 for(int j =0; j<f_vector[i].size();j++)
		 {
		 std::cout << "("<<j+1 << ")" << "edge is " <<  "("<<graph[source(f_vector[i][j],graph)].index_ << ","<<graph[target(f_vector[i][j],graph)].index_  
		 << ")" <<std::endl;
			 
		
	 
	
	     }
	 }
}
void mband::print_match(std::vector<int> vec,int ord) {
	if (vec.empty()){
		std::cout<<ord<<"empty" << std::endl;
	}
	else{
	std::cout<<"for order " << ord  << std::endl;
    for (int elem : vec) {
        std::cout << elem << " ";
    }
    std::cout << std::endl;
	}
}	
	

double mband::perturb(double number) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(0.1, 0.9);
    std::uniform_int_distribution<int> sign_dist(0, 1);
    int sign = (sign_dist(gen) == 0) ? -1 : 1;
    return number + (dist(gen) * 1.0e-7 * sign);
}


std::vector<double> mband::sumVectors(std::vector<std::vector<double>> vectors)
{
    std::vector<double> sum(vectors[0].size(), 0);
    for (auto vec : vectors)
    {
        for (int i = 0; i < vec.size(); i++)
        {   
			if (sum[i] == 0){
            sum[i] += vec[i];
			}
        }
    }
    return sum;
}
std::vector<std::complex<double>> mband::convertToComplex(const std::vector<double> vec) {
    std::vector<std::complex<double>> cplx_vec;
    cplx_vec.reserve(vec.size()); 

    for (auto elem : vec) {
        cplx_vec.emplace_back(elem, 0.0);
    }

    return cplx_vec;
}

std::vector<std::vector<double>>  mband::band_to_hab(std::vector<std::vector<int>> band) {
	std::vector<double> hab = mband::energy;
    std::vector<std::vector<double>> band_hab(band.size(), std::vector<double>(band[0].size()));
    for (int i = 0; i < band.size(); i++) {
        for (int j = 0; j < band[i].size(); j++) {
            //band_hab[i][j] =  -1*mband::perturb(hab[band[i][j] - 1]);
			band_hab[i][j] =  hab[band[i][j] - 1];
        }
    }
    std::cout << "Band index replaced by energy:\n";
    for (int i = 0; i < band_hab.size(); i++)
    {
        for (int j = 0; j < band_hab[i].size(); j++)
        {
            std::cout << band_hab[i][j] << " ";
        }
        std::cout << "\n\n";
    }
    return band_hab;
}

std::vector<std::complex<double>> mband::generate_ept(std::vector<std::vector<int>> epsilon, std::vector<double> band_value) {
	//std::sort(epsilon.begin(), epsilon.end());
    //epsilon.erase(std::unique(epsilon.begin(), epsilon.end()), epsilon.end());
    std::vector<std::vector<double>> results;
    std::vector<double> collect;
    for (int  i = 0; i < epsilon.size(); i++) {
        
        for (int j = 0; j < epsilon[i].size(); j++) {
            collect.push_back(-1*epsilon[i][j] * band_value[i]);
        }
        results.push_back(collect);
        collect.clear();;

    }
    std::vector<double> sum = mband::sumVectors(results);



    return mband::convertToComplex(sum);
    
}


double mband::Umatch(const std::vector<std::vector<int>>& int_matrix, const std::vector<double>& int_value, const std::vector<std::vector<int>>& int_species) {
    double U = 1.00;
    for (int i = 0; i < int_species.size(); i++) {
        for (int j = 0; j < int_matrix.size(); j++) {
            if ((int_species[i][1] == int_matrix[j][1] && int_species[i][2] == int_matrix[j][2] && int_species[i][3]
                == int_matrix[j][3] && int_species[i][0] == int_matrix[j][0])) {
                U = U * int_value[j];
            }
        }
       }
    return U;
}
/*
double mband::Hubbard_Energy(NewAmiCalc::ext_vars ext,std::vector<double> momenta, int species,mband::params_param param){
	if (species ==1){
		return -2*(std::cos(momenta[0]) + std::cos(momenta[1])) -ext.MU_.real()-0.5*ext.H_-4*param.tp*(std::cos(momenta[0])*std::cos(momenta[1]));	
	}
	if (species ==2){
		return -2*(std::cos(momenta[0]) + std::cos(momenta[1])) -ext.MU_.real()+0.5*ext.H_-4*param.tp*(std::cos(momenta[0])*std::cos(momenta[1]));
	}
	else{
		std::cerr<<" Species numer should be 1 or 2 for hubbard problem"<< std::endl;
		return 0.0;
	}

}

*/
double mband::Hubbard_Energy(NewAmiCalc::ext_vars ext,std::vector<double> momenta, int species,mband::params_param param){
	if (species ==1){
		return -2*(std::cos(momenta[0]) + std::cos(momenta[1])) -ext.MU_.real()-0.0*ext.H_-4*param.tp*(std::cos(momenta[0])*std::cos(momenta[1]));	
	}
	if (species ==2){
		return -2*(std::cos(momenta[0]) + std::cos(momenta[1])) -ext.MU_.real()+0.0*ext.H_-4*param.tp*(std::cos(momenta[0])*std::cos(momenta[1]));
	}
	else{
		std::cerr<<" Species numer should be 1 or 2 for hubbard problem"<< std::endl;
		return 0.0;
	}

}
/*
double mband::Bilayer_Hubbard_Energy(NewAmiCalc::ext_vars ext,std::vector<double> momenta, int species,mband::params_param param){
	if (species ==1){
		return -2*(1+param.tperp_p)*(std::cos(momenta[0]) + std::cos(momenta[1]))-param.tperp -ext.MU_.real()-0.5*ext.H_-4*param.tp*(std::cos(momenta[0])*std::cos(momenta[1]));	
	}
	if (species ==2){
		return -2*(1+param.tperp_p)*(std::cos(momenta[0]) + std::cos(momenta[1]))-param.tperp -ext.MU_.real()+0.5*ext.H_-4*param.tp*(std::cos(momenta[0])*std::cos(momenta[1]));
	}
	if (species ==3){
		return -2*(1-param.tperp_p)*(std::cos(momenta[0]) + std::cos(momenta[1]))+param.tperp -ext.MU_.real()-0.5*ext.H_-4*param.tp*(std::cos(momenta[0])*std::cos(momenta[1]));	
	}
	if (species ==4){
		return -2*(1-param.tperp_p)*(std::cos(momenta[0]) + std::cos(momenta[1]))+param.tperp -ext.MU_.real()+0.5*ext.H_-4*param.tp*(std::cos(momenta[0])*std::cos(momenta[1]));
	}
	else{
		std::cerr<<" Species numer should be 1,2,3 and 4 for bi-layer hubbard problem"<< std::endl;
		return 0.0;
	}

}*/
double mband::Bilayer_Hubbard_Energy(NewAmiCalc::ext_vars ext,std::vector<double> momenta, int species,mband::params_param param){
	if (species ==1 || species ==2){
		return -2*(1)*(std::cos(momenta[0]) + std::cos(momenta[1])) -ext.MU_.real()-4*param.tp*(std::cos(momenta[0])*std::cos(momenta[1]))-2*(param.tpp)*(std::cos(2*momenta[0]) - std::cos(2*momenta[1])) -param.tbs- param.tperp*std::pow((std::cos(momenta[0]) - std::cos(momenta[1])),2);	
	}
	
	if (species ==3 || species ==4 ){
		return -2*(1)*(std::cos(momenta[0]) + std::cos(momenta[1]))-ext.MU_.imag()-4*param.tp*(std::cos(momenta[0])*std::cos(momenta[1]))-2*(param.tpp)*(std::cos(2*momenta[0]) - std::cos(2*momenta[1]))+ param.tbs + param.tperp*std::pow((std::cos(momenta[0]) - std::cos(momenta[1])),2);
	}

	else{
		std::cerr<<" Species numer should be 1-4 for Bi-layer hubbard problem"<< std::endl;
		return 0.0;
	}
}

double mband::Trilayer_Hubbard_Energy(NewAmiCalc::ext_vars ext,std::vector<double> momenta, int species,mband::params_param param){
	if (species ==1 || species ==2){
		return -2*(1)*(std::cos(momenta[0]) + std::cos(momenta[1])) -ext.MU_.real()-4*param.tp*(std::cos(momenta[0])*std::cos(momenta[1]))-2*(param.tpp)*(std::cos(2*momenta[0]) - std::cos(2*momenta[1])) - std::sqrt(2)*(param.tbs + param.tperp*std::pow((std::cos(momenta[0]) - std::cos(momenta[1])),2));	
	}
	if (species ==3 || species ==4){
		return -2*(1)*(std::cos(momenta[0]) + std::cos(momenta[1]))-ext.MU_.imag()-4*param.tp*(std::cos(momenta[0])*std::cos(momenta[1]))-2*(param.tpp)*(std::cos(2*momenta[0]) - std::cos(2*momenta[1]));	
	}
	if (species ==5 || species ==6 ){
		return -2*(1)*(std::cos(momenta[0]) + std::cos(momenta[1]))-ext.MU_.real()-4*param.tp*(std::cos(momenta[0])*std::cos(momenta[1]))-2*(param.tpp)*(std::cos(2*momenta[0]) - std::cos(2*momenta[1])) + std::sqrt(2)*(param.tbs + param.tperp*std::pow((std::cos(momenta[0]) - std::cos(momenta[1])),2));
	}

	else{
		std::cerr<<" Species numer should be 1-6 for Tri-layer hubbard problem"<< std::endl;
		return 0.0;
	}

}
double mband::Quadlayer_Hubbard_Energy(NewAmiCalc::ext_vars ext,std::vector<double> momenta, int species,mband::params_param param){
	if (species ==1 || species ==2){
		return -2*(1)*(std::cos(momenta[0]) + std::cos(momenta[1])) -ext.MU_.real()-4*param.tp*(std::cos(momenta[0])*std::cos(momenta[1]))-2*(param.tpp)*(std::cos(2*momenta[0]) - std::cos(2*momenta[1])) - (std::sqrt(5)+1)*(param.tbs + param.tperp*std::pow((std::cos(momenta[0]) - std::cos(momenta[1])),2))/2;	
	}
	if (species ==3 || species ==4){
		return -2*(1)*(std::cos(momenta[0]) + std::cos(momenta[1]))-ext.MU_.imag()-4*param.tp*(std::cos(momenta[0])*std::cos(momenta[1]))-2*(param.tpp)*(std::cos(2*momenta[0]) - std::cos(2*momenta[1]))- (std::sqrt(5)-1)*(param.tbs + param.tperp*std::pow((std::cos(momenta[0]) - std::cos(momenta[1])),2))/2;	
	}
	if (species ==5 || species ==6 ){
		return -2*(1)*(std::cos(momenta[0]) + std::cos(momenta[1]))-ext.MU_.real()-4*param.tp*(std::cos(momenta[0])*std::cos(momenta[1]))-2*(param.tpp)*(std::cos(2*momenta[0]) - std::cos(2*momenta[1])) + (std::sqrt(5)-1)*(param.tbs + param.tperp*std::pow((std::cos(momenta[0]) - std::cos(momenta[1])),2))/2;
	}
	
	if (species ==7 || species ==8 ){
		return -2*(1)*(std::cos(momenta[0]) + std::cos(momenta[1]))-ext.MU_.real()-4*param.tp*(std::cos(momenta[0])*std::cos(momenta[1]))-2*(param.tpp)*(std::cos(2*momenta[0]) - std::cos(2*momenta[1])) + (std::sqrt(5)+1)*(param.tbs + param.tperp*std::pow((std::cos(momenta[0]) - std::cos(momenta[1])),2))/2;
	}


	else{
		std::cerr<<" Species numer should be 1-8 for Quad-layer hubbard problem"<< std::endl;
		return 0.0;
	}

}


double mband::Triangular_Hubbard_Energy(NewAmiCalc::ext_vars ext,std::vector<double> momenta, int species,mband::params_param param){
		if (species ==1){
		return -2*(std::cos(momenta[1]) + 2*std::cos(momenta[1]/2)*std::cos(std::sqrt(3)*momenta[0]/2) ) -ext.MU_.real()-ext.H_-2*param.tp*(std::cos(std::sqrt(3)*momenta[0])+ 2*cos(std::sqrt(3)*momenta[0]/2)*std::cos(3*momenta[1]/2));	
	}
	if (species ==2){
		return -2*(std::cos(momenta[1]) + 2*std::cos(momenta[1]/2)*std::cos(std::sqrt(3)*momenta[0]/2) ) -ext.MU_.real()+ext.H_-2*param.tp*(std::cos(std::sqrt(3)*momenta[0])+ 2*cos(std::sqrt(3)*momenta[0]/2)*std::cos(3*momenta[1]/2));
	}
	else{
		std::cerr<<" Species numer should be 1 or 2 for triangular hubbard problem"<< std::endl;
		return 0.0;
	}

}

double mband::non_local_U_formfactor(std::vector<double> vq) {
    return 2 * std::cos(vq[0]) + 2 * std::cos(vq[1]);
}




double mband::Stdev(double total_sq, double mean, int n) {
    double sample_var = (total_sq - (mean * mean * n)) / (n - 1.0);
    return std::sqrt(sample_var/n);
}


void mband::filter(std::vector<std::vector<int>>& possible_species, const std::vector<int>& list) {
    if (list.empty()) {

    }
    else {
        size_t x = static_cast<size_t>(list[0]);
        size_t y = static_cast<size_t>(list[1]);

        std::vector<std::vector<int>> filter;
        filter.reserve(possible_species.size());

        for (const auto& vec : possible_species) {
            if (vec[x] == vec[y]) {
                filter.emplace_back(vec);
            }
        }

        possible_species = std::move(filter);
    }

}

std::vector<int>  mband::interaction_index(const  std::vector<std::vector<int>>& int_species) {
   std::vector<int> vec;
   const std::vector<std::vector<int>> int_matrix = mband::interaction_legs;
    for (int i = 0; i < int_species.size(); i++) {
        for (int j = 0; j < int_matrix.size(); j++) {
            if ((int_species[i][1] == int_matrix[j][1] && int_species[i][2] == int_matrix[j][2] && int_species[i][3]
                == int_matrix[j][3] && int_species[i][0] == int_matrix[j][0])) {
                vec.push_back(j);
            }
        }
        
    }
   
    return vec;
   
}


std::complex<double> mband::gfunc_pp(std::vector<double> momenta, params_param &param, int bandindex){
	double kz=0;
	if (param.lattice_type ==2 && (bandindex ==3 || bandindex ==4))  {
		kz= M_PI;
	}
	
	if (param.lattice_type != 5 && param.lattice_type != 6){
		//std::cout <<" Generating it for square lattice\n";
	if (param.G_FUNC == 0){	
		return std::complex<double>(1,0);
		
	}
	else if (param.G_FUNC == 1){	
		return std::complex<double>(std::cos(momenta[0]) + std::cos(momenta[1]),0);	
	}
	else if (param.G_FUNC == 2){	
		return std::complex<double>(std::sin(momenta[0]),0);	
	}
	else if (param.G_FUNC == 3){	
		return std::complex<double>(std::cos(momenta[0])- std::cos(momenta[1]),0);	
	}
	else if (param.G_FUNC == 30){	
		return std::complex<double>( std::cos(kz)+std::cos(momenta[0])- std::cos(momenta[1]),0);	
	}
	 else if (param.G_FUNC == 4){
                return std::complex<double>((std::sin(momenta[0])* std::sin(momenta[1])),0);
      }
	else if (param.G_FUNC == 5){
                return std::complex<double>((std::sin(momenta[0])* std::sin(momenta[1]))*(std::cos(momenta[0])- std::cos(momenta[1])),0);
       }
	 else if (param.G_FUNC == 6){
		 //std::cout <<"value of kz used is " << kz <<" \n";
		 return std::complex<double>(std::cos(kz)- std::cos(momenta[0]),0);	
	   }
	 else if (param.G_FUNC == 7){
		 //std::cout <<"value of kz used is " << kz <<" \n";
		 return std::complex<double>(std::cos(kz) + std::cos(momenta[0]),0);	
	   }
	else{
		std::cerr <<" Invalid G_FUNC values"<<std::endl;
		return std::complex<double>(1,0);
		}
	}
	else {
		//std::cout <<"Applying symmetry factor for triangular lattice";
		
		if (param.G_FUNC == 0){	
		return std::complex<double>(1,0);	
	}
	if (param.G_FUNC == 1){	
	    double Setd =  (std::cos(momenta[1]) + 2*std::cos(momenta[1]/2)*std::cos(std::sqrt(3)*momenta[0]/2) );
		return std::complex<double>(Setd,0);
		
	}
	if (param.G_FUNC == 2){	
	    double p1 =  std::sqrt(3)*std::sin(std::sqrt(3)*momenta[0]/2)*std::cos(momenta[1]/2);
		double p2 =  std::sin(momenta[1])+std::cos(std::sqrt(3)*momenta[0]/2)*std::sin(momenta[1]/2);
		return std::complex<double>(p1,0);	
	}
	if (param.G_FUNC == 3){	
	    double p1 =  std::sqrt(3)*std::sin(std::sqrt(3)*momenta[0]/2)*std::cos(momenta[1]/2);
		double p2 =  std::sin(momenta[1])+std::cos(std::sqrt(3)*momenta[0]/2)*std::sin(momenta[1]/2);
		return std::complex<double>(p2,0);	
	}
	if (param.G_FUNC == 23){	
	    double p1 =  std::sqrt(3)*std::sin(std::sqrt(3)*momenta[0]/2)*std::cos(momenta[1]/2);
		double p2 =  std::sin(momenta[1])+std::cos(std::sqrt(3)*momenta[0]/2)*std::sin(momenta[1]/2);
		return std::complex<double>(p1,p2)/std::sqrt(2);	
	}
	if (param.G_FUNC == 4){	
	    double dx2y2 =  std::cos(momenta[1]) -std::cos(momenta[1]/2)*std::cos(std::sqrt(3)*momenta[0]/2);
		double dxy =  std::sqrt(3)*std::sin(momenta[1]/2)*std::sin(std::sqrt(3)*momenta[0]/2);
		return std::complex<double>(dx2y2,0);	
	}
	if (param.G_FUNC == 5){	
	    double dx2y2 =  std::cos(momenta[1]) -std::cos(momenta[1]/2)*std::cos(std::sqrt(3)*momenta[0]/2);
		double dxy =  std::sqrt(3)*std::sin(momenta[1]/2)*std::sin(std::sqrt(3)*momenta[0]/2);
		return std::complex<double>(dxy,0);	
	}
	if (param.G_FUNC == 45){	
	    double dx2y2 =  std::cos(momenta[1]) -std::cos(momenta[1]/2)*std::cos(std::sqrt(3)*momenta[0]/2);
		double dxy =  std::sqrt(3)*std::sin(momenta[1]/2)*std::sin(std::sqrt(3)*momenta[0]/2);
		return std::complex<double>(dx2y2,dxy)/std::sqrt(2);	
	}
	if (param.G_FUNC==6){
		double f = std::sin(momenta[1]/2)*(std::cos(std::sqrt(3)*momenta[0]/2)-std::cos(momenta[1]/2));
		return std::complex<double>(f,0);		
		}
	else{
		std::cerr <<" Invalid G_FUNC values"<<std::endl;
		return std::complex<double>(1,0);
		
		}
	}
	
}

std::pair<double, double> mband::generate_hex_bz() {
    std::random_device rd;
    std::default_random_engine engine(rd());
    std::uniform_real_distribution<double> distribution_y(-4 * M_PI / 3, 4 * M_PI / 3);
    std::uniform_real_distribution<double> distribution_x(-2 * M_PI / std::sqrt(3), 2 * M_PI / std::sqrt(3));
    double x;
    double y;
    bool found = false;
    double comp = 2 * M_PI / std::sqrt(3);
 
    do {
        x = distribution_x(engine);
        y = distribution_y(engine);

        if (std::fabs(y) < 2 * M_PI / 3 || (std::fabs(x) + std::sqrt(3) * (std::fabs(y) - 2 * M_PI / 3)) < comp) {
            found = true;
        }

    } while (!found);

    return std::make_pair(x, y);
}
bool mband::check_mfreq_independent(const std::vector<AmiBase::alpha_t>& matrix) {
    for (const auto& row : matrix) {
        if (!row.empty() && row.back() != 0) {
            return false;
        }
    }
    return true;
}

std::vector<AmiBase::epsilon_t> mband::updateEpsilon(const std::vector<AmiBase::epsilon_t>& epsilon, const std::vector<double>& energy) {
    std::vector<std::vector<int>> updatedEpsilon = epsilon;
    std::vector<double> uniqueEnergies;
    for (double e : energy) {
        if (std::find(uniqueEnergies.begin(), uniqueEnergies.end(), e) == uniqueEnergies.end()) {
            uniqueEnergies.push_back(e);
        }
    }

    for (int i = 0; i < energy.size(); i++) {
        double currentEnergy = energy[i];
        std::vector<int> currentEpsilonRow = updatedEpsilon[i];

        for (int j = i + 1; j < energy.size(); j++) {
            if (currentEnergy == energy[j]) {
                updatedEpsilon[j] = currentEpsilonRow;
            }
        }
    }

    return updatedEpsilon;
}


// Function to trim leading and trailing whitespaces from a string
std::string trim(const std::string& str) {
    size_t first = str.find_first_not_of(" \t\n\r\f\v");
    size_t last = str.find_last_not_of(" \t\n\r\f\v");
    if (first == std::string::npos || last == std::string::npos)
        return "";
    return str.substr(first, (last - first + 1));
}

// Function to parse the parameter name and value from a line
void parseLine(const std::string& line, std::string& paramName, std::string& paramValue) {
    size_t equalPos = line.find('=');
    if (equalPos != std::string::npos) {
        paramName = trim(line.substr(0, equalPos));
        paramValue = trim(line.substr(equalPos + 1));
    }
}

// Function to fill the struct from the provided text file
void params_loader(const std::string& filename, mband::params_param& params) {


    // Open the file and read its contents
    std::ifstream inputFile(filename);
    if (!inputFile) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }

    std::string line;
    while (std::getline(inputFile, line)) {
        std::string paramName, paramValue;
        parseLine(line, paramName, paramValue);


        // Fill the struct based on the parameter name
        if (paramName == "molecular")
            params.molecular = std::stoi(paramValue);
        else if (paramName == "E_reg")
            params.E_reg = std::stod(paramValue);
        else if (paramName == "lattice_type")
            params.lattice_type = std::stoi(paramValue);
        else if (paramName == "set_precision")
            params.set_precision = std::stod(paramValue);
        else if (paramName == "hatree_fock")
            params.hatree_fock = (paramValue == "true");
        else if (paramName == "tp")
            params.tp = std::stod(paramValue);
        else if (paramName == "tperp")
            params.tperp = std::stod(paramValue);
		 else if (paramName == "tpp")
            params.tpp = std::stod(paramValue);
		else if (paramName == "tbs")
            params.tbs = std::stod(paramValue);
		 else if (paramName == "molecular_type")
            params.molecular_type = std::stoi(paramValue);
		 else if (paramName == "max_ord")
            params.max_ord = std::stoi(paramValue);
		 else if (paramName == "min_ord")
            params.min_ord = std::stoi(paramValue);
		 else if (paramName == "MC_num")
            params.MC_num= std::stoi(paramValue);
		  else if (paramName == "in")
            params.in= std::stoi(paramValue);
		else if (paramName == "out")
            params.out= std::stoi(paramValue);
		else if (paramName == "molec_beta")
            params.molec_beta= std::stod(paramValue);
		else if (paramName == "molec_mfreq")
            params.molec_mfreq= std::stoi(paramValue);
		else if (paramName == "time")
			params.time= std::stoi(paramValue);
		else if (paramName == "V")
			params.V= std::stod(paramValue);
		else if (paramName == "cutoff")
			params.cutoff_value= std::stod(paramValue);
		else if (paramName == "mfreq_indp")
			params.mfreq_indp = std::stoi(paramValue);
		else if (paramName == "graph_type")
			params.graph_type = std::stoi(paramValue);
		else if (paramName == "graph")
                        params.graph = paramValue;
		else if (paramName == "G_FUNC")
			params.G_FUNC= std::stoi(paramValue);
		else if (paramName == "SOC")
			params.SOC= std::stod(paramValue);
		else if (paramName == "t_orb")
			params.t_orb= std::stod(paramValue);
		  
    }

    inputFile.close();
}




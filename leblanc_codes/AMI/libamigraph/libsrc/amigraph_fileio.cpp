#include "amigraph.hpp"
// #include <boost/filesystem.hpp>
// #include <boost/filesystem/operations.hpp>
// #include <boost/filesystem/path.hpp>

#include<experimental/filesystem>



void AmiGraph::graph_write(std::string filename, graph_t &g){

if( graph_type==AmiBase::density || graph_type==AmiBase::doubleocc || graph_type== AmiBase::DOS|| graph_type== AmiBase::ENERGY){

number_vertices(g);
	
}else{
systematic_vertex_label(g);	
}


std::ofstream file;
file.open(filename);
// std::cout<<"opened file"<<std::endl;

// print_all_edge_info(g);

boost::graph_traits<graph_t>::edge_iterator ei, ei_end;

for(boost::tie(ei,ei_end)=edges(g); ei!=ei_end; ++ei){
file<<g[source(*ei,g)].index_<<" "<<g[target(*ei,g)].index_<<" "<<g[*ei].g_struct_.stat_<<" "<< g[*ei].spin<<std::endl;

}	

file.close();
	
	
}

void AmiGraph::graph_read(std::string filename, graph_t &g){



std::ifstream infile_stream;

std::vector<int> svec, tvec, statvec, spinvec;

infile_stream.open(filename);

if(infile_stream.fail()) // checks to see if file opened 
    { 
	std::cout<<filename<<std::endl;
      throw std::runtime_error("Could not open input file");
    } 	
	
std::string line;
// std::getline(infile_stream,line);	

while (std::getline(infile_stream, line))
{

	std::stringstream ss(line);
	int source, target, stat,spin;
	
    std::string laststring;
	bool read = bool(ss >> laststring);
	if(read){
	source=std::stoi(laststring);
	ss >>  target >> stat>>spin;// >> kx >> ky >> realW>> imagW;
	
	svec.push_back(source);
	tvec.push_back(target);
	statvec.push_back(stat);
	spinvec.push_back(spin);
	
	// std::cout<<"s t stat spin were "<< source <<" "<<target <<" "<<stat<<" "<<spin<<std::endl;
	
	}
}

infile_stream.close();

// std::cout<<"File had num lines="<<svec.size()<<std::endl;

int val=0;

for(int i=0; i< svec.size(); i++){
if(svec[i]>val){val=svec[i];}
if(tvec[i]>val){val=tvec[i];}	
	
}

int n= val+1;
graph_t loaded(n);

boost::graph_traits<graph_t>::vertex_iterator vi, vi_end;
int id=0;
  for (boost::tie(vi, vi_end) = vertices(loaded); vi != vi_end; ++vi, ++id) {
    loaded[*vi].index_=id;
  }



for(int edge=0; edge< svec.size(); edge++){
	
boost::graph_traits<graph_t>::vertex_descriptor source,target;

boost::graph_traits<graph_t>::vertex_iterator vi, vi_end;
  for (boost::tie(vi, vi_end) = vertices(loaded); vi != vi_end; ++vi, ++id) {
    if(loaded[*vi].index_==svec[edge]){ source=*vi;}
	if(loaded[*vi].index_==tvec[edge]){ target=*vi;}
  }	
  
  add_edge(source,target, edge_info(AmiBase::stat_type(statvec[edge]), -1, spin_type(spinvec[edge])), loaded);
	
	
	
}


g=loaded;
	
	
}



void AmiGraph::pair_read(std::string filename, git_pair &p){

	std::ifstream infile_stream;

	std::vector<int> svec, tvec, statvec, spinvec;
	std::vector<int> pairID;
	std::vector<AmiBase::epsilon_t> epsvec;
	std::vector<AmiBase::alpha_t> alphavec;
	AmiBase::epsilon_t eps;
	AmiBase::alpha_t alpha;
	
	// graph_t g1, g2;
	git_perm_set_t pst;
	git_perm_t temp_pst;
	// git_pair temp_pair;

	infile_stream.open(filename);

	if(infile_stream.fail()) // checks to see if file opened 
		{ 
		std::cout<<filename<<std::endl;
		  throw std::runtime_error("Could not open input file");
		} 	
		
	std::string line;
	// std::getline(infile_stream,line);	

	while (std::getline(infile_stream, line))
	{

		std::stringstream ss(line);
		int source, target, stat,spin;
		
		std::string number;
		std::string left("{");
		std::string right("}");
		
		std::string laststring;
		bool read = bool(ss >> laststring);
		if(read){
			pairID.push_back(stoi(laststring));
		// source=std::stoi(laststring);
		// std::cout<<"Last string was "<<laststring<<std::endl;
				if(std::stoi(laststring)<2){
				ss >> source >> target >> stat>>spin;// >> kx >> ky >> realW>> imagW;
				
				svec.push_back(source);
				tvec.push_back(target);
				statvec.push_back(stat);
				spinvec.push_back(spin);
				
				int par=0;
				
				// ss>> laststring; // this is parenthesis
				
				
				
				while(ss>>number){
				
	
				if(number==left || number==right){ par++;}
				else{
				if(par<2){eps.push_back(std::stoi(number));}
				if(par>1){alpha.push_back(std::stoi(number));}		
					
				}
				
					
					// std::cout<<"par is "<<par<<" entry is "<< number<<std::endl;
				// if(is_number(number)){
				
				// if(par<2){
					// std::cout<<"par was "<<par<<" pushing back "<<number<<std::endl;
					// eps.push_back(std::stoi(number));
				// }
				// if(par>1 && par<4){
					// std::cout<<"par was "<<par<<" pushing back "<<number<<std::endl;
					// alpha.push_back(std::stoi(number));
				// }

				
				// }else{par++;}
				
				
				}
				// std::cout<<"s t stat spin were "<< source <<" "<<target <<" "<<stat<<" "<<spin<<std::endl;
				epsvec.push_back(eps);
				alphavec.push_back(alpha);
				eps.clear();
				alpha.clear();
				}
				else{
					
				// temp_pst.push_back(std::stoi(laststring));
				while(ss>>number){
					// std::cout<<number<<" "<<std::endl;
				temp_pst.push_back(std::stoi(number));	
				}	
				pst.push_back(temp_pst);
				temp_pst.clear();
				
				}
				
		}
		
	}
	
	// at this stage only the pst is loaded 
	p.pst_=pst;
	//

	infile_stream.close();

	// std::cout<<"File had num lines="<<svec.size()<<std::endl;

	int val=0;

	for(int i=0; i< svec.size(); i++){
	if(svec[i]>val){val=svec[i];}
	if(tvec[i]>val){val=tvec[i];}	
		
	}

	int n= val+1;
	graph_t loaded(n);
	graph_t g2(n);

	boost::graph_traits<graph_t>::vertex_iterator vi, vi_end;
	int id=0;
	  for (boost::tie(vi, vi_end) = vertices(loaded); vi != vi_end; ++vi, ++id) {
		loaded[*vi].index_=id;
	  }
	  
	 id=0;
	  for (boost::tie(vi, vi_end) = vertices(g2); vi != vi_end; ++vi, ++id) {
		g2[*vi].index_=id;
	  }



	for(int edge=0; edge< svec.size(); edge++){
	
	if(pairID[edge]==0){
	boost::graph_traits<graph_t>::vertex_descriptor source,target;

	boost::graph_traits<graph_t>::vertex_iterator vi, vi_end;
	  for (boost::tie(vi, vi_end) = vertices(loaded); vi != vi_end; ++vi, ++id) {
		if(loaded[*vi].index_==svec[edge]){ source=*vi;}
		if(loaded[*vi].index_==tvec[edge]){ target=*vi;}
	  }	
	  
	  add_edge(source,target, edge_info(epsvec[edge], alphavec[edge], AmiBase::stat_type(statvec[edge]), -1, spin_type(spinvec[edge])), loaded);
	}
	if(pairID[edge]==1){
	boost::graph_traits<graph_t>::vertex_descriptor source,target;

	boost::graph_traits<graph_t>::vertex_iterator vi, vi_end;
	  for (boost::tie(vi, vi_end) = vertices(g2); vi != vi_end; ++vi, ++id) {
		if(g2[*vi].index_==svec[edge]){ source=*vi;}
		if(g2[*vi].index_==tvec[edge]){ target=*vi;}
	  }	
	  
	  add_edge(source,target, edge_info(epsvec[edge], alphavec[edge], AmiBase::stat_type(statvec[edge]), -1, spin_type(spinvec[edge])), g2);
	}
		
		
	}


	p.g1_=loaded;
	p.g2_=g2;
		
		
	}

void AmiGraph::read_ggm(std::string folder, gg_matrix_t &ggm, int max_ord){

ggm.clear();
// oX_gX_nX.graph 
    std::stringstream ss;
	ss<<std::experimental::filesystem::current_path().string()<<"/"<<folder;
	std::string path=ss.str();
	
	if(!std::experimental::filesystem::is_directory(path))
	{
	throw std::runtime_error("Could not open ggm directory");
}
	
	for (const auto & entry : std::experimental::filesystem::directory_iterator(path)){
    std::cout << entry.path() << std::endl;
	// std::cout << entry.path().extension() << std::endl;
	if(entry.path().extension().string().compare(".graph")){

	std::string stem=entry.path().stem().string();
	
	
	std::size_t pos=stem.find("_");
	std::string second=stem.substr(1, pos-std::size_t(1));
	int ord=std::stoi(second);
	stem=stem.substr(pos+std::size_t(1));
	
	pos=stem.find("_");
	second=stem.substr(1,pos-std::size_t(1));
	int group=std::stoi(second);
	
	second=stem.substr(pos+std::size_t(2));
	int num=std::stoi(second);
	
	if(ord>max_ord){continue;}
	

// std::cout<<stem <<" "<<second<<" "<<pos<<std::endl;

	 // std::cout << entry.path().filename() << std::endl;
	  // std::cout<<"Ord is "<<ord<<" group is "<< group<<" num is "<< num <<std::endl;
	  
	if(ggm.size()<ord+1){ ggm.resize(ord+1);}
	if(ggm[ord].size()< group+1){ggm[ord].resize(group+1);}
	if(ggm[ord][group].graph_vec.size()<num+1){ggm[ord][group].graph_vec.resize(num+1);}
	
	graph_t temp;
	graph_read(entry.path().string(),temp);
	
	ggm[ord][group].graph_vec[num]=temp;

	
	  
	}
	
	}
	
	
// for(int ord=2; ord< ggm.size(); ord++){
	// for(int group=0; group< ggm[ord].size(); group++){
		// for(int graph=0; graph< ggm[ord][group].graph_vec.size(); graph++)
	
	
	
	
	
}


void AmiGraph::read_ggmp(std::string folder, gg_matrix_t &ggm, int max_ord){

ggm.clear();
// oX_gX_nX.graph 
    std::stringstream ss;
	ss<<std::experimental::filesystem::current_path().string()<<"/"<<folder;
	std::string path=ss.str();
	
	if(!std::experimental::filesystem::is_directory(path))
	{
	throw std::runtime_error("Could not open ggm directory");
}
	
	for (const auto & entry : std::experimental::filesystem::directory_iterator(path)){
		
		
    // std::cout << entry.path() << std::endl;
	// std::cout << entry.path().extension() <<" "<<entry.path().extension().string()<<  std::endl;
	
			
				if(entry.path().extension().string().compare(".graph")==0){
				// std::cout<<"Was true group? "<<entry.path().extension().string().compare(".group")<<std::endl;

			std::string stem=entry.path().stem().string();
			
			
			std::size_t pos=stem.find("_");
			std::string second=stem.substr(1, pos-std::size_t(1));
			int ord=std::stoi(second);
			stem=stem.substr(pos+std::size_t(1));
			
			pos=stem.find("_");
			second=stem.substr(1,pos-std::size_t(1));
			int group=std::stoi(second);
			
			second=stem.substr(pos+std::size_t(2));
			int num=std::stoi(second);
			
			if(ord>max_ord){continue;}
			

		// std::cout<<stem <<" "<<second<<" "<<pos<<std::endl;

			 // std::cout << entry.path().filename() << std::endl;
			  // std::cout<<"Ord is "<<ord<<" group is "<< group<<" num is "<< num <<std::endl;
			  
			if(ggm.size()<ord+1){ ggm.resize(ord+1);}
			if(ggm[ord].size()< group+1){ggm[ord].resize(group+1);}
			if(ggm[ord][group].graph_vec.size()<num+1){ggm[ord][group].graph_vec.resize(num+1);}
			
			graph_t temp;
			graph_read(entry.path().string(),temp);
			
			ggm[ord][group].graph_vec[num]=temp;

			
			  
			}

		// std::strcmp(entry.path().extension().string(),".group")
	// if(entry.path().extension().string().compare("pair")==0){
		// if(entry.path().extension().string()==".pair"){
		if(entry.path().extension().string().compare(".pair")==0){
		// std::cout<<"Was true pair? "<<entry.path().extension().string().compare(".group")<<std::endl;

	std::string stem=entry.path().stem().string();
	
	
	std::size_t pos=stem.find("_");
	std::string second=stem.substr(1, pos-std::size_t(1));
	int ord=std::stoi(second);
	stem=stem.substr(pos+std::size_t(1));
	
	pos=stem.find("_");
	second=stem.substr(1,pos-std::size_t(1));
	int group=std::stoi(second);
	
	second=stem.substr(pos+std::size_t(2));
	int num=std::stoi(second);
	
	if(ord>max_ord){continue;}
	

// std::cout<<stem <<" "<<second<<" "<<pos<<std::endl;

	 // std::cout << entry.path().filename() << std::endl;
	  // std::cout<<"Ord is "<<ord<<" group is "<< group<<" num is "<< num <<std::endl;
	  
	if(ggm.size()<ord+1){ ggm.resize(ord+1);}
	if(ggm[ord].size()< group+1){ggm[ord].resize(group+1);}
	if(ggm[ord][group].gp_vec.size()<num+1){ggm[ord][group].gp_vec.resize(num+1);}
	
	git_pair temppair;
	pair_read(entry.path().string(),temppair);
	
	ggm[ord][group].gp_vec[num]=temppair;

	
	  
	}
			
			
			
	
	}
	
	
// for(int ord=2; ord< ggm.size(); ord++){
	// for(int group=0; group< ggm[ord].size(); group++){
		// for(int graph=0; graph< ggm[ord][group].graph_vec.size(); graph++)
	
	
	
	
	
}


void AmiGraph::read_ggmp(std::string folder, gg_matrix_t &ggm, int min_ord, int max_ord){

ggm.clear();
// oX_gX_nX.graph 
    std::stringstream ss;
	ss<<std::experimental::filesystem::current_path().string()<<"/"<<folder;
	std::string path=ss.str();
	
	if(!std::experimental::filesystem::is_directory(path))
	{
	throw std::runtime_error("Could not open ggm directory");
}
	
	for (const auto & entry : std::experimental::filesystem::directory_iterator(path)){
		
		
    // std::cout << entry.path() << std::endl;
	// std::cout << entry.path().extension() <<" "<<entry.path().extension().string()<<  std::endl;
	
			
				if(entry.path().extension().string().compare(".graph")==0){
				// std::cout<<"Was true group? "<<entry.path().extension().string().compare(".group")<<std::endl;

			std::string stem=entry.path().stem().string();
			
			
			std::size_t pos=stem.find("_");
			std::string second=stem.substr(1, pos-std::size_t(1));
			int ord=std::stoi(second);
			stem=stem.substr(pos+std::size_t(1));
			
			pos=stem.find("_");
			second=stem.substr(1,pos-std::size_t(1));
			int group=std::stoi(second);
			
			second=stem.substr(pos+std::size_t(2));
			int num=std::stoi(second);
			
			if(ord>max_ord){continue;}
			

		// std::cout<<stem <<" "<<second<<" "<<pos<<std::endl;

			 // std::cout << entry.path().filename() << std::endl;
			  // std::cout<<"Ord is "<<ord<<" group is "<< group<<" num is "<< num <<std::endl;
			  
			if(ggm.size()<ord+1){ ggm.resize(ord+1);}
			if(ggm[ord].size()< group+1){ggm[ord].resize(group+1);}
			if(ggm[ord][group].graph_vec.size()<num+1){ggm[ord][group].graph_vec.resize(num+1);}
			
			graph_t temp;
			graph_read(entry.path().string(),temp);
			
			ggm[ord][group].graph_vec[num]=temp;

			
			  
			}

		// std::strcmp(entry.path().extension().string(),".group")
	// if(entry.path().extension().string().compare("pair")==0){
		// if(entry.path().extension().string()==".pair"){
		if(entry.path().extension().string().compare(".pair")==0){
		// std::cout<<"Was true pair? "<<entry.path().extension().string().compare(".group")<<std::endl;

	std::string stem=entry.path().stem().string();
	
	
	std::size_t pos=stem.find("_");
	std::string second=stem.substr(1, pos-std::size_t(1));
	int ord=std::stoi(second);
	stem=stem.substr(pos+std::size_t(1));
	
	pos=stem.find("_");
	second=stem.substr(1,pos-std::size_t(1));
	int group=std::stoi(second);
	
	second=stem.substr(pos+std::size_t(2));
	int num=std::stoi(second);
	
	if(ord>max_ord ){continue;}
	if(ord< min_ord){continue;}

// std::cout<<stem <<" "<<second<<" "<<pos<<std::endl;

	 // std::cout << entry.path().filename() << std::endl;
	  // std::cout<<"Ord is "<<ord<<" group is "<< group<<" num is "<< num <<std::endl;
	  
	if(ggm.size()<ord+1){ ggm.resize(ord+1);}
	if(ggm[ord].size()< group+1){ggm[ord].resize(group+1);}
	if(ggm[ord][group].gp_vec.size()<num+1){ggm[ord][group].gp_vec.resize(num+1);}
	
	git_pair temppair;
	pair_read(entry.path().string(),temppair);
	
	ggm[ord][group].gp_vec[num]=temppair;

	
	  
	}
			
			
			
	
	}
	
	
// for(int ord=2; ord< ggm.size(); ord++){
	// for(int group=0; group< ggm[ord].size(); group++){
		// for(int graph=0; graph< ggm[ord][group].graph_vec.size(); graph++)
	
	
	
	
	
}

void AmiGraph::gp_write(std::string filename, git_pair &p){

systematic_vertex_label(p.g1_);
systematic_vertex_label(p.g2_);	

std::ofstream file;
file.open(filename);

boost::graph_traits<graph_t>::edge_iterator ei, ei_end;

for(boost::tie(ei,ei_end)=edges(p.g1_); ei!=ei_end; ++ei){
file<<"0 "<<p.g1_[source(*ei,p.g1_)].index_<<" "<<p.g1_[target(*ei,p.g1_)].index_<<" "<<p.g1_[*ei].g_struct_.stat_<<" "<< p.g1_[*ei].spin<<" ";

file<<"{ ";
for(int i=0; i< p.g1_[*ei].g_struct_.eps_.size(); i++){
file<<	p.g1_[*ei].g_struct_.eps_[i]<<" ";
}
file<<"} ";

file<<"{ ";
for(int i=0; i< p.g1_[*ei].g_struct_.alpha_.size(); i++){
file<<	p.g1_[*ei].g_struct_.alpha_[i]<<" ";
}
file<<"} ";


file<<std::endl;

}	

for(boost::tie(ei,ei_end)=edges(p.g2_); ei!=ei_end; ++ei){
file<<"1 "<<p.g2_[source(*ei,p.g2_)].index_<<" "<<p.g2_[target(*ei,p.g2_)].index_<<" "<<p.g2_[*ei].g_struct_.stat_<<" "<< p.g2_[*ei].spin<<" ";

file<<"{ ";
for(int i=0; i< p.g2_[*ei].g_struct_.eps_.size(); i++){
file<<	p.g2_[*ei].g_struct_.eps_[i]<<" ";
}
file<<"} ";

file<<"{ ";
for(int i=0; i< p.g2_[*ei].g_struct_.alpha_.size(); i++){
file<<	p.g2_[*ei].g_struct_.alpha_[i]<<" ";
}
file<<"} ";


file<<std::endl;
}

for(int i=0; i< p.pst_.size(); i++){
	file<<"2 ";
for(int j=0; j<p.pst_[i].size(); j++){
 file<<p.pst_[i][j]<<" ";
}
file<<std::endl;		

}


file.close();
	
	
}

void AmiGraph::write_ggmp(std::string folder, gg_matrix_t &ggm, int min){
	
	// std::experimental::filesystem::path full_path(std::experimental::filesystem::current_path());
    // std::cout << "Current path is : " << full_path << std::endl;
	std::stringstream ss;
	ss<<std::experimental::filesystem::current_path().string()<<"/"<<folder;
	std::string path=ss.str();
	
	std::cout<< path<<std::endl;
	
	 std::experimental::filesystem::path dir(path);

    if(!(std::experimental::filesystem::exists(path))){
        std::cout<<"GGM writing folder  does not exist"<<std::endl;

		std::experimental::filesystem::create_directory(path);
        // if (std::experimental::filesystem::create_directory(dir)){
		// std::cout << "....Successfully Created Directory" << std::end;}
		
		
    }else{
		std::cout<<"GGM writing folder exists"<<std::endl;
		return;
	}


for(int ord=min; ord< ggm.size(); ord++){
	for(int group=0; group< ggm[ord].size(); group++){
	for(int graph=0; graph< ggm[ord][group].graph_vec.size(); graph++){
	
	std::ofstream file;
    std::stringstream filename;
    filename<<folder<<"/o"<<ord<<"_g"<<group<<"_n"<<graph<<".graph";

// file.open(filename.str());
graph_write(filename.str(), ggm[ord][group].graph_vec[graph]);
	
	
	
	}
	
	for(int pair=0; pair< ggm[ord][group].gp_vec.size(); pair++){
		
	std::ofstream file;
	std::stringstream filename;
	filename<<folder<<"/o"<<ord<<"_g"<<group<<"_n"<<pair<<".pair";
	
	gp_write(filename.str(), ggm[ord][group].gp_vec[pair]);	
		
	}
	
	
	}
}	
	
	
	
	
	
}



void AmiGraph::write_ggm(std::string folder, gg_matrix_t &ggm){
	
	// std::experimental::filesystem::path full_path(std::experimental::filesystem::current_path());
    // std::cout << "Current path is : " << full_path << std::endl;
	std::stringstream ss;
	ss<<std::experimental::filesystem::current_path().string()<<"/"<<folder;
	std::string path=ss.str();
	
	std::cout<< path<<std::endl;
	
	 std::experimental::filesystem::path dir(path);

    if(!(std::experimental::filesystem::exists(path))){
        std::cout<<"GGM writing folder  does not exist"<<std::endl;

		std::experimental::filesystem::create_directory(path);
        // if (std::experimental::filesystem::create_directory(dir)){
		// std::cout << "....Successfully Created Directory" << std::end;}
		
		
    }else{
		std::cout<<"GGM writing folder exists"<<std::endl;
		return;
	}

std::cout<<"Writing into directory"<<std::endl;
for(int ord=0; ord< ggm.size(); ord++){
	std::cout<<"On ord "<<ord<<std::endl;
	for(int group=0; group< ggm[ord].size(); group++){
	for(int graph=0; graph< ggm[ord][group].graph_vec.size(); graph++){
	
	std::ofstream file;
    std::stringstream filename;
    filename<<folder<<"/o"<<ord<<"_g"<<group<<"_n"<<graph<<".graph";

// file.open(filename.str());
graph_write(filename.str(), ggm[ord][group].graph_vec[graph]);
	
	
	
	}
	}
}

	
	
	
	
}



#include "mpi.h"  // technically doesn't do anything with mpi... so... unsure about this. Moved bc functions into libamigraph so it is now useful. 

#include <iostream>
#include <fstream>
#include <string> 
#include <sstream> 
#include<Eigen/Dense>
#include<Eigen/Core>

#include <complex>
#include <math.h>
#include <stdlib.h>
#include <random>
#include <chrono>
#include <vector>
#include <algorithm>
#include<experimental/filesystem>

// #include "ami.hpp"

#include "ami_base.hpp"
#include "ami_calc.hpp"
#include "ami_spec.hpp"

// Boost headers 
#include <boost/config.hpp>
#include <boost/graph/isomorphism.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/copy.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/property_map/property_map.hpp>

// #include <boost/program_options.hpp>
#include "boost/graph/graphviz.hpp"
#include <boost/graph/random.hpp>
#include <boost/random.hpp>
//#include <boost/graph/depth_first_search.hpp>
#include <boost/range/irange.hpp>
//#include <boost/pending/indirect_cmp.hpp>
//#include <boost/graph/undirected_dfs.hpp>
//#include <boost/cstdlib.hpp>
#include <boost/math/special_functions/factorials.hpp>

#include <boost/random/sobol.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/variate_generator.hpp>

// #define BOOST_NO_CXX11_SCOPED_ENUMS
// #include <boost/filesystem.hpp>
// #undef BOOST_NO_CXX11_SCOPED_ENUMS


void triangularize_matrix( std::vector< std::vector< int >> M_IN, std::vector< std::vector< int >> &M_out);






class AmiGraph
{



private:



public:

// Helper functions

bool is_number(const std::string& s);
bool bose_alphas_in_R0=false;

/// Simple factorial function - Nothing special. 
int factorial(int n);

void combinations(int max, int length, std::vector< std::vector<int>> &list);// n-choose m is max-choose length. stores indices 0 to max-1

void combinations_repetition_util(std::vector<int> &chosen, std::vector< std::vector<int>> &list, int index, int r, int start, int end);// allows repeated values in combinations 
void combinations_repetition(int n, int r, std::vector< std::vector<int>> &list );

// Instantiate an ami calculation for the graph
NewAmiCalc ami;
AmiSpec sp;

AmiBase::ami_parms ami_parameters;
NewAmiCalc::internal_state current_state;  // this isn't initialized anywhere to have dimensionality...


// Create a random number generator using std::random

    std::mt19937 rand_gen;
    std::uniform_real_distribution<double> rand_dist;
	double random_real(double max);
	double random_real(double min, double max);
    int random_int( int min, int max);  // this is a small function to get a random integer.  For real numbers use auto roll_dice example below.
	void randomize(std::vector<int> &v, int min, int max); // randomizes all elements of a vector with ints 
//    std::uniform_int_distribution<int> int_dist;
    // static int myrandom (int i); // small function for std::random_shuffle

    static const std::size_t dimension;
	 // Initialize the engine to draw randomness out of thin air
    boost::random::sobol engine; //(dimension);
	std::vector< boost::random::sobol > engine_vec;
	typedef boost::variate_generator<boost::random::sobol&, boost::uniform_01<double> > quasi_random_gen_t;

    
	
	
	// quasi-rand generators
	std::vector<double> sobol_random_real();
	std::vector<double> sobol_random_real(int size, int order);
	// double AmiGraph::sobol_random_real(double max);
	
	
    // Create a generator
   // typedef boost::variate_generator<boost::random::sobol64&, boost::uniform_01<double> > quasi_random_gen_t;

   

    // Glue the engine and the distribution together
    // quasi_random_gen_t gen(engine, boost::uniform_01<double>());

    // std::vector<double> sample(dimension);

    // At this point you can use std::generate, generate member f-n, etc.
    // std::generate(sample.begin(), sample.end(), gen);
    // engine.generate(sample.begin(), sample.end());


//    auto roll_dice=std::bind ( rand_dist, rand_gen);


enum spin_type{up,dn};
enum bubble_type{directed, both_in, both_out};
enum label_type{labelled, unlabelled};
//enum stat_type{bose,fermi}; //This may be part of the ami class

typedef int matsubara_freq_type;

enum vertex_type{three_leg, four_leg };

typedef std::vector< std::vector< std::vector <std::vector <double> > > > U_matrix_type;
void resize_U_matrix(int N, U_matrix_type &U);
void read_U_matrix(std::string filename, U_matrix_type &UM, int maxval);

typedef std::vector< std::vector< int > > Ulist_type;

typedef std::vector< std::pair<int,int> > band_map;





////////////////////////
// Vertex info structure
////////////////////////

//typedef boost::property< vertex_index_t, int> vertex_index;

struct vertex_info {
  vertex_info(vertex_type type){
    type_=type;
	visited=0;

if (type_==four_leg){  throw std::runtime_error("vertex type is 4 leg but no U_value given.");}

    }

vertex_info(vertex_type type, double U_value){
    type_=type;
    U_value_=U_value;
	visited=0;
    }
//  vertex_info(){ throw std::runtime_error("type has to be default constructible but we should never use the default constructor.");}
  vertex_info(){
type_=three_leg; // assume it is a 3-leg vertex if no type specified
visited=0;
}
  
  vertex_type type_;
  double U_value_;
  int index_;
//  size_t vertex_index;
  //int component;
  int visited;
  boost::property< boost::vertex_index_t, int> vertex_index;
  
};

/////////////////////////

////////////////////////
// Vertex info structure
////////////////////////

// This is a mirror of the g_struct from ami
struct edge_info {

AmiBase::g_struct g_struct_;
// g_struct_ includes 

// epsilon_t eps_;
// alpha_t alpha_;
// stat_type stat_;

// where these are defined as 
// typedef std::vector<int> epsilon_t;  /// the symbolic epsilon - not numeric - hence integers
// typedef std::vector<int> alpha_t;
// typedef enum {Bose,Fermi} stat_type ;



// DEBUGGING ONLY
int edge_number_;
//AmiBase::stat_type edge_stat=
int fermi_loop_id=-1;
bool bubble=false;
bool tadpole=false;
spin_type spin=up;
int band=0;
std::pair<int, int> band_element;
label_type label;

std::vector<int> fourindex;

  edge_info(){ 
label=unlabelled;
}

// todo fix this - need to remove loop id
edge_info(int bnd, AmiBase::stat_type type ){
band=bnd;
g_struct_.stat_=type;	
label=unlabelled;
}

edge_info(AmiBase::stat_type type, std::pair<int,int> be){
g_struct_.stat_=type;
band_element=be; // this is basically the original source-target labelling	
label=unlabelled;
}

edge_info(AmiBase::stat_type type, std::vector<int> four){

g_struct_.stat_=type;
fourindex=four;
label=unlabelled;	
}

edge_info(AmiBase::stat_type type, int loop_id, spin_type spin_val){ 
g_struct_.stat_=type;
fermi_loop_id=loop_id;
spin=spin_val;
label=unlabelled;

// species info to g_struct
// g_struct_.species_=(species_t)spin_val;
// in theory R0 now knows the list of species 

}

edge_info(AmiBase::epsilon_t epsilon, AmiBase::alpha_t alpha, AmiBase::stat_type type, int loop_id, spin_type spin_val){
g_struct_.eps_=epsilon;
g_struct_.alpha_=alpha;
g_struct_.stat_=type;
fermi_loop_id=loop_id;
spin=spin_val;
label=unlabelled;
  }

edge_info(AmiBase::stat_type type, int loop_id){ 
g_struct_.stat_=type;
fermi_loop_id=loop_id;
label=unlabelled;
}


edge_info(AmiBase::stat_type type){ 
g_struct_.stat_=type;
label=unlabelled;
}
  edge_info(AmiBase::epsilon_t epsilon, AmiBase::alpha_t alpha, AmiBase::stat_type type){
g_struct_.eps_=epsilon;
g_struct_.alpha_=alpha;
g_struct_.stat_=type;
label=unlabelled;
  }

 edge_info(AmiBase::epsilon_t epsilon, AmiBase::alpha_t alpha){
g_struct_.eps_=epsilon;
g_struct_.alpha_=alpha;
g_struct_.stat_=static_cast<AmiBase::stat_type>(1); // 1 here means fermi;
label=unlabelled;
  }
  	


}; // end edge_struct bracket


////////////////////////
// graph_info structure
////////////////////////
struct graph_info {
  
  int fermi_loop_id;
  int num_loops;
  int n_indep;
  int n_labelled;
  int ext_counts;
  bool is_bose;
  int ct_count=0;
	int sigma_ct_count=0;

     
};



//   A ->--- a =====b--->--B .  A and B are external vertices. a and b are the bubble vertices 

struct bubble{

std::vector<std::pair<int,int>> bubble_edges;
std::vector<std::pair<int,int>> leg_edges;	

int A,a,b,B;

	
};







// One line definition of our graph type as an undirected adjacency_list
// Open to discussion if an alternate graph type makes more sense
// Reason for this type: A directed graph cannot be traversed against its direction.
// Underlying data_structure for the graph is just a vector (vecS)
//typedef boost::adjacency_list < boost::vecS, boost::vecS, boost::undirectedS,vertex_info, edge_info, graph_info > graph_t;
//typedef boost::adjacency_list < boost::vecS, boost::vecS, boost::bidirectionalS,vertex_info, edge_info, graph_info > graph_t;

typedef boost::adjacency_list < boost::vecS, boost::listS, boost::bidirectionalS,vertex_info, edge_info, graph_info > graph_t;

 typedef boost::adjacency_list<boost::vecS, boost::listS, boost::bidirectionalS,boost::property<boost::vertex_index_t, int> > iso_graph_t;
 
 
 // Edge weight.
typedef boost::property<boost::edge_weight_t, int> EdgeWeightProperty;
//
 typedef boost::adjacency_list<boost::vecS, boost::listS, boost::bidirectionalS,boost::property<boost::vertex_index_t, int>, EdgeWeightProperty> spin_iso_graph_t;
 
 ///
 
 typedef
  boost::adjacency_list<
    boost::vecS            // edge list
  , boost::vecS            // vertex list
  , boost::undirectedS     // directedness
  , float                  // property associated with vertices
  >
cc_graph_t;
 
// vertex and edge descriptors are just indexes since we use vector (vecS)
typedef boost::graph_traits < graph_t >::vertex_descriptor vertex_t;
typedef boost::graph_traits < graph_t >::edge_descriptor edge_t;

typedef std::vector< edge_t> edge_vector_t;
typedef std::vector< vertex_t> vertex_vector_t;

typedef std::vector< edge_vector_t> edge_vector_list_t;
typedef std::vector< vertex_vector_t> vertex_vector_list_t;

using index_map_t = boost::property_map<graph_t, int vertex_info::*>::type;


// isomorphism

typedef vertex_t grandparent_vertex_t;
typedef vertex_t parent_vertex_t;
typedef std::vector< vertex_t > child_vertex_vector_t;

struct pc_vert_struct{
	pc_vert_struct(){}
	grandparent_vertex_t grandpa;
	parent_vertex_t parent;
	child_vertex_vector_t child;
};

typedef std::vector<std::vector< pc_vert_struct >> tree_t;
typedef std::vector< pc_vert_struct > tree_leg_t;

typedef std::vector<int> child_t;
typedef std::vector<int> parent_t;
typedef std::vector<int> grandparent_t;

struct pc_struct {
	pc_struct(){}
	grandparent_t grandparent;
parent_t parent;
std::vector<child_t> child;	
};

typedef std::vector< std::vector< pc_struct>> barcode_t; 
typedef std::vector< pc_struct> barcode_line_t;


typedef std::vector< std::pair<int,int>> pp_c_t;
typedef std::vector< std::pair<int,int>> pp_p_t;

struct pp_pc{
	pp_pc(){}
	pp_c_t child;
	pp_p_t parent;
};
typedef std::vector< pp_pc > pp_barcode_line_t;
typedef std::vector<pp_barcode_line_t> pp_barcode_t;

void pp_barcode(graph_t &g, pp_barcode_t &barcode);
void get_pp_line(vertex_vector_t &vv, vertex_vector_t &next_vv, graph_t &g, pp_barcode_t &barcode);
bool compare_pp_barcodes(pp_barcode_t &B1, pp_barcode_t &B2);
bool compare_bc_line(pp_barcode_line_t bl1, pp_barcode_line_t bl2);
bool pc_equal(pp_pc pc1, pp_pc pc2);
void print_bc_line(pp_barcode_line_t b);
void print_pc(pp_pc pc);

void symmetrize_spin_up(graph_t &g);
///
//typedef std::map<vertex_t, size_t> VertexDescMap; 

// We can now have as many graph objects as we need.

graph_t current_graph;
graph_t f2,f3;
graph_t proposed_graph;
// graph_t starter_graph; // Don't think we need this. but maybe we do.

// Global properties of a graph??? should be a part of the graph definition 

// Temporary stuff that will be removed later

// IO graph functions

void print_all_vertex_info();
void print_all_edge_info();
void print_all_edge_info(graph_t &g);
void print_edge_info(edge_t &e, graph_t &g);


//
////////////////////
// Constructor and member functions 
////////////////////
  AmiGraph();
  AmiGraph(AmiBase::graph_type type, int seed);
  AmiGraph(AmiBase::graph_type type, int dim, int seed);


////////////////////
// FUNCTIONS
//////////////////

void create_starter_graph(graph_t &g);
// Starter graph functions
void construct_starter_sigma(graph_t &g);
void construct_starter_phuu_bubble(graph_t &g);
void construct_starter_phud_bubble(graph_t &g);

void construct_starter_force(graph_t &g, graph_t &other, graph_t &otherp);

void construct_starter_ppuu_bubble(graph_t &g);
void construct_starter_ppud_bubble(graph_t &g);

void construct_starter_hartree(AmiGraph::graph_t &g);
void construct_starter_bare(graph_t &g);

void initialize(AmiBase::graph_type type);
void initialize(AmiBase::graph_type type, int dim);

AmiBase::graph_type graph_type;

/////////
// Update Functions
/////////
void AIL_rand(graph_t &g);
void DIL_rand(graph_t &g);

void AB_rand(graph_t &g);
void DB_rand(graph_t &g);

void AT_rand(graph_t &g);
void DT_rand(graph_t &g);

void RBI_rand(graph_t &g);
void RIB_rand(graph_t &g);

void SN_rand(graph_t &g);

/////////
// Detail Balance/acceptance Functions
/////////








// Find bosonic edges.
// helper functions

int graph_order(graph_t &g);

void number_vertices(graph_t &g);
std::string print_edge(edge_t &e, graph_t &g);

void print_edgeson_vert(int index);

void find_bosonic_edge_on_three_pointvert(graph_t &g, vertex_t vert, edge_t &bose_edge);

void find_bosonic_edges(graph_t &g, edge_vector_t &vector);
void find_non_tadpole_bosonic_edges(graph_t &g, edge_vector_t &vector); 
void find_non_entry_bosonic_edges(graph_t &g, edge_vector_t &vector);
void find_non_external_bosonic_edges(graph_t &g, edge_vector_t &vector);
void find_mutable_bosonic_edges(graph_t &g, edge_vector_t &vector);
void find_force_LR_vertices(graph_t &g, vertex_vector_t &in_vv, vertex_vector_t &out_vv);


void find_one_fermi_bose_edge(graph_t &g, edge_t &bose, edge_t &fermi);
void find_bose_fermi_edges(graph_t &g, edge_vector_t &bose, edge_vector_t &fermi);

void find_fermionic_edges(graph_t &g, edge_vector_t &vector);
void find_internal_fermionic_edges(graph_t &g, edge_vector_t &vector);
void find_internal_edges_stat(graph_t &g, edge_vector_t &vector, AmiBase::stat_type requested);
void find_internal_vertices(graph_t &g, vertex_vector_t &vector);

void find_two_unique_edges(graph_t &g, AmiBase::stat_type requested_type, edge_t &one, edge_t &two); 
void find_two_edges(graph_t &g, AmiBase::stat_type requested_type, edge_t &one, edge_t &two);

void find_in_out_edges_oftype(graph_t &g, vertex_t &v, edge_t &in_edge_output,edge_t &out_edge_output, AmiBase::stat_type stat_value);

void find_rand_internal_edge(graph_t &g, AmiBase::stat_type requested_type, edge_t &one);

// Bubble logical_and
void bubble_finder(graph_t &g, vertex_vector_list_t &bubble_vertex_list,  edge_vector_list_t &bubble_edge_list, edge_vector_list_t &legs_edge_list);
void get_bubble_type( graph_t &g, bubble_type &b_type, vertex_vector_t bubble_verts, edge_vector_t bubble_edges);

// tadpole helpers
void find_tadpoles(graph_t &g, vertex_vector_t &tp_vec, vertex_vector_t &tp_conn_vec, edge_vector_t &tp_bose_edges);
void find_tadpoles(graph_t &g, vertex_vector_t &tp_vec, vertex_vector_t &tp_conn_vec, edge_vector_t &tp_bose_edges, edge_vector_t &edge_a, edge_vector_t &edge_b );
void label_and_find_tadpoles_ami(graph_t &g, vertex_vector_t &tp_vec, vertex_vector_t &tp_conn_vec, edge_vector_t &tp_bose_edges);
void find_and_label_tadpole_strings(graph_t &g, vertex_vector_t &tp_vec, vertex_vector_t &tp_conn_vec, edge_vector_t &tp_bose_edges);
void treat_ab_edges(graph_t &g, edge_vector_t &edges_a,edge_vector_t &edges_b, edge_vector_t &new_edges_a,edge_vector_t &new_edges_b);
void get_ab_edges(graph_t &g, vertex_t &v, edge_t &edge_a, edge_t &edge_b);
int count_tadpoles(graph_t &g);
int count_bubbles(graph_t &g);
/////////
// Label graph helper functions primarily in labelling.cpp
/////////


void find_unlabelled_fermionic_edges(graph_t &g, edge_vector_t &vector);
void find_unlabelled_edges(graph_t &g, edge_vector_t &vector);
void get_edges_fermi_bose(graph_t &g, vertex_t v, edge_vector_t fermi_edges, edge_vector_t bose_edges);
void reset_g(graph_t &g);
void reset_alphas(graph_t &g);
void find_external_vertices(graph_t &g, vertex_vector_t &v, edge_vector_t &edges);
void delete_legs(graph_t &g, vertex_vector_t &v, edge_vector_t &edges);
int find_start_index(graph_t &g, vertex_t &v);

int find_end_index(graph_t &g, vertex_t &v); // only for pp bubble graphs

// needs checking but should work for any bubble 
void find_start_end_vertices( graph_t &g, vertex_t &start, vertex_t &end);

void get_kkp(std::vector<AmiBase::alpha_t> &kkp, AmiGraph::graph_t &g);
void get_pp_kkp(std::vector<AmiBase::alpha_t> &kkp, AmiGraph::graph_t &g);

void get_cond_kkp(std::vector<AmiBase::alpha_t> &kkp, AmiGraph::graph_t &g);



void get_adjacent_vertex_stat(graph_t &g, vertex_t &vin, vertex_t &next, AmiBase::stat_type stat);
void get_adjacent_vertex_stat(graph_t &g, vertex_t &vin, vertex_t &next,  AmiBase::stat_type stat, bool &success);
void get_adjacent_vertex_stat(graph_t &g, vertex_t &vin, vertex_t &next, edge_t &eout,  AmiBase::stat_type stat);

void sort_labelled_unlabelled_adjacent(graph_t &g, vertex_t &vin, vertex_vector_t &labelled_adj_vertices,vertex_vector_t &unlabelled_adj_vertices, edge_vector_t &labelled_edges, edge_vector_t &unlabelled_edges);

bool label_consv_momentum(graph_t &g, vertex_t &vin, vertex_vector_t &labelled_adj_vertices,vertex_vector_t &unlabelled_adj_vertices, edge_vector_t &labelled_edges, edge_vector_t &unlabelled_edges);

int count_labelled_edges_of_vert(vertex_t &v, graph_t &g);
void collect_vertices_to_label(graph_t &g, vertex_vector_t &vert_list);

bool take_a_step(vertex_t &current, vertex_t &next, edge_t &connect,  graph_t& g);

void label_tadpoles_ami(graph_t &g);

void label_indep_edge(edge_t &e, graph_t &g);
void label_indep_alpha(edge_t &e, graph_t &g);
void label_indep_epsilon(edge_t &e, graph_t &g);
void mark_labelled(edge_t &e, graph_t &g);

void label_extern_legs(edge_vector_t &extern_vect_list,graph_t &g);

void assign_cons_label(graph_t &g, vertex_t &vin, edge_t &fermi_one, edge_t &fermi_two, edge_t &bose);
void assign_cons_alpha(graph_t &g, vertex_t &vin, edge_t &fermi_one, edge_t &fermi_two, edge_t &bose);

// label function
void label_systematic(graph_t &g);
void label_systematic_redone(graph_t &g);
void label_half_random(graph_t &g);
void random_labelling(graph_t &g, bool &result);
void repeated_labelling(graph_t &g, bool &result);

bool is_labelled_bose_zero(edge_t &e,graph_t &g);
bool edge_labels_are_equal(edge_t &one, edge_t &two, graph_t &g);
void fix_epsilons(graph_t &g);
void reset_epsilons(graph_t &g);
void reset_epsilons(AmiBase::g_prod_t &R0);
void syk_epsilons(graph_t &g);
void fix_force_labels(graph_t &g, vertex_vector_t &in_vv, vertex_vector_t &out_vv);
void find_path_between_vertices(graph_t &g, vertex_t &v1, vertex_t &v2, bool &success, edge_vector_t &final_ev, vertex_vector_t &final_vv);
void put_back_legs(graph_t &g, vertex_vector_t &in_vv, vertex_vector_t &out_vv);

bool edge_alphas_are_equal(edge_t &one, edge_t &two, graph_t &g);
bool edge_alphas_are_negative(edge_t &one, edge_t &two, graph_t &g);

void print_label_counts(std::vector<int> in);
void  label_counts(graph_t &g, std::vector<int> &out);

// momentum conservation
void check_momentum_conservation(graph_t &g, bool &result);


/////////////////////////////////////////

// alpha operations
void add_labels(edge_t &one, edge_t &two, graph_t &g);
void subtract_labels(edge_t &one, edge_t &two, graph_t &g);


// Convert graph to ami input
AmiBase::g_prod_t graph_to_R0(graph_t &g);
void graph_to_R0(graph_t &g, AmiBase::g_prod_t &R0);

void extract_bose_alphas(graph_t g, std::vector<AmiBase::alpha_t> &bose);
// These four functions only depend on AMI content, so moved to libami
//AmiBase::ami_vars construct_ami_vars(AmiBase::g_prod_t &R0,NewAmiCalc::internal_state &state, NewAmiCalc::ext_vars &external);
//AmiCalc::energy_t construct_energy(AmiBase::g_prod_t &R0, NewAmiCalc::internal_state &state, NewAmiCalc::ext_vars &external);
//std::complex<double> eval_epsilon(NewAmiCalc::k_vector_t k, std::complex<double> mu );
//NewAmiCalc::k_vector_t construct_k(AmiBase::alpha_t alpha, NewAmiCalc::k_vect_list_t &k);

void set_amiparms(AmiBase::ami_parms &amiparms, int &graph_order, double &ereg);
void randomize_state(NewAmiCalc::internal_state &state, int order);
void randomize_state(NewAmiCalc::internal_state &state, int k_length, int gsize);
void randomize_state_spherical(NewAmiCalc::internal_state &state, int k_length, int gsize);
void shift_state(NewAmiCalc::internal_state &state, int k_index, NewAmiCalc::k_vector_t &shift);
void shift_state_list(NewAmiCalc::internal_state &state,  NewAmiCalc::k_vect_list_t &shift);
void zero_state(NewAmiCalc::internal_state &state, int k_length);
void shift_state_pi(NewAmiCalc::internal_state &state, int k_length);
void construct_state(NewAmiCalc::internal_state &state);
void print_state(NewAmiCalc::internal_state &state);
void print_fullstate(NewAmiCalc::internal_state state);
void print(AmiBase::ami_vars vars);
void print_kvector(NewAmiCalc::k_vector_t &k);

int count_fermi_loops(graph_t &g);
double get_prefactor(graph_t &g, int order);

void collect_unvisited_vertices(graph_t &g, vertex_vector_t &vertex_list);

void reset_visited(graph_t &g);

// filter functions

void get_adjacent_lists(graph_t &g, vertex_vector_t &current_list, vertex_vector_t &next_list);
void get_adjacent_lists_undirected(graph_t &g, vertex_vector_t &current_list, vertex_vector_t &next_list);
void label_visited(graph_t &g, vertex_vector_t &current_list);
bool is_connected_graph(graph_t &g);

bool is_one_particle_reducible(graph_t &g);
bool is_skeleton(graph_t &g);
bool is_hubbard(graph_t &g);
bool is_hubW_graph(graph_t g);

bool is_isomorphic(graph_t g1, graph_t g2);
bool is_AID(graph_t g1, graph_t g2);
void symmetrize_bosonic_lines(graph_t &g);
void symmetrize_internal_bosonic_lines(graph_t &g);
void symmetrize_fermionic_lines(graph_t &g);
void symmetrize_non_principle_lines(graph_t &g);
void get_principle_line(graph_t &g, edge_vector_t &ev); // TODO: does this even work?
bool is_1P_reducible_graph(graph_t g);
bool is_oneleg(graph_t g);
bool has_fock_insertion(graph_t g);
bool disconnected_loop(graph_t g);
bool is_1P_bose_reducible_graph(graph_t g);
// bool is_1P_fermi_reducible_graph(graph_t g);

void get_principle_loop(graph_t &g, vertex_vector_t &vv, edge_vector_t &ev);
void find_principle_bose_lines(graph_t &g, edge_vector_t &ev);
bool is_ladder_graph(graph_t &g);
bool is_nested(graph_t &g);
void split_graph(edge_t &e, graph_t &g);

bool is_pi_graph(graph_t &g);

bool is_ppvertex_graph(graph_t g);
bool get_next(graph_t &g, vertex_t &current, vertex_t &next, vertex_t &end_v);

// graph generation
void generate_graphs(std::vector< std::vector< AmiGraph::graph_t>> &graph_matrix, int max, int mpi_rank);
void generate_bubble_graphs(std::vector< std::vector< AmiGraph::graph_t>> &graph_matrix, int max, int mpi_rank);
void generate_sigma_graphs(std::vector< std::vector< AmiGraph::graph_t>> &graph_matrix, int max, int mpi_rank);
void generate_greens_graphs(std::vector< std::vector< AmiGraph::graph_t>> &graph_matrix, int max, int mpi_rank);
void generate_force_graphs(std::vector< std::vector< AmiGraph::graph_t>> &graph_matrix, int max, int mpi_rank);

// density functions 
void gm_close_loops(std::vector< std::vector< AmiGraph::graph_t>> &gm);
void close_loop(AmiGraph::graph_t &g);


void gm_attach_probe_lines(std::vector< std::vector< AmiGraph::graph_t>> &gm,std::vector< std::vector< AmiGraph::graph_t>> &gm_out);
void attach_probe_lines(AmiGraph::graph_t &g, std::vector< AmiGraph::graph_t> &gm_out);
void probe_lines(graph_t &g, edge_t &e1, edge_t &e2);



void label_graphs(std::vector< std::vector< AmiGraph::graph_t>> &graph_matrix, int min, int max);

void reduce_gm_rf(std::vector< std::vector< AmiGraph::graph_t>> &graph_matrix, int mpi_rank);
void reduce_gm_tp(std::vector< std::vector< AmiGraph::graph_t>> &graph_matrix, int mpi_rank, int min);
void reduce_gm_oneleg(std::vector< std::vector< AmiGraph::graph_t>> &graph_matrix, int mpi_rank);
void reduce_gm_skeleton(std::vector< std::vector< AmiGraph::graph_t>> &graph_matrix, int mpi_rank);
void reduce_gm_1PBose(std::vector< std::vector< AmiGraph::graph_t>> &graph_matrix, int mpi_rank, int min);
void reduce_gm_1PFermi(std::vector< std::vector< AmiGraph::graph_t>> &graph_matrix, int mpi_rank, int min);
void reduce_gm_fock(std::vector< std::vector< AmiGraph::graph_t>> &graph_matrix, int mpi_rank, int min);
void reduce_gm_nested(std::vector< std::vector< AmiGraph::graph_t>> &graph_matrix, int mpi_rank, int min);
void reduce_gm_RPA_chains(std::vector< std::vector< AmiGraph::graph_t>> &graph_matrix, int mpi_rank, int min);
void reduce_gm_hubbard(std::vector< std::vector< AmiGraph::graph_t>> &graph_matrix, int mpi_rank, int min);

void reduce_gm_ladder(std::vector< std::vector< AmiGraph::graph_t>> &graph_matrix, int mpi_rank, int min);

void reduce_gm_ppvertex(std::vector< std::vector< AmiGraph::graph_t>> &graph_matrix, int mpi_rank, int min);

// helper
void get_principle_pp_lines(graph_t &g, vertex_vector_list_t &vvl, edge_vector_list_t&evl);

// isomorphic functions

bool full_iso(graph_t g1, graph_t g2);
bool spin_iso(graph_t g1, graph_t g2);

edge_vector_t get_edges(int si, int ti, graph_t &g);
bool edges_equivalent( edge_vector_t ev1, edge_vector_t ev2, graph_t &g1, graph_t &g2);
// Boost connected components
int connected_components(graph_t &g);
int fermi_connected_components(graph_t &g);


bool bose_lists_equal(std::vector< std::vector<int>> b1, std::vector< std::vector<int>> b2);
void get_labelled_bose_lists(std::vector< std::vector<int>> &bose_list,  graph_t &g);
void systematic_vertex_label(graph_t &g);
void get_next_fermi_loop(vertex_vector_t &list, vertex_t &current, graph_t &g);
int get_length_principle_line(graph_t &g);
bool edges_equal(graph_t &g1, graph_t &g2);
void get_labelled_edge_list(std::vector< std::vector<int>> &bose_list,  graph_t &g);

void construct_graph_barcode(graph_t &g, barcode_t &barcode, tree_t &tree);
bool barcode_lines_equal(barcode_line_t B1_line, barcode_line_t B2_line);
bool barcodes_are_equal(barcode_t &B1, barcode_t &B2);

void childless_parent(tree_t &tree, vertex_vector_t &vertex_list);
void remove_childless(graph_t &g, vertex_vector_t &vertex_list);

void label_via_tree(tree_t &tree, graph_t &g);
void compare_trees(tree_t &T1, graph_t &g1, tree_t &T2, graph_t &g2);

void print_barcodes(barcode_t B1, barcode_t B2);

bool children_are_equal( std::vector<child_t> c1, std::vector<child_t> c2);
bool parents_are_equal( parent_t p1, parent_t p2);
bool pc_are_equal( pc_struct pc1, pc_struct pc2);

void get_ForB_line(graph_t &g, vertex_vector_t &current_list, vertex_vector_t &next_list, std::vector< std::vector<int>> &FBarray );
void get_barcode_line(graph_t &g, vertex_vector_t &current_list, vertex_vector_t &next_list, barcode_t &barcode);

void get_tree_leg(graph_t &g, vertex_vector_t &current_list, vertex_vector_t &next_list, tree_t &tree);
void convert_tree_to_barcode(tree_t &tree, barcode_t &barcode, graph_t &g);

// filter helpers
AmiGraph::spin_type return_spin_on_vertex(graph_t &g, vertex_t &v, bool &fail); // todo this one is kind of garbage...
bool spin_mismatch(vertex_t &v, graph_t &g);

// Specific defined updates for systematic generation of diagrams
void AT_defined(graph_t &g, edge_t &pickone, spin_type spin);
void AB_defined(graph_t &g, edge_t &one, edge_t &two, spin_type spin);
void AIL_defined(graph_t &g, edge_t &one, edge_t &two);
void RIB_defined(graph_t &g, edge_t &one, spin_type spin);


void graph_matrix_update(graph_t &g, std::vector< std::vector<graph_t>> &graph_matrix);

// residue and GOS
bool residue_filter(graph_t &g); // requires a labelled graph 

void reconstruct_graph_barcode(graph_t &g, barcode_t &barcode, tree_t &tree);
void reconvert_tree_to_barcode(tree_t &tree, barcode_t &barcode, graph_t &g);


// GOS structs

typedef std::vector<int> git_perm_t;
typedef std::vector<git_perm_t> git_perm_set_t;
typedef std::vector<git_perm_set_t> git_perm_list_t;

typedef std::vector<graph_t> labels_t;

	

struct git_pair{
	
git_pair(graph_t g1, graph_t g2, git_perm_set_t pst){ 

g1_=g1;
g2_=g2;
pst_=pst;

}

git_pair(){}		
	
graph_t g1_;
graph_t g2_;
git_perm_set_t pst_;
NewAmiCalc::solution_set s1_, s2_;	
};

struct graph_group {
  
  std::vector<graph_t> graph_vec;
  std::vector<git_pair> gp_vec;
  std::vector<graph_t> ct_vec; // this is not used I think due to order issues. 
  
  std::vector<double> prefactor;
  std::vector<labels_t> labels;
  std::vector<NewAmiCalc::solution_set> ss_vec; // solution set vector 
  
  std::vector<AmiSpec::spec_solution_set> sss_vec;// spectral solution set vector 
  
  // git permutation set 
  // git_perm_list_t perm_vec;
  int order_shift=0;
     
};

typedef std::vector<graph_group> gg_vec_t;
typedef std::vector< gg_vec_t > gg_matrix_t;

void close_ggm(gg_matrix_t &ggm);



// symbolic wicks determinants 

// typedef std::vector< std::vector< std::pair<int,int> > > sym_M;
// typedef std::vector< std::pair<int,int > > st_pair_list;

struct symDET_element {
std::vector< std::pair<int,int> >  plist;
int sign;
graph_t graph;	
NewAmiCalc::solution_set ss;
NewAmiCalc::solution_set rtss;// runtime solution set - get modified when epsilons are duplicated
Ulist_type Ulist;
band_map bmap;

bool use_rtss=0;
};

void print_matrix(std::vector< std::vector <int>> &M);
void extract_M(symDET_element &SE, std::vector< std::vector< int > > &M_OUT);

void extract_ABC(symDET_element &SE, Eigen::MatrixXi &A, Eigen::MatrixXi &B, Eigen::MatrixXi &C);

typedef std::vector<symDET_element> symDET;
typedef std::vector<symDET> symM;

// void construct_symwicks_matrix(int num_vert, sym_M &M);
// void convert_to_symdet(sym_M &M, symDET &SD);
void make_symDET(int L, symDET &SD);
void construct_sigma_symM(symM &SM, int min, int max);
void construct_symM_full(symM &SM, int min, int max);
void SD_multiply(symDET &SD1, symDET &SD2, symDET &SD_out, int shift);

void get_coll(std::vector< std::vector<int>> &cv, int val);

bool is_diagonal(band_map &bm, std::vector<int> &av, std::vector<int> &bv);
void assign_species(symDET_element &SD, std::vector<int> &av, std::vector<int> &bv, std::vector<int> &dup_map);

void symDET_element_construct_ami_sol(symDET_element &SE, double ereg);

void construct_rtss_ami_sol(symDET_element &SE);
void assign_solution(symDET_element &SD);

void construct_avbv_sets(std::vector< std::vector<int>> &avlist, std::vector< std::vector<int>> &bvlist, int order, int num_bands, int ext_index);
void reduce_avbv_sets(U_matrix_type &UM, std::vector< std::vector< int>> &avlist,std::vector< std::vector< int>> &bvlist);

void construct_delta_sets(std::vector< std::vector<int>> &dlist, int order, int band_num, int ext_index);
void get_avbv_from_delta( symDET_element &SD,  std::vector<int> &delta, std::vector<int> &ai, std::vector<int> &bi);

void construct_offdiag_delta_sets(std::vector< std::vector<int>> &dlist, int order, int band_num);
void get_avbv_from_offdiag_delta( symDET_element &SD,  std::vector<int> &delta, std::vector<int> &ai, std::vector<int> &bi, int in_band, int out_band);

void general_construct_avbv_sets(std::vector<std::vector<int>> &avlist,std::vector<std::vector<int>> &bvlist, int band_num, int order, int ext_index);
void general_construct_avbv_sets(std::vector<std::vector<int>> &avlist,std::vector<std::vector<int>> &bvlist, int order, int band_num, int in_index, int out_index);

void construct_minimal_avbv_sets(U_matrix_type &UM, std::vector<std::vector<int>> &avlist,std::vector<std::vector<int>> &bvlist, int order, int band_num, int ext_index);

void avbv_add_col( std::vector<std::vector<int> > &inlist, std::vector<std::vector<int> > &outlist, int maxval, int maxlength);


void construct_graphs(symDET &SD);
graph_t return_graph(symDET_element se);
std::vector< std::pair<int,int> > return_contraction(graph_t g);
void print_plist(std::vector< std::pair<int,int> > &plist);

void symM_construct_ami_sol(symM &SM, double ereg);
void extract_Ulist(graph_t &g, Ulist_type &U);

double eval_Ulist(U_matrix_type &Uvals, std::vector<int> &av, std::vector<int> &bv);
void random_avbv(std::vector<int> &av, std::vector<int> &bv, int order, int num_bands, int ext_index);

void symM_bandmaps(symM &SM);
void generate_bandmap(graph_t &g, band_map &bm);

/// probably broken due to passing struct...
void convert_symelement_to_graph(symDET_element &se);
void symM_make_graphs(symM &SM);
void symM_filter_sigma(symM &SM);

void symM_filter(symM &SM);

void symM_filter_mpi(symM &SM, int mpi_rank, int comm_size);

void symM_label(symM &SM);

//
void make_symDET_HF(int L, symDET &SD);
void construct_HF_symM(symM &SM, int min, int max);

//

void make_symDET_sigma_HF(int L, symDET &SD);
void construct_sigma_HF_symM(symM &SM, int min, int max);

// molecule stuff



// Eigen solver for GOS
void alpha_to_matrix(graph_t &g, Eigen::MatrixXd &A_out);
void full_alpha_to_matrix(graph_t &g, Eigen::MatrixXd &A_out);
void remove_simple_alphas(graph_t &g, edge_vector_t &internal);
void construct_transform(Eigen::MatrixXd &A, Eigen::MatrixXd &X, Eigen::MatrixXd &Y);
void prune(Eigen::MatrixXd &inout, float threshold);

bool abs_equal(Eigen::MatrixXd &A1, Eigen::MatrixXd &A2);
bool neg_equal(Eigen::MatrixXd &A1, Eigen::MatrixXd &A2);
bool equal(Eigen::MatrixXd &A1, Eigen::MatrixXd &A2);

// void row_col_counts(Eigen::MatrixXd &A_in, std::vector<int> &row, std::vector<int> &col);
void row_col_counts(Eigen::MatrixXd &A_in, Eigen::VectorXd &row, Eigen::VectorXd &col);
bool row_col_compare(Eigen::VectorXd r1, Eigen::VectorXd r2, Eigen::VectorXd c1, Eigen::VectorXd c2);
//

//void reorder_AID(std::vector< std::vector<graph_t>> &graph_matrix);
void reorder_AID(std::vector< std::vector<graph_t>> &graph_matrix, std::vector< std::vector<int>> &gg_lists, int mpi_rank);
//void reduce_gos( std::vector< AmiGraph::graph_t> &graph_vec, gg_vec_t &gg_vec);
void gos(std::vector< std::vector< AmiGraph::graph_t>> &graph_matrix, gg_matrix_t &ggm);
void gos(std::vector< std::vector< AmiGraph::graph_t>> &graph_matrix, gg_matrix_t &ggm, std::vector< std::vector<int>> &gg_lists);
void repeated_gos(int &max, bool &result, int &j, int &k , std::vector< AmiGraph::graph_t> &graph_vector,  graph_group &gg_temp, std::vector<int> &list);

int get_start_index(std::vector<int> list);

void gos_compare( graph_t &g1, graph_t &g2, bool &result, double &prefactor);
// void single_gos_compare(graph_t &g1, graph_t &g2, bool &result, double &prefactor);
void single_gos_compare(labels_t &L1, labels_t &L2, bool &result, double &prefactor);
void reduce_ggv_gos( gg_vec_t &ggv);
void reduce_ggm_gos( gg_matrix_t &ggm, int mpi_rank, int cutoff);

void reduce_ggm_isomorphic(gg_matrix_t &ggm);

void collapse_ggv(int i, int j, gg_vec_t &ggv, double prefactor);
void print_ggm( gg_matrix_t &ggm);
void mpi_print_ggm( gg_matrix_t &ggm, int mpi_rank);

void gg_check(gg_matrix_t &ggm);
void linearize_ggv(gg_vec_t &ggv_in, std::vector<graph_t> &gv_out);
void gm_to_ggm( std::vector< std::vector< AmiGraph::graph_t>> &gm, gg_matrix_t &ggm);
void ggm_to_gm(  gg_matrix_t &ggm, std::vector< std::vector< AmiGraph::graph_t>> &gm, int min);
void gm_to_ggm_all( std::vector< std::vector< AmiGraph::graph_t>> &gm, gg_matrix_t &ggm, int mpi_rank);
void ggm_AID_reduce(int attempts, gg_matrix_t &ggm, std::vector< std::vector<int>> &gg_lists, int mpi_rank);

void ggm_clear_labels(gg_matrix_t &ggm);
void ggm_remove_pairs(gg_matrix_t &ggm, int min);
void ggm_reset_epsilons(gg_matrix_t &ggm);
void ggm_syk_epsilons(gg_matrix_t &ggm);

void ggm_fold_in_bose_alphas(gg_matrix_t &ggm);
void fold_in_bose_alphas(graph_t &g);

void ggm_check_pair_momenta(gg_matrix_t &ggm, int mpi_rank);

void gos_sort_rows(Eigen::MatrixXd &A1, Eigen::MatrixXd &A2, Eigen::VectorXd &r1, Eigen::VectorXd &r2, Eigen::MatrixXd &perm);
void gos_sort_cols(Eigen::MatrixXd &A1, Eigen::MatrixXd &A2, Eigen::VectorXd &r1, Eigen::VectorXd &r2, Eigen::MatrixXd &perm);



double get_prefactor(Eigen::MatrixXd &A1, Eigen::MatrixXd &A2);

typedef std::vector<int> pair_t;
typedef std::vector<pair_t> set_t;
typedef std::vector<set_t> rc_set_t;
typedef std::vector<rc_set_t> permutation_set_t;




void get_rc_permutations( Eigen::VectorXd &r1, Eigen::VectorXd &c1, std::vector< std::vector<int>> &row_perm_set, std::vector< std::vector<int>> &col_perm_set);
void expand_permutations(std::vector< std::vector<int>> &row_perm_set, std::vector< std::vector<int>> &col_perm_set, permutation_set_t &pst);

bool compare_permutation_set(Eigen::MatrixXd &A1, Eigen::MatrixXd &A2, permutation_set_t &pst, double &pf);
bool compare_permutation(Eigen::MatrixXd A1, Eigen::MatrixXd &A2, rc_set_t &rcst, double &pf);
void apply_permutation(Eigen::MatrixXd &A2, rc_set_t &rcst);

void max_ext_label_counts(std::vector< std::vector< AmiGraph::graph_t>> &graph_matrix);
void min_ext_label_counts(std::vector< std::vector< AmiGraph::graph_t>> &graph_matrix);
bool in_list(int k, std::vector<int> list);


// pre labelling collection stuff
bool edge_labels_are_equal(edge_t &one,graph_t &g1,  edge_t &two, graph_t &g2);
void construct_label_sets(gg_matrix_t &ggm, int mpi_rank, int min, int size, int trim);
void construct_label_set(graph_t &g, labels_t &L, int size);
void check_label_sets(gg_matrix_t &ggm, int mpi_rank, int min);
bool equal_graph_labels(graph_t &g1, graph_t &g2);
void get_edges(edge_vector_t &e, graph_t &g);
void append_to_labels(labels_t &L, graph_t &g, bool &added);

// separation
int count_spins(graph_t &g);
void ggm_split(gg_matrix_t &ggm, gg_matrix_t &ggm_uu, gg_matrix_t &ggm_ud);
void ggm_1P_split(gg_matrix_t &ggm, gg_matrix_t &ggm_irr, gg_matrix_t &ggm_red);
void ggm_eo_split(gg_matrix_t &ggm, gg_matrix_t &ggm_even, gg_matrix_t &ggm_odd); // split based on even or odd bubble sets 
//typedef std::vector< std::pair< std::vector< std::pair<int, int > >, std::vector< std::pair<int, int > >>> permutation_set_t;


//typedef std::vector<NewAmiCalc::solution_set> ggss_t
// graph group solutions with ami
void ggm_construct_ami_sol(gg_matrix_t &ggm, double ereg, int mpi_rank);
void gg_construct_ami_sol(graph_group &gg, double ereg);
void ggm_construct_ami_term_sol(gg_matrix_t &ggm, double ereg, int mpi_rank);
void gg_construct_ami_term_sol(graph_group &gg, double ereg);

void merge_ggm(gg_matrix_t &ggm1, gg_matrix_t &ggm2);

void construct_ami_sol(graph_t &g, NewAmiCalc::solution_set &ss, double ereg);
void construct_ami_terms_sol(graph_t &g, NewAmiCalc::solution_set &ss, double ereg);
void ct_construct_ami_sol(graph_t &g, NewAmiCalc::solution_set &ss, double ereg);// differs from construct_ami_sol only in how number of integrals is determined. in this case it is hard coded to be alpha_.size()-1.  This could probably be done most of the time as a starting point. 
int alpha_size(graph_t &g);

void ggm_filter_fock(gg_matrix_t &ggm);
void ggm_filter_tp(gg_matrix_t &ggm);
void ggm_filter_hfspinsusc(gg_matrix_t &ggm);
void ggm_filter_max_eff_order(gg_matrix_t &ggm, int max, int max_bub_ct, int max_sigma_ct, int ct_size);
void ggm_group_eff_order_CT(gg_matrix_t &ggm, int max_bub_ct, int max_sigma_ct);
void ggm_reduce_ami_terms(gg_matrix_t &ggm,double ereg, int mpi_rank, int max_try);
void ggm_maximize_ami_terms(gg_matrix_t &ggm,double ereg, int mpi_rank, int max_try);
void ggm_optimize_ami(gg_matrix_t &ggm, int mpi_rank);

void ggm_opt_der(gg_matrix_t &ggm, int mpi_rank);

void ggm_to_amim(gg_matrix_t &ggm, NewAmiCalc::solution_set_matrix_t &AMI_MATRIX); // this is for non-group evaluation.
void ggm_to_ggamim(gg_matrix_t &ggm, NewAmiCalc::gg_solution_set_matrix_t &GG_AMI_MATRIX, int min); // this is for group and non-group evaluation
void ggm_to_sp_ggamim(gg_matrix_t &ggm, AmiSpec::sp_gg_solution_set_matrix_t &sp_GG_AMI_MATRIX, int min);
void ggm_generate_sp_terms(gg_matrix_t &ggm, int mpi_rank);
void ggm_generate_simple_sp_terms(gg_matrix_t &ggm, int mpi_rank);


// playing with cuba
std::vector< double > integrator(NewAmiCalc::solution_set &sol, NewAmiCalc::ext_vars &ext, int nsamples, int intseed);
// static int ami_integrand(const int *ndim, const double xx[],  const int *ncomp, double ff[], void *userdata);

// sobol test
void sobol_randomize_state(NewAmiCalc::internal_state &state, int k_length);

// ggm file io functions
void ggm_label(gg_matrix_t &ggm, int min);
void graph_write(std::string filename, graph_t &g);
void graph_read(std::string filename, graph_t &g);

void write_ggm(std::string folder, gg_matrix_t &ggm);
void read_ggm(std::string folder, gg_matrix_t &ggm, int max_order);

// ggm pair file io functions
void gp_write(std::string filename, git_pair &p);
void write_ggmp(std::string folder, gg_matrix_t &ggm, int min);
void read_ggmp(std::string folder, gg_matrix_t &ggm, int max_ord);
void read_ggmp(std::string folder, gg_matrix_t &ggm, int min_ord, int max_ord);
void pair_read(std::string filename, git_pair &p);

// void find_entry_vertices(graph_t &g, vertex_vector_t &v);
void find_entry_vertex(graph_t &g, vertex_t &v);
///
// NEW ATTEMPT AT FULL GIT
void single_git_compare(labels_t &L1, labels_t &L2, bool &result, double &prefactor, git_perm_set_t &pst, graph_t &g2_out);
void git_compare( graph_t &g1, graph_t &g2, bool &result, double &prefactor);
void git_sort_rows(Eigen::MatrixXd &A1, Eigen::MatrixXd &A2, Eigen::VectorXd &r1, Eigen::VectorXd &r2, Eigen::MatrixXd &perm);
void git_sort_cols(Eigen::MatrixXd &A1, Eigen::MatrixXd &A2, Eigen::VectorXd &r1, Eigen::VectorXd &r2, Eigen::MatrixXd &perm);

void git_swap_row_signs(int row_index, Eigen::MatrixXd &A2,Eigen::VectorXd &energy );
void git_swap_col_signs(int index, Eigen::MatrixXd &A2);
void git_swap_easy_rows(Eigen::MatrixXd &A1, Eigen::MatrixXd &A2, Eigen::VectorXd &energy );
void git_fix_ext_signs(Eigen::MatrixXd &A1, Eigen::MatrixXd &A2,Eigen::VectorXd &energy );
void git_decide_swap(Eigen::MatrixXd &A1, Eigen::MatrixXd &A2, Eigen::VectorXd &energy, Eigen::VectorXd &kshifts);

void get_git_rc_permutations( Eigen::VectorXd &r1, Eigen::VectorXd &c1, std::vector< std::vector<int>> &row_perm_set, std::vector< std::vector<int>> &col_perm_set);
void git_get_perm_set(std::vector< std::vector<int>> &row_perm_set, std::vector< std::vector<int>> &col_perm_set, git_perm_set_t &pst);

// bool git_find_permutation(Eigen::MatrixXd &A1, Eigen::MatrixXd &A2, git_perm_set_t &pst);
bool git_compare_n_permutations(int n, Eigen::MatrixXd A1, Eigen::MatrixXd &A2, git_perm_set_t &pst);

bool git_compare_rc_permutations(Eigen::MatrixXd A1, Eigen::MatrixXd &A2, git_perm_set_t &pst);

// bool git_compare_single_permutation(Eigen::MatrixXd A1, Eigen::MatrixXd &A2, permutation_set_t &pst);
void git_apply_permutation(Eigen::MatrixXd &A2, git_perm_t &gpt);

void print_permutation_set( git_perm_set_t &pst);

void reorder_matrix(Eigen::MatrixXd &A1);

void git_k_sign_change(Eigen::MatrixXd &A1, Eigen::MatrixXd &A2, Eigen::VectorXd &energy, Eigen::VectorXd &kshifts );

void git_apply_sign_permutation(Eigen::MatrixXd &A2, Eigen::VectorXd &energy, Eigen::VectorXd &kshifts, git_perm_t &gpt);
void git_shift_kp(int index, Eigen::MatrixXd &A2, Eigen::VectorXd &energy, Eigen::VectorXd &kshifts);
bool git_sign_permutations(Eigen::MatrixXd &A1, Eigen::MatrixXd &A2, Eigen::VectorXd &energy, Eigen::VectorXd &kshifts);

// git ggm
void reduce_ggm_git( gg_matrix_t &ggm, int mpi_rank, int cutoff);
void reduce_ggv_git( gg_vec_t &ggv);

void reduce_ggv_git_odd_random( gg_vec_t &ggv);
void reduce_ggm_git_odd_random( gg_matrix_t &ggm, int mpi_rank, int cutoff);
void repeated_git_compare( graph_t &g1, graph_t &g2, bool &result, double &prefactor,git_perm_set_t &pst, graph_t &g2_out);

// ugh
void replace_graph_with_matrix(Eigen::MatrixXd &A_in, graph_t &g );

void git_apply_T1(Eigen::MatrixXd &A2, Eigen::RowVectorXd &perm, int first, int second);
void git_apply_T2(Eigen::MatrixXd &A2, Eigen::RowVectorXd &perm, int first);
void git_apply_T3(int index, Eigen::MatrixXd &A2, Eigen::VectorXd &energy, Eigen::VectorXd &kshifts);
void git_compare_v2( graph_t &g1, graph_t &g2, bool &result, double &prefactor,git_perm_set_t &pst, graph_t &g2_out);
bool git_sign_permutations_v2(Eigen::MatrixXd &A1, Eigen::MatrixXd &A2, Eigen::VectorXd &energy, Eigen::VectorXd &kshifts, git_perm_set_t &pst_out);
void permute_state(NewAmiCalc::internal_state &state, NewAmiCalc::k_vector_t &k_shift, git_perm_set_t &pst);

// git eval permutations
//NOTE: these two functions didn't really work - need to make pairs only.
// void ggm_assign_git_perms(gg_matrix_t &ggm, int mpi_rank);
// void assign_git_perms(graph_t &g1, graph_t &g2, git_perm_set_t &pst);
void ggm_populate_git_pairs(gg_matrix_t &ggm, int mpi_rank, int min);
void group_reduce_to_pairs(graph_group &gg, std::vector<int> &map);
void gg_pair_remove(graph_group &gg, int i, int j);

void generate_git_map(graph_group &gg, std::vector<int> &map_out);




// systematic labelling rewrite
void sys_label_sets(labels_t &L, graph_t g, bool &result, int max);
void sys_label(graph_t &g, bool &result);
void bubble_sys_label(graph_t &g, bool &result);
void label_graphs_sys(std::vector< std::vector< AmiGraph::graph_t>> &graph_matrix,int min, int max);
void bubble_label_graphs_sys(std::vector< std::vector< AmiGraph::graph_t>> &graph_matrix, int min, int max);
void trim_label_set(labels_t &L, int trim);

// Relabelling to avoid matsubara only poles
void gg_remove_matsubara_poles(graph_group &gg);
void ggm_fix_matsubara_poles(gg_matrix_t &ggm);

bool has_matsubara_pole(AmiBase::P_t &this_P);
void fix_matsubara_pole(graph_t &g, NewAmiCalc::solution_set &sol);

// graph counter-term construction
bool not_labelled(graph_t &g);
void generate_bubble_ct( graph_t &g_in, std::vector< graph_t> &ct_vec);
void generate_sigma_ct( graph_t &g_in, std::vector< graph_t> &ct_vec, int maxdots);

void insert_chain(graph_t &g_in, graph_t &g_out,  edge_t &ei, int length);


void bg_to_ct(graph_t &g, graph_t &ctg, bubble &bub );
void bg_to_ct(graph_t gin, graph_t &g_out, std::vector< bubble > bubble_vector, std::vector<int> list);
void descriptor_to_index(graph_t &g, vertex_vector_t &bv, edge_vector_t &be, edge_vector_t &le , bubble &bub);

void save_ct_eff_orders(gg_matrix_t &ggm);
void zero_external(AmiGraph::graph_t &g, int index);
void zero_external_Q(AmiGraph::graph_t &g);
void PIPI_external_Q(AmiGraph::graph_t &g);
void swap_alphas(AmiGraph::graph_t &g, int first, int second);



// Multiband stuff
void mb_setup(graph_t &g, std::vector< std::vector< int >> &V_prod);
void get_vprod(graph_t &g, edge_vector_t &bev,std::vector< std::vector< int > > &V_prod );
std::vector< int > get_v(graph_t &g, edge_t &be);
std::complex< double > eval_symbolic_Vprod(std::vector< int > &V_values, std::vector< std::vector< int >> &V_indices,AmiGraph::U_matrix_type &UM);

// mpi broadcast functions

 //MPI specific stuff
  int mpi_rank;
  int mpi_size;

// void broadcast_ggamim( AmiCalc::gg_solution_set_matrix_t & gg,  AmiCalc::gg_solution_set_matrix_t & gg_bc);
// //// these taken from diagai_sim blindly
typedef std::vector<AmiBase::alpha_t> alpha_vec_t;
typedef std::vector< alpha_vec_t > alpha_vec_list_t;
typedef std::vector< alpha_vec_list_t> kkp_vec_t;

typedef std::vector< kkp_vec_t> kkp_matrix_t;
/////////////////////


void broadcast_stdvec_d(std::vector<double> &dub, int from);
void broadcast_stdvec_i(std::vector<int> &dub, int from);
void broadcast_stdvec_cd(std::vector<std::complex<double>> &cdub, int from);

void broadcast_pole_struct( AmiBase::pole_struct &p, int from);
void broadcast_pole_array(AmiBase::pole_array_t &p, int from);
void broadcast_g_prod_t(AmiBase::g_prod_t &gprod, int from);
void broadcast_g_struct( AmiBase::g_struct &g, int from);
void broadcast_ami_parms(AmiBase::ami_parms &p, int from);
void broadcast_k_list(NewAmiCalc::k_vect_list_t &k,  int from);
void broadcast_state(NewAmiCalc::internal_state &s, int from);
void broadcast_hopping_list(NewAmiCalc::hopping_list_t &t, int from);
void broadcast_ext_vars(NewAmiCalc::ext_vars &e, int from);
void broadcast_Ri(AmiBase::Ri_t &ri, int from);
void broadcast_Pi(AmiBase::Pi_t &pi, int from);
void broadcast_Si(AmiBase::Si_t &si, int from);

void broadcast_R_t(AmiBase::R_t &r, int from);
void broadcast_P_t(AmiBase::P_t &p, int from);
void broadcast_S_t(AmiBase::S_t &s, int from);

void broadcast_solution_set(NewAmiCalc::solution_set &s, int from);

void broadcast_solution_set_vec(NewAmiCalc::solution_set_vec_t &s, int from);
void broadcast_solution_set_matrix_t(NewAmiCalc::solution_set_matrix_t &s, int from);
void broadcast_gg_solution_set_matrix_t(NewAmiCalc::gg_solution_set_matrix_t &s, int from);

void broadcast_R_ref_t(AmiBase::R_ref_t &Rref, int from);
void broadcast_ref_v_t(AmiBase::ref_v_t &ref_v, int from);
void broadcast_pair_int_int(std::pair<int,int> &ref_t, int from);

void broadcast_bose_alphas( std::vector<AmiBase::alpha_t> &bose, int from);


void broadcast_ggv(AmiGraph::gg_vec_t &ggv, int from);
void broadcast_ggm(AmiGraph::gg_matrix_t &ggm, int from);
void broadcast_graph(AmiGraph::graph_t &g, int from); // TODO: This does not work. non-trivial to implement. avoid.

void broadcast_graph_group_pairs(AmiGraph::graph_group &gg, int from); // This only transfers AMI solutions. 
void broadcast_git_perm_set(AmiGraph::git_perm_set_t &pst, int from);
void broadcast_gp(AmiGraph::git_pair &gp, int from);
void broadcast_gp_vec(std::vector<AmiGraph::git_pair> &gp_vec, int from);

//
void broadcast_alpha_vec(alpha_vec_t &av, int from );
void broadcast_alpha_vec_list(alpha_vec_list_t &avl, int from );
void broadcast_kkpv(kkp_vec_t &kkpv, int from );
void broadcast_kkpm(kkp_matrix_t &kkpm, int from );





// end class bracket
};





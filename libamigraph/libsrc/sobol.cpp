#include <boost/random/sobol.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/variate_generator.hpp>

int main()
{
    static const std::size_t dimension = 4;

    // Create a generator
    typedef boost::variate_generator<boost::random::sobol&, boost::uniform_01<double> > quasi_random_gen_t;

    // Initialize the engine to draw randomness out of thin air
    boost::random::sobol engine(dimension);
	
	std::vector< boost::random::sobol > eng_vec;
	
	boost::random::sobol temp(dimension);
	for(int i=0; i< dimension; i++){
	
	eng_vec.push_back(temp);
	}
	
	
	eng_vec[0].seed(0);
	eng_vec[1].seed(1);
	

    // Glue the engine and the distribution together
    quasi_random_gen_t gen(engine, boost::uniform_01<double>());

    std::vector<double> sample(dimension);
	
	for(int i=0; i< sample.size(); i++){
		
		std::cout<<"Sample "<<i<<" = "<< sample[i]<<std::endl;
	}

engine.seed(12.4);
// engine.discard(1*dimension);

    // At this point you can use std::generate, generate member f-n, etc.
	
for(int j=0; j< dimension; j++){	
	
    std::generate(sample.begin(), sample.end(), gen);
	
	for(int i=0; i< sample.size(); i++){
		
		std::cout<<"on j "<<j<<" Sample "<<i<<" = "<< sample[i]<<std::endl;
	}
	

quasi_random_gen_t gen0(eng_vec[0], boost::uniform_01<double>());	
	
	
    std::generate(sample.begin(), sample.end(), gen0);
	
	for(int i=0; i< sample.size(); i++){
		
		std::cout<<"on j "<<j<<" Sample "<<i<<" = "<< sample[i]<<std::endl;
	}	

quasi_random_gen_t gen1(eng_vec[1], boost::uniform_01<double>());	
	
	
    std::generate(sample.begin(), sample.end(), gen1);
	
	for(int i=0; i< sample.size(); i++){
		
		std::cout<<"on j "<<j<<" Sample "<<i<<" = "<< sample[i]<<std::endl;
	}	
	
	
	
}
	
	
}


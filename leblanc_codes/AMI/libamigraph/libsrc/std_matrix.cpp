#include<iostream>
#include<iomanip>

#include <complex>
#include <math.h>
#include <stdlib.h>
#include <random>

#include "amigraph.hpp"
// void triangularize_matrix( std::vector< std::vector< int >> &M_IN, std::vector< std::vector< int >> &M_out);


void triangularize_matrix( std::vector< std::vector< int >> M_IN, std::vector< std::vector< int >> &M_out)
{
    int Nx,Ny;
	
	if(M_IN.size()==0){return;}
	
	Nx=M_IN.size();
	Ny=M_IN[0].size();
	
	// Triangularization
for (int i = 0; i < Nx - 1; i++){
    for (int h = i + 1; h < Nx; h++)
    {
        double t = M_IN[h][i] / M_IN[i][i];
        for (int j = 0; j <= Nx; j++)
        {
            M_IN[h][j] = M_IN[h][j] - t * M_IN[i][j];
        }
    }
}	
	
M_out=M_IN;    
		  
		  
    return;
}
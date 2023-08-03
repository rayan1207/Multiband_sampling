cmake_minimum_required(VERSION 3.1)
project(mini_ami CXX)

find_package(Boost 1.71.0)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
endif()


find_package( Eigen3 3.3 REQUIRED NO_MODULE)
include_directories(${EIGEN3_INCLUDE_DIR})

#add_executable(${PROJECT_NAME} main.cpp ising.cpp)
link_directories(${Boost_LIBRARY_DIRS})


add_executable(${PROJECT_NAME} 
./src/mini_main.cpp 
./src/sparse_sampler_sigma.cpp
./src/helper.cpp
./src/molecule_solver.cpp
./src/mini_ami.hpp
./src/sparse_pp_sampler.cpp
)


SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -DUSE_MPI -lstdc++fs")# -msse4.2")# -ffast-math -funsafe-math-optimizations -msse4.2")# -ffast-math")


link_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../libami/install/lib)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../libami/install/include)



include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../libamicalc/install/include)
link_directories(${CMAKE_CURRENT_SOURCE_DIR}/../libamicalc/install/lib)
link_directories(${CMAKE_CURRENT_SOURCE_DIR}/../libamigraph/install/lib)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../libamigraph/install/include)


#set_target_properties(${EXEC_NAME} PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
# Use ALPSCore_LIBRARIES variable to link to ALPSCore 
#target_link_libraries(${PROJECT_NAME} ${ALPSCore_LIBRARIES})

#NOTE: The relative paths in libami and libamigraph compilation don't allow cmake to install. fix later
target_link_libraries(${PROJECT_NAME} ${ALPSCore_LIBRARIES} ${CMAKE_CURRENT_SOURCE_DIR}/../../../libami/install/lib/libami.so ${CMAKE_CURRENT_SOURCE_DIR}/../libamicalc/install/lib/libamicalc.so
${CMAKE_CURRENT_SOURCE_DIR}/../libamigraph/install/lib/libamigraph.so m stdc++fs boost_filesystem boost_system)

install (TARGETS ${EXEC_NAME} DESTINATION ${CMAKE_INSTALL_PREFIX}/bin )
#install(TARGETS dmft RUNTIME DESTINATION "bin/")

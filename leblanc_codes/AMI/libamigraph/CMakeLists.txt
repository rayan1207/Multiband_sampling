cmake_minimum_required(VERSION 2.8)

project (graph)
link_directories("/usr/lib")

find_package(Boost)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
endif()

find_package( Eigen3 3.3 REQUIRED NO_MODULE)
include_directories(${EIGEN3_INCLUDE_DIR})

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set(MACOSX TRUE)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -stdlib=libc++")
else()
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -lstdc++fs")
endif()


# include_directories(${AMI_INCLUDE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../libami/install/include)
link_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../libami/install/lib)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../libamicalc/install/include)
link_directories(${CMAKE_CURRENT_SOURCE_DIR}/../libamicalc/install/lib)
# include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../Cuba-3.0/install/include)
link_directories(${Boost_LIBRARY_DIRS})

set(LIB_FILES 
 #	./libsrc/graph_io.cpp 
#	./libsrc/insert_remove.cpp 
#	./libsrc/graph_helper.cpp 
#	./libsrc/graph_constructs.cpp
	./libsrc/amigraph_constructs.cpp
	./libsrc/amigraph_helper.cpp
	./libsrc/amigraph_updates.cpp
	./libsrc/amigraph_updates_ABDB.cpp
	./libsrc/amigraph_io.cpp
	./libsrc/amigraph_labelling.cpp
	./libsrc/amigraph_amiconvert.cpp
	./libsrc/amigraph_filter.cpp
	./libsrc/amigraph_isomorphism.cpp
	./libsrc/amigraph_graph_generator.cpp
	./libsrc/amigraph_rf.cpp
	./libsrc/amigraph_gos.cpp
	./libsrc/amigraph_eigen.cpp
	./libsrc/amigraph_gg.cpp
	./libsrc/amigraph_fileio.cpp
	./libsrc/amigraph_git.cpp
	./libsrc/amigraph_git_permutations.cpp
	./libsrc/amigraph_syslabel.cpp
	./libsrc/amigraph_git_eval.cpp
	./libsrc/amigraph_ct.cpp
	# ./libsrc/amigraph_wicks.cpp
	./libsrc/amigraph_wicks.cpp
	./libsrc/amigraph_bc.cpp
	# ./libsrc/amigraph_ct.cpp
	./libsrc/std_matrix.cpp
	./libsrc/wicks_label.cpp
	./libsrc/wicks_eigen.cpp
	
  )


ADD_LIBRARY(amigraph SHARED ${LIB_FILES})


if(MACOSX)
	target_link_libraries(amigraph ${CMAKE_CURRENT_SOURCE_DIR}/../../../libami/install/lib/libami.dylib ${CMAKE_CURRENT_SOURCE_DIR}/../libamicalc/install/lib/libamicalc.dylib)
endif()

INSTALL(TARGETS amigraph 
   ARCHIVE DESTINATION "lib/"
   LIBRARY DESTINATION "lib/"
   COMPONENT library)

#FILE(GLOB files "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
INSTALL(DIRECTORY "${CMAKE_SOURCE_DIR}/libsrc/." # source directory
        DESTINATION include # target directory
        FILES_MATCHING PATTERN "*.hpp" # select header files
)

# # add_executable(${PROJECT_NAME}_test 
	# # ./libsrc/main.cpp 
	# # ./libsrc/amigraph_constructs.cpp 
	# # ./libsrc/amigraph_io.cpp 
	# # ./libsrc/amigraph_updates.cpp 
	# # ./libsrc/amigraph_updates_ABDB.cpp
	# # ./libsrc/amigraph_helper.cpp
	# # ./libsrc/amigraph_labelling.cpp
	# # ./libsrc/amigraph_amiconvert.cpp
	# # ./libsrc/amigraph_filter.cpp
	# # ./libsrc/amigraph_isomorphism.cpp
	# # ./libsrc/amigraph_graph_generator.cpp
	# # ./libsrc/amigraph_rf.cpp
	# # ./libsrc/amigraph_gos.cpp
	# # ./libsrc/amigraph_eigen.cpp
	# # ./libsrc/amigraph_gg.cpp
	# # ./libsrc/amigraph_fileio.cpp
	# # ./libsrc/amigraph_git.cpp
	# # ./libsrc/amigraph_git_permutations.cpp
	# # ./libsrc/amigraph_syslabel.cpp
	# # ./libsrc/amigraph_git_eval.cpp
	# # ./libsrc/amigraph_ct.cpp
	# # # ./libsrc/amigraph_wicks.cpp
	# # )
	




# # target_link_libraries(${PROJECT_NAME}_test ${CMAKE_CURRENT_SOURCE_DIR}/../libami/install/lib/libamibase.so m stdc++fs boost_system boost_filesystem )
# # install (TARGETS ${EXEC_NAME} DESTINATION bin )

#if(!MACOSX)
add_executable(bubble_gen 
	./libsrc/bubble_main.cpp 
	./libsrc/amigraph_constructs.cpp 
	./libsrc/amigraph_io.cpp 
	./libsrc/amigraph_updates.cpp 
	./libsrc/amigraph_updates_ABDB.cpp
	./libsrc/amigraph_helper.cpp
	./libsrc/amigraph_labelling.cpp
	./libsrc/amigraph_amiconvert.cpp
	./libsrc/amigraph_filter.cpp
	./libsrc/amigraph_isomorphism.cpp
	./libsrc/amigraph_graph_generator.cpp
	./libsrc/amigraph_rf.cpp
	./libsrc/amigraph_gos.cpp
	./libsrc/amigraph_eigen.cpp
	./libsrc/amigraph_gg.cpp
	./libsrc/amigraph_fileio.cpp
	./libsrc/amigraph_git.cpp
	./libsrc/amigraph_git_permutations.cpp
	./libsrc/amigraph_syslabel.cpp
	./libsrc/amigraph_git_eval.cpp
	./libsrc/amigraph_ct.cpp
	
	# ./libsrc/amigraph_wicks.cpp
	)
	




target_link_libraries(bubble_gen ${CMAKE_CURRENT_SOURCE_DIR}/../../../libami/install/lib/libami.so ${CMAKE_CURRENT_SOURCE_DIR}/../libamicalc/install/lib/libamicalc.so m stdc++fs boost_system boost_filesystem )
install (TARGETS ${EXEC_NAME} DESTINATION bin )

#endif()

# add_executable(wickstest
# ./libsrc/wicksmain.cpp 
	# ./libsrc/amigraph_constructs.cpp 
	# ./libsrc/amigraph_io.cpp 
	# ./libsrc/amigraph_updates.cpp 
	# ./libsrc/amigraph_updates_ABDB.cpp
	# ./libsrc/amigraph_helper.cpp
	# ./libsrc/amigraph_labelling.cpp
	# ./libsrc/amigraph_amiconvert.cpp
	# ./libsrc/amigraph_filter.cpp
	# ./libsrc/amigraph_isomorphism.cpp
	# ./libsrc/amigraph_graph_generator.cpp
	# ./libsrc/amigraph_rf.cpp
	# ./libsrc/amigraph_gos.cpp
	# ./libsrc/amigraph_eigen.cpp
	# ./libsrc/amigraph_gg.cpp
	# ./libsrc/amigraph_fileio.cpp
	# ./libsrc/amigraph_git.cpp
	# ./libsrc/amigraph_git_permutations.cpp
	# ./libsrc/amigraph_syslabel.cpp
	# ./libsrc/amigraph_git_eval.cpp
	# ./libsrc/amigraph_wicks.cpp
	# ./libsrc/amigraph_ct.cpp
	# ./libsrc/std_matrix.cpp
	# ./libsrc/wicks_label.cpp
	# ./libsrc/wicks_eigen.cpp
# )
# target_link_libraries(wickstest ${CMAKE_CURRENT_SOURCE_DIR}/../libami/install/lib/libamibase.so m stdc++fs boost_system boost_filesystem Eigen3::Eigen  )
# install (TARGETS ${EXEC_NAME} DESTINATION bin )


# #########

# add_executable(wicks_sc
# ./libsrc/wicks_sc.cpp 
	# ./libsrc/amigraph_constructs.cpp 
	# ./libsrc/amigraph_io.cpp 
	# ./libsrc/amigraph_updates.cpp 
	# ./libsrc/amigraph_updates_ABDB.cpp
	# ./libsrc/amigraph_helper.cpp
	# ./libsrc/amigraph_labelling.cpp
	# ./libsrc/amigraph_amiconvert.cpp
	# ./libsrc/amigraph_filter.cpp
	# ./libsrc/amigraph_isomorphism.cpp
	# ./libsrc/amigraph_graph_generator.cpp
	# ./libsrc/amigraph_rf.cpp
	# ./libsrc/amigraph_gos.cpp
	# ./libsrc/amigraph_eigen.cpp
	# ./libsrc/amigraph_gg.cpp
	# ./libsrc/amigraph_fileio.cpp
	# ./libsrc/amigraph_git.cpp
	# ./libsrc/amigraph_git_permutations.cpp
	# ./libsrc/amigraph_syslabel.cpp
	# ./libsrc/amigraph_git_eval.cpp
	# ./libsrc/amigraph_wicks.cpp
	# ./libsrc/amigraph_ct.cpp
	# ./libsrc/std_matrix.cpp
	# ./libsrc/wicks_label.cpp
	# ./libsrc/wicks_eigen.cpp
# )
# target_link_libraries(wicks_sc ${CMAKE_CURRENT_SOURCE_DIR}/../libami/install/lib/libamibase.so m stdc++fs boost_system boost_filesystem Eigen3::Eigen  )
# install (TARGETS ${EXEC_NAME} DESTINATION bin )


# ##--------------------------------##
# SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3 -std=c++11 -lstdc++fs -fopenmp")#  -ffast-math -funsafe-math-optimizations -msse4.2")
# # add_executable(ami_cuba_mpi 
	# # ./libsrc/ami_cuba_mpi.cpp 
	# # ./libsrc/amigraph_constructs.cpp 
	# # ./libsrc/amigraph_io.cpp 
	# # ./libsrc/amigraph_updates.cpp 
	# # ./libsrc/amigraph_updates_ABDB.cpp
	# # ./libsrc/amigraph_helper.cpp
	# # ./libsrc/amigraph_labelling.cpp
	# # ./libsrc/amigraph_amiconvert.cpp
	# # ./libsrc/amigraph_filter.cpp
	# # ./libsrc/amigraph_isomorphism.cpp
	# # ./libsrc/amigraph_graph_generator.cpp
	# # ./libsrc/amigraph_rf.cpp
	# # ./libsrc/amigraph_gos.cpp
	# # ./libsrc/amigraph_eigen.cpp
	# # ./libsrc/amigraph_gg.cpp
	# # ./libsrc/amigraph_fileio.cpp
	# # ./libsrc/amigraph_git.cpp
	# # ./libsrc/amigraph_git_permutations.cpp
	# # ./libsrc/amigraph_syslabel.cpp
	# # ./libsrc/amigraph_git_eval.cpp
	# # ./libsrc/amigraph_ct.cpp
	# # )
	

	
	
# # target_link_libraries(ami_cuba_mpi ${CMAKE_CURRENT_SOURCE_DIR}/../libami/install/lib/libamibase.so m stdc++fs Eigen3::Eigen ${CMAKE_CURRENT_SOURCE_DIR}/../../../Cuba-3.0/install/lib/libcuba.a boost_system boost_filesystem)
# # install(TARGETS ${EXEC_NAME} DESTINATION ${CMAKE_INSTALL_PREFIX}/bin )


# # add_executable(newiso ./libsrc/iso.cpp )

# # add_executable(soboltest ./libsrc/sobol.cpp)

# ##--------------------------------##


	# # add_executable(git_compare 
	# # ./libsrc/git_main.cpp 
	# # ./libsrc/amigraph_constructs.cpp 
	# # ./libsrc/amigraph_io.cpp 
	# # ./libsrc/amigraph_updates.cpp 
	# # ./libsrc/amigraph_updates_ABDB.cpp
	# # ./libsrc/amigraph_helper.cpp
	# # ./libsrc/amigraph_labelling.cpp
	# # ./libsrc/amigraph_amiconvert.cpp
	# # ./libsrc/amigraph_filter.cpp
	# # ./libsrc/amigraph_isomorphism.cpp
	# # ./libsrc/amigraph_graph_generator.cpp
	# # ./libsrc/amigraph_rf.cpp
	# # ./libsrc/amigraph_gos.cpp
	# # ./libsrc/amigraph_eigen.cpp
	# # ./libsrc/amigraph_gg.cpp
	# # ./libsrc/amigraph_fileio.cpp
	# # ./libsrc/amigraph_git.cpp
	# # ./libsrc/amigraph_git_permutations.cpp
	# # ./libsrc/amigraph_syslabel.cpp
	# # ./libsrc/amigraph_git_eval.cpp
	# # )
	

# # target_link_libraries(git_compare ${CMAKE_CURRENT_SOURCE_DIR}/../libami/install/lib/libami.so m stdc++fs Eigen3::Eigen  boost_system boost_filesystem)
# # install(TARGETS ${EXEC_NAME} DESTINATION ${CMAKE_INSTALL_PREFIX}/bin )


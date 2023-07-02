cmake_minimum_required(VERSION 3.18...3.22.1)

project (amicalc CXX)
#find_package(ALPSCore REQUIRED)
#For efficiency, compile with optimization by default
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RelWithDebInfo")
endif()

if (NOT BUILD_DOC)
	set(BUILD_DOC OFF)
endif()

find_package(BOOST)
if(Boost_FOUND)
	include_directories(${Boost_INCLUDE_DIRS})
endif()

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set(MACOSX TRUE)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -stdlib=libc++")
else()
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -lstdc++fs")
endif()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
link_directories("/usr/lib")





set(LIB_FILES 
  ./src/ami_calc_io.cpp
  ./src/ami_calc_evaluate.cpp
	./src/ami_calc_spectral.cpp
	./src/ami_calc_terms_spectral.cpp
	./src/ami_spec.cpp
	./src/ami_spec_io.cpp
	./src/ami_spec_sigma.cpp
	./src/ami_spec_eval.cpp
	./src/ami_spec_opt.cpp
  )

ADD_LIBRARY(amicalc SHARED ${LIB_FILES})


#target_link_libraries(${PROJECT_NAME} ${ALPSCore_LIBRARIES} m)
if(MACOSX)
target_link_libraries(${PROJECT_NAME}  m  ${CMAKE_CURRENT_SOURCE_DIR}/../../../libami/install/lib/libami.dylib)
else()
	target_link_libraries(${PROJECT_NAME} m stdc++fs)
endif(MACOSX)


include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../libami/install/include)
link_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../libami/install/lib)

   
   INSTALL(TARGETS amicalc 
   ARCHIVE DESTINATION "lib/"
   LIBRARY DESTINATION "lib/"
   COMPONENT library)

INSTALL(DIRECTORY "${CMAKE_SOURCE_DIR}/src/." # source directory
        DESTINATION include # target directory
        FILES_MATCHING 
		PATTERN "*.hpp"
		PATTERN	"*.h" # select header files
)







if(BUILD_DOC)

add_subdirectory(test)
enable_testing()
add_test(fb test/fb)
add_test(helper_functions test/helper_functions)
add_test(num_test test/num_test)
add_test(construct_tests test/construct_tests)

message("Adding doc subdirectory")
add_subdirectory("doc")

#########################################################
#Documentation
find_package(Doxygen)
if(DOXYGEN_FOUND)
  # exclude sqlite code
  set(DOXYGEN_EXCLUDE_PATTERNS
        */sqlite3/*
  )
  message("run 'make DOCS' to generate documentation")
  # doxygen settings can be set here, prefixed with "DOXYGEN_"
  set(DOXYGEN_EXTRACT_PRIVATE YES)
  set(DOXYGEN_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/DOCS")
  set(DOXYGEN_GENERATE_LATEX YES)
  set(DOXYGEN_GENERATE_XML YES)
  # this target will only be built if specifically asked to.
  # run "make LibAMI-docs" to create the doxygen documentation
  # run "make DOCS" to create the doxygen documentation
  doxygen_add_docs(
    DOCS
    ${PROJECT_SOURCE_DIR}/src
    COMMENT "Generate API-documents for LibAMI."
  )
endif(DOXYGEN_FOUND)


endif(BUILD_DOC)






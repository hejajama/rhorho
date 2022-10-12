configure_file("${CMAKE_CURRENT_SOURCE_DIR}/gitsha1.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/gitsha1.cpp" @ONLY)

add_library(
	rholib
	STATIC
	functions.cpp
	vector.cpp
	proton.cpp
	diagram_integrator.cpp
	diagram_sum.cpp
	interpolation.cpp
)

add_executable(
	rhorho
	main.cpp
	${CMAKE_CURRENT_BINARY_DIR}/gitsha1.cpp
)

add_executable(
	rhotest
	function_tests.cpp
	vector_tests.cpp
	test_main.cpp
	proton_tests.cpp
)

target_link_libraries(
	rhorho
	PRIVATE
	rholib
	GSL::gsl
	GSL::gslcblas
)
target_link_libraries(
	rhotest
	PRIVATE
	rholib
	${Boost_LIBRARIES}
	GSL::gsl
	GSL::gslcblas
)

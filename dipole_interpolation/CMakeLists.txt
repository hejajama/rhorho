configure_file("${CMAKE_CURRENT_SOURCE_DIR}/gitsha1.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/gitsha1.cpp" @ONLY)

add_library(
	dipolelib
	STATIC
	dipoleamplitude.cpp
	interpolation2d.cpp
)

add_executable(
	dipole
	main.cpp
	${CMAKE_CURRENT_BINARY_DIR}/gitsha1.cpp
)

target_link_libraries(
	dipole
	PRIVATE
	dipolelib
	GSL::gsl
	GSL::gslcblas
)



enable_language(Fortran)

FILE(GLOB slsrc "*.f")
add_library(depblas ${slsrc})

INSTALL(TARGETS depblas DESTINATION lib)


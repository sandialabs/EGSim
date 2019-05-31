
enable_language(Fortran)

FILE(GLOB slsrc "*.f")
add_library(deplapack ${slsrc})

INSTALL(TARGETS deplapack DESTINATION lib)

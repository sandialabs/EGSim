FILE(GLOB slsrc "*.cpp")

# egsim
add_executable (egsim ${slsrc})

target_link_libraries (egsim deplapack)
target_link_libraries (egsim depblas  )

# Link fortran libraries
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  target_link_libraries (egsim gfortran)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel
  target_link_libraries (egsim ifcore)
  target_link_libraries (egsim ifport)
endif()

add_definitions(-DIEEE_STRICT)
add_definitions(-g)

include_directories(../pg_dist/app/egsim)

INSTALL(TARGETS egsim DESTINATION bin)

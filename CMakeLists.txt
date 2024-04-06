cmake_minimum_required(VERSION 3.15)
project(structblas Fortran)

# ------------------------------------------------------------------------------ #
# Find FFTW                                                                      #
# ------------------------------------------------------------------------------ #
find_package(OpenMP REQUIRED)

find_package(PkgConfig REQUIRED)     
pkg_check_modules(FFTW IMPORTED_TARGET REQUIRED fftw3)
          
if(NOT FFTW_ROOT AND DEFINED ENV{FFTWDIR})
    set(FFTW_ROOT $ENV{FFTWDIR})
endif()
          
find_library(
    FFTW_DOUBLE_THREADS_LIB
    NAMES "fftw3_threads"
    PATHS ${PKG_FFTW_LIBRARY_DIRS} ${LIB_INSTALL_DIR}
)
          
if (FFTW_DOUBLE_THREADS_LIB)
    set(FFTW_DOUBLE_THREADS_LIB_FOUND TRUE)
    set(FFTW_LIBRARIES ${FFTW_LIBRARIES} ${FFTW_DOUBLE_THREADS_LIB})
    add_library(FFTW::DoubleThreads INTERFACE IMPORTED)
    set_target_properties(FFTW::DoubleThreads
        PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${FFTW_INCLUDE_DIRS}"
        INTERFACE_LINK_LIBRARIES "${FFTW_DOUBLE_THREADS_LIB}"
    )
else()
    set(FFTW_DOUBLE_THREADS_LIB_FOUND FALSE)
endif()
     
# ------------------------------------------------------------------------------ #
# Fortuno                                                                        #
# This is used for the unit tests                                                #
# ------------------------------------------------------------------------------ #
include(FetchContent)
FetchContent_Declare(
  Fortuno
  GIT_REPOSITORY "https://github.com/fortuno-repos/fortuno"
  GIT_TAG "main"
)
FetchContent_MakeAvailable(Fortuno)

# ------------------------------------------------------------------------------ #
# Enable Optimization                                                            #
# Add option to enable optimization -------------------------------------------- #
option(ENABLE_OPTIMIZATION "Enable optimization" ON)

# Set compiler flags
if(ENABLE_OPTIMIZATION)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3")
else()
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O0 -g -ggdb")
endif()

# ------------------------------------------------------------------------------ #
# Link the library with FFTW                                                     #
# ------------------------------------------------------------------------------ #

# Add all the Fortran source files in the src directory
file(GLOB SOURCES "src/*.f90")
# Add the library target
add_library(structblas ${SOURCES})
target_link_libraries(structblas PUBLIC PkgConfig::FFTW ${FFTW_DOUBLE_THREADS_LIB})

# ------------------------------------------------------------------------------- #
# Tests                                                                           #
# ------------------------------------------------------------------------------- #

enable_testing()
file(GLOB TEST_SOURCES "test/*.f90")

# Add the executable target for each test file
foreach(TEST_SOURCE ${TEST_SOURCES})
    get_filename_component(TEST_NAME ${TEST_SOURCE} NAME_WE)
    add_executable(${TEST_NAME} ${TEST_SOURCE} ${LIB_SOURCES})
    target_link_libraries(${TEST_NAME} PRIVATE structblas 
        Fortuno::fortuno_serial 
        PkgConfig::FFTW ${FFTW_DOUBLE_THREADS_LIB})
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
endforeach()


# Set any compiler flags you need
# For example, you might want to set optimization flags or warn flags.
# set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -Wall")

# ------------------------------------------------------------------------------- #
# Add option to generate documentation                                            #
# ------------------------------------------------------------------------------- #
option(GENERATE_DOCUMENTATION "Generate documentation using Ford" ON)
if(GENERATE_DOCUMENTATION)
    find_program(FORD_EXECUTABLE ford)
    if(NOT FORD_EXECUTABLE)
        message(FATAL_ERROR "Ford executable not found. Make sure it is installed and accessible in your PATH.")
    endif()

    add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/ford_output
                       COMMAND ${FORD_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/ford.md
                       WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                       COMMENT "Generating documentation with Ford"
                       VERBATIM)

    add_custom_target(documentation DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/ford_output)

    install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/ford_output
            DESTINATION ${CMAKE_INSTALL_DOCDIR}/ford_output)
endif()


# ------------------------------------------------------------------------------- #
# Set where to install the library                                                #
# ------------------------------------------------------------------------------- #
# Set where to install the library
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/install CACHE PATH "Installation directory" FORCE)
endif()

install(TARGETS structblas DESTINATION lib)
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/
         DESTINATION include
         FILES_MATCHING PATTERN "*.mod")
foreach(TEST_SOURCE ${TEST_SOURCES})
    get_filename_component(TEST_NAME ${TEST_SOURCE} NAME_WE)
    install(TARGETS ${TEST_NAME} DESTINATION test) 
endforeach()
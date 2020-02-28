
cmake_minimum_required(VERSION 3.5)

set(PROJECT_NAME FisherSolver)
set(EXECUTABLE_OUTPUT_PATH ../../bin)
project(${PROJECT_NAME})

# Set CMake behavior
cmake_policy(SET CMP0004 NEW)

# inlcude boost
FIND_PACKAGE(Boost)
IF (Boost_FOUND)
    INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
    ADD_DEFINITIONS( "-DHAS_BOOST" )
ENDIF()

# Get DOLFIN configuration data (DOLFINConfig.cmake must be in
# DOLFIN_CMAKE_CONFIG_PATH)
find_package(DOLFIN REQUIRED)

if (EXISTS ${DOLFIN_USE_FILE})
  include(${DOLFIN_USE_FILE})

  # Default build type (can be overridden by user)
  if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING
      "Choose the type of build, options are: Debug MinSizeRel Release RelWithDebInfo." FORCE)
  endif()

  # Do not throw error for 'multi-line comments' (these are typical in
  # rst which includes LaTeX)
  include(CheckCXXCompilerFlag)
  CHECK_CXX_COMPILER_FLAG("-Wno-comment" HAVE_NO_MULTLINE)
  if (HAVE_NO_MULTLINE)
    set(CMAKE_CXX_FLAGS "-Wno-comment ${CMAKE_CXX_FLAGS}")
  endif()

else()
  # Compiler definitions
  add_definitions(${DOLFIN_CXX_DEFINITIONS})

  # Compiler flags
  set(CMAKE_CXX_FLAGS "${DOLFIN_CXX_FLAGS} ${CMAKE_CXX_FLAGS}")

  # Do not throw error for 'multi-line comments' (these are typical in
  # rst which includes LaTeX)
  include(CheckCXXCompilerFlag)
  CHECK_CXX_COMPILER_FLAG("-Wno-comment" HAVE_NO_MULTLINE)
  if (HAVE_NO_MULTLINE)
    set(CMAKE_CXX_FLAGS "-Wno-comment ${CMAKE_CXX_FLAGS}")
  endif()

  # Include directories
  include_directories(${DOLFIN_INCLUDE_DIRS})
  include_directories(SYSTEM ${DOLFIN_3RD_PARTY_INCLUDE_DIRS})
endif()



# Executable
add_executable(${PROJECT_NAME} main.cpp FisherNewtonContainer.cpp RuntimeTracker.cpp TimeStepper.cpp FisherProblem.cpp Initializers.cpp Tensors.cpp ReaderWriter.cpp)

# Target libraries
target_link_libraries(${PROJECT_NAME} ${DOLFIN_LIBRARIES})

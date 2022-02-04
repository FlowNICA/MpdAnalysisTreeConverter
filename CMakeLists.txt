# CMakeLists.txt for FAS package. It creates a library with dictionary and a main program
cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(MpdAnalysisTreeConverter)

find_package(Git)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Net MathMore EG)

#---Locate AnalysisTree library
find_library(AnalysisTree NAMES AnalysisTreeBase PATHS $ENV{ANALYSISTREE_LIB})
if (AnalysisTree)
  message(STATUS "AnalysisTree is found: ${AnalysisTree}")
endif()
if (NOT AnalysisTree)
  message(FATAL_ERROR "Error: AnalysisTree library wasn't found!")
endif()

#---Locate FairSoft/FairRoot package
set(FAIRSOFT_PATH $ENV{FAIRROOTPATH}/share/fairbase/cmake/modules)
list(APPEND CMAKE_MODULE_PATH ${FAIRSOFT_PATH})
find_package(FairRoot QUIET COMPONENTS Base)
if(NOT FAIRROOT_FOUND)
  message(FATAL_ERROR "FairRoot package NOT found.")
endif()
set(FAIRROOT_INCLUDE_DIR $ENV{FAIRROOTPATH}/include)
set(FAIRROOT_LIBRARY_DIR $ENV{FAIRROOTPATH}/lib)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

add_definitions(${ROOT_CXX_FLAGS})

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -Wall -g")

if(NOT DEFINED ENV{VMCWORKDIR})
  message(FATAL_ERROR "MpdRoot environment is NOT found. Skip corresponding step.")
endif()

set(MPDROOTDIR $ENV{VMCWORKDIR})
message(STATUS "MpdRoot environment is found: ${MPDROOTDIR}")

set(MPD_AT_CONVERTER_INCLUDE_DIRECTORIES
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${ROOT_INLCUDE_DIRS}
  $ENV{ANALYSISTREE_INC}
  ${FAIRROOT_INCLUDE_DIR}
  ${MPDROOTDIR}/include
)

include_directories(${MPD_AT_CONVERTER_INCLUDE_DIRECTORIES})

#---Add Needed MpdRoot/FairRoot libraries
find_library(BASE NAMES Base PATHS ${FAIRROOT_LIBRARY_DIR})
find_library(MCSTACK NAMES MCStack PATHS $ENV{FAIRLIBDIR})
if (MCSTACK)
  message(STATUS "MCStack library is found")
  add_definitions("-D_MCSTACK_")
endif()
if (NOT MCSTACK)
  find_library(MPDMCSTACK NAMES MpdMCStack PATHS $ENV{FAIRLIBDIR})
  message(STATUS "MpdMCStack library is found")
  add_definitions("-D_MPDMCSTACK_")
endif()
find_library(PASSIVE NAMES Passive PATHS $ENV{FAIRLIBDIR})
find_library(MPDFIELD NAMES MpdField PATHS $ENV{FAIRLIBDIR})
find_library(MPDBASE NAMES MpdBase PATHS $ENV{FAIRLIBDIR})
find_library(MPDDST NAMES MpdDst PATHS $ENV{FAIRLIBDIR})
find_library(MPDGENERALGENERATOR NAMES MpdGeneralGenerator PATHS $ENV{FAIRLIBDIR})
find_library(KALMAN NAMES Kalman PATHS $ENV{FAIRLIBDIR})
find_library(CLUSTER NAMES Cluster PATHS $ENV{FAIRLIBDIR})
find_library(TPC NAMES tpc PATHS $ENV{FAIRLIBDIR})
find_library(TOF NAMES Tof PATHS $ENV{FAIRLIBDIR})
find_library(ETOF NAMES Etof PATHS $ENV{FAIRLIBDIR})
# find_library(STRAWENDCAP NAMES Strawendcap PATHS $ENV{FAIRLIBDIR})
find_library(ZDC NAMES Zdc PATHS $ENV{FAIRLIBDIR})
find_library(FFD NAMES Ffd PATHS $ENV{FAIRLIBDIR})
find_library(CPC NAMES Cpc PATHS $ENV{FAIRLIBDIR})
# find_library(STRAWECT NAMES StrawECT PATHS $ENV{FAIRLIBDIR})
# find_library(LHETRACK NAMES LHETrack PATHS $ENV{FAIRLIBDIR})
find_library(STS NAMES Sts PATHS $ENV{FAIRLIBDIR})
find_library(MPDPID NAMES MpdPid PATHS $ENV{FAIRLIBDIR})

#---Get all needed libraries together
if(MCSTACK)
set(MPD_AT_CONVERTER_INCLUDE_LIBRARIES
  ${AnalysisTree} ${ROOT_LIBRARIES} ${BASE}
  ${MCSTACK} ${PASSIVE} ${KALMAN} ${CLUSTER} ${TPC} ${STS} ${ZDC} ${FFD}
  # ${STRAWECT} ${STRAWENDCAP} 
  # ${LHETRACK}
  ${TOF}
  ${TOF}
  ${MPDFIELD} 
  ${MPDBASE} 
  # ${MPDDST} 
  ${MPDGENERALGENERATOR} 
  ${MPDPID}
)
elseif(MPDMCSTACK)
set(MPD_AT_CONVERTER_INCLUDE_LIBRARIES
  ${AnalysisTree} ${ROOT_LIBRARIES} ${BASE}
  ${MPDMCSTACK} ${PASSIVE} ${KALMAN} ${CLUSTER} ${TPC} ${STS} ${ZDC} ${FFD}
  # # ${STRAWECT} ${STRAWENDCAP} 
  # ${LHETRACK}
  ${TOF}
  ${MPDFIELD} 
  ${MPDBASE} 
  # ${MPDDST} 
  ${MPDGENERALGENERATOR} 
  ${MPDPID}
)
elseif()
  message(FATAL_ERROR "No MCSTACK nor MPDMCSTACK libraries was found!")
endif()

#---Compile converter's executable file
add_executable(MpdDst2AnalysisTree "${CMAKE_CURRENT_SOURCE_DIR}/MpdDstAnalysisTree.cpp")
#target_include_directories(MpdDst2AnalysisTree PUBLIC ${MPD_AT_CONVERTER_INCLUDE_DIRECTORIES})
target_link_libraries(MpdDst2AnalysisTree ${MPD_AT_CONVERTER_INCLUDE_LIBRARIES})

# #---Compile dca getter's executable file
# add_executable(get-dca "${CMAKE_CURRENT_SOURCE_DIR}/get-dca.cpp")
# target_link_libraries(get-dca ${AnalysisTree} ${ROOT_LIBRARIES})

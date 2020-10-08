CMAKE_MINIMUM_REQUIRED(VERSION 3.1 FATAL_ERROR)
if (POLICY CMP0025)
  cmake_policy(SET CMP0025 NEW)
endif ()
PROJECT(VEGAS1D CXX)

SET(VEGAS1D_VERSION "0.1.")

SET(ROOTINCLUDEDIR "/usr/local/Cellar/root/6.20.04_1/include/root")
SET(ROOTLIBDIR "/usr/local/Cellar/root/6.20.04_1/lib/root")

INCLUDE_DIRECTORIES(${ROOTINCLUDEDIR})
LINK_DIRECTORIES(${ROOTLIBDIR})
# LIST(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
# FIND_PACKAGE(ROOT REQUIRED)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pthread -m64 -std=c++11 -stdlib=libc++")

ADD_EXECUTABLE(RandomTest RandomGeneratorTest.cpp)
ADD_EXECUTABLE(ImproveGridTest ImproveGridTest.cpp OneDIntegralSimple.cpp)
ADD_EXECUTABLE(OneDTest OneDTest.cpp OneDIntegralSimple.cpp)
TARGET_LINK_LIBRARIES(RandomTest Core Imt RIO Net Hist Graf Graf3d Gpad ROOTVecOps Tree TreePlayer Rint Postscript Matrix Physics MathCore Thread MultiProc ROOTDataFrame m dl pthread)
TARGET_LINK_LIBRARIES(ImproveGridTest Core Imt RIO Net Hist Graf Graf3d Gpad ROOTVecOps Tree TreePlayer Rint Postscript Matrix Physics MathCore Thread MultiProc ROOTDataFrame m dl pthread)
TARGET_LINK_LIBRARIES(OneDTest Core Imt RIO Net Hist Graf Graf3d Gpad ROOTVecOps Tree TreePlayer Rint Postscript Matrix Physics MathCore Thread MultiProc ROOTDataFrame m dl pthread)
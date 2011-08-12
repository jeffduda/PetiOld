##############################################################################
#
# Project:   Petiole
#
# Copyright (c) 2011
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
# 1. Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.
# 3. Neither the name of the consortium nor the names of its contributors
#   may be used to endorse or promote products derived from this software
#   without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE CONSORTIUM AND CONTRIBUTORS ``AS IS'' AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
# OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
# OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
# SUCH DAMAGE.
#
##############################################################################

PROJECT( PETIOLE )

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
if(COMMAND CMAKE_POLICY)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND CMAKE_POLICY)

include_regular_expression("^([Nn]rrd|itk|vtk|vnl|vcl|vxl|v3p|f2c|netlib|ce|itpack|DICOM|meta|png|dbh|tif|jpeg|zlib|gdcm).*$")
source_group("XML Files" REGULAR_EXPRESSION "[.]xml$")

# On Visual Studio 8 MS deprecated C. This removes all 1.276E1265 security
# warnings
if(WIN32)
  if(not CYGWIN)
    if(not MINGW)
      if(not ITK_ENABLE_VISUAL_STUDIO_DEPRECATED_C_WARNINGS)
        add_definitions(
          -D_CRT_FAR_MAPPINGS_NO_DEPRECATE
          -D_CRT_IS_WCTYPE_NO_DEPRECATE
          -D_CRT_MANAGED_FP_NO_DEPRECATE
          -D_CRT_NONSTDC_NO_DEPRECATE
          -D_CRT_SECURE_NO_DEPRECATE
          -D_CRT_SECURE_NO_DEPRECATE_GLOBALS
          -D_CRT_SETERRORMODE_BEEP_SLEEP_NO_DEPRECATE
          -D_CRT_TIME_FUNCTIONS_NO_DEPRECATE
          -D_CRT_VCCLRIT_NO_DEPRECATE
          -D_SCL_SECURE_NO_DEPRECATE
          )
      endif(not ITK_ENABLE_VISUAL_STUDIO_DEPRECATED_C_WARNINGS)
    endif(not MINGW)
  endif(not CYGWIN)
endif(WIN32)

#-----------------------------------------------------------------------------
# Petiole version number.
set(PETIOLE_VERSION_MAJOR "0")
set(PETIOLE_VERSION_MINOR "1")
set(PETIOLE_VERSION_PATCH "0")

# Version string should not include patch level.  The major.minor is
# enough to distinguish available features of the toolkit.
set(PETIOLE_VERSION_STRING "${PETIOLE_VERSION_MAJOR}.${PETIOLE_VERSION_MINOR}")

#-----------------------------------------------------------------------------
# Dashboard setup

INCLUDE( CTest )
ENABLE_TESTING()
IF( BUILD_TESTING )
  SET( BUILDNAME "${BUILDNAME}" CACHE STRING "Name of build on the dashboard" )
  MARK_AS_ADVANCED( BUILDNAME )
ENDIF( BUILD_TESTING )
SET( BUILD_WARNING_REPORT_LIMIT -1)
SET( BUILD_ERROR_REPORT_LIMIT -1)

INCLUDE(Dart)

#-----------------------------------------------------------------------------
# Optional libraries
#-----------------------------------------------------------------------------
OPTION( USE_VTK "Use VTK Libraries" OFF )

# Set up VTK
IF(USE_VTK)
  FIND_PACKAGE(VTK)
  IF(VTK_FOUND)
    INCLUDE(${VTK_USE_FILE})
  ELSE(VTK_FOUND)
     MESSAGE("Cannot build some programs without VTK.  Please set VTK_DIR if you need these programs.")
  ENDIF(VTK_FOUND)
ENDIF(USE_VTK)


OPTION( USE_BOOST "Use Boost Libraries" ON )
if( USE_BOOST )

  SET(Boost_USE_STATIC_LIBS ON)
  SET(Boost_ADDITIONAL_VERSIONS "1.36.0" "1.37.0")  
  SET(Boost_USE_MULTITHREAD OFF)
  
  FIND_PACKAGE( Boost 1.34.0)

  IF(Boost_FOUND)

    INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})

  ELSE(Boost_FOUND)

    MESSAGE(FATAL_ERROR "Cannot build some programs without Boost.  Please set BOOST_INCLUDEDIR if you need these programs.")

  ENDIF(Boost_FOUND)
ENDIF (USE_BOOST)

include(CTest)
include(CPack)

FIND_PACKAGE( ITK REQUIRED )
INCLUDE( ${ITK_USE_FILE} )

#FIND_PACKAGE( ANTS REQUIRED )
#INCLUDE( ${ANTS_USE_FILE} )

add_subdirectory( UserTools )
#add_subdirectory( Code/Testing )

# ENH: adding explicit bin/lib output
SET(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib/${CMAKE_BUILD_TYPE} CACHE PATH "Library output path")
SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin/${CMAKE_BUILD_TYPE} CACHE PATH "Executable output path")

 
# First Test - Setup
# This must be run first. It sets up the Documentation dir in source dir, and
# copies static documentation source files there.
#ADD_TEST(PETIOLE_SETUP_DOCS ${CMAKE_CURRENT_SOURCE_DIR}/DocumentationStatic/ctestSetupDocs.sh ${CMAKE_CURRENT_SOURCE_DIR} ${EXECUTABLE_OUTPUT_PATH})

# Example Test
# This is a simple example of a test. Define the name of the test here
# and pass the test script to run, followed by any command line params the
# script expects.
#ADD_TEST(PETIOLE_EXAMPLE_TEST ${CMAKE_CURRENT_SOURCE_DIR}/Testing/petioleExampleTest ${EXECUTABLE_OUTPUT_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/ExampleData/data.nii.gz)

# Other Tests
#
#ADD_TEST(PETIOLE_TEST_1  ${CMAKE_CURRENT_SOURCE_DIR}/Testing/petioleCodeTest ${EXECUTABLE_OUTPUT_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/ExampleData/data.nii.gz)

# This compiles results from previous tests so should be after 
# all tests other than web page compilation
#ADD_TEST(PETIOLE_TESTS  ${CMAKE_CURRENT_SOURCE_DIR}/Testing/petioleTests ${EXECUTABLE_OUTPUT_PATH})
  
# Last Test - Build documentation
#
# This must be the last test, to be run after other tests have run and generated 
# the source for the documentation. This builds the documentation into the 
# build dir.
#ADD_TEST(PETIOLE_BUILD_DOCS ${EXECUTABLE_OUTPUT_PATH}/Documentation/ctestBuildDocs.sh)




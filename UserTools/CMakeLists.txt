##############################################################################
#
# Project:   Petiole-UserTools
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

PROJECT( PETIOLE-USERTOOLS )

FIND_PACKAGE( ITK REQUIRED )
INCLUDE( ${USE_ITK_FILE} )

FIND_PACKAGE( VTK REQUIRED )
INCLUDE( ${USE_VTK_FILE} )

#-----------------------------------------------------------------------------#
# Add ANTS command line library
#


#-----------------------------------------------------------------------------#
# Compile ANTS user tools
#
FIND_PATH( ANTS_SOURCE_DIR   /usr/pkg )
FIND_PATH( ANTS_LIBRARY_DIR  /usr/pkg )

SET( PETIOLE_INCLUDE_DIRS
   ${CMAKE_SOURCE_DIR}/Code/Utilities/
   ${CMAKE_SOURCE_DIR}/Code/Graph/
   ${ANTS_SOURCE_DIR}/Code/Utilities/
   ${ANTS_SOURCE_DIR}/Code/Tensor/
   ${ANTS_SOURCE_DIR}/UserTools/
   )

INCLUDE_DIRECTORIES( ${PETIOLE_INCLUDE_DIRS} )

LINK_DIRECTORIES( ${ITK_LIBRARY_PATH} ${VTK_LIBRARY_PATH}
${ANTS_LIBRARY_DIR} )

#-----------------------------------------------------------------------------#
# Add ANTS command line library
#
ADD_LIBRARY( ANTSCommandLine
../Code/Utilities/antsCommandLineParser.cxx
../Code/Utilities/antsCommandLineOption.cxx )

ADD_EXECUTABLE( petioleConnectivity petioleConnectivity.cxx )
TARGET_LINK_LIBRARIES( petioleConnectivity ANTSCommandLine ${ITK_LIBRARIES}  )

ADD_EXECUTABLE( petioleFibers petioleFibers.cxx )
TARGET_LINK_LIBRARIES( petioleFibers ANTSCommandLine ${ITK_LIBRARIES} )

ADD_EXECUTABLE( petioleUnwrapImage petioleUnwrapImage.cxx )
TARGET_LINK_LIBRARIES( petioleUnwrapImage ANTSCommandLine ${ITK_LIBRARIES} )

ADD_EXECUTABLE( petioleWrapImage petioleWrapImage.cxx )
TARGET_LINK_LIBRARIES( petioleWrapImage ANTSCommandLine ${ITK_LIBRARIES} )

# -----------------------------------------------------------------------------#
# Perform corresponding testing
#

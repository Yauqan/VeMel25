CMAKE_MINIMUM_REQUIRED(VERSION 3.13.4)
SET (CMAKE_CXX_STANDARD 20)

# set ( Aspect_DIR ~/bin/aspect-3.0.0/ )
# set ( deal.II_DIR ~/bin/dealii-9.6.2/ )

FIND_PACKAGE(Aspect 3.0.0 QUIET HINTS ${Aspect_DIR} ../ ../../ $ENV{ASPECT_DIR})

IF (NOT Aspect_FOUND)
  MESSAGE(FATAL_ERROR "\n"
	"Could not find a valid ASPECT build/installation directory. "
	"Please specify the directory where you are building ASPECT by passing\n"
	"   -D Aspect_DIR=<path to ASPECT>\n"
	"to cmake or by setting the environment variable ASPECT_DIR in your shell "
	"before calling cmake. See the section 'How to write a plugin' in the "
        "manual for more information.")
ENDIF ()

DEAL_II_INITIALIZE_CACHED_VARIABLES()
SET(TARGET "VeMel25")
PROJECT(${TARGET})
ADD_LIBRARY(${TARGET} SHARED 	src/MM/MM.cc src/HM/HM.cc src/PP/PP.cc src/better-dynamic-topography/better-dynamic-topography.cc src/prescribed-velocity/prescribed-velocity.cc )
target_compile_features(${TARGET} PUBLIC cxx_std_20)
# target_compile_features(${TARGET} PRIVATE -fconcepts)
ASPECT_SETUP_PLUGIN(${TARGET})














































# Copyright (C) 2011 - 2016 by the authors of the ASPECT code.
#
# This file is part of ASPECT.
#
# ASPECT is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# ASPECT is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with ASPECT; see the file LICENSE.  If not see
# <http://www.gnu.org/licenses/>.


# -----------------------------------------------------------------
# This file is a sample CMakeLists.txt file that plugins can use to
# configure and build a shared library that can then be loaded via
# the input file. For more information, see the manual.
#
# For actual use, you may want to use a different name than
# <my_plugin> for the project, and for the names of the source files
# below.
# -----------------------------------------------------------------


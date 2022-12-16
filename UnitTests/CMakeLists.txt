# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
cmake_minimum_required( VERSION 3.20 FATAL_ERROR )

project( FrostTest )

find_package( thinkboxcmlibrary REQUIRED )
include( ThinkboxCMLibrary )
include( PrecompiledHeader )

add_executable( test_frost "" )
target_include_directories( test_frost PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} )

find_package( GTest REQUIRED )
target_include_directories( test_frost PRIVATE ${GTest_INCLUDE_DIRS} )
target_link_libraries( test_frost PRIVATE GTest::GTest )

file( GLOB_RECURSE H_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
      "*.h"
      "*.hpp"
)

file( GLOB_RECURSE CXX_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
      "*.cpp"
      "*.cc"
)

target_sources( test_frost PRIVATE
	stdafx.cpp
	stdafx.h
	${H_FILES}
	${CXX_FILES}
)

frantic_common_platform_setup( test_frost )
frantic_default_source_groups( test_frost )
frantic_link_apple_core_libraries( test_frost )

add_precompiled_header( test_frost stdafx.h SOURCE_CXX stdafx.cpp )

# Disable optimization for the RelWithDebInfo configuration on Windows.
# This allows breakpoints to be hit reliably when debugging in Visual Studio.
if( WIN32 )
	target_compile_options( test_frost PRIVATE "$<$<CONFIG:Release>:/O2>$<$<CONFIG:RelWithDebInfo>:/Od>" )
endif()

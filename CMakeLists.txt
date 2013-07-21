cmake_minimum_required (VERSION 2.6.0)
project (MTMCSim)


function(info var)
    message("${var}: ${${var}}")
endfunction()

file (GLOB found ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
set (srcfiles ${found})

file (GLOB found ${CMAKE_CURRENT_SOURCE_DIR}/include/*.h)
set (includefiles ${found})

set(BOOST_USE_STATIC_LIBS OFF)

find_package(Boost REQUIRED COMPONENTS thread)


info(Boost_INCLUDE_DIRS)
info(Boost_FOUND)
info(Boost_LIBRARY_DIRS)
include_directories (${CMAKE_CURRENT_SOURCE_DIR}/include/ ${Boost_INCLUDE_DIRS})
link_directories( ${Boost_LIBRARY_DIRS})

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/simulations) # could be bin  

add_executable (MTMCSim ${srcfiles})

#add_library(boost_thread SHARED IMPORTED)

info(Boost_THREAD_LIBRARY)

#The ${Boost_THREAD_LIBRARY} has format optimized;filename;debug;filename, we only want file name

if ( CMAKE_SYSTEM_NAME STREQUAL "Windows")
if ( CMAKE_BUILD_TYPE STREQUAL "Debug")
    set (Boost_THREAD_LIBFILE ${Boost_THREAD_LIBRARY})
    list(REMOVE_AT Boost_THREAD_LIBFILE 0 1 2)
else()
    set (Boost_THREAD_LIBFILE ${Boost_THREAD_LIBRARY})
    list(REMOVE_AT Boost_THREAD_LIBFILE 0 2 3)        
endif()
endif()


if (CMAKE_SYSTEM_NAME STREQUAL "Windows")
    string(REPLACE .lib .dll dll_loc ${Boost_THREAD_LIBFILE})
    set(BOOST_FILES ${dll_loc})
elseif (CMAKE_SYSTEM_NAME STREQUAL "Linux" )
    set(BOOST_FILES ${Boost_THREAD_LIBFILE})
endif()

info(BOOST_FILES)
# Need to copy DLL file only if on Windows. On Linux, set LD_LIBRARY_PATH
# Note: Depending on different versions of boost, one may need more libaries, the libs that boost thread depends on
if ( CMAKE_SYSTEM_NAME STREQUAL "Windows")
add_custom_command(TARGET MTMCSim
                    POST_BUILD
                    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${BOOST_FILES} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
endif()


#target_include_directories( MTMCSim PRIVATE ${Boost_INCLUDE_DIRS} )

target_link_libraries( MTMCSim ${Boost_THREAD_LIBRARY} )





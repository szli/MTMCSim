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

include_directories (${CMAKE_CURRENT_SOURCE_DIR}/include/ ${Boost_INCLUDE_DIRS})
link_directories( ${BOOST_ROOT}/lib)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)

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
add_custom_command(TARGET MTMCSim
                    POST_BUILD
                    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${BOOST_FILES} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})


include_directories( ${BOOST_INCLUDEDIR} )

target_link_libraries( MTMCSim ${Boost_THREAD_LIBRARY} )





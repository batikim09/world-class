cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(world_vocoder)

#set(FFTW3_LIBRARY /usr/lib/x86_64-linux-gnu/libfftw3.so.3)
#set(SNDFIEL_LIBRARY /usr/lib/x86_64-linux-gnu/libsndfile.so.1)

find_library(FFTW3 fftw3)
find_library(SNDFILE sndfile)

include_directories(src include tools sptk sptk/bin/fft sptk/bin/fftr sptk/bin/ifftr sptk/bin/freqt /usr/local/include/eigen3/Eigen)

file(GLOB SPTK_HEADER_FILES
            "sptk/*.h"
    )

#We avoid generating executable files.
file(GLOB SPTK_SOURCE_FILES
            "sptk/*.c"
            "sptk/bin/fft/_*.c"
            "sptk/bin/fftr/_*.c"
            "sptk/bin/ifftr/_*.c"
            "sptk/bin/freqt/_*.c"
    )


file(GLOB TOOL_HEADER_FILES
            "tools/*.hpp"
            "tools/*.cpp"
    )
file(GLOB HEADER_FILES
            "include/*.hpp"

    )
file(GLOB SOURCE_FILES
            "src/*.cpp"
    )


add_library(world SHARED ${SOURCE_FILES} ${HEADER_FILES} ${TOOL_HEADER_FILES} ${SPTK_SOURCE_FILES} ${SPTK_HEADER_FILES})

add_executable(test test/test.cpp ${HEADER_FILES} ${TOOL_HEADER_FILES})
target_link_libraries(test ${FFTW3} ${SNDFILE} world)
set_property(TARGET test PROPERTY CXX_STANDARD 11)

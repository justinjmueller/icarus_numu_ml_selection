cmake_minimum_required(VERSION 3.12)
project(run_systematics)

#set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
#set(THREADS_PREFER_PTHREAD_FLAG TRUE)
#find_package(Threads REQUIRED)

add_compile_options(-Wall -Werror -Wno-error=maybe-uninitialized -Ofast)

# Find the ROOT package
find_package(ROOT REQUIRED)

find_package(sbnanaobj)
set(SBNANAOBJ_INCLUDE_DIRS "$ENV{MRB_SOURCE}/sbnanaobj/")

file(GLOB SYSINC "include/*.h")

# Add the executable target
add_executable(run_systematics src/main.cc ${SYSINC})

# Link the ROOT libraries to the target
target_link_libraries(run_systematics ${ROOT_LIBRARIES} ${sbnanaobj_LIBRARY_DIRS}/libsbnanaobj_StandardRecord.so)# Threads::Threads)

# Include the ROOT headers
include_directories(${ROOT_INCLUDE_DIRS} ${SBNANAOBJ_INCLUDE_DIRS} include/)

# Add ROOT definitions
add_definitions(${ROOT_CXX_FLAGS})
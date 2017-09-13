#include experimental::filesystem
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lstdc++fs")

seqan3_test(matching_statistics_test.cpp)

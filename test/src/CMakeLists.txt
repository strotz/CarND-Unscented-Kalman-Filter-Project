set(test_sources
        main.cpp
        assignments.cpp)

add_executable(testUKF ${test_sources})

include_directories(../../src)

target_link_libraries(testUKF
        libgtest
        libgmock
        )

add_test(NAME testUKF
        COMMAND testUKF)
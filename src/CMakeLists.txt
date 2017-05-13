set(sources
    ukf.cpp
    main.cpp
    tools.cpp)

add_executable(UnscentedKF ${sources})

install(TARGETS UnscentedKF DESTINATION bin)

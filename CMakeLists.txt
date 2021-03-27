cmake_minimum_required(VERSION 3.17)
project(DeformableCellModel)

set(CMAKE_CXX_STANDARD 17)

include_directories(/usr/local/include/)
link_directories(/usr/local/lib/)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -flto")

add_executable(DeformableCellModel main.cpp
        Engine.cpp Engine.hpp
        CellMesh.cpp CellMesh.hpp
        Definitions.hpp
        Observer.cpp Observer.hpp
        EquationsOfMotion.cpp EquationsOfMotion.hpp
        Parameters.cpp Parameters.hpp
        StlObserver.cpp StlObserver.hpp
        Stepper.cpp Stepper.hpp)

//
// Created by Nikita Kruk on 15.03.21.
//

#include "Engine.hpp"
#include "EquationsOfMotion.hpp"
#include "Observer.hpp"
#include "StlObserver.hpp"

#include <iostream>
#include <string>

Engine::Engine() :
    parameters_()
{
  std::cout << "Start of the simulation" << std::endl;

  std::string
      cell_mesh_file_name("/Users/nikita/CLionProjects/cgal_sphere_mesh_generation/cmake-build-debug/sphere.off");
  cell_mesh_ = CellMesh(cell_mesh_file_name, parameters_);
}

Engine::~Engine()
{
  std::cout << "End of the simulation" << std::endl;
}

void Engine::RunSimulation()
{
  EquationsOfMotion equations_of_motion(parameters_);
  Observer observer;

  observer.SaveNodes(cell_mesh_);
  observer.SaveNormalsForNodes(cell_mesh_);
  Observer::SaveFaces(cell_mesh_);

  const int T = 1000;
  for (int n = 0; n < T; ++n)
  {
    std::cout << "n = " << n << std::endl;
    equations_of_motion.DoStep(cell_mesh_);
    observer.SaveNodes(cell_mesh_);
    observer.SaveNormalsForNodes(cell_mesh_);
//    cell_mesh_.SetInitialVolume(cell_mesh_.GetInitialVolume() + 1e-15 * parameters_.GetDt());
  } // n

  StlObserver stl_observer;
  stl_observer.SaveCellMesh(cell_mesh_);
}
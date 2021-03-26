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
    parameters_(),
    cell_meshes_()
{
  std::cout << "Start of the simulation" << std::endl;

  std::string
      cell_mesh_file_name("/Users/nikita/CLionProjects/cgal_sphere_mesh_generation/cmake-build-debug/sphere.off");
//  cell_mesh_ = CellMesh(cell_mesh_file_name, parameters_);
  cell_meshes_.push_back(CellMesh(cell_mesh_file_name, parameters_));
}

Engine::~Engine()
{
  std::cout << "End of the simulation" << std::endl;
}

void Engine::RunSimulation()
{
  EquationsOfMotion equations_of_motion(parameters_);
  Observer observer;

  observer.SaveNodes(cell_meshes_);
  observer.SaveNormalsForNodes(cell_meshes_);
  Observer::SaveFaces(cell_meshes_);

  const int T = 1000;
  for (int n = 0; n < T; ++n)
  {
    std::cout << "n = " << n << std::endl;
    equations_of_motion.DoStep(cell_meshes_);
    observer.SaveNodes(cell_meshes_);
    observer.SaveNormalsForNodes(cell_meshes_);
//    cell_mesh_.SetInitialVolume(cell_mesh_.GetInitialVolume() + 1e-15 * parameters_.GetDt());
  } // n

  StlObserver stl_observer;
  stl_observer.SaveCellMeshes(cell_meshes_);
}
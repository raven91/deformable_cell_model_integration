//
// Created by Nikita Kruk on 15.03.21.
//

#ifndef DEFORMABLECELLMODEL_ENGINE_HPP
#define DEFORMABLECELLMODEL_ENGINE_HPP

#include "CellMesh.hpp"
#include "Parameters.hpp"

#include <vector>

class Engine
{
 public:

  Engine();
  ~Engine();

  void RunSimulation();

 private:

  Parameters parameters_;
//  CellMesh cell_mesh_;
  std::vector<CellMesh> cell_meshes_;

};

#endif //DEFORMABLECELLMODEL_ENGINE_HPP

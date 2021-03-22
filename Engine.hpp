//
// Created by Nikita Kruk on 15.03.21.
//

#ifndef DEFORMABLECELLMODEL_ENGINE_HPP
#define DEFORMABLECELLMODEL_ENGINE_HPP

#include "CellMesh.hpp"
#include "Parameters.hpp"

class Engine
{
 public:

  Engine();
  ~Engine();

  void RunSimulation();

 private:

  CellMesh cell_mesh_;
  Parameters parameters_;

};

#endif //DEFORMABLECELLMODEL_ENGINE_HPP

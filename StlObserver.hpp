//
// Created by Nikita Kruk on 18.03.21.
//

#ifndef DEFORMABLECELLMODEL_STLOBSERVER_HPP
#define DEFORMABLECELLMODEL_STLOBSERVER_HPP

#include "CellMesh.hpp"

#include <vector>

/*
 * Observer that saves a cell mesh in a binary STL format
 * (https://en.wikipedia.org/wiki/STL_(file_format))
 */
class StlObserver
{
 public:

  StlObserver();
  ~StlObserver();

  void SaveCellMeshes(const std::vector<CellMesh> &cell_meshes) const;

 private:

  void SaveCellMesh(const CellMesh &cell_mesh) const;

};

#endif //DEFORMABLECELLMODEL_STLOBSERVER_HPP

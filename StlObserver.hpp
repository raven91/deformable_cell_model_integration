//
// Created by Nikita Kruk on 18.03.21.
//

#ifndef DEFORMABLECELLMODEL_STLOBSERVER_HPP
#define DEFORMABLECELLMODEL_STLOBSERVER_HPP

#include "CellMesh.hpp"

class StlObserver
{
 public:

  StlObserver();
  ~StlObserver();

  void SaveCellMesh(const CellMesh &cell_mesh) const;

};

#endif //DEFORMABLECELLMODEL_STLOBSERVER_HPP

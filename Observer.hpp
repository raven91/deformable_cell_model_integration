//
// Created by Nikita Kruk on 15.03.21.
//

#ifndef DEFORMABLECELLMODEL_OBSERVER_HPP
#define DEFORMABLECELLMODEL_OBSERVER_HPP

#include "Definitions.hpp"
#include "CellMesh.hpp"

#include <string>
#include <fstream>

class Observer
{
 public:

  Observer();
  ~Observer();

  void SaveNodes(const CellMesh &cell_mesh);
  void SaveNormalsForNodes(const CellMesh &cell_mesh);
  static void SaveFaces(const CellMesh &cell_mesh);

 private:

  std::string node_file_name_;
  std::ofstream node_file_;
  std::string normal_file_name_;
  std::ofstream normal_file_;

};

#endif //DEFORMABLECELLMODEL_OBSERVER_HPP

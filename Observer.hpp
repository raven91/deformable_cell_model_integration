//
// Created by Nikita Kruk on 15.03.21.
//

#ifndef DEFORMABLECELLMODEL_OBSERVER_HPP
#define DEFORMABLECELLMODEL_OBSERVER_HPP

#include "Definitions.hpp"
#include "CellMesh.hpp"

#include <string>
#include <fstream>
#include <vector>

class Observer
{
 public:

  Observer();
  ~Observer();

  void SaveNodes(const std::vector<CellMesh> &cell_meshes);
  void SaveNormalsForNodes(const std::vector<CellMesh> &cell_meshes);
  static void SaveFaces(const std::vector<CellMesh> &cell_meshes);

 private:

  std::string node_file_name_;
  std::ofstream node_file_;
  std::string normal_file_name_;
  std::ofstream normal_file_;

};

#endif //DEFORMABLECELLMODEL_OBSERVER_HPP

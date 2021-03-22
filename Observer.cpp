//
// Created by Nikita Kruk on 15.03.21.
//

#include "Observer.hpp"

#include <iostream>
#include <sstream>

Observer::Observer()
{
  node_file_name_ = "/Users/nikita/Documents/Projects/DeformableCellModel/nodes.txt";
  node_file_.open(node_file_name_, std::ios::out | std::ios::trunc);
  if (!node_file_.is_open())
  {
    std::cerr << "Cannot write to a file" << std::endl;
    exit(-1);
  }

  normal_file_name_ = "/Users/nikita/Documents/Projects/DeformableCellModel/normals.txt";
  normal_file_.open(normal_file_name_, std::ios::out | std::ios::trunc);
  if (!normal_file_.is_open())
  {
    std::cerr << "Cannot write to a file" << std::endl;
    exit(-1);
  }
}

Observer::~Observer()
{
  if (node_file_.is_open())
  {
    node_file_.close();
  }
  if (normal_file_.is_open())
  {
    normal_file_.close();
  }
}

void Observer::SaveNodes(const CellMesh &cell_mesh)
{
  std::ostringstream output_buffer;
  for (const std::array<double, kDim> &node : cell_mesh.GetNodes())
  {
    output_buffer << node[0] << '\t' << node[1] << '\t' << node[2] << '\t';
  } // node
  output_buffer << std::endl;

  node_file_ << output_buffer.str();
}

void Observer::SaveNormalsForNodes(const CellMesh &cell_mesh)
{
  std::ostringstream output_buffer;
  for (const std::array<double, kDim> &normal : cell_mesh.GetNormalsForNodes())
  {
    output_buffer << normal[0] << '\t' << normal[1] << '\t' << normal[2] << '\t';
  } // normal
  output_buffer << std::endl;

  normal_file_ << output_buffer.str();
}

void Observer::SaveFaces(const CellMesh &cell_mesh)
{
  std::ostringstream output_buffer;
  for (const std::array<int, kFaceDim> &face : cell_mesh.GetFaces())
  {
    for (int node_index : face)
    {
      output_buffer << node_index << '\t';
    } // node_index
  } // face
  output_buffer << std::endl;

  std::ofstream
      faces_file("/Users/nikita/Documents/Projects/DeformableCellModel/faces.txt", std::ios::out | std::ios::trunc);
  if (!faces_file.is_open())
  {
    std::cerr << "Cannot write to a file" << std::endl;
    exit(-1);
  } else
  {
    faces_file << output_buffer.str();
  }
  faces_file.close();
}
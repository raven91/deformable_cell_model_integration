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

void Observer::SaveNodes(const std::vector<CellMesh> &cell_meshes)
{
  for (int cell_idx = 0; cell_idx < cell_meshes.size(); ++cell_idx)
  {
    std::ostringstream output_buffer;

    for (const VectorType &node : cell_meshes[cell_idx].GetNodes())
    {
      output_buffer << node[0] << '\t' << node[1] << '\t' << node[2] << '\t';
    } // node
    output_buffer << std::endl;

    node_file_ << output_buffer.str();
  } // cell_idx
}

void Observer::SaveNormalsForNodes(const std::vector<CellMesh> &cell_meshes)
{
  for (int cell_idx = 0; cell_idx < cell_meshes.size(); ++cell_idx)
  {
    std::ostringstream output_buffer;
    for (const VectorType &normal : cell_meshes[cell_idx].GetNormalsForNodes())
    {
      output_buffer << normal[0] << '\t' << normal[1] << '\t' << normal[2] << '\t';
    } // normal
    output_buffer << std::endl;

    normal_file_ << output_buffer.str();
  } // cell_idx
}

void Observer::SaveFaces(const std::vector<CellMesh> &cell_meshes)
{
  for (int cell_idx = 0; cell_idx < cell_meshes.size(); ++cell_idx)
  {
    std::ostringstream output_buffer;
    for (const FaceType &face : cell_meshes[cell_idx].GetFaces())
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
  } // cell_idx
}
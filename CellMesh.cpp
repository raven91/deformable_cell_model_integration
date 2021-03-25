//
// Created by Nikita Kruk on 15.03.21.
//

#include "CellMesh.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <numeric> // std::accumulate
#include <eigen3/Eigen/Dense>

CellMesh::CellMesh() :
    nodes_(),
    faces_(),
    edges_(),
    adjacent_faces_for_nodes_(),
    surface_areas_for_nodes_(),
    normals_for_faces_(),
    normals_for_nodes_(),
    initial_cell_volume_(0.0)
{

}

CellMesh::CellMesh(int n_nodes, int n_faces, int n_edges) :
    nodes_(n_nodes, VectorType{0.0}),
    faces_(n_faces, FaceType{-1}),
    edges_(n_edges, EdgeType{-1}),
    adjacent_faces_for_nodes_(),
    surface_areas_for_nodes_(),
    normals_for_faces_(),
    normals_for_nodes_(),
    initial_cell_volume_(0.0)
{

}

CellMesh::CellMesh(const std::string &off_fine_name, const Parameters &parameters) :
    surface_areas_for_nodes_(),
    normals_for_faces_(),
    normals_for_nodes_(),
    initial_cell_volume_(0.0)
{
  std::string file_name("/Users/nikita/CLionProjects/cgal_sphere_mesh_generation/cmake-build-debug/sphere.off");
  std::ifstream file(file_name, std::ios::in);
  if (file.is_open())
  {
    std::string header;
    std::getline(file, header);

    std::string parameters_line;
    std::getline(file, parameters_line);
    std::istringstream parameters_line_buffer(parameters_line);
    int n_nodes = 0, n_faces = 0, n_edges = 0;
    parameters_line_buffer >> n_nodes >> n_faces >> n_edges;

    double x, y, z;
    double norm, scaling;
    for (int i = 0; i < n_nodes; ++i)
    {
      file >> x >> y >> z;
      norm = std::hypot(x, y, z);
      scaling = parameters.GetRadius() / norm;
      nodes_.push_back({x * scaling, y * scaling, z * scaling});
    } // i

    int n_nodes_per_face, n_0, n_1, n_2;
    for (int i = 0; i < n_faces; ++i)
    {
      file >> n_nodes_per_face >> n_0 >> n_1 >> n_2;
      faces_.push_back({n_0, n_1, n_2});
    } // i

    n_edges = 0;
    // TODO: construct a list of edges

    adjacent_faces_for_nodes_.resize(n_nodes);
    for (int f = 0; f < faces_.size(); ++f)
    {
      for (int vertex : faces_[f])
      {
        adjacent_faces_for_nodes_[vertex].insert(f);
      } // vertex
    } // f

    CalculateFaceSurfaceAreas();
    CalculateNodeSurfaceAreas();
    initial_cell_surface_area_ = CalculateCellSurfaceArea();
    initial_cell_volume_ = CalculateCellVolume();
    MakeFacesOriented();
    CalculateNodeNormals();

    file.close();
  } else
  {
    std::cerr << "Cannot open off-file" << std::endl;
    exit(-1);
  }
}

CellMesh::~CellMesh()
{
  nodes_.clear();
  faces_.clear();
  edges_.clear();
}

const std::vector<VectorType> &CellMesh::GetNodes() const
{
  return nodes_;
}

std::vector<VectorType> &CellMesh::GetNodes()
{
  return nodes_;
}

const std::vector<FaceType> &CellMesh::GetFaces() const
{
  return faces_;
}

const std::vector<std::set<int>> &CellMesh::GetAdjacentFacesForNodes() const
{
  return adjacent_faces_for_nodes_;
}

const std::vector<VectorType> &CellMesh::GetNormalsForNodes() const
{
  return normals_for_nodes_;
}

const std::vector<VectorType> &CellMesh::GetNormalsForFaces() const
{
  return normals_for_faces_;
}

int CellMesh::GetNumNodes() const
{
  return nodes_.size();
}

int CellMesh::GetNumFaces() const
{
  return faces_.size();
}

void CellMesh::CalculateFaceSurfaceAreas() const
{
  surface_areas_for_faces_ = std::vector<double>(faces_.size(), 0.0);
  for (int f = 0; f < faces_.size(); ++f)
  {
    surface_areas_for_faces_[f] = FaceArea(f);
  } // face
}

/*
 * Requires surface areas of each face
 * Calculates the surface area of a node as the average of surface areas of adjacent faces
 */
const std::vector<double> &CellMesh::CalculateNodeSurfaceAreas() const
{
  surface_areas_for_nodes_ = std::vector<double>(nodes_.size(), 0.0);
  for (int i = 0; i < nodes_.size(); ++i)
  {
    double total_area = 0.0;
    for (int face_index : adjacent_faces_for_nodes_[i])
    {
      total_area += surface_areas_for_faces_[face_index];
    } // face
    surface_areas_for_nodes_[i] = total_area / adjacent_faces_for_nodes_[i].size();
  } // i
  return surface_areas_for_nodes_;
}

double CellMesh::FaceArea(int face_index) const
{
  int n_0 = faces_[face_index][0], n_1 = faces_[face_index][1], n_2 = faces_[face_index][2];
  const VectorType &p_0 = nodes_[n_0], &p_1 = nodes_[n_1], &p_2 = nodes_[n_2];
  double a = std::hypot(p_0[0] - p_1[0], p_0[1] - p_1[1], p_0[2] - p_1[2]);
  double b = std::hypot(p_1[0] - p_2[0], p_1[1] - p_2[1], p_1[2] - p_2[2]);
  double c = std::hypot(p_2[0] - p_0[0], p_2[1] - p_0[1], p_2[2] - p_0[2]);
  double s = (a + b + c) / 2.0;
  double area = std::sqrt(s * (s - a) * (s - b) * (s - c));
  return area;
}

double CellMesh::GetInitialSurfaceArea() const
{
  return initial_cell_surface_area_;
}

double CellMesh::GetInitialVolume() const
{
  return initial_cell_volume_;
}

void CellMesh::SetInitialVolume(double new_volume)
{
  initial_cell_volume_ = new_volume;
}

/*
 * Requires surface areas of each face
 */
double CellMesh::CalculateCellSurfaceArea() const
{
  return std::accumulate(surface_areas_for_faces_.begin(), surface_areas_for_faces_.end(), 0.0);
}

double CellMesh::CalculateCellVolume() const
{
  std::array<double, kDim> center_of_mass{0.0, 0.0, 0.0};
  for (const VectorType &node : nodes_)
  {
    center_of_mass[0] += node[0];
    center_of_mass[1] += node[1];
    center_of_mass[2] += node[2];
  } // node
  center_of_mass[0] /= nodes_.size();
  center_of_mass[1] /= nodes_.size();
  center_of_mass[2] /= nodes_.size();

  double volume = 0.0;
  Eigen::Vector3d p_0{center_of_mass[0], center_of_mass[1], center_of_mass[2]}, p_1, p_2, p_3;
  int n_1, n_2, n_3;
  for (const FaceType &face : faces_)
  {
    n_1 = face[0];
    n_2 = face[1];
    n_3 = face[2];
    p_1 = Eigen::Vector3d(nodes_[n_1][0], nodes_[n_1][1], nodes_[n_1][2]);
    p_2 = Eigen::Vector3d(nodes_[n_2][0], nodes_[n_2][1], nodes_[n_2][2]);
    p_3 = Eigen::Vector3d(nodes_[n_3][0], nodes_[n_3][1], nodes_[n_3][2]);
    volume += (p_1 - p_0).dot((p_2 - p_0).cross(p_3 - p_0)) / 6.0;
  } // face
  return volume;
}

void CellMesh::MakeFacesOriented()
{
  VectorType center_of_mass{0.0, 0.0, 0.0};
  for (const VectorType &node : nodes_)
  {
    center_of_mass[0] += node[0];
    center_of_mass[1] += node[1];
    center_of_mass[2] += node[2];
  } // node
  center_of_mass[0] /= nodes_.size();
  center_of_mass[1] /= nodes_.size();
  center_of_mass[2] /= nodes_.size();

  Eigen::Vector3d p_0{center_of_mass[0], center_of_mass[1], center_of_mass[2]}, p_1, p_2, p_3, p_12, p_23;
  int n_1, n_2, n_3;
  for (FaceType &face : faces_)
  {
    n_1 = face[0];
    n_2 = face[1];
    n_3 = face[2];
    p_1 = Eigen::Vector3d(nodes_[n_1][0], nodes_[n_1][1], nodes_[n_1][2]);
    p_2 = Eigen::Vector3d(nodes_[n_2][0], nodes_[n_2][1], nodes_[n_2][2]);
    p_3 = Eigen::Vector3d(nodes_[n_3][0], nodes_[n_3][1], nodes_[n_3][2]);
    p_12 = p_2 - p_1;
    p_23 = p_3 - p_2;
    if ((p_1 - p_0).dot(p_12.cross(p_23)) < 0.0)
    {
      face = FaceType{n_1, n_3, n_2};
    }
  } // face
}

/*
 * Requires faces to be CCW
 */
const std::vector<VectorType> &CellMesh::CalculateFaceNormals() const
{
  normals_for_faces_ = std::vector<VectorType>(faces_.size(), {0.0});
  Eigen::Vector3d p_1, p_2, p_3, p_12, p_23;
  Eigen::Vector3d normal;
  int n_1, n_2, n_3;
  for (int f = 0; f < faces_.size(); ++f)
  {
    n_1 = faces_[f][0];
    n_2 = faces_[f][1];
    n_3 = faces_[f][2];
    p_1 = Eigen::Vector3d(nodes_[n_1][0], nodes_[n_1][1], nodes_[n_1][2]);
    p_2 = Eigen::Vector3d(nodes_[n_2][0], nodes_[n_2][1], nodes_[n_2][2]);
    p_3 = Eigen::Vector3d(nodes_[n_3][0], nodes_[n_3][1], nodes_[n_3][2]);
    normal = (p_2 - p_1).cross(p_3 - p_2).normalized();
    normals_for_faces_[f][0] = normal[0];
    normals_for_faces_[f][1] = normal[1];
    normals_for_faces_[f][2] = normal[2];
  } // f
  return normals_for_faces_;
}

/*
 * Requires face normals
 */
const std::vector<VectorType> &CellMesh::CalculateNodeNormals() const
{
  CalculateFaceNormals();
  normals_for_nodes_ = std::vector<VectorType>(nodes_.size(), {0.0});
  for (int i = 0; i < nodes_.size(); ++i)
  {
    VectorType node_normal{0.0};
    for (int face_index : adjacent_faces_for_nodes_[i])
    {
      node_normal[0] += normals_for_faces_[face_index][0];
      node_normal[1] += normals_for_faces_[face_index][1];
      node_normal[2] += normals_for_faces_[face_index][2];
    } // face_index
    double norm = std::hypot(node_normal[0], node_normal[1], node_normal[2]);
    node_normal[0] /= norm;
    node_normal[1] /= norm;
    node_normal[2] /= norm;
    normals_for_nodes_[i] = node_normal;
  } // i
  return normals_for_nodes_;
}
//
// Created by Nikita Kruk on 15.03.21.
//

#include "CellMesh.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <numeric> // std::accumulate
#include <algorithm> // std::min, std::max, std::find
#include <iterator> // std::distance

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

CellMesh::CellMesh(const std::string &off_fine_name, const Parameters &parameters) :
    surface_areas_for_nodes_(),
    normals_for_faces_(),
    normals_for_nodes_(),
    initial_cell_surface_area_(0.0),
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
      nodes_.emplace_back(x * scaling, y * scaling, z * scaling);
    } // i

    int n_nodes_per_face, n_0, n_1, n_2;
    for (int i = 0; i < n_faces; ++i)
    {
      file >> n_nodes_per_face >> n_0 >> n_1 >> n_2;
      faces_.push_back({n_0, n_1, n_2});
    } // i

    // construct the list of edges
    // and simultaneously their adjacent faces
    EdgeType e_1, e_2, e_3;
    std::vector<EdgeType>::iterator it;
    int edge_idx;
    for (int f = 0; f < faces_.size(); ++f)
    {
      n_0 = faces_[f][0];
      n_1 = faces_[f][1];
      n_2 = faces_[f][2];

      e_1 = std::make_pair(std::min(n_0, n_1), std::max(n_0, n_1));
      e_2 = std::make_pair(std::min(n_1, n_2), std::max(n_1, n_2));
      e_3 = std::make_pair(std::min(n_0, n_2), std::max(n_0, n_2));
      const std::set<EdgeType> edges_per_face{e_1, e_2, e_3};

      for (const EdgeType &edge : edges_per_face)
      {
        it = std::find(edges_.begin(), edges_.end(), edge);
        if (it != edges_.end()) // edge already exists
        {
          edge_idx = std::distance(edges_.begin(), it);
          adjacent_faces_for_edges_[edge_idx][1] = f;
        } else // edge is new
        {
          edges_.emplace_back(edge);
          adjacent_faces_for_edges_.emplace_back(std::array<int, 2>{f, -1});
        }
      } // edge
    } // f
    n_edges = edges_.size();

    VectorType initial_separation;
    double initial_distance = 0.0;
    adjacent_nodes_for_nodes_.resize(n_nodes);
    initial_lenghts_between_adjacent_nodes_.resize(n_nodes);
    for (const EdgeType &edge : edges_)
    {
      adjacent_nodes_for_nodes_[edge.first].push_back(edge.second);
      adjacent_nodes_for_nodes_[edge.second].push_back(edge.first);

      initial_separation = nodes_[edge.first] - nodes_[edge.second];
      initial_distance = initial_separation.norm();
      initial_lenghts_between_adjacent_nodes_[edge.first].push_back(initial_distance);
      initial_lenghts_between_adjacent_nodes_[edge.second].push_back(initial_distance);
    } // edge

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
    CalculateInitialCurvatureAngleForEdges();

    file.close();
  } else
  {
    std::cerr << "Cannot open off-file" << std::endl;
    exit(-1);
  }
}

CellMesh::CellMesh(const std::string &off_plane_file_name) :
    initial_cell_surface_area_(0.0),
    initial_cell_volume_(0.0)
{
  std::string file_name("/Users/nikita/CLionProjects/cgal_sphere_mesh_generation/cmake-build-debug/plane.off");
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
    const double norm = 1.0, scaling = 10e-6; // todo: remove hardcored constants
    const VectorType downward_shift(0.0, 0.0, -2.0 * scaling);
    for (int i = 0; i < n_nodes; ++i)
    {
      file >> x >> y >> z;
      nodes_.emplace_back(VectorType(x * scaling, y * scaling, z * scaling) + downward_shift);
    } // i

    int n_nodes_per_face, n_0, n_1, n_2;
    for (int i = 0; i < n_faces; ++i)
    {
      file >> n_nodes_per_face >> n_0 >> n_1 >> n_2;
      faces_.push_back({n_0, n_1, n_2});
    } // i

    adjacent_faces_for_nodes_.resize(n_nodes);
    for (int f = 0; f < faces_.size(); ++f)
    {
      for (int vertex : faces_[f])
      {
        adjacent_faces_for_nodes_[vertex].insert(f);
      } // vertex
    } // f

    MakeFacesOriented();
    CalculateNodeNormals();
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

const std::vector<EdgeType> &CellMesh::GetEdges() const
{
  return edges_;
}

const std::vector<IndexSet> &CellMesh::GetAdjacentFacesForNodes() const
{
  return adjacent_faces_for_nodes_;
}

const std::vector<std::vector<int>> &CellMesh::GetAdjacentNodesForNodes() const
{
  return adjacent_nodes_for_nodes_;
}

const std::vector<std::array<int, 2>> &CellMesh::GetAdjacentFacesForEdges() const
{
  return adjacent_faces_for_edges_;
}

const std::vector<VectorType> &CellMesh::GetNormalsForNodes() const
{
  return normals_for_nodes_;
}

const std::vector<VectorType> &CellMesh::GetNormalsForFaces() const
{
  return normals_for_faces_;
}

const std::vector<std::vector<double>> &CellMesh::GetInitialLengthsBetweenAdjacentNodes() const
{
  return initial_lenghts_between_adjacent_nodes_;
}

const std::vector<double> &CellMesh::GetInitialCurvatureAngleForEdges() const
{
  return initial_curvature_angle_for_edges_;
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
    surface_areas_for_nodes_[i] = total_area / 3.0;
  } // i
  return surface_areas_for_nodes_;
}

/*
 * Calculate face area using Heron's formula
 */
double CellMesh::FaceArea(int face_index) const
{
  int n_0 = faces_[face_index][0], n_1 = faces_[face_index][1], n_2 = faces_[face_index][2];
  const VectorType &p_0 = nodes_[n_0], &p_1 = nodes_[n_1], &p_2 = nodes_[n_2];
  double a = (p_0 - p_1).norm();
  double b = (p_1 - p_2).norm();
  double c = (p_2 - p_0).norm();
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
  VectorType center_of_mass{0.0, 0.0, 0.0};
  for (const VectorType &node : nodes_)
  {
    center_of_mass += node;
  } // node
  center_of_mass /= nodes_.size();

  double volume = 0.0;
  VectorType p_0(center_of_mass), p_1, p_2, p_3;
  int n_1, n_2, n_3;
  for (const FaceType &face : faces_)
  {
    n_1 = face[0];
    n_2 = face[1];
    n_3 = face[2];
    p_1 = nodes_[n_1];
    p_2 = nodes_[n_2];
    p_3 = nodes_[n_3];
    volume += (p_1 - p_0).dot((p_2 - p_0).cross(p_3 - p_0)) / 6.0;
  } // face
  return volume;
}

void CellMesh::MakeFacesOriented()
{
  VectorType center_of_mass = VectorType::Zero();
  for (const VectorType &node : nodes_)
  {
    center_of_mass += node;
  } // node
  center_of_mass /= nodes_.size();

  Eigen::Vector3d p_0(center_of_mass), p_1, p_2, p_3, p_12, p_23;
  int n_1, n_2, n_3;
  for (FaceType &face : faces_)
  {
    n_1 = face[0];
    n_2 = face[1];
    n_3 = face[2];
    p_1 = nodes_[n_1];
    p_2 = nodes_[n_2];
    p_3 = nodes_[n_3];
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
  normals_for_faces_ = std::vector<VectorType>(faces_.size(), VectorType::Zero());
  VectorType p_1, p_2, p_3, p_12, p_23;
  VectorType normal;
  int n_1, n_2, n_3;
  for (int f = 0; f < faces_.size(); ++f)
  {
    n_1 = faces_[f][0];
    n_2 = faces_[f][1];
    n_3 = faces_[f][2];
    p_1 = nodes_[n_1];
    p_2 = nodes_[n_2];
    p_3 = nodes_[n_3];
    normal = (p_2 - p_1).cross(p_3 - p_2).normalized();
    normals_for_faces_[f] = normal;
  } // f
  return normals_for_faces_;
}

/*
 * Requires face normals
 */
const std::vector<VectorType> &CellMesh::CalculateNodeNormals() const
{
  CalculateFaceNormals();
  normals_for_nodes_ = std::vector<VectorType>(nodes_.size(), VectorType::Zero());
  for (int i = 0; i < nodes_.size(); ++i)
  {
    VectorType node_normal = VectorType::Zero();
    for (int face_index : adjacent_faces_for_nodes_[i])
    {
      node_normal += normals_for_faces_[face_index];
    } // face_index
    // todo: rewrite using Eigen's functionality
    normals_for_nodes_[i] = node_normal.normalized();
  } // i
  return normals_for_nodes_;
}

/*
 * Requires face normals
 */
void CellMesh::CalculateInitialCurvatureAngleForEdges()
{
  CalculateFaceNormals();
  initial_curvature_angle_for_edges_ = std::vector<double>(edges_.size(), 0.0);
  int i, j, k, l;
  int face_1, face_2;
  for (int edge_idx = 0; edge_idx < edges_.size(); ++edge_idx)
  {
    k = edges_[edge_idx].first;
    j = edges_[edge_idx].second;

    face_1 = adjacent_faces_for_edges_[edge_idx][0];
    face_2 = adjacent_faces_for_edges_[edge_idx][1];

    const std::set<int> edge_nodes{k, j};
    std::set_difference(faces_[face_1].begin(), faces_[face_1].end(), edge_nodes.begin(), edge_nodes.end(), &i);
    std::set_difference(faces_[face_2].begin(), faces_[face_2].end(), edge_nodes.begin(), edge_nodes.end(), &l);

    VectorType p_kj = nodes_[j] - nodes_[k], p_ki = nodes_[i] - nodes_[k], p_kl = nodes_[l] - nodes_[k];
    VectorType proj_i = p_ki - p_ki.dot(p_kj) / p_kj.squaredNorm() * p_kj;
    VectorType proj_l = p_kl - p_kl.dot(p_kj) / p_kj.squaredNorm() * p_kj;

    const VectorType &n_1 = normals_for_faces_[face_1], &n_2 = normals_for_faces_[face_2];
    double theta_0 = std::acos(n_1.dot(n_2));
    initial_curvature_angle_for_edges_[edge_idx] = theta_0;
  } // edge_idx
}
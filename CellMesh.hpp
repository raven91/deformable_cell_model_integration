//
// Created by Nikita Kruk on 15.03.21.
//

#ifndef DEFORMABLECELLMODEL_CELLMESH_HPP
#define DEFORMABLECELLMODEL_CELLMESH_HPP

#include "Definitions.hpp"
#include "Parameters.hpp"

#include <vector>
#include <array>
#include <string>
#include <set>

class CellMesh
{
 public:

  CellMesh();
  CellMesh(int n_nodes, int n_faces, int n_edges);
  CellMesh(const std::string &off_fine_name, const Parameters &parameters);
  ~CellMesh();

  const std::vector<std::array<double, kDim>> &GetNodes() const;
  std::vector<std::array<double, kDim>> &GetNodes();
  const std::vector<std::array<int, kFaceDim>> &GetFaces() const;
  const std::vector<std::set<int>> &GetAdjacentFacesForNodes() const;
  const std::vector<std::array<double, kDim>> &GetNormalsForNodes() const;
  const std::vector<std::array<double, kDim>> &GetNormalsForFaces() const;

  int GetNumNodes() const;
  int GetNumFaces() const;
  void CalculateFaceSurfaceAreas() const;
  const std::vector<double> &CalculateNodeSurfaceAreas() const;
  double GetInitialSurfaceArea() const;
  double GetInitialVolume() const;
  void SetInitialVolume(double new_volume);
  double CalculateCellSurfaceArea() const;
  double CalculateCellVolume() const;
  const std::vector<std::array<double, kDim>> &CalculateFaceNormals() const;
  const std::vector<std::array<double, kDim>> &CalculateNodeNormals() const;

 private:

  std::vector<std::array<double, kDim>> nodes_;
  std::vector<std::array<int, kFaceDim>> faces_;
  std::vector<std::array<int, kEdgeDim>> edges_;
  std::vector<std::set<int>> adjacent_faces_for_nodes_;

  mutable std::vector<double> surface_areas_for_faces_;
  mutable std::vector<double> surface_areas_for_nodes_;
  double initial_cell_surface_area_;
  double initial_cell_volume_;
  mutable std::vector<std::array<double, kDim>> normals_for_faces_;
  mutable std::vector<std::array<double, kDim>> normals_for_nodes_;

  double FaceArea(int face_index) const;
  void MakeFacesOriented();

};

#endif //DEFORMABLECELLMODEL_CELLMESH_HPP

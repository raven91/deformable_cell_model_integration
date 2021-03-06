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
  CellMesh(const std::string &off_fine_name, const Parameters &parameters);
  CellMesh(const std::string &off_plane_file_name); // todo: refactor with the hierarchy of mesh classes
  ~CellMesh();

  const std::vector<VectorType> &GetNodes() const;
  std::vector<VectorType> &GetNodes();
  const std::vector<FaceType> &GetFaces() const;
  const std::vector<EdgeType> &GetEdges() const;
  const std::vector<IndexSet> &GetAdjacentFacesForNodes() const;
  const std::vector<std::vector<int>> &GetAdjacentNodesForNodes() const;
  const std::vector<std::array<int, 2>> &GetAdjacentFacesForEdges() const;
  const std::vector<VectorType> &GetNormalsForNodes() const;
  const std::vector<VectorType> &GetNormalsForFaces() const;
  const std::vector<std::vector<double>> &GetInitialLengthsBetweenAdjacentNodes() const;
  const std::vector<double> &GetInitialCurvatureAngleForEdges() const;

  int GetNumNodes() const;
  int GetNumFaces() const;
  void CalculateFaceSurfaceAreas() const;
  const std::vector<double> &CalculateNodeSurfaceAreas() const;
  double GetInitialSurfaceArea() const;
  double GetInitialVolume() const;
  void SetInitialVolume(double new_volume);
  double CalculateCellSurfaceArea() const;
  double CalculateCellVolume() const;
  const std::vector<VectorType> &CalculateFaceNormals() const;
  const std::vector<VectorType> &CalculateNodeNormals() const;
  const std::vector<VectorType> &CalculateNodeCurvatures() const;

 private:

  std::vector<VectorType> nodes_;
  std::vector<FaceType> faces_;
  std::vector<EdgeType> edges_;
  std::vector<IndexSet> adjacent_faces_for_nodes_;
  std::vector<std::vector<int>> adjacent_nodes_for_nodes_;
  std::vector<std::array<int, 2>> adjacent_faces_for_edges_;

  mutable std::vector<double> surface_areas_for_faces_;
  mutable std::vector<double> surface_areas_for_nodes_;
  double initial_cell_surface_area_;
  double initial_cell_volume_;
  mutable std::vector<VectorType> normals_for_faces_;
  mutable std::vector<VectorType> normals_for_nodes_;
  mutable std::vector<std::vector<double>>
      initial_lenghts_between_adjacent_nodes_; // has same lengths as adjacent_nodes_for_nodes_
  mutable std::vector<double> initial_curvature_angle_for_edges_;
  mutable std::vector<VectorType> curvatures_for_nodes_;

  double FaceArea(int face_index) const;
  void MakeFacesOriented();
  void CalculateInitialCurvatureAngleForEdges();

};

#endif //DEFORMABLECELLMODEL_CELLMESH_HPP

//
// Created by Nikita Kruk on 15.03.21.
//

#include "EquationsOfMotion.hpp"

#include <unordered_set>
#include <cmath>
#include <iostream>
#include <algorithm> // std::set_difference, std::sort
#include <iterator> // std::back_inserter
#include <unordered_map>

EquationsOfMotion::EquationsOfMotion(const Parameters &parameters) :
    parameters_(parameters),
    mersenne_twister_generator_(std::random_device{}())
{

}

EquationsOfMotion::~EquationsOfMotion() = default;

/*
 * Formulate the model in the form A v = b, and solve for v,
 * where v stands for velocity
 */
void EquationsOfMotion::ComputeForces(const CellMesh &cell_mesh, Eigen::SparseMatrix<double> &A, Eigen::VectorXd &b)
{
  // fill A
  // we will construct a sparse matrix A from a list of triplets
  // note: triplet_list must be the same for all forces
  // otherwise each call to setFromTriplets destroys the previous contents from the matrix A
  std::vector<Eigen::Triplet<double>> triplet_list;

  // \Gamma_{nn}^{cc}
  NodeToNodeFrictionSameCell(cell_mesh, triplet_list);

  // \Gamma_{ns}
  NodeToExtracellularMatrixFriction(cell_mesh, triplet_list);

  A.setFromTriplets(triplet_list.begin(), triplet_list.end());

  // fill b
  // F_{m}
  CytoskeletonBendingElasticityForce(cell_mesh, b);

  // F_{vol}
  VolumePreservationForce(cell_mesh, b);

  // F_{T}
  AreaConservationForce(cell_mesh, b);

  // F_{migration} = F_{random} + F_{morphogen}
  CellMigrationForce(cell_mesh, b);
}

void EquationsOfMotion::NodeToNodeFrictionSameCell(const CellMesh &cell_mesh,
                                                   std::vector<Eigen::Triplet<double>> &triplet_list)
{
  const std::vector<FaceType> &faces = cell_mesh.GetFaces();
  const std::vector<IndexSet> &adjacent_faces = cell_mesh.GetAdjacentFacesForNodes();
  for (int i = 0; i < cell_mesh.GetNumNodes(); ++i)
  {
    // find adjacent neighbors
    std::unordered_set<int> adjacent_neighbors;
    for (int face_index : adjacent_faces[i])
    {
      auto neighbors = faces[face_index];
      adjacent_neighbors.insert(neighbors.begin(), neighbors.end());
    } // face_index
    adjacent_neighbors.erase(i); // remove the current point from the list of its neighbors

    const int node_i_idx = i * kDim;
    for (int j : adjacent_neighbors)
    {
      const int node_j_idx = j * kDim;
      triplet_list.emplace_back(node_i_idx, node_i_idx, parameters_.GetNodalFriction());
      triplet_list.emplace_back(node_i_idx + 1, node_i_idx + 1, parameters_.GetNodalFriction());
      triplet_list.emplace_back(node_i_idx + 2, node_i_idx + 2, parameters_.GetNodalFriction());

      triplet_list.emplace_back(node_i_idx, node_j_idx, -parameters_.GetNodalFriction());
      triplet_list.emplace_back(node_i_idx + 1, node_j_idx + 1, -parameters_.GetNodalFriction());
      triplet_list.emplace_back(node_i_idx + 2, node_j_idx + 2, -parameters_.GetNodalFriction());
    } // j
  } // i
}

void EquationsOfMotion::NodeToExtracellularMatrixFriction(const CellMesh &cell_mesh,
                                                          std::vector<Eigen::Triplet<double>> &triplet_list)
{
  cell_mesh.CalculateFaceSurfaceAreas();
  const std::vector<double> &surface_areas = cell_mesh.CalculateNodeSurfaceAreas();
  for (int i = 0; i < cell_mesh.GetNumNodes(); ++i)
  {
    const int node_i_idx = i * kDim;

    double surface_area = surface_areas[i];
    double cell_ecm_friction = parameters_.GetCellEcmFrictionPerSurfaceArea() * surface_area;

    triplet_list.emplace_back(node_i_idx, node_i_idx, cell_ecm_friction);
    triplet_list.emplace_back(node_i_idx + 1, node_i_idx + 1, cell_ecm_friction);
    triplet_list.emplace_back(node_i_idx + 2, node_i_idx + 2, cell_ecm_friction);
  } // i
}

void EquationsOfMotion::CytoskeletonBendingElasticityForce(const CellMesh &cell_mesh, Eigen::VectorXd &b)
{
  const std::vector<VectorType> &nodes = cell_mesh.GetNodes();
  const std::vector<EdgeType> &edges = cell_mesh.GetEdges();
  const std::vector<std::array<int, 2>> &adjacent_faces_for_edges = cell_mesh.GetAdjacentFacesForEdges();
  const std::vector<FaceType> &faces = cell_mesh.GetFaces();
  const std::vector<VectorType>
      &face_normals = cell_mesh.CalculateFaceNormals(); // todo: optimaze calls to CalculateFaceNormals
  const std::vector<double> &theta_0 = cell_mesh.GetInitialCurvatureAngleForEdges();
  int face_1, face_2;
  int i, j, k, l;
  VectorType force_i, force_j, force_k, force_l;
  const double bending_constant = parameters_.GetCortexYoungsModulus() * std::pow(parameters_.GetCortexThickness(), 3.0)
      / (12.0 * (1.0 - std::pow(parameters_.GetCortexPoissonRatio(), 2.0)));
  double bending_moment = 0.0;
  for (int edge_idx = 0; edge_idx < edges.size(); ++edge_idx)
  {
    k = edges[edge_idx].first;
    j = edges[edge_idx].second;

    face_1 = adjacent_faces_for_edges[edge_idx][0];
    face_2 = adjacent_faces_for_edges[edge_idx][1];

    const std::set<int> edge_nodes{k, j};
    std::vector<int> i_ptr, l_ptr;
    // note: both containers must be sorted
    FaceType sorted_face_1 = faces[face_1], sorted_face_2 = faces[face_2];
    std::sort(sorted_face_1.begin(), sorted_face_1.end());
    std::sort(sorted_face_2.begin(), sorted_face_2.end());
    std::set_difference(sorted_face_1.begin(),
                        sorted_face_1.end(),
                        edge_nodes.begin(),
                        edge_nodes.end(),
                        std::back_inserter(i_ptr));
    std::set_difference(sorted_face_2.begin(),
                        sorted_face_2.end(),
                        edge_nodes.begin(),
                        edge_nodes.end(),
                        std::back_inserter(l_ptr));
    i = i_ptr.front();
    l = l_ptr.front();

    VectorType p_kj = nodes[j] - nodes[k], p_ki = nodes[i] - nodes[k], p_kl = nodes[l] - nodes[k];
    VectorType proj_i = p_ki - p_ki.dot(p_kj) / p_kj.squaredNorm() * p_kj;
    VectorType proj_l = p_kl - p_kl.dot(p_kj) / p_kj.squaredNorm() * p_kj;

    const VectorType &n_1 = face_normals[face_1], &n_2 = face_normals[face_2];
    double theta = std::acos(-n_1.dot(n_2));
    bending_moment = bending_constant * std::sin(theta - theta_0[edge_idx]);
    force_i = bending_moment / proj_i.norm() * n_1;
    force_l = bending_moment / proj_l.norm() * n_2;
    force_j = force_k = -0.5 * (force_i + force_l);

    const std::unordered_map<int, VectorType> force_system{{i, force_i}, {j, force_j}, {k, force_k}, {l, force_l}};
    for (const std::pair<const int, VectorType> &force : force_system)
    {
      const int node_idx = force.first * kDim; // beginning of z-th node values in b
      b[node_idx] += force.second[0];
      b[node_idx + 1] += force.second[1];
      b[node_idx + 2] += force.second[2];
    } // force
  } // edge_idx
}

void EquationsOfMotion::VolumePreservationForce(const CellMesh &cell_mesh, Eigen::VectorXd &b)
{
//  double p = -parameters_.GetCellBulkModulus() * std::log(cell_mesh.CalculateCellVolume() / cell_mesh.GetInitialVolume());
  double p = -parameters_.GetCellBulkModulus() * (cell_mesh.CalculateCellVolume() - cell_mesh.GetInitialVolume())
      / cell_mesh.GetInitialVolume();
  const std::vector<double> &surface_areas = cell_mesh.CalculateNodeSurfaceAreas();
  const std::vector<VectorType> &node_normals = cell_mesh.CalculateNodeNormals();
  for (int i = 0; i < cell_mesh.GetNumNodes(); ++i)
  {
    const int node_idx = i * kDim; // beginning of i-th node values in b
    b[node_idx] += p * surface_areas[i] * node_normals[i][0];
    b[node_idx + 1] += p * surface_areas[i] * node_normals[i][1];
    b[node_idx + 2] += p * surface_areas[i] * node_normals[i][2];
  } // i
}

void EquationsOfMotion::AreaConservationForce(const CellMesh &cell_mesh, Eigen::VectorXd &b)
{
  double area_compression_stiffness = parameters_.GetMembraneAreaCompression() * parameters_.GetRestLength();
  double area_conservation_force =
      -area_compression_stiffness * (cell_mesh.CalculateCellSurfaceArea() - cell_mesh.GetInitialSurfaceArea())
          / cell_mesh.GetInitialSurfaceArea();
  const std::vector<IndexSet> &adjacent_faces = cell_mesh.GetAdjacentFacesForNodes();
  const std::vector<VectorType> &face_normals = cell_mesh.GetNormalsForFaces();
  for (int i = 0; i < cell_mesh.GetNumNodes(); ++i)
  {
    VectorType force = VectorType::Zero();
    for (int face_index : adjacent_faces[i])
    {
      force += face_normals[face_index];
    } // face_index
    force *= area_conservation_force;

    const int node_idx = i * kDim; // beginning of i-th node values in b
    b[node_idx] += force[0];
    b[node_idx + 1] += force[1];
    b[node_idx + 2] += force[2];
  } // i
}

void EquationsOfMotion::CellMigrationForce(const CellMesh &cell_mesh, Eigen::VectorXd &b)
{
  // F_{random}
  StochasticForce(cell_mesh, b);

  // F_{morphogen}
  DirectedForce(cell_mesh, b);
}

void EquationsOfMotion::StochasticForce(const CellMesh &cell_mesh, Eigen::VectorXd &b)
{
  const double mean = 0.0;
  const double spectral_density = parameters_.GetMotility() * parameters_.GetMotility();
  std::normal_distribution<double> normal_distribution(mean, std::sqrt(spectral_density));
  for (int i = 0; i < cell_mesh.GetNumNodes(); ++i)
  {
    const int node_idx = i * kDim; // beginning of i-th node values in b
    b[node_idx] += normal_distribution(mersenne_twister_generator_);
    b[node_idx + 1] += normal_distribution(mersenne_twister_generator_);
    b[node_idx + 2] += normal_distribution(mersenne_twister_generator_);
  } // i
}

void EquationsOfMotion::DirectedForce(const CellMesh &cell_mesh, Eigen::VectorXd &b)
{
//  const double morphogen_strength = parameters_.GetMorphogenStrength();
//  const VectorType morphogen_direction(parameters_.GetMorphogenDirection());
  const std::vector<VectorType> morphogen_directions =
      {{1.0, 0.0, 0.0}};//, {-0.5, 0.86, 0.0}};//, {0.0, -1.0, 0.0}, {0.0, 0.0, 1.0}, {0.0, 0.0, -1.0}};
  const std::vector<double> morphogen_strengths = {0e-11};//, 8e-11, 2e-11, 2e-11, 2e-11};
  double dot_product = 0.0, nonlinear_impact = 0.0;
  const std::vector<VectorType> &node_normals = cell_mesh.CalculateNodeNormals();
  for (int i = 0; i < cell_mesh.GetNumNodes(); ++i)
  {
    const int node_idx = i * kDim; // beginning of i-th node values in b
//    for (const VectorType &morphogen_direction : morphogen_directions)
    for (int m = 0; m < morphogen_directions.size(); ++m)
    {
      const VectorType &morphogen_direction = morphogen_directions[m];
      const double morphogen_strength = morphogen_strengths[m];
      dot_product = morphogen_direction.dot(node_normals[i]);
      if (dot_product < 0.0)
      {
//        nonlinear_impact = std::abs(dot_product);
//        nonlinear_impact = std::exp(-8.0 * (1.0 - std::abs(dot_product)) * (1.0 - std::abs(dot_product)));
        nonlinear_impact = std::exp(-5.0 * (1.0 - std::abs(dot_product)));
        b[node_idx] += morphogen_strength * morphogen_direction[0] * nonlinear_impact;
        b[node_idx + 1] += morphogen_strength * morphogen_direction[1] * nonlinear_impact;
        b[node_idx + 2] += morphogen_strength * morphogen_direction[2] * nonlinear_impact;
      }
    } // m
  } // i
}
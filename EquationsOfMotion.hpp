//
// Created by Nikita Kruk on 15.03.21.
//

#ifndef DEFORMABLECELLMODEL_EQUATIONSOFMOTION_HPP
#define DEFORMABLECELLMODEL_EQUATIONSOFMOTION_HPP

#include "Definitions.hpp"
#include "CellMesh.hpp"
#include "Parameters.hpp"

#include <vector>
#include <random>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/IterativeLinearSolvers>

class EquationsOfMotion
{
 public:

  explicit EquationsOfMotion(const Parameters &parameters);
  ~EquationsOfMotion();

  void ComputeForces(const CellMesh &cell_mesh, Eigen::SparseMatrix<double> &A, Eigen::VectorXd &b);

 private:

  const Parameters &parameters_;
  std::mt19937 mersenne_twister_generator_;
  CellMesh plane_;

  void NodeToNodeFrictionSameCell(const CellMesh &cell_mesh, std::vector<Eigen::Triplet<double>> &triplet_list);
//  void NodeToNodeFrictionDifferentCells(const CellMesh &cell_mesh, std::vector<Eigen::Triplet<double>> &triplet_list);
  void NodeToExtracellularMatrixFriction(const CellMesh &cell_mesh, std::vector<Eigen::Triplet<double>> &triplet_list);

  void CytoskeletonInPlaneElasticityForce(const CellMesh &cell_mesh, Eigen::VectorXd &b);
  void CytoskeletonBendingElasticityForce(const CellMesh &cell_mesh, Eigen::VectorXd &b);
  void VolumePreservationForce(const CellMesh &cell_mesh, Eigen::VectorXd &b);
  void AreaConservationForce(const CellMesh &cell_mesh, Eigen::VectorXd &b);
  void CellToPlaneContactForce(const CellMesh &cell_mesh, Eigen::VectorXd &b);
//  void CellToCellContactForce(const CellMesh &cell_mesh, Eigen::VectorXd &b);
  void CellMigrationForce(const CellMesh &cell_mesh, Eigen::VectorXd &b);
  void StochasticForce(const CellMesh &cell_mesh, Eigen::VectorXd &b);
  void DirectedForce(const CellMesh &cell_mesh, Eigen::VectorXd &b);
  void GravitationalForce(const CellMesh &cell_mesh, Eigen::VectorXd &b);

};

#endif //DEFORMABLECELLMODEL_EQUATIONSOFMOTION_HPP

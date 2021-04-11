//
// Created by Nikita Kruk on 26.03.21.
//

#include "Stepper.hpp"
#include "EquationsOfMotion.hpp"

#include <eigen3/Eigen/IterativeLinearSolvers>

Stepper::Stepper(const Parameters &parameters) :
    parameters_(parameters)
{

}

Stepper::~Stepper() = default;

void Stepper::DoStep(std::vector<CellMesh> &cell_meshes)
{
  // todo: process the colleciton of cells
  CellMesh &cell_mesh = cell_meshes.front();

  const int n = cell_mesh.GetNumNodes() * kDim;
  Eigen::VectorXd b = Eigen::VectorXd::Zero(n);
  Eigen::SparseMatrix<double> A(n, n);

  EquationsOfMotion eom(parameters_);
  eom.ComputeForces(cell_mesh, A, b);

  UpdatePositionsAndVelocities(cell_mesh, A, b);
}

/*
 * Update v by solving A v = b
 * Update x by x += v dt
 */
void Stepper::UpdatePositionsAndVelocities(CellMesh &cell_mesh,
                                           const Eigen::SparseMatrix<double> &A,
                                           const Eigen::VectorXd &b)
{
  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper> conjugate_solver;
  conjugate_solver.compute(A);
  Eigen::VectorXd velocities = conjugate_solver.solve(b);

  std::vector<VectorType> &nodes = cell_mesh.GetNodes();
  VectorType velocity;
  for (int i = 0; i < nodes.size(); ++i)
  {
    const int node_idx = i * kDim; // beginning of i-th node values in v
    velocity = VectorType{velocities[node_idx], velocities[node_idx + 1], velocities[node_idx + 2]};
    nodes[i] += velocity * parameters_.GetDt();
  } // i
}
//
// Created by Nikita Kruk on 15.03.21.
//

#ifndef DEFORMABLECELLMODEL_EQUATIONSOFMOTION_HPP
#define DEFORMABLECELLMODEL_EQUATIONSOFMOTION_HPP

#include "Definitions.hpp"
#include "CellMesh.hpp"
#include "Parameters.hpp"

#include <random>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/IterativeLinearSolvers>

class EquationsOfMotion
{
 public:

  explicit EquationsOfMotion(const Parameters &parameters);
  ~EquationsOfMotion();

  void DoStep(CellMesh &cell_mesh);

 private:

  const Parameters &parameters_;
  std::mt19937 mersenne_twister_generator_;

  void ComputeForces(const CellMesh &cell_mesh, Eigen::SparseMatrix<double> &A, Eigen::VectorXd &b);
  void UpdatePositionsAndVelocities(CellMesh &cell_mesh,
                                    const Eigen::SparseMatrix<double> &A,
                                    const Eigen::VectorXd &b);

};

#endif //DEFORMABLECELLMODEL_EQUATIONSOFMOTION_HPP

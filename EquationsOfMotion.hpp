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

};

#endif //DEFORMABLECELLMODEL_EQUATIONSOFMOTION_HPP

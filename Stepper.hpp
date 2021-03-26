//
// Created by Nikita Kruk on 26.03.21.
//

#ifndef DEFORMABLECELLMODEL_STEPPER_HPP
#define DEFORMABLECELLMODEL_STEPPER_HPP

#include "Definitions.hpp"
#include "Parameters.hpp"
#include "CellMesh.hpp"

#include <vector>
#include <eigen3/Eigen/Sparse>

class Stepper
{
 public:

  explicit Stepper(const Parameters &parameters);
  ~Stepper();

  void DoStep(std::vector<CellMesh> &cell_meshes);

 private:

  const Parameters &parameters_;

  void UpdatePositionsAndVelocities(CellMesh &cell_mesh,
                                    const Eigen::SparseMatrix<double> &A,
                                    const Eigen::VectorXd &b);

};

#endif //DEFORMABLECELLMODEL_STEPPER_HPP

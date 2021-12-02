//
// Created by Nikita Kruk on 15.03.21.
//

#ifndef DEFORMABLECELLMODEL_DEFINITIONS_HPP
#define DEFORMABLECELLMODEL_DEFINITIONS_HPP

#include <array>
#include <set>
#include <utility> // std::pair

#include <eigen3/Eigen/Dense>

const int kDim = 3;
const int kFaceDim = 3;
const int kEdgeDim = 2;

//using VectorType = std::array<double, kDim>; // as an element of a vector space
using VectorType = Eigen::Vector3d; // as an element of a vector space
using FaceType = std::array<int, kFaceDim>;
using EdgeType = std::pair<int, int>;
using IndexSet = std::set<int>;

#endif //DEFORMABLECELLMODEL_DEFINITIONS_HPP

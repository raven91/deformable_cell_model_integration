//
// Created by Nikita Kruk on 15.03.21.
//

#ifndef DEFORMABLECELLMODEL_DEFINITIONS_HPP
#define DEFORMABLECELLMODEL_DEFINITIONS_HPP

#include <array>
#include <set>

const int kDim = 3;
const int kFaceDim = 3;
const int kEdgeDim = 2;

typedef std::array<double, kDim> VectorType; // as an element of a vector space
typedef std::array<int, kFaceDim> FaceType;
typedef std::array<int, kEdgeDim> EdgeType;
typedef std::set<int> IndexSet;

#endif //DEFORMABLECELLMODEL_DEFINITIONS_HPP

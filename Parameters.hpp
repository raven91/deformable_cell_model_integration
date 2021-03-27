//
// Created by Nikita Kruk on 16.03.21.
//

#ifndef DEFORMABLECELLMODEL_PARAMETERS_HPP
#define DEFORMABLECELLMODEL_PARAMETERS_HPP

#include "Definitions.hpp"

#include <unordered_map>
#include <string>
#include <array>

class Parameters
{
 public:

  Parameters();
  ~Parameters();

  double GetDt() const;
  double GetRadius() const;
  double GetCortexYoungsModulus() const;
  double GetCortexThickness() const;
  double GetCortexPoissonRatio() const;
  double GetMembraneAreaCompression() const;
  double GetCellBulkModulus() const;
  double GetNodalFriction() const;
  double GetCellEcmFrictionPerSurfaceArea() const;
  double GetMotility() const;
  double GetMorphogenStrength() const;
  VectorType GetMorphogenDirection() const;
  double GetRestLength() const;

 private:

  std::unordered_map<std::string, double> parameter_dictionary_;

};

#endif //DEFORMABLECELLMODEL_PARAMETERS_HPP

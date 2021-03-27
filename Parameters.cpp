//
// Created by Nikita Kruk on 16.03.21.
//

#include "Parameters.hpp"

#include <fstream>

Parameters::Parameters()
{
  std::string file_name("/Users/nikita/CLionProjects/DeformableCellModel/config.cfg");
  std::ifstream file(file_name, std::ios::in);
  std::string key;
  double value;
  while (file >> key >> value)
  {
    parameter_dictionary_[key] = value;
  }
}

Parameters::~Parameters()
{
  parameter_dictionary_.clear();
}

double Parameters::GetDt() const
{
  return parameter_dictionary_.at("dt");
}

double Parameters::GetRadius() const
{
  return parameter_dictionary_.at("radius");
}

double Parameters::GetCortexYoungsModulus() const
{
  return parameter_dictionary_.at("cortex_youngs_modulus");
}

double Parameters::GetCortexThickness() const
{
  return parameter_dictionary_.at("cortex_thickness");
}

double Parameters::GetCortexPoissonRatio() const
{
  return parameter_dictionary_.at("cortex_poisson_ratio");
}

double Parameters::GetMembraneAreaCompression() const
{
  return parameter_dictionary_.at("membrane_area_compression");
}

double Parameters::GetCellBulkModulus() const
{
  return parameter_dictionary_.at("cell_bulk_modulus");
}

double Parameters::GetNodalFriction() const
{
  return parameter_dictionary_.at("nodal_friction");
}

double Parameters::GetCellEcmFrictionPerSurfaceArea() const
{
  return parameter_dictionary_.at("cell_ecm_friction_per_surface_area");
}

double Parameters::GetMotility() const
{
  return parameter_dictionary_.at("motility");
}

double Parameters::GetMorphogenStrength() const
{
  return parameter_dictionary_.at("morphogen_strength");
}

VectorType Parameters::GetMorphogenDirection() const
{
  return VectorType{parameter_dictionary_.at("morphogen_direction_x"),
                    parameter_dictionary_.at("morphogen_direction_y"),
                    parameter_dictionary_.at("morphogen_direction_z")};
}

double Parameters::GetRestLength() const
{
  return GetRadius();
}
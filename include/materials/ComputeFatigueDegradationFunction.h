//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "ADMaterial.h"
#include "BaseNameInterface.h"

class ComputeFatigueDegradationFunction : public ADMaterial,
                                          public BaseNameInterface
{
public:
  static InputParameters validParams();

  ComputeFatigueDegradationFunction(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;

  /// active part of elastic energy from coupled variable
  const VariableValue * _sigma_1;
  const VariableValue * _sigma_1_old;

  /// fatigue history variable
  ADMaterialProperty<Real> & _alpha_bar;

  /// old fatigue history variable
  const MaterialProperty<Real> & _alpha_bar_old;

  /// fatigue degradation function
  ADMaterialProperty<Real> & _f_alpha;

  /// fatigue degradation function type
  const MaterialPropertyName _f_alpha_type;

  /// parameter in logaritmic fatigue degradation function
  const Real & _kappa;
  const Real & _k;
  const Real & _p;

  /// fatigue threshold
  const Real & _alpha_T;

  /// initial fatigue history variable (alpha_bar)
  const VariableValue * _alpha_bar_init;

  /// nucleation energy
  const ADMaterialProperty<Real> & _sigma_c;

  /// Fatigue flag to check if fatigue has occured
  MaterialProperty<bool> & _fatigue_flag;
};

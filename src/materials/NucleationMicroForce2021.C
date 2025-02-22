//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "Function.h"
#include "NucleationMicroForce2021.h"

registerADMooseObject("raccoonApp", NucleationMicroForce2021);

InputParameters
NucleationMicroForce2021::validParams()
{
  InputParameters params = Material::validParams();
  params += BaseNameInterface::validParams();

  params.addClassDescription(
      "This class computes the external driving force for nucleation given "
      "a Drucker-Prager strength envelope. The implementation follows Kumar et. al. (2021).");

  params.addParam<MaterialPropertyName>(
      "fracture_toughness", "Gc", "energy release rate or fracture toughness");
  params.addParam<MaterialPropertyName>(
      "normalization_constant", "c0", "The normalization constant $c_0$");
  params.addParam<MaterialPropertyName>(
      "regularization_length", "l", "the phase field regularization length");

  params.addParam<MaterialPropertyName>("lambda", "lambda", "Lame's first parameter lambda");
  params.addParam<MaterialPropertyName>("shear_modulus", "G", "shear modulus mu or G");

  params.addRequiredParam<Real>(
      "tensile_strength", "The tensile strength of the material beyond which the material fails.");
  params.addRequiredParam<Real>(
      "compressive_strength",
      "The compressive strength of the material beyond which the material fails.");

  params.addRequiredParam<Real>("delta", "delta");
  params.addParam<MaterialPropertyName>(
      "external_driving_force_name",
      "ex_driving",
      "Name of the material that holds the external_driving_force");
  params.addParam<MaterialPropertyName>(
      "stress_balance_name", "stress_balance", "Name of the stress balance");
  params.addParam<MaterialPropertyName>("stress_name", "stress", "Name of the stress tensor");
  return params;
}

NucleationMicroForce2021::NucleationMicroForce2021(const InputParameters & parameters)
  : Material(parameters),
    BaseNameInterface(parameters),
    _ex_driving(declareADProperty<Real>(prependBaseName("external_driving_force_name", true))),
    _Gc(getADMaterialProperty<Real>(prependBaseName("fracture_toughness", true))),
    _c0(getADMaterialProperty<Real>(prependBaseName("normalization_constant", true))),
    _L(getADMaterialProperty<Real>(prependBaseName("regularization_length", true))),
    _lambda(getADMaterialProperty<Real>(prependBaseName("lambda", true))),
    _mu(getADMaterialProperty<Real>(prependBaseName("shear_modulus", true))),
    _sigma_ts(getParam<Real>("tensile_strength")),
    _sigma_cs(getParam<Real>("compressive_strength")),
    _delta(getParam<Real>("delta")),
    _stress(getADMaterialProperty<RankTwoTensor>(prependBaseName("stress_name", true))),
    _stress_balance(declareADProperty<Real>(prependBaseName("stress_balance_name", true)))
{
}

void
NucleationMicroForce2021::computeQpProperties()
{
  // The bulk modulus
  ADReal K = _lambda[_qp] + 2 * _mu[_qp] / 3;

  // Parameters in the strength surface
  ADReal gamma_0 =
      _sigma_ts / 6.0 / (3.0 * _lambda[_qp] + 2.0 * _mu[_qp]) + _sigma_ts / 6.0 / _mu[_qp];
  ADReal gamma_1 = (1.0 + _delta) / (2.0 * _sigma_ts * _sigma_cs);
  // ADReal gamma_2 = (8 * _mu[_qp] + 24 * K - 27 * _sigma_ts) / 144 / _mu[_qp] / K;

  // The mobility
  ADReal M = _Gc[_qp] / _L[_qp] / _c0[_qp];

  // Invariants of the stress
  ADReal I1 = _stress[_qp].trace();
  ADRankTwoTensor stress_dev = _stress[_qp].deviatoric();
  ADReal J2 = 0.5 * stress_dev.doubleContraction(stress_dev);
  // ADReal M = _Gc[_qp] * 3.0 / _L[_qp] / 8.0;
  // ADReal I1 = _stress[_qp](0, 0) + _stress[_qp](1, 1) + _stress[_qp](2, 2);
  // ADReal J2 = (pow(_stress[_qp](0, 0) - _stress[_qp](1, 1), 2) +
  //            pow(_stress[_qp](1, 1) - _stress[_qp](2, 2), 2) +
  //            pow(_stress[_qp](0, 0) - _stress[_qp](2, 2), 2)) /
  //               6.0 +
  //           pow(_stress[_qp](0, 1), 2) + pow(_stress[_qp](0, 2), 2) +
  //           pow(_stress[_qp](2, 1), 2);

  // Just to be extra careful... J2 is for sure non-negative.
  mooseAssert(J2 >= 0, "Negative J2");

  // define zero J2's derivative
  if (MooseUtils::absoluteFuzzyEqual(J2, 0))
    J2.value() = libMesh::TOLERANCE * libMesh::TOLERANCE;

  // Compute the external driving force required to recover the desired strength envelope.
  ADReal beta_0 = _delta * M;
  ADReal beta_1 = -gamma_1 * M * (_sigma_cs - _sigma_ts) - gamma_0;
  ADReal beta_2 = std::sqrt(3.0) * (-gamma_1 * M * (_sigma_cs + _sigma_ts) + gamma_0);
  ADReal beta_3 = _L[_qp] * _sigma_ts / _mu[_qp] / K / _Gc[_qp];
  _ex_driving[_qp] =
      (beta_2 * std::sqrt(J2) + beta_1 * I1 + beta_0) +
      (1.0 - std::sqrt(I1 * I1) / I1) *
          (J2 / 2.0 / _mu[_qp] + I1 * I1 / 6.0 / (3.0 * _lambda[_qp] + 2.0 * _mu[_qp]));
  _stress_balance[_qp] = J2 / _mu[_qp] + pow(I1, 2) / 9.0 / K - _ex_driving[_qp] - M;
}

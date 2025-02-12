//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "ComputeFatigueDegradationFunction.h"
#include "Function.h"
#include "MooseMesh.h"
#include <cmath>

registerMooseObject("raccoonApp", ComputeFatigueDegradationFunction);

InputParameters
ComputeFatigueDegradationFunction::validParams()
{
  InputParameters params = ADMaterial::validParams();
  params += BaseNameInterface::validParams();

  params.addClassDescription("computes fatigue degradation function");
  params.addCoupledVar("sigma_1",
                       "principalStress_max");
  params.addParam<MaterialPropertyName>(
      "alpha_bar_name",
      "alpha_bar",
      "name of the parameter used in fatigue degradation function equation");
  params.addParam<MaterialPropertyName>(
      "f_alpha_name", "f_alpha", "name of the fatigue degradation function");
  params.addParam<MaterialPropertyName>(
      "f_alpha_type",
      "f_alpha_type",
      "type of the fatigue degradation function (logarithmic or asymptotic");
  params.addParam<Real>("kappa", 0, "parameter in logaritmic fatigue degaradtion function");
  params.addParam<Real>("k", 0, "residual in asmpytotic fatigue degaradtion function");
  params.addParam<Real>(
      "p", 1, "constribution parameter in asmpytotic fatigue degaradtion function");
  params.addRequiredParam<Real>("alpha_T", "fatigue threshold");
  params.addCoupledVar("initial_alpha_bar", "initial value for alpha_bar");
  // params.addParam<MaterialPropertyName>(
  //     "degradation_mat", "g", "name of the material that holds the degradation");
  params.addRequiredParam<MaterialPropertyName>(
      "sigma_c", "The strength of the material beyond which the fatigue effect accumulates.");
  return params;
}

ComputeFatigueDegradationFunction::ComputeFatigueDegradationFunction(
    const InputParameters & parameters)
  : ADMaterial(parameters),
    BaseNameInterface(parameters),
    _sigma_1(isParamValid("sigma_1") ? &coupledValue("sigma_1")
                                                        : nullptr),
    _sigma_1_old(isParamValid("sigma_1") ? &coupledValueOld("sigma_1")
                                                            : nullptr),
    _alpha_bar(declareADProperty<Real>(getParam<MaterialPropertyName>("alpha_bar_name"))),
    _alpha_bar_old(getMaterialPropertyOld<Real>("alpha_bar_name")),
    _f_alpha(declareADProperty<Real>(getParam<MaterialPropertyName>("f_alpha_name"))),
    _f_alpha_type(getParam<MaterialPropertyName>("f_alpha_type")),
    _kappa(getParam<Real>("kappa")),
    _k(getParam<Real>("k")),
    _p(getParam<Real>("p")),
    _alpha_T(getParam<Real>("alpha_T")),
    _alpha_bar_init(isParamValid("initial_alpha_bar") ? &coupledValue("initial_alpha_bar")
                                                      : nullptr),
    _sigma_c(getADMaterialProperty<Real>(prependBaseName("sigma_c", true))),
    _fatigue_flag(declareProperty<bool>("fatigue_flag"))
{
}

void
ComputeFatigueDegradationFunction::initQpStatefulProperties()
{
  _alpha_bar[_qp] = _alpha_bar_init ? (*_alpha_bar_init)[_qp] : 0.0;
  _fatigue_flag[_qp] = false;
}

void
ComputeFatigueDegradationFunction::computeQpProperties()
{
  ADReal sigma_1 = (*_sigma_1)[_qp];
  ADReal sigma_1_old = (*_sigma_1_old)[_qp];

  // calculate degraded active elastic energy
  ADReal alpha = ((sigma_1 - _sigma_c[_qp]) > 0) ? ((sigma_1- _sigma_c[_qp]) / _sigma_c[_qp]) : 0;
  ADReal alpha_old =
      ((sigma_1_old - _sigma_c[_qp]) > 0) ? ((sigma_1_old - _sigma_c[_qp]) / _sigma_c[_qp]) : 0;

  // update alpha_bar
  if (alpha > alpha_old)
    _alpha_bar[_qp] = _alpha_bar_old[_qp] + alpha * _dt;
  else
    _alpha_bar[_qp] = _alpha_bar_old[_qp];

  if (_alpha_bar[_qp] > _alpha_T)
    _fatigue_flag[_qp] = true;

  // calculate f_alpha
  // asymptotic
  if (_f_alpha_type == "asymptotic")
  {

    if (_alpha_bar[_qp] < _alpha_T)
      _f_alpha[_qp] = 1.0;
    else
      _f_alpha[_qp] =
          (1 - _k) * std::pow(2 * _alpha_T / ((2 - _p) * _alpha_T + _p * _alpha_bar[_qp]), 2.0) +
          _k;

    // logarithmic
  }
  else if (_f_alpha_type == "logarithmic")
  {

    if (_alpha_bar[_qp] < _alpha_T)
      _f_alpha[_qp] = 1.0;
    else if ((_alpha_bar[_qp] > _alpha_T) &&
             (_alpha_bar[_qp] < _alpha_T * std::pow(10.0, 1.0 / _kappa)))
      _f_alpha[_qp] = std::pow(1 - _kappa * std::log10(_alpha_bar[_qp] / _alpha_T), 2.0);
    else
      _f_alpha[_qp] = 0.0;
  }
  else
    mooseError("Invalid type of fatigue degaradtion fucntion");
}

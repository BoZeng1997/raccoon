//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "ADPenaltyVariableDirichletBC.h"

registerMooseObject("raccoonApp", ADPenaltyVariableDirichletBC);

InputParameters
ADPenaltyVariableDirichletBC::validParams()
{
  InputParameters params = ADIntegratedBC::validParams();
  params.addRequiredParam<Real>("penalty", "Penalty scalar");
  params.addRequiredCoupledVar("v", "The variable whose value we are to match.");
  params.addParam<Real>("u_coeff", 1.0, " A coefficient for primary variable u");
  params.addParam<Real>("v_coeff", 1.0, " A coefficient for coupled variable v");
  //   params.addParam<Real>("value", 0.0, "Boundary value of the variable");
  //   params.declareControllable("value");
  params.addClassDescription("Enforces a Dirichlet boundary condition "
                             "in a weak sense by penalizing differences between the current "
                             "solution and the Dirichlet data.");
  return params;
}

ADPenaltyVariableDirichletBC::ADPenaltyVariableDirichletBC(const InputParameters & parameters)
  : ADIntegratedBC(parameters),
    _p(getParam<Real>("penalty")),
    //   _v(getParam<Real>("value"))
    _v(adCoupledValue("v")),
    _u_coeff(getParam<Real>("u_coeff")),
    _v_coeff(getParam<Real>("v_coeff"))
{
}

ADReal
ADPenaltyVariableDirichletBC::computeQpResidual()
{
  return _p * _test[_i][_qp] * (-_v_coeff * _v[_qp] + _u_coeff * _u[_qp]);
}

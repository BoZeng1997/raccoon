//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "ADIntegratedBC.h"

/**
 * Weakly enforce a Dirichlet BC using a penalty term. This class is
 * an alternative to the DirichletBC that maintains the symmetry (if
 * any) present in the original problem, and does not involve
 * explicitly zeroing matrix rows for its implementation. The main
 * drawback of this approach is that the penalty parameter must tend
 * to infinity in order for the constraint to be satisfied in the
 * limit as h->0, and this causes the Jacobian matrix to be
 * ill-conditioned.
 *
 * The weak form contribution for this term is:
 *
 * \f$ (p (u - g), \psi)_{\Gamma} \f$,
 *
 * where:
 * p = penalty parameter (user-selectable)
 * u = the unknown
 * g = Dirichlet data (given)
 * Gamma = the part of the boundary where the penalty BC is applied.
 */
class ADPenaltyVariableDirichletBC : public ADIntegratedBC
{
public:
  static InputParameters validParams();

  ADPenaltyVariableDirichletBC(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

private:
  /// Penalty value
  const Real _p;

  /// Value of u on the boundary
  //   const Real & _v;
  const ADVariableValue & _v;

  /// Coefficient for primary variable
  const Real _u_coeff;
  /// Coefficient for coupled variable
  const Real _v_coeff;
};
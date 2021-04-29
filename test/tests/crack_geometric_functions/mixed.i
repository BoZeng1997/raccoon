[Problem]
  solve = false
[]

[Mesh]
  [one]
    type = GeneratedMeshGenerator
    dim = 1
  []
[]

[AuxVariables]
  [d]
  []
[]

[AuxKernels]
  [d]
    type = FunctionAux
    variable = d
    function = 't'
  []
[]

[Materials]
  [crack_geometric]
    type = CrackGeometricFunction
    f_name = alpha
    function = 'xi*d+(1-xi)*d^2'
    d = d
    parameter_names = 'xi'
    parameter_values = '2'
    tolerance = 1e-8
  []
[]

[Postprocessors]
  [d]
    type = ElementAverageValue
    variable = d
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [alpha]
    type = ADElementAverageMaterialProperty
    mat_prop = alpha
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [dalpha_dd]
    type = ADElementAverageMaterialProperty
    mat_prop = dalpha/dd
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [dalpha_dd_0]
    type = ADElementAverageMaterialProperty
    mat_prop = dalpha/dd(0)
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [c0]
    type = ADElementAverageMaterialProperty
    mat_prop = c0
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Executioner]
  type = Transient
  dt = 0.1
  end_time = 1
[]

[Outputs]
  exodus = true
[]
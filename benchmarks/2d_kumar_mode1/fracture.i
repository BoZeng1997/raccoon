[Mesh]
  [top_half]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 50
    ny = 25
    ymin = 0
    ymax = 0.5
    boundary_id_offset = 0
    boundary_name_prefix = top_half
  []
  [top_stitch] #top_stitch
    type = BoundingBoxNodeSetGenerator
    input = top_half
    new_boundary = top_stitch
    bottom_left = '0.5 0 0'
    top_right = '1 0 0'
  []
  [bottom_half]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 50
    ny = 25
    ymin = -0.5
    ymax = 0
    boundary_id_offset = 5
    boundary_name_prefix = bottom_half
  []
  [bottom_stitch]
    type = BoundingBoxNodeSetGenerator
    input = bottom_half
    new_boundary = bottom_stitch
    bottom_left = '0.5 0 0'
    top_right = '1 0 0'
  []
  [stitch]
    type = StitchedMeshGenerator
    inputs = 'top_stitch bottom_stitch'
    stitch_boundaries_pairs = 'top_stitch bottom_stitch'
  []
  construct_side_list_from_node_list = true
[]

# [Adaptivity]
#   marker = marker
#   initial_marker = marker
#   initial_steps = 2
#   stop_time = 0
#   max_h_level = 2
#   [Markers]
#     [marker]
#       type = BoxMarker
#       bottom_left = '0.4 -0.15 0'
#       top_right = '1 0.15 0'
#       outside = DO_NOTHING
#       inside = REFINE
#     []
#   []
# []

[Variables]
  [d]
  []
[]

[AuxVariables]
  [bounds_dummy]
  []
  [psie_active]
    order = CONSTANT
    family = MONOMIAL
  []
  [ce]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[Bounds]
  [irreversibility]
    type = VariableOldValueBoundsAux
    variable = bounds_dummy
    bounded_variable = d
    bound_type = lower
  []
  [upper]
    type = ConstantBoundsAux
    variable = bounds_dummy
    bounded_variable = d
    bound_type = upper
    bound_value = 1
  []
[]

[Kernels]
  [diff]
    type = ADPFFDiffusion
    variable = d
    fracture_toughness = Gc
    regularization_length = l
    normalization_constant = c0
  []
  [source]
    type = ADPFFSource
    variable = d
    free_energy = psi
  []
[]

[Materials]
  [fracture_properties]
    type = ADGenericConstantMaterial
    prop_names = 'Gc l Lambda G'
    prop_values = '${Gc} ${l} ${Lambda} ${G}'
  []
  [degradation]
    type = PowerDegradationFunction
    f_name = g
    function = (1-d)^p*(1-eta)+eta
    phase_field = d
    parameter_names = 'p eta '
    parameter_values = '2 0'
  []
  [crack_geometric]
    type = CrackGeometricFunction
    f_name = alpha
    function = 'd'
    phase_field = d
  []
  [psi]
    type = ADDerivativeParsedMaterial
    f_name = psi
    function = '2*d/c0*(4.0/3.0*d*psie_active-8.0/3.0*psie_active+4.0/3.0*ce+Gc/2.0/l)'
    args = 'd psie_active'
    material_property_names = 'Gc c0 l ce' #alpha(d) g(d)
    derivative_order = 1
  []
  # [kumar_material]
  #   type = GeneralizedExternalDrivingForce
  #   first_invariant = invar_1
  #   second_invariant = invar_2
  #   tensile_strength = '${sigma_ts}' #27MPa
  #   compressive_strength = '${sigma_cs}' #77MPa
  #   delta = '${delta}'
  #   energy_release_rate = '${Gc}'
  #   phase_field_regularization_length = '${l}'
  #   Lame_first_parameter = '${Lambda}'
  #   shear_modulus = '${G}'
  #   external_driving_force_name = ce
  # []

[]

[Postprocessors]
  # [extdriving]
  #   type = ADElementAverageMaterialProperty
  #   mat_prop = 'ce'
  # []
  [extdriving_v]
    type = ElementAverageValue
    variable = 'ce'
  []
  [d_avg]
    type = AverageNodalVariableValue
    variable = d
  []
  # [invar_1]
  #   type = AverageNodalVariableValue
  #   variable = invar_1
  # []
  # [invar_2]
  #   type = AverageNodalVariableValue
  #   variable = invar_2
  # []
[]

[VectorPostprocessors]
  [damage]
    type = NodalValueSampler
    variable = 'd'
    sort_by = id
  []
  [ext]
    type = ElementValueSampler
    variable = 'ce'
    sort_by = id
  []
[]

[Executioner]
  type = Transient

  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -snes_type'
  petsc_options_value = 'lu       superlu_dist                  vinewtonrsls'
  automatic_scaling = true

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
[]

[Outputs]
  [csv_]
    type = CSV
    file_base = kumar_mode1_2d_Gc4L0.35del4.41_frac
    append_date = true
    execute_vector_postprocessors_on = final
  []
  print_linear_residuals = false
[]
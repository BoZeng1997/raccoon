[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = ${nx}
    ny = ${ny}
    xmax = 30
    ymin = -5
    ymax = 5
  []
[]

[Adaptivity]
  marker = marker
  initial_marker = marker
  initial_steps = ${refine}
  stop_time = 0
  max_h_level = ${refine}
  [Markers]
    [marker]
      type = BoxMarker
      bottom_left = '0 -0.7 0'
      top_right = '10 0.7 0'
      outside = DO_NOTHING
      inside = REFINE
    []
  []
[]

[Variables]
  [d]
    [InitialCondition]
      type = FunctionIC
      function = 'if(y=0&x>=0&x<=5,1,0)'
    []
  []
[]

[AuxVariables]
  [bounds_dummy]
  []
  [disp_x]
  []
  [disp_y]
  []
  [strain_zz]
  []
  [psie_active]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[Bounds]
  [conditional]
    type = ConditionalBoundsAux
    variable = bounds_dummy
    bounded_variable = d
    fixed_bound_value = 0.0
    threshold_value = 0.0
  []
  [upper]
    type = ConstantBounds
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
    fracture_toughness = Gc_delta
    regularization_length = l
    normalization_constant = c0
  []
  [source]
    type = ADPFFSource
    variable = d
    free_energy = psi
  []
  [nuc_force]
    type = ADCoefMatSource
    variable = d
    # prop_names = 'ce_deg'
    prop_names = 'ce'
    # coefficient = -1.0
  []
  # [dummy]
  #   type = NullKernel
  #   variable = d
  # []
[]

[Materials]
  [fracture_properties]
    type = ADGenericConstantMaterial
    prop_names = 'E K G lambda Gc l sigma_ts sigma_hs'
    prop_values = '${E} ${K} ${G} ${Lambda} ${Gc} ${l} ${sigma_ts} ${sigma_hs}'
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
    output_properties = c0
    outputs = exodus
  []
  [psi]
    type = ADDerivativeParsedMaterial
    f_name = psi
    function = 'g*psie_active + (Gc*delta/c0/l)*alpha ' #+ (1-d)^3/3*ce
    args = 'd psie_active'
    material_property_names = 'delta alpha(d) g(d) Gc c0 l'
    derivative_order = 1
  []
  [Gc_delta]
    type = ADParsedMaterial
    property_name = Gc_delta
    expression = 'Gc*delta'
    material_property_names = 'Gc delta'
    output_properties = 'Gc_delta'
    outputs = exodus
  []
  # [kumar_material]
  #   type = KLRNucleationMicroForce
  #   phase_field = d
  #   # type = KLBFNucleationMicroForce
  #   normalization_constant = c0
  #   tensile_strength = sigma_ts
  #   compressive_strength = sigma_cs
  #   delta = delta
  #   external_driving_force_name = ce
  #   output_properties = 'ce'
  # []
  [nucleation_micro_force]
    type = LDLNucleationMicroForce
    regularization_length = l
    normalization_constant = c0
    fracture_toughness = Gc
    tensile_strength = sigma_ts
    hydrostatic_strength = sigma_hs
    delta = delta
    external_driving_force_name = ce
    stress_balance_name = stress_balance
    output_properties = 'ce delta stress_balance'
    outputs = exodus
    h = 0.0
  []
  # [ce_deg]
  #   type = ADParsedMaterial
  #   property_name = ce_deg
  #   expression = 'ce*(1-d)^2'
  #   material_property_names = 'ce'
  #   coupled_variables = 'd'
  #   output_properties = 'ce_deg'
  #   outputs = exodus
  # []
  [strain]
    type = ADComputePlaneSmallStrain
    out_of_plane_strain = 'strain_zz'
    displacements = 'disp_x disp_y'
  []
  [elasticity]
    type = SmallDeformationIsotropicElasticity
    bulk_modulus = K
    shear_modulus = G
    phase_field = d
    degradation_function = g
    decomposition = NONE
  []
  [stress]
    type = ComputeSmallDeformationStress
    elasticity_model = elasticity
    output_properties = 'stress'
  []
[]

[Executioner]
  type = Transient

  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -snes_type'
  petsc_options_value = 'lu       superlu_dist                  vinewtonrsls'
  automatic_scaling = true
  line_search = bt
  compute_scaling_once = true
  nl_max_its = 100

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
[]

[Outputs]
  [exodus]
    type = Exodus
  []
  print_linear_residuals = false
  file_base = '${filebase}_frac'
[]

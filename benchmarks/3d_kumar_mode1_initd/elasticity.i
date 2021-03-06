E = 9.8e3 #9.8Gpa
nu = 0.13
K = '${fparse E/3/(1-2*nu)}'
G = '${fparse E/2/(1+nu)}'
Lambda = '${fparse E*nu/(1+nu)/(1-2*nu)}'

Gc = 9.1e-4 # 91N/m
l = 0.20 #0.5
sigma_ts = 27
sigma_cs = 77
delta = 3.22 #13

[MultiApps]
  [fracture]
    type = TransientMultiApp
    input_files = fracture.i
    cli_args = 'Gc=${Gc};l=${l};G=${G};Lambda=${Lambda};sigma_ts=${sigma_ts};sigma_cs=${sigma_cs};delta=${delta}'
    execute_on = 'TIMESTEP_END'
  []
[]

[Transfers]
  [from_d]
    type = MultiAppCopyTransfer
    multi_app = fracture
    direction = from_multiapp
    variable = d
    source_variable = d
  []
  [to_psie_active]
    type = MultiAppCopyTransfer
    multi_app = fracture
    direction = to_multiapp
    variable = 'psie_active invar_1 invar_2'
    source_variable ='psie_active invar_1 invar_2'
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  [top_half]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 10
    ny = 100
    nz = 20
    ymin = -5
    ymax = 5
    zmax = 2
    boundary_id_offset = 0
    boundary_name_prefix = top_half
  []
  # [noncrack]
  #   type = BoundingBoxNodeSetGenerator
  #   input = top_half
  #   new_boundary = noncrack
  #   bottom_left = '0.5 0 0'
  #   top_right = '1 0 0.2'
  # []
  # construct_side_list_from_node_list = true
[]

[Adaptivity]
  marker = marker
  initial_marker = marker
  initial_steps = 2
  stop_time = 0
  max_h_level = 2
  [Markers]
    [marker]
      type = BoxMarker
      bottom_left = '0 -0.4 0'
      top_right = '1 0.4 0.2'
      outside = DO_NOTHING
      inside = REFINE
    []
  []
[]
[Variables]
  [disp_x]
  []
  [disp_y]
  []
  [disp_z]
  []
[]

[AuxVariables]
  [fy]
  []
  [d]
    [InitialCondition]
     type = FunctionIC
     function = 'if(x<=0.5&y=0,0.5,0)'
    []
  []
[]


[Kernels]
  [solid_x]
    type = ADStressDivergenceTensors
    variable = disp_x
    component = 0
    displacements = 'disp_x disp_y disp_z'
  []
  [solid_y]
    type = ADStressDivergenceTensors
    variable = disp_y
    component = 1
    displacements = 'disp_x disp_y disp_z'
    save_in = fy
  []
  [solid_z]
    type = ADStressDivergenceTensors
    variable = disp_z
    component = 2
    displacements = 'disp_x disp_y disp_z'
  []
[]

[BCs]
  # [left_x]
  #   type = DirichletBC
  #   variable = disp_x
  #   boundary = left
  #   value = 0
  # []
  # [left_y]
  #   type = DirichletBC
  #   variable = disp_y
  #   boundary = left
  #   value = 0
  # []
  # [right_x]
  #   type = DirichletBC
  #     variable = disp_x
  #     boundary = 'top_half_right bottom_half_right'
  #     value = 0
  # []
  # [right_y]
  #   type = DirichletBC
  #   variable = disp_y
  #   boundary = 'top_half_right bottom_half_right'
  #   value = 0
  # []
  # [right_z]
  #   type = DirichletBC
  #   variable = disp_z
  #   boundary = 'top_half_right bottom_half_right'
  #   value = 0
  # []
  [bottom_x]
    type = DirichletBC
    variable = disp_x
    boundary = top_half_bottom
    value = 0
  []
  [bottom_y]
    type = DirichletBC
    variable = disp_y
    boundary = top_half_bottom
    value = 0
  []
  [bottom_z]
    type = DirichletBC
    variable = disp_z
    boundary = top_half_bottom
    value = 0
  []
  # [top_x]
  #   type = DirichletBC
  #   variable = disp_x
  #   boundary = top
  #   value = 0
  # []
  [top_y]
    type = FunctionDirichletBC
      variable = disp_y
      boundary = top_half_top
      function = 't'
  []
  # [back_x]
  #   type = DirichletBC
  #   variable = disp_x
  #   boundary = noncrack
  #   value = 0
  # []
  # [back_y]
  #   type = DirichletBC
  #   variable = disp_y
  #   boundary = noncrack
  #   value = 0
  # []
  # [back_z]
  #   type = DirichletBC
  #   variable = disp_z
  #   boundary = back
  #   value = 0
  # []
  # [front_x]
  #   type = DirichletBC
  #   variable = disp_x
  #   boundary = front
  #   value = 0
  # []
  # [front_y]
  #   type = DirichletBC
  #   variable = disp_y
  #   boundary = front
  #   value = 0
  # []
  # [front_z]
  #   type = FunctionDirichletBC
  #     variable = disp_z
  #     boundary = front
  #     function = 't*0.03'
  # []
  # [front_z]
  #   type = FunctionDirichletBC
  #   variable = disp_z
  #   boundary = front
  #   function = 't/100'
  # []
[]

[Materials]
  [bulk]
    type = ADGenericConstantMaterial
    prop_names = 'K G'
    prop_values = '${K} ${G}'
  []
  [degradation]
    type = PowerDegradationFunction
    f_name = g
    function = (1-d)^p*(1-eta)+eta
    phase_field = d
    parameter_names = 'p eta '
    parameter_values = '2 0'
  []
  [strain]
    type = ADComputeSmallStrain
  []
  [elasticity]
    type = SmallDeformationIsotropicElasticity
    bulk_modulus = K
    shear_modulus = G
    phase_field = d
    degradation_function = g
    decomposition = NONE
    output_properties = 'elastic_strain psie_active'
    outputs = exodus
  []
  [stress]
    type = ComputeSmallDeformationStress
    elasticity_model = elasticity
    output_properties = 'stress'
    outputs = exodus
  []
  [I1]
    type = ADRankTwoInvariant
    property_name = invar_1
    rank_two_tensor = 'stress'
    invariant = FirstInvariant
    output_properties = 'invar_1'
    outputs = exodus
  []
  [I2]
    type = ADRankTwoInvariant
    property_name = invar_2
    rank_two_tensor = 'stress'
    invariant = SecondInvariant
    output_properties = 'invar_2'
    outputs = exodus
  []
[]

[Postprocessors]
  [psie_active]
    type = ADElementAverageMaterialProperty
    mat_prop = psie_active
  []
  [Fy]
    type = NodalSum
    variable = fy
    boundary = top_half_top
  []
  # [sigma_00_0]
  #   type = ElementalVariableValue
  #   elementid = 0
  #   variable = sigma_00
  # []
[]

[VectorPostprocessors]
  [nodal]
    type = NodalValueSampler
    variable = 'd disp_x disp_y disp_z' # sigma_00 sigma_11 sigma_22 invar_1 invar_2'
    sort_by = id
  []
[]

[Executioner]
  type = Transient

  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu       superlu_dist                 '
  automatic_scaling = true

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10

  dt = 1e-5
  end_time = 1e-2

  picard_max_its = 20
  accept_on_max_picard_iteration = true
  picard_rel_tol = 1e-8
  picard_abs_tol = 1e-10
[]

[Outputs]
  [csv_]
type = CSV
file_base = kumar_mode1_initd_Gc4L0.2del3.22_ela
append_date = true
#show = 'var_u'
execute_vector_postprocessors_on = final
[]
  exodus = true
  file_base = kumar_mode1_initd_Gc4L0.2del3.22
  append_date = true
  print_linear_residuals = false
[]

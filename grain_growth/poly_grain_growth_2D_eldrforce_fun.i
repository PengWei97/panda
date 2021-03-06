[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 20
  ny = 20
  nz = 0
  xmax = 1000
  ymax = 1000
  zmax = 0
  elem_type = QUAD4
  uniform_refine = 2
[]

[GlobalParams]
  op_num = 8
  var_name_base = gr
  grain_num = 36
[]

[Variables]
  [./PolycrystalVariables]
  # 根据op_num自动创建序参数的数目
  # 本例中为8个变量，都使用一阶拉格朗日形函数
  [../]
  [./disp_x] # 这是位移变量
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[UserObjects]
  [./euler_angle_file]
    type = EulerAngleFileReader
    # Read Euler angle data from a file and provide it to other objects.
    file_name = grn_36_rand_2D.tex
  [../]
  [./voronoi]
    type = PolycrystalVoronoi
    # 根据晶粒数目来产生晶体结构
    coloring_algorithm = bt
  [../]
  [./grain_tracker]
    type = GrainTrackerElasticity
    # 用于无晶粒粗化来减少序参数的数目
    # Grain Tracker object for running reduced order parameter simulations without grain coalescence.
    # 将初始弹性模量根据每个晶粒的欧拉角来赋予给每个晶粒
    threshold = 0.2
    # The threshold value for which a new feature may be started
    compute_var_to_feature_map = true
    # Instruct the Postprocessor to compute the active vars to features map
    execute_on = 'initial timestep_begin'
    flood_entity_type = ELEMENTAL
    # Determines whether the flood algorithm runs on nodes or elements

    C_ijkl = '1.27e5 0.708e5 0.708e5 1.27e5 0.708e5 1.27e5 0.7355e5 0.7355e5 0.7355e5' # copper
    # C_ijkl = '1.94e5 0.655e5 0.698e5 1.94e5 0.698e5 1.98e5 0.4627e5 0.4627e5 0.6435e5' # Titanium,2,0Pa，可行
    # C_ijkl = '1.94305e5 0.65597e5 0.69870e5 1.94305e5 0.69870e5 1.98907e5 0.46270e5 0.46270e5 0.64354e5' # Titanium,2,0Pa，不可行
    fill_method = symmetric9
    euler_angle_provider = euler_angle_file
    # <--UserObjects/euler_angle_file
  [../]
[]

[ICs]
  [./PolycrystalICs]
    [./PolycrystalColoringIC]
      polycrystal_ic_uo = voronoi
    [../]
  [../]
[]

[AuxVariables]
  [./bnds]
    order = FIRST
    family = LAGRANGE
  [../]
  [./elastic_strain11]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./elastic_strain22]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./elastic_strain12]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./unique_grains]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./var_indices]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./vonmises_stress]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./C1111]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./euler_angle]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Kernels]
  [./PolycrystalKernel]
  # ACGrGrPoly<--ACGrGrBase<--ACBulk
  [../]
  [./PolycrystalElasticDrivingForce]
  # PolycrystalElasticDrivingForceAction：Action that addes the elastic driving force for each order parameter
  # ACGrGrElasticDrivingForce，Adds elastic energy contribution to the Allen-Cahn equation
    # D_tensor_name：The elastic tensor derivative ，Stiffness Tensor
    # 多少个序参数,就创建多少个kernel
    # ACBulk<--ACGrGrElasticDrivingForce
  [../]
  [./TensorMechanics]
  # Set up stress divergence kernels with coordinate system aware logic
  # Optional: Add AuxVariables and AuxKernels for various tensor components and quantity outputs
  # 
  # StressDivergenceTensors:Calculate stress divergence equilibrium for the given coordinate system
    use_displaced_mesh = true
    displacements = 'disp_x disp_y'
  [../]
[]

[AuxKernels]
  [./BndsCalc]
    type = BndsCalcAux
    variable = bnds
    execute_on = timestep_end
  [../]
  [./elastic_strain11]
    type = RankTwoAux
    variable = elastic_strain11
    rank_two_tensor = elastic_strain
    index_i = 0
    index_j = 0
    execute_on = timestep_end
  [../]
  [./elastic_strain22]
    type = RankTwoAux
    variable = elastic_strain22
    rank_two_tensor = elastic_strain
    index_i = 1
    index_j = 1
    execute_on = timestep_end
  [../]
  [./elastic_strain12]
    type = RankTwoAux
    variable = elastic_strain12
    rank_two_tensor = elastic_strain
    index_i = 0
    index_j = 1
    execute_on = timestep_end
  [../]
  [./unique_grains]
    type = FeatureFloodCountAux
    variable = unique_grains
    execute_on = timestep_end
    flood_counter = grain_tracker
    field_display = UNIQUE_REGION
  [../]
  [./var_indices]
    type = FeatureFloodCountAux
    variable = var_indices
    execute_on = timestep_end
    flood_counter = grain_tracker
    field_display = VARIABLE_COLORING
  [../]
  [./C1111]
    type = RankFourAux
    variable = C1111
    rank_four_tensor = elasticity_tensor
    index_l = 0
    index_j = 0
    index_k = 0
    index_i = 0
    execute_on = timestep_end
  [../]
  [./vonmises_stress]
    type = RankTwoScalarAux
    variable = vonmises_stress
    rank_two_tensor = stress
    scalar_type = VonMisesStress
    execute_on = timestep_end
  [../]
  [./euler_angle]
    type = OutputEulerAngles
    variable = euler_angle
    euler_angle_provider = euler_angle_file
    grain_tracker = grain_tracker
    output_euler_angle = 'phi1'
    execute_on = 'initial timestep_end'
  [../]
[]

[BCs]
  [./Periodic]
    [./All]
      auto_direction = 'x'
      variable = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7'
    [../]
  [../]
  [./top_displacement]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = top
    function = 'if(t<20,2.5*t,50)'
  [../]
  [./x_anchor]
    type = DirichletBC
    variable = disp_x
    boundary = 'left right'
    value = 0.0
  [../]
  [./y_anchor]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0.0
  [../]
[]

[Materials]
  [./Copper]
    type = GBEvolution
    block = 0
    T = 500 # K
    wGB = 15 # nm
    GBmob0 = 2.5e-6 # m^4/(Js) from Schoenfelder 1997
    Q = 0.23 # Migration energy in eV
    GBenergy = 0.708 # GB energy in J/m^2
  [../]
  [./ElasticityTensor]
    type = ComputePolycrystalElasticityTensor

    grain_tracker = grain_tracker
  [../]
  [./strain]
    type = ComputeSmallStrain
    block = 0
    displacements = 'disp_x disp_y'
  [../]
  [./stress]
    type = ComputeLinearElasticStress
    block = 0
  [../]
[]

[Postprocessors]
  [./ngrains]
    type = FeatureFloodCount
    variable = bnds
    threshold = 0.7
  [../]
  [./dofs]
    type = NumDOFs
  [../]
  [./dt]
    type = TimestepSize
  [../]
  [./run_time]
    type = PerfGraphData
    section_name = "Root"
    data_type = total
  [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    coupled_groups = 'disp_x,disp_y'
  [../]
[]

[Executioner]
  type = Transient
  scheme = bdf2
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -pc_hypre_boomeramg_strong_threshold'
  petsc_options_value = 'hypre boomeramg 31 0.7'
  l_tol = 1.0e-4
  l_max_its = 30
  nl_max_its = 25
  nl_rel_tol = 1.0e-7
  
  start_time = 0.0
  num_steps = 50
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 1.5
    growth_factor = 1.2
    cutback_factor = 0.8
    optimal_iterations = 8
  [../]
  [./Adaptivity]
    initial_adaptivity = 2
    refine_fraction = 0.8
    coarsen_fraction = 0.05
    max_h_level = 3
  [../]
[]

[Outputs]
  file_base = poly36_grtracker
  exodus = true
[]

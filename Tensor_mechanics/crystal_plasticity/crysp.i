[Mesh]
  type = GeneratedMesh
  dim = 3
  elem_type = HEX8
  displacements = 'ux uy uz'
  # The nonlinear displacement variables for the problem
[]

[Variables]
  [./ux]
    block = 0
  [../]
  [./uy]
    block = 0
  [../]
  [./uz]
    block = 0
  [../]
[]

[AuxVariables]
  [./stress_zz]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./fp_zz] # Plastic deformation gradient of previous increment
  # <--materials\FiniteStrainCrystalPlasticity.C
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./rotout]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./e_zz]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./gss1]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
[]

[Functions]
  [./tdisp]
    type = ParsedFunction
    value = 0.01*t
  [../]
[]

[Kernels]
  [./TensorMechanics]
    displacements = 'ux uy uz'
    use_displaced_mesh = true
    # Stress divergence kernel for the Cartesian coordinate system
    # This kernel can be automatically created with the TensorMechanics Master Action.
  [../]
[]

[AuxKernels]
  [./stress_zz]
    type = RankTwoAux
    variable = stress_zz
    rank_two_tensor = stress
    index_j = 2
    index_i = 2
    execute_on = timestep_end
    block = 0
  [../]
  [./fp_zz] # Plastic deformation gradient of previous increment
    type = RankTwoAux
    variable = fp_zz
    rank_two_tensor = fp
    index_j = 2
    index_i = 2
    execute_on = timestep_end
    block = 0
  [../]
  [./e_zz]
    type = RankTwoAux
    variable = e_zz
    rank_two_tensor = lage # Lagrangian strain
    index_j = 2
    index_i = 2
    execute_on = timestep_end
    block = 0
  [../]
  [./gss1]
    type = MaterialStdVectorAux
    variable = gss1 # Slip system resistances of previous increment
    property = gss
    index = 0
    execute_on = timestep_end
    block = 0
  [../]
[]

[BCs]
  [./symmy]
    type = DirichletBC
    variable = uy
    boundary = bottom
    value = 0
  [../]
  [./symmx]
    type = DirichletBC
    variable = ux
    boundary = left
    value = 0
  [../]
  [./symmz]
    type = DirichletBC
    variable = uz
    boundary = back
    value = 0
  [../]
  [./tdisp]
    type = FunctionDirichletBC
    variable = uz
    boundary = front
    function = tdisp
  [../]
[]

[Materials]
  [./crysp]
    # Calculated stress
    type = FiniteStrainCrystalPlasticity
    # /**
    # * FiniteStrainCrystalPlasticity uses the multiplicative decomposition of deformation gradient
    # * and solves the PK2 stress residual equation at the intermediate configuration to evolve the
    # * material state.
    # * The internal variables are updated using an interative predictor-corrector algorithm.
    # * Backward Euler integration rule is used for the rate equations.
    # */
    # power law flow rule

    block = 0
    gtol = 1e-2
    slip_sys_file_name = input_slip_sys.txt
    # 12 slip system for fcc
    nss = 12
    # Number of slip systems
    num_slip_sys_flowrate_props = 2 # Number of properties in a slip system
    # Used for reading flow rate parameters
    flowprops = '1 4 0.001 0.1 5 8 0.001 0.1 9 12 0.001 0.1'
    # 1-4,a0,xm
    # Parameters used in slip rate equations
    hprops = '1.0 541.5 60.8 109.8 2.5'
    # Hardening properties
    gprops = '1 4 60.8 5 8 60.8 9 12 60.8'
    # Initial values of slip system resistances
    tan_mod_type = exact
    # Type of tangent moduli for preconditioner: default elastic
  [../]
  [./elasticity_tensor]
    type = ComputeElasticityTensorCP
    block = 0
    C_ijkl = '1.684e5 1.214e5 1.214e5 1.684e5 1.214e5 1.684e5 0.754e5 0.754e5 0.754e5'
    fill_method = symmetric9
  [../]
  [./strain]
    type = ComputeFiniteStrain
    block = 0
    displacements = 'ux uy uz'
  [../]
[]

[Postprocessors]
  [./stress_zz]
    type = ElementAverageValue
    variable = stress_zz
    block = 'ANY_BLOCK_ID 0'
  [../]
  [./fp_zz]
    type = ElementAverageValue
    variable = fp_zz
    block = 'ANY_BLOCK_ID 0'
  [../]
  [./e_zz]
    type = ElementAverageValue
    variable = e_zz
    block = 'ANY_BLOCK_ID 0'
  [../]
  [./gss1]
    type = ElementAverageValue
    variable = gss1
    block = 'ANY_BLOCK_ID 0'
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient

  #Preconditioned JFNK (default)
  solve_type = 'PJFNK'

  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-10

  dt = 0.05
  dtmax = 10.0
  dtmin = 0.05

  num_steps = 10
[]

[Outputs]
  file_base = out
  exodus = true
[]

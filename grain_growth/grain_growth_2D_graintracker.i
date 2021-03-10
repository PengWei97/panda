# This simulation predicts GB migration of a 2D copper polycrystal with 100 grains represented with 18 order parameters
# Mesh adaptivity and time step adaptivity are used
# An AuxVariable is used to calculate the grain boundary locations
# Postprocessors are used to record time step and the number of grains

[Mesh]
  # Mesh block.  Meshes can be read in or automatically generated
  type = GeneratedMesh
  dim = 2 # Problem dimension
  nx = 11 # Number of elements in the x-direction
  ny = 11 # Number of elements in the y-direction
  xmin = 0    # minimum x-coordinate of the mesh
  xmax = 1000 # maximum x-coordinate of the mesh
  ymin = 0    # minimum y-coordinate of the mesh
  ymax = 1000 # maximum y-coordinate of the mesh
  elem_type = QUAD4  # Type of elements used in the mesh
  uniform_refine = 3 # Initial uniform refinement of the mesh

  parallel_type = replicated # Periodic BCs
[]

[GlobalParams]
  # Parameters used by several kernels that are defined globally to simplify input file
  op_num = 8 # Number of order parameters used
  var_name_base = gr # Base name of grains
[]

[Variables]
  # Variable block, where all variables in the simulation are declared
  [./PolycrystalVariables]
  [../]
[]

[UserObjects]
  [./voronoi]
    type = PolycrystalVoronoi
    # The PolycrystalVoronoi UserObject either generates a set of random points or reads a set of 
    # grain centroids from a file and performs a Voronoi tesslation to produce a grain structure.
    # 2D:8,3D:25
    grain_num = 100 # Number of grains
    rand_seed = 10
    # The random seed
    # file_name = 
    # File containing grain centroids, if file_name is provided, the centroids from the file will be used.
    # int_width = 5.0 # test,Width of diffuse interfaces
  [../]
  [./grain_tracker]
    type = GrainTracker
    # The Grain Tracker is a utility that may be used in phase-field 
    # simulations to reduce the number of order parameters needed to model a large polycrystal system. 
    threshold = 0.2
    connecting_threshold = 0.08
    compute_halo_maps = true # Only necessary for displaying HALOS
  [../]
[]

[ICs]
  [./PolycrystalICs]
    [./PolycrystalColoringIC]
      polycrystal_ic_uo = voronoi
      # <--UserObjects\voronoi
    [../]
  [../]
[]

[AuxVariables]
  # Dependent variables
  [./bnds]
    # Variable used to visualize the grain boundaries in the simulation
  [../]
  [./unique_grains]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./var_indices]
    order = CONSTANT
    family = MONOMIAL
  [../]
  # [./ghost_regions]
  #   order = CONSTANT
  #   family = MONOMIAL
  # [../]
  # [./halos]
  #   order = CONSTANT
  #   family = MONOMIAL
  # [../]
[]

[Kernels]
  # Kernel block, where the kernels defining the residual equations are set up.
  [./PolycrystalKernel]
    # Custom action creating all necessary kernels for grain growth.  All input parameters are up in GlobalParams
    # <--PolycrystalKernelAction,Set up ACGrGrPoly, ACInterface, TimeDerivative, and ACGBPoly kernels
    # ACGrGrPoly,see PolycrystalKernelAction.JPG
    # ACInterface
    # TimeDerivative
    # ACGBPoly,'ACGBPoly' is used only when there are bubbles in the grain-growth simulations. If there are no bubbles, this kernel is not set by the action.
  [../]
[]

[AuxKernels]
  # AuxKernel block, defining the equations used to calculate the auxvars
  [./bnds_aux]
    # AuxKernel that calculates the GB term
    type = BndsCalcAux
    variable = bnds
    execute_on = 'initial timestep_end'
  [../]
  [./unique_grains]
    type = FeatureFloodCountAux
    # Feature detection by connectivity analysis
    variable = unique_grains
    # <--AuxVariables\unique_grains
    flood_counter = grain_tracker
    # <--UserObjects\grain_tracker
    field_display = UNIQUE_REGION
    # UNIQUE_REGION VARIABLE_COLORING GHOSTED_ENTITIES HALOS CENTROID ACTIVE_BOUNDS
    execute_on = 'initial timestep_end'
  [../]
  [./var_indices]
    type = FeatureFloodCountAux
    variable = var_indices
    flood_counter = grain_tracker
    field_display = VARIABLE_COLORING
    execute_on = 'initial timestep_end'
  [../]
  # [./ghosted_entities]
  #   type = FeatureFloodCountAux
  #   variable = ghost_regions
  #   flood_counter = grain_tracker
  #   field_display = GHOSTED_ENTITIES
  #   execute_on = 'initial timestep_end'
  # [../]
  # [./halos]
  #   type = FeatureFloodCountAux
  #   variable = halos
  #   flood_counter = grain_tracker
  #   field_display = HALOS
  #   execute_on = 'initial timestep_end'
  # [../]
[]

[BCs]
  # Boundary Condition block
  [./Periodic]
    [./top_bottom]
      auto_direction = 'x y' # Makes problem periodic in the x and y directions
    [../]
  [../]
[]

[Materials]
  [./CuGrGr]
    # Material properties
    type = GBEvolution
    T = 450 # Constant temperature of the simulation (for mobility calculation)
    # GBEvolutionBase:Computes necessary material properties for the isotropic grain growth model
    wGB = 14 # Width of the diffuse GB
    GBmob0 = 2.5e-6 # m^4(Js) for copper from Schoenfelder1997,Grain boundary mobility prefactor
    Q = 0.23 #eV for copper from Schoenfelder1997,Grain boundary migration activation energy
    GBenergy = 0.708 #J/m^2 from Schoenfelder1997
    # length_scale = 1.0e-9 # Length scale in m, where default is nm
    # time_scale = 1.0e-9 # Time scale in s, where default is ns"
    # GBMobility = -1 # GB mobility input in m^4/(J*s), that overrides the temperature dependent calculation
    # molar_volume = 24.62e-6 # m^3/mol, needed for temperature gradient driving force
  [../]
[]

[Postprocessors]
  # Scalar postprocessors
  [./dt]
    # Outputs the current time step
    type = TimestepSize
  [../]
[]

[Executioner]
  type = Transient # Type of executioner, here it is transient with an adaptive time step
  scheme = bdf2 # Type of time integration (2nd order backward euler), defaults to 1st order backward euler

  #Preconditioned JFNK (default)
  solve_type = 'PJFNK'

  # Uses newton iteration to solve the problem.
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -mat_mffd_type'
  petsc_options_value = 'hypre boomeramg 101 ds'

  l_max_its = 30 # Max number of linear iterations
  l_tol = 1e-4 # Relative tolerance for linear solves
  nl_max_its = 40 # Max number of nonlinear iterations
  nl_rel_tol = 1e-10 # Absolute tolerance for nonlienar solves

  start_time = 0.0
  end_time = 4000

  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 25 # Initial time step.  In this simulation it changes.
    optimal_iterations = 6 # Time step will adapt to maintain this number of nonlinear iterations
  [../]

  [./Adaptivity]
    # Block that turns on mesh adaptivity. Note that mesh will never coarsen beyond initial mesh (before uniform refinement)
    initial_adaptivity = 2 # Number of times mesh is adapted to initial condition
    refine_fraction = 0.7 # Fraction of high error that will be refined
    coarsen_fraction = 0.1 # Fraction of low error that will coarsened
    max_h_level = 4 # Max number of refinements used, starting from initial mesh (before uniform refinement)
  [../]
[]

[Outputs]
  file_base = poly1000_graintracker
  [./csv]
    type = CSV
    interval = 8
  [../]
  [./exodus]
    type = Exodus
    interval = 16
  [../]
  [./console]
    type = Console
    max_rows = 20 # Will print the 20 most recent postprocessor values to the screen
  [../]
  [out]
    type = Checkpoint
    interval = 10
    num_files = 6
  []
  [pgraph]
    type = PerfGraphOutput
    execute_on = 'initial final'  # Default is "final"
    level = 2                     # Default is 1
    heaviest_branch = true        # Default is false
    heaviest_sections = 7         # Default is 0
  []
[]

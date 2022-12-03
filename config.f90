module config

    integer, parameter :: N = 4000 ! number of bodies
    integer, parameter :: steps = 10 ! simulation time in steps
    real, parameter    :: dt = 0.1 ! time step in seconds
    logical            :: new_initial_conditions = .true. ! set to .true. to generate new initial conditions
    logical            :: verbose = .false.
    logical            :: save_results = .true. ! set to .false. to disable saving results (for bulk runs)

end module config
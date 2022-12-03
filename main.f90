program nbody

    implicit none

    ! Safe to edit
    integer, parameter :: N = 1000 ! number of bodies
    integer, parameter :: steps = 100 ! simulation time in steps
    real, parameter    :: dt = 0.1 ! time step in seconds
    logical            :: new_initial_conditions = .true. ! set to .true. to generate new initial conditions
    logical            :: verbose = .false.
    logical            :: save_results = .true. ! set to .false. to disable saving results (for bulk runs)

    ! Don't touch to that, you fool!
    ! integer                    :: MAX_TENTATIVES = 1000
    ! INTEGER                    :: tentatives = 0
    INTEGER                    :: i, t, tmp = -1
    REAL                       :: PI
    REAL                       :: x, y, z
    real, dimension(3,N,steps) :: p=0, v=0, a=0
    real                       :: m
    real, dimension(steps)     :: Ep=0, Ec=0, Et=0
    real, dimension(3,steps)   :: b ! barycenter

    PI = acos(-1.0)

    ! print *, "NBody simulation..."

    !--------------------
    ! Initial conditions
    !--------------------

    m = 1./N

    ! Position
    if (new_initial_conditions .eqv. .false.) then
        if (verbose .eqv. .true.) then
            print *, "Getting initial conditions"
        end if
        open(42, file='data/initial_positions.txt', status='old')
        do i= 1, N
            read(42, *) p(:, i, 1)
        end do
        close(42)
    
    else
        do i= 1, N

            if (verbose .eqv. .true.) then
                if (i * 100 / N .ne. tmp) then
                    tmp = i * 100 / N
                    print *, "Generating initial conditions ", tmp, "%"
                end if
            end if

            ! generating random coordinate
            call random_coord_in_sphere(x,y,z)

            ! Storing the data
            p(1, i, 1) = x
            p(2, i, 1) = y
            p(3, i, 1) = z

        end do

        if (verbose .eqv. .true.) then
            print *, "Saving initial conditions"
        end if

        open(42, file='data/initial_positions.txt', status='replace')
        do i= 1, N
            write(42, *) p(:, i, 1)
        end do
        close(42)
    end if

    ! Velocities
    
    do i=1,N
        ! Highly complex computation of object velocity
        ! Assumption to compute it: angular velocity around z axes with constant value = 1
        v(1, i, 1) = -y
        v(2, i, 1) = x
        v(3, i, 1) = 0
        
    end do

    ! Accelerations

    ! Computing initial acceleration
    call acceleration(p(:, :, 1), m, N, a(:, :, 1), Ep(1))

    ! Barycenter

    call barycenter(p(:,:,1), m, N, b(:,1))

    ! Energies
    
    ! do i=1,N
    !     Ec(1) = Ec(1) + 0.5 * m * (v(1, i, 1) ** 2 + v(2, i, 1) ** 2 + v(3, i, 1) ** 2)
    ! end do
    ! Et(1) = Ep(1) + Ec(1)

    !--------------------
    ! Time evolution
    !--------------------

    do t=1,steps-1
    
        ! ! Print progress
        if (verbose .eqv. .true.) then
            if (t * 100 / steps .ne. tmp) then
                tmp = t * 100 / steps
                print *, "Computing evolution ", tmp, "%"
            end if
        end if

        ! Compute particle position, velocity and acceleration at time t+1
        call next_state(p, v, a, t, m, dt, N, steps, Ep(t+1), Ec(t+1), Et(t+1))
        call barycenter(p(:,:,t+1), m, N, b(:,t+1))

    end do

    !--------------------
    ! Export results
    !--------------------

    if (save_results .eqv. .true.) then

        ! Writing bodies position, velocity and acceleration
        if (verbose .eqv. .true.) then
            print *, "Saving bodies data"
        end if
        open(42, file='data/nbody.txt', status='replace')
        write(42, *) "# p_x     p_y     p_z     v_x     v_y     v_z     a_x     a_y     a_z"
        do t = 1, steps
            write(42, *) " "
            write(42, *) "# t = ", t*dt, " step = ", t
            do i = 1, N
                write(42, *) p(:,i,t), v(:,i,t), a(:,i,t)
            end do
        end do
        close(42)
        
        ! Writing system energies
        if (verbose .eqv. .true.) then
            print *, "Saving energies"
        end if
        open(42, file='data/energies.txt', status='replace')
        write(42, *) "# t     Ep     Ec     Et"
        do t = 1, steps
            write(42, *) t*dt, Ep(t), Ec(t), Et(t)
        end do
        close(42)

        ! Writing barycenter
        if (verbose .eqv. .true.) then
            print *, "Saving energies"
        end if
        open(42, file='data/barycenter.txt', status='replace')
        write(42, *) "# t, barycenter"
        do t = 1, steps
            write(42, *) t*dt, b(:,t)
        end do
        close(42)
        
        ! Writing simulationn parameters
        if (verbose .eqv. .true.) then
            print *, "Saving simulation parameters"
        end if
        open(42, file='data/parameters.txt', status='replace')
        write(42, *) "# N     steps     dt"
        write(42, *) N, steps, dt
        close(42)
    
    end if

    ! print *, "Done!"

end program nbody
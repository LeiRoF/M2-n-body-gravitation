program nbody

    ! use settings

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
    real, dimension(N)         :: m
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
    Ec(1) = Ec(1) + 0.5 * m(i) * (v(1, i, 1) ** 2 + v(2, i, 1) ** 2 + v(3, i, 1) ** 2)
    Et(1) = Ep(1) + Ec(1)

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





    ! ====================================================================================================
    ! SUBROUTINES & FUNCTIONS
    ! ====================================================================================================




    

    contains

        ! ----------------------------------------------------------------------------------------------------
        ! generate a random cartesian coordinate in a sphere of radius 1

        subroutine random_coord_in_sphere(x,y,z)
            real, intent(  out) :: x, y, z
            logical             :: lock

            lock = .true.
            do while (lock .eqv. .true.)
                CALL RANDOM_NUMBER(x)
                x = x * 2.0 - 1.0
                CALL RANDOM_NUMBER(y)
                y = y * 2.0 - 1.0
                CALL RANDOM_NUMBER(z)
                z = z * 2.0 - 1.0

                if (x**2 + y**2 + z**2 < 1) then
                    lock = .false.
                end if
            end do

        end subroutine random_coord_in_sphere

        !----------------------------------------------------------------------------------------------------
        ! Computing next position, speed and acceleration using leap frog algorithm

        subroutine next_state(p, v, a, t, m, dt, N, steps, Ep, Ec, Et)
            real,    intent(inout), dimension(3, N, steps) :: p, v, a ! time t and t+1
            integer, intent(in   )                         :: t ! step
            real,    intent(in   ), dimension(N)           :: m 
            integer, intent(in   )                         :: N, steps
            real,    intent(in   )                         :: dt
            real,    intent(  out)                         :: Ep, Ec, Et
            real,                   dimension(3,N)         :: pi,vi,ai ! initial state at t
            real,                   dimension(3,N)         :: p2, a2   ! intermediate state at t+1/2 
            real,                   dimension(3,N)         :: pf,vf,af ! final state at t+1
            integer                                        :: i
            ! integer                                        :: omp_get_num_threads

            pi = p(:,:,t)
            vi = v(:,:,t)   
            ai = a(:,:,t)
            pf = 0
            vf = 0
            af = 0
            p2 = 0
            a2 = 0
            
            do i=1,N

                ! Compute p at t+1/2
                p2(:,i) = pi(:,i) + vi(:,i) * dt/2.

            end do

            ! Compute a at t+1/2
            call acceleration(p2(:, :), m, N, a2(:,:))

            do i=1,N

                ! Compute velocity at t+1
                vf(:, i) = vi(:, i) + ai(:, i) * dt/2 + a2(:, i) * dt/2

                Ec = Ec + 0.5 * m(i) * (vf(1, i)**2 + vf(2, i)**2 + vf(3, i)**2)

                ! Position at t+1
                pf(:, i) = p2(:, i) + vi(:, i) * dt/2

            end do

            ! Compute a at i+1
            call acceleration(pf(:, :), m, N, af(:,:), Ep)

            ! omp end do
            ! omp end parallel

            p(:, :, t+1) = pf
            v(:, :, t+1) = vf
            a(:, :, t+1) = af
            Et = Ep + Ec

        end subroutine next_state

        !----------------------------------------------------------------------------------------------------
        ! Compute the acceleration of a body considering the position of the others

        subroutine acceleration(p, m, N, a, Ep)
            real,    intent(in   ), dimension(3, N) :: p
            real,    intent(in   ), dimension(N)    :: m
            integer, intent(in   )                  :: N
            real,    intent(  out), dimension(3,N)  :: a
            real,    intent(  out), optional        :: Ep
            real                                    :: x, y, z
            real                                    :: G = 1, eps = 0.05
            real                                    :: r
            integer                                 :: i, j

            ! Compute the acceleration of the i-th body
            !private(dist, r, eps, G, m) reduction(+:a)
            ! schedule(dynamic,N)
            !$omp parallel reduction(+:Ep)
            !$omp do
            do i=1,N
                do j=1,N

                    if (j==i) then
                        cycle
                    end if
                    
                    x = p(1,j) - p(1,i)
                    y = p(2,j) - p(2,i)
                    z = p(3,j) - p(3,i)

                    r = sqrt(x*x + y*y + z*z + eps*eps)

                    a(1, i) = a(1, i) + x * G * m(j) / r**3
                    a(2, i) = a(2, i) + y * G * m(j) / r**3
                    a(3, i) = a(3, i) + z * G * m(j) / r**3

                    Ep = Ep - 0.5 * G * m(j)**2 / r
                end do
            end do
            !$omp end do
            !$omp end parallel

        end subroutine acceleration

        !----------------------------------------------------------------------------------------------------
        ! Compute the energy of the system

        ! function potential_energy(p, m, N)
        
        !     real,    intent(in   ), dimension(3, N) :: p
        !     real,    intent(in   ), dimension(N)    :: m
        !     integer, intent(in   )                  :: N
        !     real,                   dimension(3)    :: a, b, dist
        !     real                                    :: potential_energy
        !     real                                    :: G=1, eps = 0.05
        !     real                                    :: r
        !     integer                                 :: i, j, k

        !     potential_energy = 0

        !     !$omp parallel
        !     !$omp do
        !     do i=1,N
        !         do j=1,N
        !             if (j==i) then
        !                 cycle
        !             end if

        !             dist(:) = p(:,j) - p(:,j)

        !             r = sqrt(sum(dist**2) + eps**2)

        !         end do
        !     end do
        !     !$omp parallel
        !     !$omp do

        ! end function potential_energy

        !----------------------------------------------------------------------------------------------------
        ! Compute the kinetic energy of the system

        ! function kinetic_energy(v, m, N)
        
        !     real,    intent(in   ), dimension(3, N) :: v
        !     real,    intent(in   ), dimension(N)    :: m
        !     integer, intent(in   )                  :: N
        !     real                                    :: kinetic_energy
        !     integer                                 :: i

        !     kinetic_energy = 0

        !     ! Ec = sum on i of 1/2 * m_i * v_i^2
        !     do i=1,N
        !         kinetic_energy = kinetic_energy + 0.5 * m(i) * sum(v(:, i)**2)
        !     end do

        ! end function kinetic_energy

        !----------------------------------------------------------------------------------------------------
        ! Get the barycenter

        subroutine barycenter(p, m, N, b)
        
            real,    intent(in   ), dimension(3, N) :: p
            real,    intent(in   ), dimension(N)    :: m
            integer, intent(in   )                  :: N
            real,    intent(  out), dimension(3)    :: b
            integer                                 :: i

            b = 0

            ! Ec = sum on i of 1/2 * m_i * v_i^2
            do i=1,N
                b = b + m(i) * p(:, i)
            end do

            b = b / sum(m)

        end subroutine barycenter

end program nbody
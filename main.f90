program nbody
    implicit none

    ! Safe to edit
    INTEGER, PARAMETER         :: N = 1000 ! number of bodies
    integer, parameter         :: steps = 1000 ! simulation time in steps
    real, parameter            :: dt = 0.1 ! time step in seconds
    logical, parameter         :: new_initial_conditions = .false. ! set to .true. to generate new initial conditions

    ! Don't touch to that, you fool!
    integer                    :: MAX_TENTATIVES = 1000
    INTEGER                    :: tentatives = 0
    INTEGER                    :: i,j, t, tmp = -1
    REAL                       :: PI, tmp2, dist
    REAL                       :: x, y, z
    real                       :: vx, vy, vz
    real, dimension(steps,N,3) :: p=0, v=0, a=0
    real, dimension(N)         :: m
    real, dimension(steps)     :: Ep, Ec, Et
    real, dimension(steps,3)   :: b ! barycenter

    PI = 4*ATAN(1.)

    !--------------------
    ! Initial conditions
    !--------------------

    m = 1./N

    ! Position

    if (new_initial_conditions .eqv. .false.) then
        print *, "Getting initial conditions"
        open(42, file='initial_positions.txt', status='old')
        do i= 1, N
            read(42, *) p(1, i, :)
        end do
        close(42)
    else
        do i= 1, N

            if (i * 100 / N .ne. tmp) then
                tmp = i * 100 / N
                print *, "Generating initial conditions ", tmp, "%"
            end if

            ! generating random coordinate
            call random_coord_in_sphere(x,y,z)

            if (i > 1) then

                ! if the coordinate is to near from another one, retry
                do while ((is_too_near_from_another(x, y, z, p(1, :, 1), p(1, :, 2), p(1, :, 3), i) .eqv. .true.) &
                .and. (tentatives < MAX_TENTATIVES))

                    tentatives = tentatives + 1
                    call random_coord_in_sphere(x,y,z)

                end do
            end if

            ! Storing the data
            p(1, i, 1) = x
            p(1, i, 2) = y
            p(1, i, 3) = z

        end do
        
        print *, "Saving initial conditions"
        open(42, file='initial_positions.txt', status='replace')
        do i= 1, N
            write(42, *) p(1, i, :)
        end do
        close(42)
    end if

    ! Velocities
    
    do i=1,N
        ! Highly complex computation of object velocity
        ! Assumption to compute it: angular velocity around z axes with constant value = 1
        v(1, i, 1) = -y
        v(1, i, 2) = x
        v(1, i, 3) = 0
    end do

    ! Accelerations

    do i = 1,N
        ! Computing initial acceleration
        call acceleration(p(1, :, :), m, i, N, a(1, i, :))
    end do

    ! Barycenter

    call barycenter(p(1,:,:), m, N, b(1,:))

    ! Energies

    Ep(1) = potential_energy(p(1,:,:), m, N)
    Ec(1) = kinetic_energy(v(1,:,:), m, N)
    Et(1) = Ep(1) + Ec(1)

    !--------------------
    ! Time evolution
    !--------------------

    do t=1,steps-1
    
        ! Print progress
        if (t * 100 / steps .ne. tmp) then
            tmp = t * 100 / steps
            print *, "Computing evolution ", tmp, "%"
        end if

        ! Compute particle position, velocity and acceleration at time t+1
        call next_state(p, v, a, t, m, dt, N, steps)
        call barycenter(p(t+1,:,:), m, N, b(t+1,:))

        ! Compute energies at time t+1
        Ep(t+1) = potential_energy(p(t+1,:,:), m, N)
        Ec(t+1) = kinetic_energy(v(t+1,:,:), m, N)
        Et(t+1) = Ep(t+1) + Ec(t+1)
    end do

    !--------------------
    ! Export results
    !--------------------

    ! Writing bodies position, velocity and acceleration
    print *, "Saving bodies data"
    open(42, file='nbody.txt', status='replace')
    write(42, *) "# p_x     p_y     p_z     v_x     v_y     v_z     a_x     a_y     a_z"
    do t = 1, steps
        write(42, *) " "
        write(42, *) "# t = ", t*dt, " step = ", t
        do i = 1, N
            write(42, *) p(t,i,:), v(t,i,:), a(t,i,:)
        end do
    end do
    close(42)
    
    ! Writing system energies
    print *, "Saving energies"
    open(42, file='energies.txt', status='replace')
    write(42, *) "# t     Ep     Ec     Et"
    do t = 1, steps
        write(42, *) t*dt, Ep(t), Ec(t), Et(t)
    end do
    close(42)

    ! Writing barycenter
    print *, "Saving energies"
    open(42, file='barycenter.txt', status='replace')
    write(42, *) "# t, barycenter"
    do t = 1, steps
        write(42, *) t*dt, b(t,:)
    end do
    close(42)
    
    ! Writing simulationn parameters
    print *, "Saving simulation parameters"
    open(42, file='parameters.txt', status='replace')
    write(42, *) "# N     steps     dt"
    write(42, *) N, steps, dt
    close(42)






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

        ! ----------------------------------------------------------------------------------------------------
        ! Function that ensure the new object will not appear too near from another one

        function is_too_near_from_another(x, y, z, x_L, y_L, z_L, N)
            integer, intent(in   )               :: N
            real,    intent(in   )               :: x,y,z
            real,    intent(in   ), dimension(N) :: x_L, y_L, z_L
            logical                              :: is_too_near_from_another
            real                                 :: volume, distance, eps
            real                                 :: PI = 3.14159

            is_too_near_from_another = .false.
            return

            volume = 4./3. * PI

            eps = volume / N ! get the volume occupied by one object
            eps = eps**(1./3.) ! get the size of the volume occupied by object
            eps = eps / 5 ! devide it by an arbitrary number
            ! here 5 in order that local density cannot exceed 5 times homogenous density

            ! looking at all existing objects
            do i= 1,N

                ! distance between the two objects
                distance = sqrt((x-x_L(i))**2 + (y-y_L(i))**2 + (z-z_L(i))**2)  

                ! if eps > distance > -eps -> reject (return true)
                if (distance < eps .and. distance > -eps) then
                    is_too_near_from_another = .true.
                    return
                end if
            end do
            is_too_near_from_another = .false.

        end function is_too_near_from_another

        !----------------------------------------------------------------------------------------------------
        ! Computing next position, speed and acceleration using leap frog algorithm

        subroutine next_state(p, v, a, t, m, dt, N, steps)
            real,    intent(inout), dimension(steps, N, 3) :: p, v, a ! time t and t+1
            integer, intent(in   )                         :: t ! step
            real,    intent(in   ), dimension(N)           :: m 
            integer, intent(in   )                         :: N, steps
            real,    intent(in   )                         :: dt
            real,                   dimension(N,3)         :: tmp_p, tmp_v, tmp_a ! time t+1/2
            real,                   dimension(3)           :: dist
            integer,                dimension(3)           :: sign
            integer                                        :: i, j, k
            
            ! ! Compute p at t+1/2
            ! !OMP PARALLEL DO SHARE(tmp_p)
            ! do i=1,N
            !     tmp_p(i,:) = p(t, i,:) + v(t, i,:) * dt/2.
            ! end do
            ! !OMP END PARALLEL DO

            ! ! Compute a at t+1/2
            ! !OMP PARALLEL DO private(dist,sign,j,m,tmp_p) shared(tmp_a)
            ! do i=1,N
            !     call acceleration(tmp_p(:, :), m, i, N, tmp_a(i, :))
            ! end do
            ! !OMP END PARALLEL DO

            ! Compute velocity at t+1
            !OMP PARALLEL DO SHARE(v) PRIVATE(a, dt)
            do i=1,N
                v(t+1, i, :) = v(t, i, :) + a(t, i, :) * dt
            end do
            !OMP END PARALLEL DO

            ! Position at t+1
            !OMP PARALLEL DO SHARE(p) PRIVATE(v, dt)
            do i=1,N
                p(t+1, i, :) = p(t,i, :) + v(t, i, :) * dt
            end do
            !OMP END PARALLEL DO

            ! Compute a at i+1
            !OMP PARALLEL DO private(dist,sign,j,m, p) shared(a)
            do i=1,N
                call acceleration(p(t+1, :, :), m, i, N, a(t+1, i, :))
            end do
            !OMP END PARALLEL DO

        end subroutine next_state

        !----------------------------------------------------------------------------------------------------
        ! Compute the acceleration of a body considering the position of the others

        subroutine acceleration(p, m, i, N, a)
            real,    intent(in   ), dimension(N, 3) :: p
            real,    intent(in   ), dimension(N)    :: m
            integer, intent(in   )                  :: i, N
            real,    intent(  out), dimension(3)    :: a
            real,                   dimension(3)    :: dist
            integer,                dimension(3)    :: s
            real                                    :: G = 1, eps = 0.05, r
            integer                                 :: j, k

            ! print *, "--------------------------------------------------"
            ! print *, "p_i = ", i, " | ", p(i,:)
            ! print *, " "

            a = 0

            ! Compute a at i+1
            do j=1,N
                if (j==i) then
                    cycle
                end if
                
                dist(:) = p(j,:) - p(i,:)

                r = sqrt(sum(dist**2) + eps**2)

                a(:) = a(:) + dist(:) * G * m(j) / r**3

                ! Compute the potential energy
                ! ep = ep - 0.5*gg*m**2/r
            end do

        end subroutine acceleration

        !----------------------------------------------------------------------------------------------------
        ! Compute the energy of the system

        function potential_energy(p, m, N)
        
            real,    intent(in   ), dimension(N, 3) :: p
            real,    intent(in   ), dimension(N)    :: m
            integer, intent(in   )                  :: N
            real,                   dimension(3)    :: a, b, dist
            real                                    :: potential_energy, G=1, eps = 0.05
            real                                    :: r
            integer                                 :: i, j, k

            potential_energy = 0

            do i=1,N
                do j=1,N
                    if (j==i) then
                        cycle
                    end if

                    dist(:) = p(j,:) - p(i,:)

                    r = sqrt(sum(dist**2) + eps**2)

                    potential_energy = potential_energy - 0.5 * G * m(j)**2 / r
                end do
            end do

        end function potential_energy

        !----------------------------------------------------------------------------------------------------
        ! Compute the kinetic energy of the system

        function kinetic_energy(v, m, N)
        
            real,    intent(in   ), dimension(N, 3) :: v
            real,    intent(in   ), dimension(N)    :: m
            integer, intent(in   )                  :: N
            real                                    :: kinetic_energy
            integer                                 :: i

            kinetic_energy = 0

            ! Ec = sum on i of 1/2 * m_i * v_i^2
            do i=1,N
                kinetic_energy = kinetic_energy + 0.5 * m(i) * sum(v(i, :)**2)
            end do

        end function kinetic_energy

        !----------------------------------------------------------------------------------------------------
        ! Get the barycenter

        subroutine barycenter(p, m, N, b)
        
            real,    intent(in   ), dimension(N, 3) :: p
            real,    intent(in   ), dimension(N)    :: m
            integer, intent(in   )                  :: N
            real,    intent(  out), dimension(3)    :: b
            integer                                 :: i

            b = 0

            ! Ec = sum on i of 1/2 * m_i * v_i^2
            do i=1,N
                b = b + m(i) * p(i, :)
            end do

            b = b / sum(m)

        end subroutine barycenter

end program nbody
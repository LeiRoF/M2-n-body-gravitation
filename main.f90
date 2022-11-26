program nbody
    implicit none

    ! Safe to edit
    INTEGER, PARAMETER         :: N = 100 ! number of bodies
    integer, parameter         :: steps = 1000 ! simulation time in steps
    real, parameter            :: dt = 0.001 ! time step in seconds

    ! Don't touch to that, you fool!
    integer                    :: MAX_TENTATIVES = 1000
    INTEGER                    :: tentatives = 0
    INTEGER                    :: i, t, tmp = -1
    REAL                       :: PI
    REAL                       :: x, y, z
    real                       :: vx, vy, vz
    real, dimension(steps,N,3) :: p, v, a
    real, dimension(N)         :: m

    PI = 4*ATAN(1.)

    !--------------------
    ! Initial conditions
    !--------------------

    m = 1

    ! generating n objects
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

        ! storing the data
        p(1, i, 1) = x  
        p(1, i, 2) = y
        p(1, i, 3) = z

        ! highly complex computation of object velocity
        ! assumption to compute it: angular velocity around z axes with constant value = 1
        v(1, i, 1) = -y
        v(1, i, 2) = x
        v(1, i, 3) = 0

        ! setting default acceleration to 0
        a(1, i, 1) = 0
        a(1, i, 2) = 0
        a(1, i, 3) = 0

    end do

    !--------------------
    ! Time evolution
    !--------------------

    do t=1,steps-1
    
        if (t * 100 / steps .ne. tmp) then
            tmp = t * 100 / steps
            print *, "Computing evolution ", tmp, "%"
        end if

        call next_state(p, v, a, t, m, dt, N, steps)
    end do

    !--------------------
    ! Export results
    !--------------------

    open(42, file='nbody.txt', status='replace')

    write(42, *) "# p_x     p_y     p_z     v_x     v_y     v_z     a_x     a_y     a_z"

    do t = 1, steps
        write(42, *) " "
        write(42, *) "# t = ", t*dt

        do i = 1, N
            write(42, *) p(t,i,:), v(t,i,:), a(t,i,:)
        end do
    end do

    close(42)

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
        ! Subroutine to convert from spherical to cartesian coordinates

        subroutine spherical_to_cartesian(r, phi, theta, x, y, z)
            real, intent(in   ) :: r, phi, theta
            real, intent(  out) :: x, y, z

            x = r * SIN(phi) * COS(theta)
            y = r * SIN(phi) * SIN(theta)
            z = r * COS(phi)

        end subroutine spherical_to_cartesian

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
            real,    intent(inout), dimension(steps, N, 3) :: p, v, a ! time t
            integer, intent(in   )                         :: t ! step
            real,    intent(in   ), dimension(N)           :: m 
            integer, intent(in   )                         :: N, steps
            real,    intent(in   )                         :: dt
            real,                   dimension(N,3)         :: tmp_p, tmp_v, tmp_a ! time t+1/2
            real,                   dimension(3)           :: dist
            integer,                dimension(3)           :: sign
            integer                                        :: i, j, k

            ! print *, " "
            ! print *, " "
            ! print *, " "
            ! print *, "--------------------------------------------------"

            ! print *, "dt = ", dt, "dt/2 = ", dt/2.
            ! print *, "m = ", m(:)

            ! print *, " "
            ! print *, "p(t) = "
            ! do i=1,n
            !     print *, p(t, i, :)
            ! end do
            ! print *, "v(t) = "
            ! do i=1,n
            !     print *, v(t, i, :)
            ! end do

            
            !fromage SHARE(tmp_p)

            ! Compute p at i + 1/2
            do i=1,N
                tmp_p(i,:) = p(t, i,:) + v(t, i,:) * dt/2.
            end do

            !OMP END PARALLEL DO

            ! print *, " "
            ! print *, "p(t+1/2) = "
            ! do i=1,n
            !     print *, tmp_p(i, :)
            ! end do

            !fromage private(dist,sign,j,m,tmp_p) shared(tmp_a)

            ! Compute a at i+1/2
            do i=1,N
                call acceleration(tmp_p(:, :), m, i, tmp_a(i, :), N)
            end do

            !OMP END PARALLEL DO

            ! print *, " "
            ! print *, "a(t+1/2) = "
            ! do i=1,n
            !     print *, tmp_a(i, :)
            ! end do

            !fromage SHARE(v)

            ! Compute v at i+1
            do i=1,N
                v(t+1, i, :) = v(t, i, :) + tmp_a(i, :) * dt
            end do

            !OMP END PARALLEL DO

            ! print *, " "
            ! print *, "v(t+1) = "
            ! do i=1,n
            !     print *, v(t+1, i, :)
            ! end do

            !fromage SHARE(p)

            ! Compute p at i+1
            do i=1,N
                p(t+1, i, :) = tmp_p(i, :) + tmp_v(i, :) * dt/2.
            end do

            !OMP END PARALLEL DO

            ! print *, " "
            ! print *, "p(t+1) = "
            ! do i=1,n
            !     print *, p(t+1, i, :)
            ! end do

            !fromage private(dist,sign,j,m, p) shared(a)

            ! Compute a at i+1
            do i=1,N
                call acceleration(p(t+1, :, :), m, i, a(t+1, i, :), N)
            end do

            !OMP END PARALLEL DO

            ! print *, " "
            ! print *, "a(t+1) = "
            ! do i=1,n
            !     print *, a(t+1, i, :)
            ! end do

        end subroutine next_state

        !----------------------------------------------------------------------------------------------------
        ! Compute the acceleration of a body considering the position of the others

        subroutine acceleration(p, m, i, a, N)
            real,    intent(in   ), dimension(N, 3) :: p
            real,    intent(in   ), dimension(N)    :: m
            integer, intent(in   )                  :: i, N
            real,    intent(  out), dimension(3)    :: a
            real,                   dimension(3)    :: dist
            integer,                dimension(3)    :: s
            real                                    :: G = 1.0, eps = 0.05
            integer                                 :: j, k

            ! print *, "--------------------------------------------------"
            ! print *, "p_i = ", i, " | ", p(i,:)
            ! print *, " "

            a(:) = 0

            ! Compute a at i+1
            do j=1,N
                if (i .ne. j) then
                    ! print *, "p_j = ", j, " | ", p(j,:)

                    ! Computing distance between i and j
                    dist = p(j, :) - p(i, :)

                    ! print *, "d_j = ", j, " | ", dist(:)

                    s = 1
                    do k=1,3
                        if (dist(k) < 0) then
                            s(k) = -1
                        end if
                    end do

                    ! print *, "s_j = ", j, " | ", s(:)
                    ! print *, "a_j = ", j, " | ", s(:) * G * m(j) / (dist(:) * dist(:) + eps)
                    ! print *, " "

                    ! Computing acceleration
                    a(:) = a(:) + s(:) * G * m(j) / (dist(:) * dist(:) + eps)                    

                end if
            end do
            
            ! print *, "a = ", a(:)

        end subroutine acceleration

end program nbody
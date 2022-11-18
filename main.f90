program nbody
    implicit none

    INTEGER, PARAMETER :: N = 1000, MAX_TENTATIVES = 1000
    INTEGER            :: tentatives = 0
    INTEGER            :: i
    REAL               :: PI = 3.14159
    REAL               :: x, y, z, r, phi, theta
    REAL, DIMENSION(N) :: x_L, y_L, z_L, r_L, phi_L, theta_L
    real               :: vx, vy, vz, vr, vphi, vtheta
    real, DIMENSION(N) :: vx_L, vy_L, vz_L, vr_L, vphi_L, vtheta_L

    ! generating n objects
    do i= 1, N

        ! generating random coordinate
        call random_coord_in_sphere(x,y,z)

        if (i > 1) then

            ! if the coordinate is to near from another one, retry
            do while ((is_too_near_from_another(x, y, z, x_L, y_L, z_L, i) .eqv. .true.) .and. (tentatives < MAX_TENTATIVES))

                tentatives = tentatives + 1
                call random_coord_in_sphere(x,y,z)

            end do
        end if

        ! highly complex computation of object velocity
        ! assumption to compute it: angular velocity around z axes with constant value = 1

        vx = -y
        vy = x
        vz = 0

        ! storing the data

        x_L(i) = x  
        y_L(i) = y
        z_L(i) = z

        vx_L(i) = vx
        vy_L(i) = vy
        vz_L(i) = vz
    end do

    ! saving the data in a file
    open(42, file='nbody.txt', status='replace')
    do i= 1, N
        write(42, *) x_L(i), y_L(i), z_L(i), vx_L(i), vy_L(i), vz_L(i)
    end do
    close(42)

    ! ====================================================================================================
    ! SUBROUTINES & FUNCTIONS
    ! ====================================================================================================

    contains

        ! ----------------------------------------------------------------------------------------------------
        ! generate a random cartesian coordinate in a sphere of radius 1

        subroutine random_coord_in_sphere(x,y,z)
            real    :: x, y, z
            logical :: lock

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
            REAL, INTENT(IN)    :: r, phi, theta
            REAL, INTENT(OUT)   :: x, y, z

            x = r * SIN(phi) * COS(theta)
            y = r * SIN(phi) * SIN(theta)
            z = r * COS(phi)

        end subroutine spherical_to_cartesian

        ! ----------------------------------------------------------------------------------------------------
        ! Function that ensure the new object will not appear too near from another one

        function is_too_near_from_another(x, y, z, x_L, y_L, z_L, N)
            integer             :: N
            real                :: volume, distance, eps
            real                :: PI = 3.14159
            real                :: x,y,z
            real, dimension(N)  :: x_L, y_L, z_L
            logical             :: is_too_near_from_another

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

end program nbody
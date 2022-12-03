! ----------------------------------------------------------------------------------------------------
! generate a random cartesian coordinate in a sphere of radius 1

subroutine random_coord_in_sphere(x,y,z)
    implicit none
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
    implicit none
    real,    intent(inout), dimension(3, N, steps) :: p, v, a ! time t and t+1
    integer, intent(in   )                         :: t ! step
    real,    intent(in   )                         :: m 
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
    call acceleration(p2(:, :), m, N, a2(:,:), 0.)

    do i=1,N

        ! Compute velocity at t+1
        vf(:, i) = vi(:, i) + ai(:, i) * dt/2 + a2(:, i) * dt/2

        ! Ec = Ec + 0.5 * m * (vf(1, i)**2 + vf(2, i)**2 + vf(3, i)**2)

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
    ! Et = Ep + Ec

end subroutine next_state

!----------------------------------------------------------------------------------------------------
! Compute the acceleration of a body considering the position of the others

subroutine acceleration(p, m, N, a, Ep)
    implicit none
    real,    intent(in   ), dimension(3, N) :: p
    real,    intent(in   )                  :: m
    integer, intent(in   )                  :: N
    real,    intent(  out), dimension(3,N)  :: a
    real        :: Ep
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

            a(1, i) = a(1, i) + x * G * m / r**3
            a(2, i) = a(2, i) + y * G * m / r**3
            a(3, i) = a(3, i) + z * G * m / r**3

            ! Ep = Ep - 0.5 * G * m**2 / r
        end do
    end do
    !$omp end do
    !$omp end parallel

end subroutine acceleration

!----------------------------------------------------------------------------------------------------
! Compute the energy of the system

! function potential_energy(p, m, N)

!     real,    intent(in   ), dimension(3, N) :: p
!     real,    intent(in   )                  :: m
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
!     real,    intent(in   )                  :: m
!     integer, intent(in   )                  :: N
!     real                                    :: kinetic_energy
!     integer                                 :: i

!     kinetic_energy = 0

!     ! Ec = sum on i of 1/2 * m_i * v_i^2
!     do i=1,N
!         kinetic_energy = kinetic_energy + 0.5 * m * sum(v(:, i)**2)
!     end do

! end function kinetic_energy

!----------------------------------------------------------------------------------------------------
! Get the barycenter

subroutine barycenter(p, m, N, b)
    implicit none
    real,    intent(in   ), dimension(3, N) :: p
    real,    intent(in   )                  :: m
    integer, intent(in   )                  :: N
    real,    intent(  out), dimension(3)    :: b
    integer                                 :: i

    b = 0

    ! Ec = sum on i of 1/2 * m_i * v_i^2
    do i=1,N
        b = b + m * p(:, i)
    end do

    b = b / N

end subroutine barycenter
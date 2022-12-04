program nbody

    implicit none

    !----------------------------------------------------------------------------------------------------
    ! Config

    ! Safe to edit
    integer, parameter :: N = 1000 ! number of bodies
    integer, parameter :: steps = 1000 ! simulation time in steps
    real, parameter    :: dt = 0.1 ! time step in seconds
    ! logical            :: use_initial_conditions = .false. ! set to .true. to generate new initial conditions
    ! logical            :: verbose = .false.
    ! logical            :: save_results = .true. ! set to .false. to disable saving results (for bulk runs)

    !----------------------------------------------------------------------------------------------------
    ! Other declarations

    ! Don't touch to that, you fool!
    integer                    :: i, t
    real                       :: PI
    real, dimension(N)         :: x, y, z
    real, dimension(N)         :: x2, y2, z2
    real, dimension(N)         :: vx, vy, vz
    real, dimension(N)         :: ax, ay, az
    real, dimension(N)         :: ax2, ay2, az2
    real                       :: Ep, Ec, Et
    real                       :: bx, by, bz ! barycenter
    logical                    :: lock

    PI = acos(-1.0)

    !----------------------------------------------------------------------------------------------------
    ! Output files

    open(10, file='data/positions.txt', status='replace')
    write(10, '(a)') '# x y z'
    open(11, file='data/velocities.txt', status='replace')
    write(11, '(a)') '# vx vy vz'
    open(12, file='data/accelerations.txt', status='replace')
    write(12, '(a)') '# ax ay az'
    open(13, file='data/energies.txt', status='replace')
    write(13, '(a)') '# Ep Ec Et'
    open(14, file='data/barycenter.txt', status='replace')
    write(14, '(a)') '# bx by bz'
    open(15, file='data/parameters.txt', status='replace')
    write(15, '(a)') '# N steps dt'
    write(15, *) N, steps, dt

    !----------------------------------------------------------------------------------------------------
    ! Initial positions

    do i=1,N
        lock = .true.
        do while (lock .eqv. .true.)
            CALL RANDOM_NUMBER(x(i))
            x(i) = x(i) * 2.0 - 1.0
            CALL RANDOM_NUMBER(y(i))
            y(i) = y(i) * 2.0 - 1.0
            CALL RANDOM_NUMBER(z(i))
            z(i) = z(i) * 2.0 - 1.0

            if (x(i)**2 + y(i)**2 + z(i)**2 < 1) then
                lock = .false.
            end if
        end do
    end do

    !----------------------------------------------------------------------------------------------------
    ! Initial speeds

    vx = -y
    vy = x
    vz = 0.0

    !----------------------------------------------------------------------------------------------------
    ! Initial accelerations

    call acceleration(x, y, z, N, ax, ay, az, Ep)

    !----------------------------------------------------------------------------------------------------
    ! Initial energies

    Ec = 0.0
    do i=1,N
        Ec = Ec + 0.5 * (vx(i)**2 + vy(i)**2 + vz(i)**2)
    end do
    Et = Ec + Ep

    !----------------------------------------------------------------------------------------------------
    ! Initial barycenter

    bx = 0.0
    by = 0.0
    bz = 0.0
    do i=1,N
        bx = bx + x(i)
        by = by + y(i)
        bz = bz + z(i)
    end do
    bx = bx / N
    by = by / N
    bz = bz / N

    !----------------------------------------------------------------------------------------------------
    ! Main loop

    do t=1,steps

        !------------------------------
        ! Position at t+1/2
        do i=1,N
            x2(i) = x(i) + vx(i) * dt/2. 
            y2(i) = y(i) + vy(i) * dt/2.
            z2(i) = z(i) + vz(i) * dt/2.
        end do

        !------------------------------
        ! Acceleration at t+1/2
        call acceleration(x2, y2, z2, N, ax2, ay2, az2)

        !------------------------------
        ! Compute velocity at t+1
        do i=1,N
            vx(i) = vx(i) + ax(i) * dt/2 + ax2(i) * dt/2
            vy(i) = vy(i) + ay(i) * dt/2 + ay2(i) * dt/2
            vz(i) = vz(i) + az(i) * dt/2 + az2(i) * dt/2
        end do

        !------------------------------
        ! Position at t+1
        do i=1,N
            x(i) = x2(i) + vx(i) * dt/2
            y(i) = y2(i) + vy(i) * dt/2
            z(i) = z2(i) + vz(i) * dt/2
        end do

        !------------------------------
        ! Acceleration at i+1
        call acceleration(x, y, z, N, ax, ay, az, Ep)

        !------------------------------
        ! Energy at t+1
        Ec = 0.0
        do i=1,N
            Ec = Ec + 0.5 * (vx(i)**2 + vy(i)**2 + vz(i)**2)
        end do
        Et = Ec + Ep

        !------------------------------
        ! Barycenter
        bx = 0.0
        by = 0.0
        bz = 0.0
        do i=1,N
            bx = bx + x(i)
            by = by + y(i)
            bz = bz + z(i)
        end do
        bx = bx / N
        by = by / N
        bz = bz / N

        !------------------------------
        ! Print results

        write(10, *) x, y, z
        write(11, *) vx, vy, vz
        write(12, *) ax, ay, az
        write(13, *) Ep, Ec, Et
        write(14, *) bx, by, bz

    end do

    !----------------------------------------------------------------------------------------------------
    ! Close files

    close(10)
    close(11)
    close(12)
    close(13)
    close(14)
    close(15)

    !====================================================================================================
    ! Subroutines & Functions
    !====================================================================================================

    contains

        !----------------------------------------------------------------------------------------------------   
        ! Compute accelerations

        subroutine acceleration(x, y, z, N, ax, ay, az, Ep)
            real,    intent(in   ), dimension(N) :: x, y, z
            real,    intent(out  ), dimension(N) :: ax, ay, az
            real,    intent(out  ), optional     :: Ep
            integer, intent(in   )               :: N
            real                                 :: m, G=1.0, eps=0.05
            real                                 :: dx, dy, dz, r
            integer                              :: i, j

            ax = 0
            ay = 0
            az = 0
            Ep = 0
            m = 1.0/N

            ! omp parallel reduction(+:Ep)
            ! omp do
            do i=1,N
                do j=1,N
                    if (j==i) then
                        cycle
                    end if
                    
                    dx = x(j) - x(i)
                    dy = y(j) - y(i)
                    dz = z(j) - z(i)

                    r = sqrt(dx*dx + dy*dy + dz*dz + eps*eps)

                    ax = ax + dx * G * m / (r**3)

                    Ep = ep - 0.5*G*m**2/r
                end do
            end do
            ! omp parallel
            ! omp do

        end subroutine acceleration


end program nbody
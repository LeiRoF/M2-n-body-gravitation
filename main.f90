program nbody

    implicit none

    !----------------------------------------------------------------------------------------------------
    ! Config

    ! Safe to edit
    integer, parameter :: N = 2000 ! number of bodies
    integer, parameter :: steps = 1000 ! simulation time in steps
    real, parameter    :: dt = 0.01 ! time step in seconds
    ! logical            :: use_initial_conditions = .false. ! set to .true. to generate new initial conditions
    logical            :: verbose = .false.
    ! logical            :: save_results = .true. ! set to .false. to disable saving results (for bulk runs)
    integer            :: method = 1 ! Method to compute acceleration. Complexity 1: N^2 , 2: N^2 /2

    !----------------------------------------------------------------------------------------------------
    ! Other declarations

    ! Don't touch to that, you fool!
    integer                    :: i, t
    real                       :: PI
    real, dimension(N)         :: x, y, z
    real, dimension(N)         :: vx, vy, vz
    real, dimension(N)         :: ax, ay, az
    real                       :: Ep, Ec, Et
    real                       :: bx, by, bz ! barycenter
    logical                    :: lock
    integer                    :: omp_get_num_threads

    PI = acos(-1.0)

    ! Verification that the program well run on desired number of threads
    !$omp parallel
    if (i .eq. 1 .and. verbose .eqv. .true.) then
        print *, 'Running on ', omp_get_num_threads(), ' threads'
    end if
    !$omp end parallel

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

    if (method .eq. 1) then
        call acceleration(x, y, z, N, ax, ay, az, Ep)
    else if (method .eq. 2) then
        call acceleration2(x, y, z, N, ax, ay, az, Ep)
    end if

    !----------------------------------------------------------------------------------------------------
    ! Initial energies

    Ec = 0.0
    do i=1,N
        Ec = Ec + 0.5 * (vx(i)**2 + vy(i)**2 + vz(i)**2) / N
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
    ! Saving initial conditions

    do i=1,N
        write(10, *) x(i), y(i), z(i)
        write(11, *) vx(i), vy(i), vz(i)
        write(12, *) ax(i), ay(i), az(i)
    end do
    write(13, *) Ep, Ec, Et
    write(14, *) bx, by, bz

    
    if (verbose .eqv. .true.) then
        print *, 'Step ', 1 , ' of ', steps, ' completed'
    end if

    !----------------------------------------------------------------------------------------------------
    ! Main loop

    do t=2,steps

        !------------------------------
        ! Position at t+1
        do i=1,N
            x(i) = x(i) + vx(i) * dt + 0.5 * ax(i) * dt**2
            y(i) = y(i) + vy(i) * dt + 0.5 * ay(i) * dt**2
            z(i) = z(i) + vz(i) * dt + 0.5 * az(i) * dt**2
        end do

        !------------------------------
        ! Compute velocity at t+1
        do i=1,N
            vx(i) = vx(i) + ax(i) * dt * 0.5
            vy(i) = vy(i) + ay(i) * dt * 0.5
            vz(i) = vz(i) + az(i) * dt * 0.5
        end do

        !------------------------------
        ! Acceleration at i+1

        if (method .eq. 1) then
            call acceleration(x, y, z, N, ax, ay, az, Ep)
        else if (method .eq. 2) then
            call acceleration2(x, y, z, N, ax, ay, az, Ep)
        end if

        !------------------------------
        ! Compute velocity at t+1
        do i=1,N
            vx(i) = vx(i) + ax(i) * dt * 0.5
            vy(i) = vy(i) + ay(i) * dt * 0.5
            vz(i) = vz(i) + az(i) * dt * 0.5
        end do

        !------------------------------
        ! Energy at t+1
        Ec = 0.0
        do i=1,N
            Ec = Ec + 0.5 * (vx(i)**2 + vy(i)**2 + vz(i)**2) / N
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
        ! Saving results

        do i=1,N
            write(10, *) x(i), y(i), z(i)
            write(11, *) vx(i), vy(i), vz(i)
            write(12, *) ax(i), ay(i), az(i)
        end do
        write(13, *) Ep, Ec, Et
        write(14, *) bx, by, bz

        if (verbose .eqv. .true.) then
            print *, 'Step ', t , ' of ', steps, ' completed'
        end if

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
            integer                              :: omp_get_num_threads

            Ep = 0
            m = 1.0/N
            ax = 0
            ay = 0
            az = 0

            !$omp parallel reduction(+:Ep,ax,ay,az) private(dx,dy,dz,r,G,m,eps,i,j) firstprivate(x,y,z,N)
            !$omp do schedule(auto)
            do i=1,N
                G=1.0
                eps=0.05
                m = 1./N

                do j=1,N

                    if (i .eq. j) cycle
                    
                    dx = x(j) - x(i)
                    dy = y(j) - y(i)
                    dz = z(j) - z(i)

                    r = sqrt(dx*dx + dy*dy + dz*dz + eps*eps)

                    ax(i) = ax(i) + dx * G * m / (r**3)
                    ay(i) = ay(i) + dy * G * m / (r**3)
                    az(i) = az(i) + dz * G * m / (r**3)

                    Ep = Ep - 0.5 * G * m*m / r
                end do

            end do
            !$omp end do
            !$omp end parallel

            
        end subroutine acceleration

        !----------------------------------------------------------------------------------------------------   
        ! Compute accelerations

        subroutine acceleration2(x, y, z, N, ax, ay, az, Ep)
            real,    intent(in   ), dimension(N) :: x, y, z
            real,    intent(out  ), dimension(N) :: ax, ay, az
            real,    intent(out  ), optional     :: Ep
            integer, intent(in   )               :: N
            real                                 :: m, G=1.0, eps=0.05
            real                                 :: dx, dy, dz, r
            integer                              :: i, j
            integer                              :: omp_get_num_threads

            Ep = 0
            m = 1.0/N
            ax = 0
            ay = 0
            az = 0

            !$omp parallel reduction(+:Ep) private(dx,dy,dz,r,G,m,eps,i,j) firstprivate(x,y,z,N) shared(ax,ay,az)
            !$omp do schedule(auto)
            do i=1,N
                G=1.0
                eps=0.05
                m = 1./N

                do j=i+1,N
                    
                    dx = x(j) - x(i)
                    dy = y(j) - y(i)
                    dz = z(j) - z(i)

                    r = sqrt(dx*dx + dy*dy + dz*dz + eps*eps)

                    ax(i) = ax(i) + dx * G * m / (r**3)
                    ay(i) = ay(i) + dy * G * m / (r**3)
                    az(i) = az(i) + dz * G * m / (r**3)

                    ax(j) = ax(j) - dx * G * m / (r**3)
                    ay(j) = ay(j) - dy * G * m / (r**3)
                    az(j) = az(j) - dz * G * m / (r**3)

                    Ep = Ep - G * m*m / r
                end do

            end do
            !$omp end do
            !$omp end parallel

            
        end subroutine acceleration2


end program nbody
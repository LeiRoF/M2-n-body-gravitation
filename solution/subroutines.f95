module subroutines

    use settings

    implicit none

    contains

    !##########################################
    ! Initialisation of the particles (random)
    !##########################################
    subroutine init_part(pos, vel)
        real(kind=dp), intent(inout) :: pos(3, n), vel(3, n)
        integer :: i
        real(kind=dp) :: r

        do i=1, n

            ! Initialize positions
            r = 2

            do while (r > 1.0)

                call random_number(pos(:, i))
                pos(:, i) = 1-2*pos(:, i)
                r = (sum(pos(:, i)**2))**0.5

            end do

            ! Initialize velocities based on solid rotation
            vel(1, i) = -pos(2, i)
            vel(2, i) = pos(1, i)
            vel(3, i) = 0

        end do

        ! Uncomment to save new initial positions 

        ! open(unit=300, file='pos.dat', action='write')
        ! do i=1, n
        !     write(300, *) pos(1, i), pos(2, i), pos(3, i)
        ! end do
        ! close(300)

    end subroutine init_part

    !###############################################
    ! Initialisation of the particles (from a file)
    !###############################################
    subroutine init_part_file(pos, vel)
        real(kind=dp), intent(inout) :: pos(3, n), vel(3, n)
        integer :: i

        open(200, file = init_file, action = 'read')

        do i=1, n
            read(200, *) pos(:, i)
        
            vel(1, i) = -pos(2, i)
            vel(2, i) = pos(1, i)
            vel(3, i) = 0
        end do

        close(200)

    end subroutine init_part_file


    !###########################################
    !           Compute accelerations
    !###########################################
    subroutine accel(pos, acc, ep)
        real(kind=dp), intent(in) :: pos(3, n)
        real(kind=dp), intent(out) :: acc(3, n), ep
        real(kind=dp) :: dpos(3), r
        integer :: i, j

        ! Reset potential energy and accelerations
        ep = 0.
        acc = 0.

        do i=1, n
            do j=1, n
                if (j==i) then
                    cycle
                end if
                
                dpos(:) = pos(:, j) - pos(:, i)

                r = sqrt(sum(dpos**2) + epsilon**2)

                acc(:, i) = acc(:, i) + (gg*m/r**3)*dpos(:)

                ! Compute the potential energy
                ep = ep - 0.5*gg*m**2/r

            end do
        end do

    end subroutine accel


    !###########################################
    !            Leapfrog integrator
    !###########################################
    subroutine leapfrog(pos, vel, acc)
        real(kind=dp), intent(inout) :: pos(3, n), vel(3, n), acc(3, n)
        real (kind=dp) :: acc_new(3 ,n), ec, ep
        integer :: i

        ! Reset kinetic energy
        ec = 0

        ! Update positions
        do i=1, n
            pos(:, i) = pos(:, i) + vel(:, i)*dt + 0.5*acc(:, i)*dt**2
        end do

        ! Compute new accelerations
        call accel(pos, acc_new, ep)

        ! Update velocities
        do i=1, n
            vel(:, i) = vel(:, i) + 0.5*(acc(:, i) + acc_new(:, i))*dt
        end do

        ! Update accelerations and compute kinetic energy 
        do i=1, n
            acc(:, i) = acc_new(:, i)
            ec = ec + sum(vel(:, i)**2)
        end do

        ! Save kinetic and potantial energies
        write(100, *) ep, ec*0.5*m

    end subroutine leapfrog

    !###########################################
    !         Evolution of the system
    !###########################################
    subroutine evolve(pos, vel, acc)
        real(kind=dp), intent(inout) :: pos(3, n), vel(3, n), acc(3, n)
        real(kind=dp) :: ep
        integer :: i

        ! Open files to save the positions, accelerations and the energies at each time step
        open(unit = 10, file = outpos_name, form = 'unformatted', access = 'stream', action = 'write')
        open(unit = 100, file = energ_name, action="write")

        ! Compute accelerations at t=0
        call accel(pos, acc, ep)

        ! Save positions and accelerations at t=0
        write(10) pos
        
        ! Main loop
        do i=1, nt

            call leapfrog(pos, vel, acc)

            write(10) pos

            if (verbose) then
                print*, "step ", i, " done"
            end if

        end do

        close(10)
        close(100)

    end subroutine evolve

end module subroutines

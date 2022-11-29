program main

    use settings
    use subroutines

    implicit NONE

    real(kind=dp), allocatable :: pos(:,:), vel(:,:), acc(:,:)

    ! Read some settings in 'input.ini'
    call readsettings()

    ! Allocate the main arrays
    allocate(pos(3, n))
    allocate(vel(3, n))
    allocate(acc(3, n))

    ! Initalisation of the particles
    if (use_init_file) then
        call init_part_file(pos, vel)
    else
        call init_part(pos, vel)
    end if
    
    ! Evolution of the system
    call evolve(pos, vel, acc)

    ! Deallocate the main arrays
    deallocate(pos)
    deallocate(vel)
    deallocate(acc)
    
end program main
module settings
    implicit none

    ! Precision selector
    integer, public, parameter :: dp = SELECTED_REAL_KIND(P=8)
    
    ! Number of particles
    integer, public :: n

    ! Number of time step
    integer, public :: nt

    ! Time step
    real(kind=dp), public :: dt

    ! Mass of each particle
    real(kind=dp), public :: m

    ! Name of the file to save the positions
    character(len=*), public, parameter :: outpos_name = 'position_fort_seq.bin'

    ! Name of the file to save the energies
    character(len=*), public, parameter :: energ_name = 'energy_fort_seq.dat'

    ! Choice between initial positions from a file or computed randomly
    logical, public :: use_init_file

    ! Name of the file to load the initial positions (not necessary if use_init_file is set to .false.)
    character(len=*), public, parameter :: init_file = 'pos.dat'

    ! Minimal distance between two particles (smoothing lenght)
    real(kind=dp), public :: epsilon

    ! Pi
    real(kind=dp), public, parameter :: pi = acos(-1.0)

    ! Gravitational constant
    real(kind=dp), public, parameter :: gg = 1.

    ! Show simulation evolution
    logical, public :: verbose

    contains

    subroutine readsettings()

        integer :: fileend, l
        character(len=20) :: settingname
        real(kind=dp) :: value
    
        fileend = 0
        l= 0
        open(unit = 101, file = 'input.ini', action = 'read')

        do while (fileend >= 0)
            read(101, *, IOSTAT=fileend) settingname, value
            
            if (settingname == "npart") then
                n = int(value)
                m = 1/value
            end if

            if (settingname == "nstep") then
                nt = int(value)
            end if

            if (settingname == "dt") then
                dt = value
            end if

            if (settingname == "epsilon") then
                epsilon = value
            end if

            if (settingname == 'init') then
                if (int(value) == 0) then
                    use_init_file = .false.
                else
                    use_init_file = .true.
                    n = 2000
                    m = 1.0/2000.0
                end if
            end if

            if (settingname == 'verbose') then
                if (int(value) == 0) then
                    verbose = .false.
                else
                    verbose = .true.
                end if
            end if
        end do
        
        close(101)
    end subroutine readsettings

end module settings

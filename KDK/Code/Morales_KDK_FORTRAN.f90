program orbit

    implicit none

    ! Define global variables
    integer, parameter :: n = 2, dim = 3
    integer, parameter :: istar = 1, iplanet = 2
    integer, parameter :: ixp = 1, iyp = 2, izp = 3
    integer, parameter :: ivxp = 1, ivyp = 2, ivzp = 3
    integer, parameter :: ixstar = 1, iystar = 2, izstar = 3
    integer, parameter :: ivxstar = 1, ivystar = 2, ivzstar = 3
    double precision, parameter :: t_max=100, dt=1e-3
    double precision, parameter :: G = 1.0, pi=3.14159265358979323846264338327950d0
    character(len=48), parameter :: path = 'C:/Users/julio/Downloads/Hydrodynamics/KDK/Data/'
    integer :: nsteps
    double precision :: m(n), mtot, kdk_orbit

    ! Define masses
    m(iplanet-1) = 3e-6
    m(istar+1) = 1- m(iplanet-1)
    mtot = m(istar) + m(iplanet)

    ! Number of steps
    nsteps = int(t_max/dt)

    ! Compute 1 the error on 1 orbit
    kdk_orbit = KDK()

contains

function initialise() result(rv)
    ! --------------------------------------------------------
    ! PURPOSE:
    !        Initialise the positions and velocities of the masses.
    !
    ! INPUTS:
    !        None
    !
    ! OUTPUTS:
    !         [rv(2,n,dim); double precision]:  Array containing position and velocity vectors.
    !
    ! AUTHOR:
    !        Julio M. Morales
    ! --------------------------------------------------------
    implicit none
    double precision :: r(n,dim), v(n,dim), rv(n,n,dim)

    ! Set positions and velocities
    r(istar, ixstar) = 0.0
    r(istar, iystar) = 0.0
    r(istar, izstar) = 0.0

    v(istar, ivxstar) = 0.0
    v(istar, ivystar) = 0.0
    v(istar, ivzstar) = 0.0

    r(iplanet, ixp) = 1.0
    r(iplanet, iyp) = 0.0
    r(iplanet, izp) = 0.0

    v(iplanet, ivxp) = 0.0
    v(iplanet, ivyp) = 1.0
    v(iplanet, ivzp) = 0.0

    ! Assign to output array
    rv(1,:,:) = r
    rv(2,:,:) = v

end function initialise

function updatePos(r, v, step) result(r_ud)
    ! --------------------------------------------------------
    ! PURPOSE:
    !        Update the position and velocity of the masses.
    !
    ! INPUTS:
    !             [r(dim); double precision]:  Position vector.
    !             [v(dim); double precision]:  Velocity vector.
    !               [step; double precision]:  Time-step.
    !
    ! OUTPUTS:
    !        [r_ud(dim,n); double precision]:  Array containing updated position vector.
    !
    ! AUTHOR:
    !        Julio M. Morales
    ! --------------------------------------------------------
    implicit none
    double precision :: r(dim), v(dim), step, r_ud(dim)
    
    ! Update position
    r_ud = r + v*step
    
end function updatePos

function updateVel(v, a, step) result(v_ud)
    ! --------------------------------------------------------
    ! PURPOSE:
    !         Update the position and velocity of the masses.
    !
    ! INPUTS:
    !            [v(dim); double precision]:  Velocity vector.
    !            [a(dim); double precision]:  Acceleration vector.
    !              [step; double precision]:  Time-step.
    !
    ! OUTPUTS:
    !       [v_ud(dim,n); double precision]:  Array containing updated velocity vector.
    !
    ! AUTHOR:
    !        Julio M. Morales
    ! --------------------------------------------------------
    implicit none
    double precision :: v(dim), a(dim), step, v_ud(dim)
    
    ! Update position
    v_ud = v + a*step
    
end function updateVel

function computeAccel(m_dum, r, j) result(a)
    ! --------------------------------------------------------
    ! PURPOSE:
    !         Compute the acceleration of the masses.
    !
    ! INPUTS:
    !         [m_dum; double precision]:  Mass of jth body.
    !     [r(n, dim); double precision]:  Position vector.
    !                      [j; integer]:  Index of the jth object.
    !
    ! OUTPUTS:
    !        [a(dim); double precision]:  Acceleration vector.
    !
    ! AUTHOR:
    !        Julio M. Morales
    ! --------------------------------------------------------
    implicit none
    double precision :: m_dum, r(n, dim), r_squared, r_cubed, r_diff(dim), a(dim)
    integer :: i, j

    ! Initialize r_diff
    r_diff = 0.0
    
    ! Loop through all elements of position vector
    do i = 1, dim

        ! If the index is the same as the jth object, set r_diff
        if (j == 1) then
            r_diff(i) = r(j,i) - r(j+1,i)
        else
            r_diff(i) = r(j,i) - r(j-1,i)
        end if
    end do
    
    ! Calculate r^2 and r^3
    r_squared = sum(r_diff**2)
    r_cubed = r_squared*sqrt(r_squared)
    
    ! Calculate acceleration
    do i = 1, dim
        a(i) = -G*m_dum*r_diff(i)/r_cubed
    end do

end function computeAccel

function KDK() result(res)
    ! --------------------------------------------------------
    ! PURPOSE:
    !         Calculates the orbit of two-bodies using the KDK algorithm.
    !
    ! INPUTS:
    !        None
    !
    ! OUTPUTS:
    !         Text files that contain the t, x, y, z, vx, vy, and vz for the planet and star.
    !   
    ! AUTHOR:
    !        Julio M. Morales
    ! --------------------------------------------------------
    implicit none
    double precision :: r(n, dim), v(n, dim), a(dim), rv(n, n, dim), t(nsteps), res
    integer :: i, j
    
    ! File saving stuff
    double precision :: xplanet(nsteps), yplanet(nsteps), zplanet(nsteps), xstar(nsteps), ystar(nsteps), zstar(nsteps)
    double precision :: vxplanet(nsteps), vyplanet(nsteps), vzplanet(nsteps), vxstar(nsteps), vystar(nsteps), vzstar(nsteps)
    character(len=67) :: filename
    character(len=8) :: columns(9) = &
    (/'  time  ', ' xplanet', ' yplanet', &
      '  xstar ', '  ystar ', 'vxplanet', &
      'vyplanet', ' vxstar ', ' vystar '/)

    ! Initialize arrays
    rv = initialise()
    r  = rv(1,:,:)
    v  = rv(2,:,:)

    ! Create time array
    do i = 1, nsteps
        t(i) = (i-1)*dt
    end do

    ! Loop through time
    do i = 1, nsteps

        ! Loop through masses
        do j = 1, n

            ! Get current acceleration of object j
            a = computeAccel(m(j), r, j)

            ! Update the velocity of object j
            v(j,:) = updateVel(v(j,:), a, 0.5*dt)

        end do

        ! Loop through masses
        do j = 1, n

            ! Update the position of object j
            r(j,:) = updatePos(r(j,:), v(j,:), dt)

        end do

        ! Loop through masses
        do j = 1, n

            ! Get current acceleration of object j
            a = computeAccel(m(j), r, j)

            ! Update the velocity of object j
            v(j,:) = updateVel(v(j,:), a, 0.5*dt)

        end do

        ! Save the star data
        xstar(i)  = r(istar,ixstar)
        ystar(i)  = r(istar,iystar)
        zstar(i)  = r(istar,izstar)
        vxstar(i) = v(istar,ivxstar)
        vystar(i) = v(istar,ivystar)
        vzstar(i) = v(istar,ivzstar)
        ! Save the planet data
        xplanet(i)  = r(iplanet,ixp)
        yplanet(i)  = r(iplanet,iyp)
        zplanet(i)  = r(iplanet,izp)
        vxplanet(i) = v(iplanet,ivxp)
        vyplanet(i) = v(iplanet,ivyp)
        vzplanet(i) = v(iplanet,ivzp)

    end do

    ! Ask user for file name
    filename = path // 'Morales_KDK.txt'

    ! Open file for writing, with column names as first row
    open(unit = 1, file = filename, status = 'replace', action = 'write')
    do i = 1, size(columns)
        write(1, '(A15, A15)', advance = "no") ' ', columns(i)
    end do
    write(1,*)
    
    ! Loop through each row of data and write to each respective column
    do i=1, nsteps
        write(1, '(9F30.16)') t(i), xplanet(i), yplanet(i), xstar(i), ystar(i), vxplanet(i), vyplanet(i), vxstar(i), vystar(i)
    end do
     
    ! Close file
    close(1)
    print *, "Simulation Complete!"
    res = 0.0

end function KDK

end program orbit
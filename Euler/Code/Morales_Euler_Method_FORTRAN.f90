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
    double precision :: m(n), mtot, euler_orbit
    character(len=45), parameter :: path = 'C:/Users/julio/Downloads/Hydrodynamics/Euler/'
    integer :: nsteps
    
    ! Define masses
    m(iplanet-1) = 3e-6
    m(istar+1) = 1- m(iplanet-1)
    mtot = m(istar) + m(iplanet)

    ! Number of steps
    nsteps = int(t_max/dt)

    ! Compute 1 the error on 1 orbit
    euler_orbit = euler() 

contains

function initialise() result(rv)
    ! --------------------------------------------------------
    ! PURPOSE:
    !   Initialise the positions and velocities of the masses.
    !
    ! INPUTS:
    !   None
    !
    ! OUTPUTS:
    !   [rv(2,n,dim); double precision]:  Array containing position and velocity vectors.
    !
    ! AUTHOR:
    !   Julio M. Morales
    ! --------------------------------------------------------
    implicit none
    double precision :: r(n,dim), v(n,dim), rv(2,n,dim)

    ! Set initial positions and velocities
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

function updatePosVel(r, v, a) result(rv_ud)
    ! --------------------------------------------------------
    ! PURPOSE:
    !   Update the position and velocity of the masses.
    !
    ! INPUTS:
    !            [r(dim); double precision]:  Position vector.
    !            [v(dim); double precision]:  Velocity vector.
    !            [a(dim); double precision]:  Acceleration vector.
    !
    ! OUTPUTS:
    !      [rv_ud(dim,n); double precision]:  Array containing position and velocity vectors.
    !
    ! AUTHOR:
    !   Julio M. Morales
    ! --------------------------------------------------------
    implicit none
    double precision :: r(dim), v(dim), a(dim), rv_ud(dim,n)
    
    ! Update position and velocity
    r = r + v*dt
    v = v + a*dt
    
    ! Store updated position and velocity vectors in the result list
    rv_ud(:,1) = r
    rv_ud(:,2) = v
    
end function updatePosVel

function computeAccel(m_dum, r, j) result(a)
    ! --------------------------------------------------------
    ! PURPOSE:
    !   Compute the acceleration of the masses.
    !
    ! INPUTS:
    !         [m_dum; double precision]:  Mass of jth body.
    !     [r(n, dim); double precision]:  Position vector.
    !             [j; integer]:  Index of the jth object.
    !
    ! OUTPUTS:
    !        [a(dim); double precision]:  Acceleration vector.
    !
    ! AUTHOR:
    !   Julio M. Morales
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

function euler() result(res)

    ! --------------------------------------------------------
    ! PURPOSE:
    !   Computes the relative error in energy for a given time step and returns the result.
    !
    ! INPUTS:
    !   None
    !
    ! OUTPUTS:
    !         Text file containing position and velocity data of both bodies.
    !
    ! AUTHOR:
    !   Julio M. Morales
    ! --------------------------------------------------------
    implicit none
    integer :: i, j
    double precision ::  r(n,dim), v(n,dim), a(dim), rv(n,n,dim), rv_ud(dim, n), t(nsteps), res

    ! File saving stuff
    double precision  :: xplanet(nsteps), yplanet(nsteps), zplanet(nsteps), xstar(nsteps), ystar(nsteps), zstar(nsteps)
    double precision  :: vxplanet(nsteps), vyplanet(nsteps), vzplanet(nsteps), vxstar(nsteps), vystar(nsteps), vzstar(nsteps)
    character(len=67) :: filename
    character(len=8) :: columns(9) = &
    (/'  time  ', ' xplanet', ' yplanet', &
      '  xstar ', '  ystar ', 'vxplanet', &
      'vyplanet', ' vxstar ', ' vystar '/)


    ! Declare the initial position and velocity for both the Sun and the Earth
    rv = initialise()
    r  = rv(1,:,:)
    v  = rv(2,:,:)

    ! Initialize the time array
    do i = 1, nsteps
        t(i) = (i-1)*dt
    end do

    ! Loop over time and amount of objects
    do i = 1, nsteps
        do j = 1, n

            ! Store the position and velocity of the star
            if (j == istar) then
                xstar(i) = r(j,ixstar)
                ystar(i) = r(j,iystar)
                zstar(i) = r(j,izstar)

                vxstar(i) = v(j,ivxstar)
                vystar(i) = v(j,ivystar)
                vzstar(i) = v(j,ivzstar)
            end if

            ! Store the position and velocity of the planet
            if (j == iplanet) then
                xplanet(i) = r(j,ixp)
                yplanet(i) = r(j,iyp)
                zplanet(i) = r(j,izp)

                vxplanet(i) = v(j,ivxp)
                vyplanet(i) = v(j,ivyp)
                vzplanet(i) = v(j,ivzp)
            end if

            ! Get the current acceleration of object j 
            a = computeAccel(m(j), r, j)

            ! Update the position and velocity
            rv_ud  = updatePosVel(r(j,:), v(j,:), a)
            r(j,:) = rv_ud(:,1)
            v(j,:) = rv_ud(:,2)

        end do

    end do

    ! Ask user for file name
    filename = path // 'Morales_Euler.txt'

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

end function euler

end program orbit
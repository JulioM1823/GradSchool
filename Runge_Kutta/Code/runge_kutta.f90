program orbit

    implicit none

    ! Define global variables
    integer, parameter :: n = 2, dim = 3
    integer, parameter :: istar = 1, iplanet = 2
    integer, parameter :: ixp = 1, iyp = 2, izp = 3
    integer, parameter :: ivxp = 1, ivyp = 2, ivzp = 3
    integer, parameter :: ixstar = 1, iystar = 2, izstar = 3
    integer, parameter :: ivxstar = 1, ivystar = 2, ivzstar = 3
    real*10, parameter :: t_max=100, dt=1e-3
    real*10, parameter :: G = 1.0, pi=3.14159265358979323846264338327950d0
    character(len=51), parameter :: path = 'C:/Users/julio/Downloads/Hydrodynamics/Runge_Kutta/'
    integer :: nsteps
    real*10 :: m(n), mtot, runge_orbit

    ! Define masses
    m(iplanet-1) = 3e-6
    m(istar+1) = 1- m(iplanet-1)
    mtot = m(istar) + m(iplanet)

    ! Number of steps
    nsteps = int(t_max/dt)

    ! Compute 1 the error on 1 orbit
    runge_orbit = rungeKutta()

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
    !   rv(2,n,dim):  Array containing position and velocity vectors.
    !
    ! AUTHOR:
    !   Julio M. Morales
    ! --------------------------------------------------------
    implicit none
    real*10 :: r(n,dim), v(n,dim), rv(n,n,dim)

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

function updatePosVel(r, v, a) result(rv_ud)
    ! --------------------------------------------------------
    ! PURPOSE:
    !   Update the position and velocity of the masses.
    !
    ! INPUTS:
    !   [r(3); real*10]:  Position vector.
    !   [v(3); real*10]:  Velocity vector.
    !   [a(3); real*10]:  Acceleration vector.
    !
    ! OUTPUTS:
    !      [rv_ud(3,2)]:  Array containing position and velocity vectors.
    !
    ! AUTHOR:
    !   Julio M. Morales
    ! --------------------------------------------------------
    implicit none
    real*10 :: r(dim), v(dim), a(dim), rv_ud(dim,n)
    
    ! Update position and velocity
    r = r + v*dt
    v = v + a*dt
    
    ! Store updated position and velocity vectors in the result list
    rv_ud(:,1) = r
    rv_ud(:,2) = v
    
end function updatePosVel

function computeAccel(r) result(a)
    ! --------------------------------------------------------
    ! PURPOSE:
    !   Compute the acceleration of the masses.
    !
    ! INPUTS:
    !       [m_dum; real*10]:  Mass of jth body.
    !   [r(n, dim); real*10]:  Position vector.
    !           [j; integer]:  Index of the jth object.
    !
    ! OUTPUTS:
    !        [a(n, dim); real*10]:  Acceleration vector.
    !
    ! AUTHOR:
    !   Julio M. Morales
    ! --------------------------------------------------------
    implicit none
    real*10 :: m_dum, r(n, dim), r_squared, r_cubed, r_diff(dim), a(n, dim)
    integer :: i, j

    ! Initialize r_diff
    r_diff = 0.0

    ! Loop over all objects
    do j = 1, n

        ! Set mass of jth object
        m_dum = m(j)
        
        ! Loop through all elements of position vector
        do i = 1, dim

            ! If the index is the same as the jth object, set r_diff
            if (j == 1) then
                r_diff(i) = r(j,i) - r(j+1,i)
            else
                r_diff(i) = r(j,i) - r(j-1,i)
            end if
                    
            ! Calculate r^2 and r^3
            r_squared = sum(r_diff**2)
            r_cubed = r_squared*sqrt(r_squared)

            ! Calculate acceleration
            a(j, i) = -G*m_dum*r_diff(i)/r_cubed

        end do

    end do

end function computeAccel

function computeEnergy(r, v) result(E)
    ! --------------------------------------------------------
    ! PURPOSE:
    !   Compute the energy of the system.
    !
    ! INPUTS:
    !       [r(n, dim); real*10]:  Position vector.
    !       [v(n, dim); real*10]:  Velocity vector.
    !
    ! OUTPUTS:
    !       [E; real*10]:  Energy of the system.
    !
    ! AUTHOR:
    !   Julio M. Morales
    ! --------------------------------------------------------
    implicit none
    real*10 :: r(n, dim), v(n, dim), E

    ! Compute energy
    E = 0.5*(m(1)*dot_product(v(1,:), v(1,:)) + m(2)*dot_product(v(2,:), v(2,:))) - m(1)*m(2)/norm2(r(2,:) - r(1,:))

end function computeEnergy

function computeAngularMomentum(r, v) result(L)
    ! --------------------------------------------------------
    ! PURPOSE:
    !   Compute the angular momentum of the system.
    !
    ! INPUTS:
    !   [r(n, dim); real*10]:  Position vector.
    !   [v(n, dim); real*10]:  Velocity vector.
    !
    ! OUTPUTS:
    !      [L(dim); real*10]:  Angular momentum of the system.
    !
    ! AUTHOR:
    !   Julio M. Morales
    ! --------------------------------------------------------
    implicit none
    real*10 :: r(n, dim), v(n, dim), rDiff(dim), vDiff(dim), rr(dim), phi(dim), vrad(dim), vphi(dim), phidot(dim), L(dim)
    
    ! Define rDiff and vDiff
    rDiff = r(iplanet,:) - r(istar,:)
    vDiff = v(iplanet,:) - v(istar,:)

    ! Converting the coordinates from Cartesian to Cylindrical ####
    rr  = sqrt(rDiff(istar)**2 + rDiff(iplanet)**2)
    phi = atan2(rDiff(iplanet), rDiff(istar))

    ! Velocity
    vrad = vDiff(istar)*cos(phi) + vDiff(iplanet)*sin(phi)
    vphi = -vDiff(istar)*sin(phi) + vDiff(iplanet)*cos(phi)

    ! Time derivative of phi
    phidot = vphi/rr

    ! Calculate the angular momentum
    L = (rr**2)*(phidot)

end function computeAngularMomentum

function rungeKutta() result(sig_E)
    ! --------------------------------------------------------
    ! PURPOSE:
    !        Calculates the relative error in the energy evolution for a given timestep
    !
    ! INPUTS:
    !        None
    !
    ! OUTPUTS:
    !        Text files that contain the t, x, y, z, vx, vy, vz, and energy error for the planet and star.
    !
    ! AUTHOR:
    !        Julio M. Morales
    ! --------------------------------------------------------
    implicit none
    real*10 :: r(n, dim), v(n, dim), E0, coeff(4), kx(5,n,dim) = 0.0d0, kv(5,n,dim) = 0.0d0, rv(n, n, dim), t(nsteps), E, sig_E
    integer :: i, k
    
    ! File saving stuff
    real*10 :: xplanet(nsteps), yplanet(nsteps), zplanet(nsteps), xstar(nsteps), ystar(nsteps), zstar(nsteps)
    real*10 :: vxplanet(nsteps), vyplanet(nsteps), vzplanet(nsteps), vxstar(nsteps), vystar(nsteps), vzstar(nsteps)
    real*10 :: En_error(nsteps)
    character(len=76) :: filename1, filename2, filename5
    character(len=79) :: filename3, filename4
    character(len=7)  :: col_names1(4) = (/'time(n)', 'xplanet', 'yplanet', 'zplanet'/)
    character(len=5)  :: col_names2(3) = (/'xstar', 'ystar', 'zstar' /)
    character(len=8)  :: col_names3(3) = (/'vxplanet', 'vyplanet', 'vzplanet'/)
    character(len=6)  :: col_names4(3) = (/'vxstar', 'vystar', 'vzstar' /)
    character(len=10) :: col_names5(1) = (/'energy_err' /)

    ! Initialize the position and velocity vectors
    rv = initialise()
    r  = rv(1,:,:)
    v  = rv(2,:,:)

    ! Compute initial energy
    E0 = computeEnergy(r, v)

    ! Create time array
    do i = 1, nsteps
        t(i) = (i-1)*dt
    end do

    ! Coefficients for Runge-Kutta
    coeff(1) = 1.0
    coeff(2) = 0.5*dt
    coeff(3) = 0.5*dt
    coeff(4) = dt

    ! Loop over time
    do i = 1, nsteps

        ! Compute kth element of kx and kv
        do k = 1, size(coeff)

            ! Calculate kx and kv
            kx(k,:,:) = v
            kv(k,:,:) = computeAccel(r)

        end do

        ! Calculate the new position and velocity
        r = r + dt*(kx(1,:,:) + 2*kx(2,:,:) + 2*kx(3,:,:) + kx(4,:,:))/6
        v = v + dt*(kv(1,:,:) + 2*kv(2,:,:) + 2*kv(3,:,:) + kv(4,:,:))/6

        ! Compute the energy
        E = computeEnergy(r, v)

        ! Compute relative error in energy
        sig_E = abs((E - E0)/E0)

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

        ! Save the energy error
        En_error(i) = sig_E

    end do
    
    ! Ask user for file name
    filename1 = path // 'Morales_RungeKutta_rp.txt'
    filename2 = path // 'Morales_RungeKutta_vp.txt'
    filename3 = path // 'Morales_RungeKutta_rstar.txt'
    filename4 = path // 'Morales_RungeKutta_vstar.txt'
    filename5 = path // 'Morales_RungeKutta_En.txt'

    ! Open file for writing, with column names as first row
    open(unit = 1, file = filename1, status = 'replace', action = 'write')
    do i = 1, size(col_names1)
        write(1, '(A15, A15)', advance = "no") ' ', col_names1(i)
    end do
    write(1,*)

    ! Open file for writing, with column names as first row
    open(unit = 2, file = filename2, status = 'replace', action = 'write')
    do i = 1, size(col_names2)
        write(2, '(A15, A15)', advance = "no") ' ', col_names2(i)
    end do
    write(2,*)

    ! Open file for writing, with column names as first row
    open(unit = 3, file = filename3, status = 'replace', action = 'write')
    do i = 1, size(col_names3)
        write(3, '(A15, A15)', advance = "no") ' ', col_names3(i)
    end do
    write(3,*)

    ! Open file for writing, with column names as first row
    open(unit = 4, file = filename4, status = 'replace', action = 'write')
    do i = 1, size(col_names4)
        write(4, '(A15, A15)', advance = "no") ' ', col_names4(i)
    end do
    write(4,*)

    ! Open file for writing, with column names as first row
    open(unit = 5, file = filename5, status = 'replace', action = 'write')
    do i = 1, size(col_names5)
        write(5, '(A15, A15)', advance = "no") ' ', col_names5(i)
    end do
    write(5,*)
    
    ! Loop through each row of data and write to each respective column
    do i=1, nsteps
        write(1, '(4F30.8)') t(i), xplanet(i), yplanet(i), zplanet(i)
        write(2, '(3F30.8)') xstar(i), ystar(i), zstar(i)
        write(3, '(3F30.8)') vxplanet(i), vyplanet(i), vzplanet(i)
        write(4, '(3F30.8)') vxstar(i), vystar(i), vzstar(i)
        write(5, '(1F30.8)') En_error(i)
    end do
     
    ! Close file
    close(1)
    close(2)
    close(3)
    close(4)
    close(5)
    print *, "Simulation Complete!"

end function rungeKutta

end program orbit
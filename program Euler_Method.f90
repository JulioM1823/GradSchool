program Euler_Method

    ! Declare that no variables be implicit
    implicit none
    
    ! Initialize mass vector and number of orbits to be calculated
    real(kind = 8), parameter :: m(1:2) = [1.0d0, 1.0d-3]

    ! Function to initialise the position and velocity vectors.
    function initialize() result (vectors)

        ! Define type vectors
        type vectors
            real(kind = 8), dimension(2, 3) :: r
            real(kind = 8), dimension(2, 3) :: v
        end type vectors

        ! Define initial positions of Sun and Earth
        real(kind = 8), parameter :: r_S = [0.0d0, 0.0d0, 0.0d0]
        real(kind = 8), parameter :: r_E = [1.0d0, 0.0d0, 0.0d0]
        real(kind = 8), parameter :: v_S = [0.0d0, 0.0d0, 0.0d0]
        real(kind = 8), parameter :: v_E = [0.0d0, 1.0d0, 0.0d0]

        ! Initialize the position and velocity vectors
        vectors%r = [r_S, r_E]
        vectors%v = [v_S, v_E]

    end function initialize
     
    ! Function to accept the mass, the current position index, and the position array (for all the objects in the system) as arguments and return the acceleration.
    function computeAcc(m, r, j, unitVec)

        ! Declare variables
        real(kind = 8), parameter :: m(1:2)
        real(kind = 8), parameter :: r(1:3)
        integer, parameter :: j
        logical :: unitVec
        real(kind = 8), parameter :: rDiff_norm
        real(kind = 8), parameter :: rUnit
        real(kind = 8), parameter :: a

        ! Compute the difference in position
        r_new(1:j) = r(1:j)
        r_new(j+1:j+n-1) = r(j+1:n)
        rDiff = r(j) - r_new

        ! Define positional difference vector
        rDiff_norm = (abs(rDiff(1)**2 + rDiff(2)**2))**0.5

        ! Check which method to use to calculate acceleration
        if (unitVec==.true.) then

            ! Define the unit vector
            rUnit = rDiff/(rDiff_norm)
            
            ! Compute acceleration
            a = -m*rUnit/(rDiff_norm)**2

        else
            ! Calculate acceleration
            a = -m*rDiff/(rDiff_norm)**3

        return a

    end function computeAcc

    ! Function to update the position and velocity in a timestep 'dt' and returns the updated position and velocity.
    function updatePosVel(r, v, a, dt), result (updatedVectors)
        
        ! Define type vectors
        type updatedVectors
            real(kind = 8), dimension(2, 3) :: r
            real(kind = 8), dimension(2, 3) :: v
        end type updatedVectors

        ! Updating the position
        r = r + v*dt

        ! Updating the velocity
        v = v + a*dt

    end function updatePosVel

    function computeEnergy(m, r, v), result (E)

        ! Declare variables
        real(kind = 8), parameter :: m(1:2)
        real(kind = 8), parameter :: r(1:3)
        real(kind = 8), parameter :: v(1:3)
        real(kind = 8), parameter :: r21_norm
        real(kind = 8), parameter :: v2_norm
        real(kind = 8), parameter :: v1_norm
        real(kind = 8), parameter :: E

        ! Calculate the normalized vectors
        r21_norm = (abs(( r(2)(1) - r(1)(1) )**2 + ( r(2)(2) - r(1)(2) )**2 + ( r(2)(3) - r(1)(3) )**2))**0.5
        v2_norm  = (abs(v(2)(1)**2 + v(2)(2)**2 + v(2)(3)**2))**0.5
        v1_norm  = (abs(v(1)(1)**2 + v(1)(2)**2 + v(1)(3)**2))**0.5

        ! Calculate the total mechanical energy
        E = 0.5*m(1)*v2_norm**2 + 0.5*m(2)*v1_norm**2 - m(1)*m(2)/r21_norm

    end function computeEnergy

    function computeAngularMomentum(r, v), result (h)

        real(kind = 8), parameter :: r(1:3)
        real(kind = 8), parameter :: v(1:3)
        real(kind = 8), parameter :: rDiff
        real(kind = 8), parameter :: vDiff
        real(kind = 8), parameter :: rr
        real(kind = 8), parameter :: phi
        real(kind = 8), parameter :: vrad
        real(kind = 8), parameter :: vphi
        real(kind = 8), parameter :: phidot
        real(kind = 8), parameter :: h

        ! Define the difference in positions and velocities
        rDiff = r(2) - r(1)
        vDiff = v(2) - v(1)

        ! Converting the coordinates from Cartesian to Cylindrical
        rr  = (rDiff(1)**2 + rDiff(2)**2)**0.5
        phi = atan2(rDiff(1), rDiff(2))

        ! Velocity
        vrad =  vDiff(1)*cos(phi) + vDiff(2)*sin(phi)
        vphi = -vDiff(1)*sin(phi) + vDiff(2)*cos(phi)

        ! Time derivative off phi
        phidot = vphi/rr

        ! Calculate the angular momentum
        h = (rr**2)*(phidot)

    end function computeAngularMomentum

    ! Computes the relative error in energy for a given time step.
    function computesigE(dt, m), result (ET_vec)

        real(kind = 8), parameter :: dt
        real(kind = 8), parameter :: m(1:2)
        real(kind = 8), parameter :: pos_vec
        real(kind = 8), parameter :: r
        real(kind = 8), parameter :: v
        real(kind = 8), parameter :: t
        real(kind = 8), parameter :: a
        real(kind = 8), parameter :: E
        real(kind = 8), parameter :: sig_E
        real, allocatable :: E(:)
        integer, :: i
        integer, :: j
                
        ! Define type vectors
        type ET_vec
            real(kind = 8), dimension(2, 3) :: t
            real(kind = 8), dimension(2, 3) :: sig_E
        end type ET_vec

        ! Declare the initial position and velocity for both the Sun and the Earth
        pos_vec = call initialise()
        r = pos_vec(1)
        v = pos_vec(2)

        ! Loop range
        real :: t(0:100/dt), dt
        t = spread(0, dt, 100)

        i = 100/dt

        do 1, i

            do j = 1, 2

                ! Get the current acceleration of object j 
                a = call computeAcc(m(j), r, j)

                ! Update the position and velocity
                r(j), v(j) = call updatePosVel(r(j), v(j), a, dt)

            end do

            ! Compute the total energy of the system
            allocate(E(i))
            E(i) = call computeEnergy(m, r, v)
            E = reshape([E, E(i)], [size(E) + 1])

        end do
            
        ! Compute the relative errors in the energy
        sig_E = abs((E - E(1))/E(1))

    end function computesigE

    ! Computes the relative error in energy for a given time step.
    function sigE_1orbit(dt, m), result (sig_E)
        
        real(kind = 8), parameter :: dt
        real(kind = 8), parameter :: m(1:2)
        real(kind = 8), parameter :: pos_vec
        real(kind = 8), parameter :: r
        real(kind = 8), parameter :: v
        real(kind = 8), parameter :: t
        real(kind = 8), parameter :: i
        real(kind = 8), parameter :: a
        real(kind = 8), parameter :: E
        real(kind = 8), parameter :: E0
        real(kind = 8), parameter :: sig_E
        integer, :: j
        real, allocatable :: E(:)

        ! Declare the initial position and velocity for both the Sun and the Earth
        pos_vec = call initialise()
        r = pos_vec(1)
        v = pos_vec(2)

        ! Compute the initial energy
        E0 = call computeEnergy(m, r, v)

        ! Loop range
        real :: t(0:100/dt), dt
        t = spread(0, dt, 100)
        
        integer :: t
        integer :: j

        ! Loop through the time steps and the masses (Earth and Sun)
        do i = t

            do j = 1, 2

                ! Get the current acceleration of object j 
                a = call computeAcc(m(j), r, j)

                ! Update the position and velocity
                r(j), v(j) = call updatePosVel(r(j), v(j), a, dt)

            end do

            if i > 2*pi then

            ! Compute the total energy of the system
            E = call computeEnergy(m, r, v)
            
            ! Compute the relative errors in the energy
            sig_E = abs((E - E0)/E0)

        end do

    end function sigE_1orbit

end program Euler_Method

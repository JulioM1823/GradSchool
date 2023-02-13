program Euler_Method

    ! Declare that no variables be implicit
    implicit none

    ! Initialize mass vector and number of orbits to be calculated
    real(kind = 8), parameter :: m(1:2) = [1.0d0, 1.0d-3]

    !Function to initialise the position and velocity vectors.
    function initialize() 

        ! Define initial positions of Sun and Earth
        real(kind = 8), parameter :: r_S, v_S, r_E, v_E
        real(kind = 8), dimension(2, 3) :: r, v
                
        ! Initialize the position and velocity vectors
        r_E = [1.0d0, 0.0d0, 0.0d0]
        r_S = [0.0d0, 0.0d0, 0.0d0]
        v_S = [0.0d0, 0.0d0, 0.0d0]
        v_E = [0.0d0, 1.0d0, 0.0d0]
        r = [r_S, r_E]
        v = [v_S, v_E]

        return r, v

    end function initialize
     
    ! Function to accept the mass, the current position index, and the position array (for all the objects in the system) as arguments and return the acceleration.
    function computeAcc(m, r, j, unitVec)

        ! Declare variables
        real(kind = 8), parameter :: m(1:2)
        real(kind = 8), parameter :: r(1:3)
        integer, parameter :: j
        logical :: unitVec

        ! Compute the difference in position    
        rDiff = r(j) - np.append(r(:j), r(j+1:))
        !! Needs to be properly translated
        r_new(1:j) = r(1:j)
        r_new(j+1:j+n-1) = r(j+1:n)

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

end program Euler_Method

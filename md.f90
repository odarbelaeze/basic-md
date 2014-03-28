program md

    implicit none
    
    ! parameters
    double precision, parameter :: pi = acos( - 1.0)

    ! Variables
    double precision, allocatable, dimension(:,:) :: r, rph
    double precision, allocatable, dimension(:,:) :: f, fph
    double precision, allocatable, dimension(:,:) :: v, vph
    double precision, allocatable, dimension(:,:) :: pbc

    integer :: l, n
    double precision :: a, h

    ! Counters and auxiliars
    integer :: counter, i, j, k

    ! Settings
    h = 0.004
    a = 1.2d0
    l = 3
    n = 4 * l ** 3

    ! Allocate required memory
    allocate(r(3, n))
    allocate(rph(3, n))
    allocate(f(3, n))
    allocate(fph(3, n))
    allocate(v(3, n))
    allocate(vph(3, n))
    allocate(pbc(3, 3 ** 3))

    ! Build the sample

    counter = 1
    do i = 0, l - 1, 1
        do j = 0, l - 1, 1
            do k = 0, l - 1, 1
                ! FCC lattice basis (0,0,0) (0.5, 0, 0.5)
                r(:, counter) = (/ i, j, k /) + (/ 0.0, 0.0, 0.0 /); counter = counter + 1
                r(:, counter) = (/ i, j, k /) + (/ 0.5, 0.5, 0.0 /); counter = counter + 1
                r(:, counter) = (/ i, j, k /) + (/ 0.0, 0.5, 0.5 /); counter = counter + 1
                r(:, counter) = (/ i, j, k /) + (/ 0.5, 0.0, 0.5 /); counter = counter + 1
            end do
        end do
    end do

    r = r * a

    ! Set the periodic boundary conditions

    counter = 1
    do i = - 1, 1, 1
        do j = - 1, 1, 1
            do k = - 1, 1, 1
               pbc(:, counter) = (/ i, j, k /); counter = counter + 1
            end do
        end do
    end do

    pbc = pbc * l * a

    ! Set the speeds

    call random_number(v)
    ! We use vph as an auxiliar but it will be set to its proper value later
    call random_number(vph)
    v = ( -2.0 * log(v)) ** 0.5 * sin(2 * pi * vph)

    ! We need now to scale speeds according to the temperature
    ! and keep the particles within the domain
    ! * pacmanEffect
    ! * centerSpeeds
    ! * temperatureControl

    call centerSpeeds(v)
    call temperatureControl(v, 1.0d0)

    ! We need now to compute the forces given the positions and periodic boundary 
    ! conditions
    ! * computeForces

    call computeForces(f, r, pbc)

    ! Let's initialize the vectors fph, rph and vph

    rph = r + h * v + 0.5 * h ** 2 * f
    call pacmanEffect(rph, l * a)
    call computeForces(fph, rph, pbc)
    vph = v + 0.5 * h * (fph + f)

    ! Now the updating loop can be made
    
    ! Relaxing
    do i = 1, 2500, 1
        r = rph
        v = vph
        f = fph

        rph = r + h * v + 0.5 * h ** 2 * f
        call pacmanEffect(rph, l * a)
        call computeForces(fph, rph, pbc)
        vph = v + 0.5 * h * (fph + f)

        if ( mod(i, 20) == 0 ) then
            call centerSpeeds(vph)
            call temperatureControl(vph, 1.0d0)
        end if

    end do

    ! Sampling
    do i = 1, 5000, 1
        r = rph
        v = vph
        f = fph

        rph = r + h * v + 0.5 * h ** 2 * f
        call pacmanEffect(rph, l * a)
        call computeForces(fph, rph, pbc)
        vph = v + 0.5 * h * (fph + f)

        if ( mod(i, 20) == 0 ) then
            call centerSpeeds(vph)
            call temperatureControl(vph, 1.0d0)
        end if

        print*, i * h, kineticEnergy(v), potentialEnergy(r, pbc)

    end do

    ! Subprograms

    contains

    subroutine pacmanEffect(rc, lac)
        double precision, allocatable, dimension(:,:), intent(inout) :: rc
        double precision, intent(in) :: lac

        integer :: ii, jj

        do ii = 1, size(rc, dim = 1), 1
            do jj = 1, size(rc, dim = 2), 1
                do while ( rc(ii, jj) < 0 )
                    rc(ii, jj) = rc(ii, jj) + lac
                end do
                do while ( rc(ii, jj) >= lac )
                    rc(ii, jj) = rc(ii, jj) - lac
                end do
            end do
        end do

    end subroutine pacmanEffect

    subroutine centerSpeeds(vc)
        double precision, allocatable, dimension(:,:), intent(inout) :: vc

        double precision, dimension(size(vc, dim = 1)) :: vmean
        integer :: ii

        vmean = sum(vc, dim = 2) / size(vc, dim = 2)

        do ii = 1, size(vc, dim = 2), 1
            vc(:,ii) = vc(:,ii) - vmean
        end do

    end subroutine centerSpeeds

    subroutine temperatureControl(vc, tdc)
        double precision, allocatable, dimension(:,:), intent(inout) :: vc
        double precision, intent(in) :: tdc

        double precision :: lambda

        lambda = sqrt((size(vc, dim = 2) - 1) * tdc / sum(vc ** 2))
        vc = lambda * vc
        
    end subroutine temperatureControl

    subroutine computeForces(fc, rc, pbcc)
        double precision, allocatable, dimension(:,:), intent(inout) :: fc
        double precision, allocatable, dimension(:,:), intent(in) :: rc
        double precision, allocatable, dimension(:,:), intent(in) :: pbcc

        double precision :: dd
        integer :: ii, jj, mu

        fc = 0.0

        do ii = 1, size(rc, dim = 2), 1
            do jj = 1, size(rc, dim = 2), 1
                do mu = 1, size(pbcc, dim = 2), 1
                    dd = sqrt(sum((r(:, ii) - (r(:, jj) + pbcc(:, mu))) ** 2))
                    if ( dd > 0.0 ) then
                        fc(:, ii) = fc(:, ii) +  12.0 * (r(:, ii) - (r(:, jj) + pbcc(:, mu))) &
                                    * (dd ** ( - 14) - dd ** ( - 8))
                    end if
                end do
            end do
        end do
        
    end subroutine computeForces

    function potentialEnergy(rc, pbcc)
        double precision, allocatable, dimension(:,:), intent(in) :: rc
        double precision, allocatable, dimension(:,:), intent(in) :: pbcc
        double precision :: potentialEnergy

        double precision :: dd
        integer :: ii, jj, mu

        potentialEnergy = 0.0

        do ii = 1, size(rc, dim = 2), 1
            do jj = ii, size(rc, dim = 2), 1
                do mu = 1, size(pbcc, dim = 2), 1
                    dd = sqrt(sum((r(:, ii) - (r(:, jj) + pbcc(:, mu))) ** 2))
                    if ( dd > 0.0 ) then
                        potentialEnergy = potentialEnergy + (dd ** ( - 12) - 2.0 * dd ** ( - 6))
                    end if
                end do
            end do
        end do

    end function potentialEnergy

    function kineticEnergy(vc)
        double precision, allocatable, dimension(:,:), intent(in) :: vc 
        double precision :: kineticEnergy

        kineticEnergy = 0.5 * sum(vc ** 2)

    end function kineticEnergy

end program md

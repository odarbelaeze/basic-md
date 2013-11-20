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
    double precision :: a

    ! Counters and auxiliars
    integer :: counter, i, j, k

    ! Settings
    a = 1.56d0
    l = 3
    n = 4 * l ** 3

    ! Allocate required memory
    allocate(r(3, n))
    allocate(rph(3, n))
    allocate(f(3, n))
    allocate(fph(3, n))
    allocate(v(3, n))
    allocate(vph(3, n))

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

    ! checks
    
    do i = 1, n
        print*, r(:,i), v(:,i)
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

end program md

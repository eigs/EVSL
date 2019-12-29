program driver

    implicit none

    ! Variable declarations
    ! n - The number of rows in the matrix
    ! nx, ny, nz - The dimension of the mesh to generate the laplacean from
    ! nslices - The number of slices to divide the spectrum into
    ! Mdeg - Polynomial degree
    ! nvec - Number of sample vectors
    ! ev_int - Number of eigenvalues per slice
    ! mlan - Dimension of krylov subspace
    ! nev - Approximate number of eigenvalues wanted
    ! max_it - Maximum number of iterations before a restart occurs
    integer :: n, nx, ny, nz, nslices, Mdeg, nvec, ev_int, mlan, nev, max_it

    ! Find the eigenvalues in the interval [a, b]
    ! a - lower bound of the interval
    ! b - upper bound of the interval
    ! lmax - largest eigenvalue
    ! lmin - smallest eigenvalue
    ! thresh_int, thresh_ext - Polynomial variables
    double precision :: a, b, lmax, lmin, tol, thresh_int, thresh_ext
    ! sli - The spectrum slices
    double precision, dimension(:), pointer :: sli
    double precision, dimension(4) :: xintv

    double precision, dimension(:), pointer :: eigval, eigvec
    ! Matrix peices
    double precision, dimension(:), pointer :: vals
    integer, dimension(:), pointer :: ia
    integer, dimension(:), pointer :: ja

    integer*8 :: pol, csr

    ! Loop varialbe declarations
    integer :: i, j, k

    ! DEBUG Variables
    integer :: s, e

    ! Variables for reading command line arguments
    character (len=32) :: arg
    integer :: readerr

    ! Variables for using five point gen from SPARSKIT
    double precision, dimension(:), pointer :: rhs
    double precision, dimension(6) :: al
    integer, dimension(:), pointer :: iau
    integer :: mode

    ! Read in command line arguments to set important values.
    ! The user can pass the phrase help to get the program to
    ! print a usage statement then terminate

    !Set default values
    nx = 16
    ny = 16
    nz = 20
    a = .40D0
    b = .80D0
    nslices = 4

    !Crude but works.
    ! This loop and if statements process the command line arguments.
    ! Type ./LapPLanN.out help for usage information
    do i = 1, iargc()
        call getarg(i, arg)
        arg = trim(arg)

        if(arg(1:2) == 'nx') then
            read(arg(4:), *, iostat = readerr) nx
        elseif(arg(1:2) == 'ny') then
            read(arg(4:), *, iostat = readerr) ny
        elseif(arg(1:2) == 'nz') then
            read(arg(4:), *, iostat = readerr) nz
        elseif(arg(1:1) == 'a') then
            read(arg(3:), *, iostat = readerr) a
        elseif(arg(1:1) == 'b') then
            read(arg(3:), *, iostat = readerr) b
        elseif(arg(:7) == 'nslices') then
            read(arg(9:), *, iostat = readerr) nslices
        elseif(arg == 'help') then
            write(*,*) 'Usage: ./testL.ex nx=[int] ny=[int] nz=[int] a=[double] b=[double] nslices=[int]'
            stop
        endif
        if(readerr /= 0) then
            write(*,*) 'There was an error while reading argument: ', arg
            stop
        endif
    enddo

    ! Initialize eigenvalue bounds set by hand
    lmin = 0.0D0
    if(nz == 1) then
        lmax = 8.0D0
    else
        lmax = 12.0D0
    endif
    xintv(1) = a
    xintv(2) = b
    xintv(3) = lmin
    xintv(4) = lmax
    tol = 1e-08

    ! Use SPARSKIT to generate the 2D/3D Laplacean matrix
    ! We change the grid size to account for the boundaries that
    ! SPARSKIT uses but not used by the LapGen tests in EVSL
    nx = nx+2
    if(ny > 1) then
        ny = ny+2
        if(nz > 1) then
            nz = nz+2
        endif
    endif
    n = nx*ny*nz

    ! allocate our csr matrix
    allocate(vals(n*7)) !Size of number of nonzeros
    allocate(ja(n*7)) !Size of number of nonzeros
    allocate(ia(n+1)) !Size of number of rows + 1

    ! Allocate sparskit things
    allocate(rhs(n*7)) ! Righthand side of size n
    allocate(iau(n)) ! iau of size n
    al(1) = 0.0D0; al(2) = 0.0D0;
    al(3) = 0.0D0; al(4) = 0.0D0;
    al(5) = 0.0D0; al(6) = 0.0D0;
    mode = 0
    call gen57pt(nx,ny,nz,al,mode,n,vals,ja,ia,iau,rhs)

    ! Since we're using this array with C we need to accomodate for C indexing
    ia = ia - 1
    ja = ja - 1
    ! Cleanup extra sparskit information
    deallocate(rhs)
    deallocate(iau)

    ! This section of the code will run the EVSL code.
    ! This file is not utilizing the matrix free format and we'll pass
    ! a CSR matrix in

    ! Initialize the EVSL global data
    call EVSL_START_F90()

    ! Set the A matrix in EVSL global data to point to the arrays built here
    call EVSL_ARR2DEVICECSR_F90(n, ia, ja, vals, csr)

    call EVSL_SETA_DEVICECSR_F90(csr)

    ! kmpdos in EVSL for the DOS for dividing the spectrum
    ! Set up necessary variabls for kpmdos
    Mdeg = 300;
    nvec = 60;
    allocate(sli(nslices+1))
    ! Call EVSL kpmdos and spslicer
    call EVSL_KPM_SPSLICER_F90(Mdeg, nvec, xintv, nslices, sli, ev_int)

    ! For each slice call ChebLatr
    do i = 1, nslices
        ! Prepare parameters for this slice
        xintv(1) = sli(i)
        xintv(2) = sli(i+1)
        thresh_int = .5
        thresh_ext = .15

        ! Call the EVSL function to create the polynomial
        call EVSL_FIND_POL_F90(xintv, thresh_int, thresh_ext, pol)

        ! Necessary paramters
        nev = ev_int + 2
        mlan = max(4*nev, 100)
        mlan = min(mlan, n)
        max_it = 3*mlan

        ! Call the EVSL cheblannr function to find the eigenvalues in the slice
        call EVSL_CHEBLANTR_F90(mlan, nev, xintv, max_it, tol, pol)

        ! Extract the number of eigenvalues found from the EVSL global data
        call EVSL_GET_NEV_F90(nev)

        ! Allocate storage for the eigenvalue and vectors found from cheblannr
        allocate(eigval(nev))
        allocate(eigvec(nev*size(ia))) ! number of eigen values * number of rows

        ! Extract the arrays of eigenvalues and eigenvectors from the EVSL global data
        call EVSL_COPY_RESULT_F90(eigval, eigvec)
        write(*,*) nev, ' Eigs in this slice'

        ! Here you can do something with the eigenvalues that were found
        ! The eigenvalues are stored in eigval and eigenvectors are in eigvec

        ! Be sure to deallocate the polynomial stored by EVSL
        call EVSL_FREE_POL_F90(pol)
        deallocate(eigval)
        deallocate(eigvec)
    enddo
    deallocate(sli)

    call EVSL_FREE_DEVICECSR_F90(csr)

    call EVSL_FINISH_F90()

    ! Necessary Cleanup
    deallocate(vals)
    deallocate(ja)
    deallocate(ia)
end program driver

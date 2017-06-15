module csrmatrix

    ! This module contains a simple CSR matrix storage type and an associated
    ! matvec routine to accompany it.

    ! This is used to demonstrate usage of the matrix free version of EVSL from
    ! Fortran

    type csrmat
        integer :: nrows ! Number of rows in the matrix
        integer :: ncols ! Number of cols in the matrix
        integer :: nnz ! Number of nonzero entries in the matrix

        integer, dimension(:), pointer :: ia ! The row map
        integer, dimension(:), pointer :: ja ! The column indices
        double precision, dimension(:), pointer :: a ! The values in the matrix
    end type csrmat

contains
    ! This will perform a matrix vector multiplication Ax = y
    subroutine csrmatvec(x, y, mat)
        ! Determine the inputs
        type(csrmat) :: mat
        double precision, dimension(mat%nrows), intent(in) :: x
        double precision, dimension(mat%nrows), intent(out) :: y
        
        ! Information to parse the matrix
        integer :: row, col, i, j, s, e

        ! For safety, parse the y vector and set entries to 0
        do i = 1, mat % nrows
            y(i) = 0
        enddo

        ! For each row perform the matvec
        do i = 1, mat % nrows
            s = mat%ia(i)
            e = mat%ia(i+1)
            do j = s, e-1
                col = mat%ja(j)
                y(i) = y(i) + mat%a(j) * x(col)
            enddo
        enddo
    end subroutine csrmatvec
end module csrmatrix

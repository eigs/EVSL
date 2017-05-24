module csrmatrix

    type csrmat
        integer :: nrows
        integer :: ncols
        integer :: nnz

        integer, dimension(:), pointer :: ia
        integer, dimension(:), pointer :: ja
        double precision, dimension(:), pointer :: a
    end type csrmat

contains
    subroutine csrmatvec(x, y, mat)
        type(csrmat) :: mat
        double precision, dimension(mat%nrows), intent(in) :: x
        double precision, dimension(mat%nrows), intent(out) :: y

        integer :: row, col, i, j, s, e
        do i = 1, mat % nrows
            y(i) = 0
        enddo

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

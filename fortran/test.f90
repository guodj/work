program test
    implicit none ! this means variables must be declared
    ! variables declaration
    integer :: i, a(10), b(10)  ! learn to use array
    character(10) :: c  ! learn to use string
    a = (/(i, i=1,10,1)/)
    b = (/(i, i=1,10,1)/)
    c = 'Fortran'
    write(*, '(10i3, /, 10i3, /, A10)') a, b, c
end program test

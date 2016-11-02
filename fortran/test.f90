program test
    implicit none ! this means variables must be declared
    ! variables declaration
    integer, dimension(2) :: a, b
    common a ! a is a global variable
    external dat  ! some variables will be initialized in dat block data
    data b/3,4/  ! initialize variables (f77)
    write(*, *) a
    write(*, *) b
end program test

block data dat
    implicit none
    integer, dimension(2) :: a
    common a
    data a/1,2/
end

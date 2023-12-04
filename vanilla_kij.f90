program vanilla_jki
    use matrixop
    use timings
    implicit none

    integer, parameter                 :: N = 2000
    real(kind=dp)                      :: A(N,N), B(N,N), C(N,N)
    integer, dimension(:), allocatable :: seed
    integer                            :: k

    call random_seed(size=k)
    allocate(seed(k))
    seed = N
    call random_seed(put=seed)

    call random_number(A)
    call random_number(B)
    call tic()
    call mm_kij(A,B,C)
    call toc()

    deallocate(seed)

end program

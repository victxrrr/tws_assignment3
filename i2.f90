program i2
    use matrixop
    use timings
    implicit none

    abstract interface
        subroutine mm_interface(a, b, c)
            import dp
            real(kind=dp), dimension(:,:), intent(in)  :: a, b
            real(kind=dp), dimension(:,:), intent(out) :: c
        end subroutine mm_interface
        subroutine mm_blocks_interface(a,b,c,blocksize)
            import dp
            real(kind=dp), dimension(:,:), intent(in)  :: a, b
            real(kind=dp), dimension(:,:), intent(out) :: c
            integer, intent(in) :: blocksize
        end subroutine mm_blocks_interface
    end interface

    integer, parameter :: N_warmup = 2000, blocksize = 100
    integer, dimension(:), allocatable :: seed
    real(kind=dp), allocatable, dimension(:,:) :: A_warmup,&
         B_warmup,C_warmup, A, B, C
    integer :: statusAllocate
    integer :: N, k
    character(len=32) :: FIG
    
    ! Compute the matrix matrix multiplication of two large matrices to avoid warmup effects
    allocate(A_warmup(N_warmup,N_warmup),B_warmup(N_warmup,N_warmup),C_warmup(N_warmup,N_warmup),STAT=statusAllocate)
    if (statusAllocate>0) then
        print *, "Failed to allocate warmup matrices. Terminating program."
        call exit;
    endif
    call random_number(A_warmup)
    call random_number(B_warmup)
    C_warmup = matmul(A_warmup,B_warmup)
    C_warmup = matmul(A_warmup,B_warmup)
    deallocate(A_warmup,B_warmup,C_warmup)

    ! Add your implementation here. 
    ! Do not forget to call random_number for each matrix size.

    call random_seed(size=k)
    allocate(seed(k))
    seed = N
    call random_seed(put=seed)

    if (command_argument_count() /= 1) then 
        print *, 'Wrong number of arguments'
        print *, 'USAGE : ./main [f1|f2]'
        stop
    endif
    call get_command_argument(1, FIG)

    if (trim(FIG) == 'f1') then
        do N = 10, 100, 10
            allocate(A(N,N), B(N,N), C(N,N), STAT=statusAllocate)
            if (statusAllocate > 0) then
                print *, "Failed to allocated matrices. Terminating program"
                stop statusAllocate
            endif
            call random_number(A)
            call random_number(B)
            call measure_flops( 'flops_f1.csv' , slow_method=mm_kij , fast_method=mm_jki )
            deallocate(A,B,C)
        enddo

        do N = 200, 1600, 100
            allocate(A(N,N), B(N,N), C(N,N), STAT=statusAllocate)
            if (statusAllocate > 0) then
                print *, "Failed to allocated matrices. Terminating program"
                stop statusAllocate
            endif
            call random_number(A)
            call random_number(B)
            call measure_flops( 'flops_f1.csv' , slow_method=mm_kij , fast_method=mm_jki )
            deallocate(A,B,C)
        enddo
    else if (trim(FIG) == 'f2') then
        do N = 100, 1600, 100
            allocate(A(N,N), B(N,N), C(N,N), STAT=statusAllocate)
            if (statusAllocate > 0) then
                print *, "Failed to allocated matrices. Terminating program"
                stop statusAllocate
            endif
            call random_number(A)
            call random_number(B)
            call measure_flops( 'flops_f2.csv' , block_method_a=mm_blocks_a , block_method_b=mm_blocks_b )
            deallocate(A,B,C)
        enddo 
    else 
        print *, 'Unknown command line argument' 
        stop
    end if

    deallocate(seed)

contains

    subroutine measure_flops(fname, slow_method, fast_method, block_method_a, block_method_b)
        procedure(mm_interface), optional :: slow_method
        procedure(mm_interface), optional :: fast_method
        procedure(mm_blocks_interface), optional :: block_method_a
        procedure(mm_blocks_interface), optional :: block_method_b
        character(*), intent(in) :: fname
        real(kind(0.e0)) :: slow_elapsedTime, fast_elapsedTime
        real(kind(0.e0)) :: rN
        rN = real(N, kind(0.e0))

        ! Do the timing
        if (present(slow_method) .and. present(fast_method)) then
            call tic()
            call slow_method(A, B, C)
            call toc(slow_elapsedTime)
            call tic()
            call fast_method(A, B, C)
            call toc(fast_elapsedTime)
        else if (present(block_method_a) .and. present(block_method_b)) then
            call tic()
            call block_method_a(A, B, C, blocksize)
            call toc(slow_elapsedTime)
            call tic()
            call block_method_b(A, B, C, blocksize)
            call toc(fast_elapsedTime)
        else
            print *, 'Invalid arguments for measure_flops'
            stop
        endif

        open(7, file=trim(fname), access='append')
        write(*, '(i4, f12.5, f14.6)') N, 2*rN**3/(slow_elapsedTime*10**6), 2*rN**3/(fast_elapsedTime*10**6)
        write(7, '(i4, A, f12.5, A, f14.6)') N, ',', 2*rN**3/(slow_elapsedTime*10**6), ',', 2*rN**3/(fast_elapsedTime*10**6)
        close(7)
      
    end subroutine measure_flops

        
end program i2




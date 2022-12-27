PROGRAM main
    implicit none
 
    double precision, dimension (:, :), allocatable :: A
    double precision, dimension (:), allocatable :: tau, work, l_work
    integer :: n, i, j, info, LDA, lwork
    real ::  start, end
    character(len=32) :: arg


    n = 1000
    LDA = n
    lwork = -1
    
    CALL get_command_argument(1, arg)
    IF (LEN_TRIM(arg) /= 0) read(arg, *) n

    allocate (A(n, n), tau(n), l_work(n))
    call random_number(A)

    call dgeqrf(n, n, A, LDA, tau, l_work, lwork, info)
    lwork = l_work(1)

    allocate(work(lwork))
    call cpu_time(start)    

    call dgeqrf(n, n, A, LDA, tau, work, lwork, info)

    call cpu_time(end)
    write (*,*) "Duration (sec): "
    print *, end - start

    deallocate (A, tau, work, l_work)

END PROGRAM main
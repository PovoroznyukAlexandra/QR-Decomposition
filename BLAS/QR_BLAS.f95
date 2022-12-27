PROGRAM main
    implicit none
    integer :: n, i, j
    DOUBLE PRECISION, dimension (:, :), allocatable :: A, Q, R
    DOUBLE PRECISION :: Mod, c, s
    
    write(*, *) 'Enter N: '
    read(*, *) n

    allocate (A(n, n), Q(n, n), R(n, n))
    
    ! Исходная матрица
    call random_number(A)
    !write(*, *) A
    
    
    ! Унитарная матрица Q
    do i = 1, n
        do j = 1, n
            Q(j, i) = 0
        end do
        Q(i, i) = 1
    end do
    !write(*, *) Q
    
    ! Правая треугольная матрца R
    do i = 1, n
        do j = 1, n
            R(j, i) = A(j, i)
        end do
    end do
    !write(*, *) R
    
    ! Вращение матриц
    do i = 1, (n - 1)
        do j  = (i + 1), n
        Mod = sqrt(R(i, i)*R(i, i) + R(j, i)*R(j, i))
        c = R(i, i) / Mod
        s = -R(j, i) / Mod
        call drot(n, Q(:, i), 1, Q(:, j), 1, c, -s)
        call drot(n, R(:, i), 1, R(:, j), 1, c, -s)
    end do
    end do
    
    call Transpose(n, Q)
    
    ! Проверка того, что R - правая треугольная матрица
    !Check_Right_Triang(n, R);

    ! Проверка унитарности Q
    !Check_Unitary(n, Q);

    ! Подсчет ошибки вычисления
    call Check_Accuracy(n, A, Q, R);

    deallocate (A, Q, R)

END PROGRAM main


! Вращение матрицы
subroutine Rotation(n, i, j, c, s, A)
    implicit none
    integer :: k
    integer, INTENT(IN) :: n, i, j
    DOUBLE PRECISION, INTENT(IN) :: c, s
    DOUBLE PRECISION, INTENT(OUT), dimension(n, n) :: A
    DOUBLE PRECISION :: Aik, Ajk
    
    do k = 1, n
        Aik = c * A(i, k) - s * A(j, k)
        Ajk = s * A(i, k) + c * A(j, k)
        A(i, k) = Aik
        A(j, k) = Ajk
    end do
end subroutine Rotation


! Транспонирование матрицы Q
subroutine Transpose(n, Q)
    implicit none
    integer :: i, j, k
    double precision t
    integer, INTENT(IN) :: n
    DOUBLE PRECISION, INTENT(OUT), dimension(n, n) :: Q
    
    do i = 1, n
        do j = 1, n
            t = Q(i, j);
            Q(i, j) = Q(j, i);
            Q(j, i) = t;
        end do
    end do
end subroutine Transpose


subroutine Check_Accuracy(n, A, Q, R)
    implicit none
    integer :: i, j, k
    double precision tmp
    integer, INTENT(IN) :: n
    DOUBLE PRECISION, INTENT(OUT), dimension(n, n) :: A, Q, R
    
    do i = 1, n
        do j = 1, n
            tmp = A(i, j)
            do k = 1, n
                tmp = tmp - Q(i, k) * R(k, j);
            end do
            if (abs(tmp) > 1e-10) EXIT
        end do
    end do
    write (*, *) 'QR decomposition is correct'
end subroutine Check_Accuracy



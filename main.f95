program main
  implicit none
  integer i, j
  double precision x

  ! Матрица
  double precision, dimension(7, 7) :: A, U, D, V, AInv, UDV1, UDV

  ! Интерфейсы процедур
  interface
    subroutine SVD(A, U, D, V)
        double precision, dimension(:,:) :: A, U, D, V
    end subroutine SVD
    subroutine InvSVD(AInv, U, D, VT)
        double precision, dimension(:,:) :: AInv, U, D, VT
    end subroutine InvSVD
  end interface

  ! Матрица для расчета
  a(1,1) = 1.0; a(1,2) = 7.0; a(1,3) = 5.0; a(1,4) = 4.0; a(1,5) = 6.0; a(1,6) = 3.0; a(1,7) = 2.0;
  a(2,1) = 8.0; a(2,2) = 9.0; a(2,3) = 2.0; a(2,4) = 3.0; a(2,5) = 4.0; a(2,6) = 5.0; a(2,7) = 7.0;
  a(3,1) = 5.0; a(3,2) = 4.0; a(3,3) = 7.0; a(3,4) = 3.0; a(3,5) = 5.0; a(3,6) = 6.0; a(3,7) = 2.0;
  a(4,1) = 6.0; a(4,2) = 4.0; a(4,3) = 7.0; a(4,4) = 8.0; a(4,5) = 2.0; a(4,6) = 3.0; a(4,7) = 9.0;
  a(5,1) = 4.0; a(5,2) = 2.0; a(5,3) = 9.0; a(5,4) = 3.0; a(5,5) = 6.0; a(5,6) = 4.0; a(5,7) = 7.0;
  a(6,1) = 3.0; a(6,2) = 2.0; a(6,3) = 4.0; a(6,4) = 8.0; a(6,5) = 6.0; a(6,6) = 2.0; a(6,7) = 6.0;
  a(7,1) = 8.0; a(7,2) = 1.0; a(7,3) = 6.0; a(7,4) = 4.0; a(7,5) = 5.0; a(7,6) = 3.0; a(7,7) = 3.0;

  ! Сингулярное разложение
  call SVD(A, U, D, V)

  write(*,*) "A"
  write(*,'(7(F3.1,X))') transpose(A)

  write(*,*) "D"
  write(*,'(7(F7.1,X))') transpose(D)
  write(*,*) "U"
  write(*,'(7(F7.1,X))') transpose(U)
  write(*,*) "V"
  write(*,'(7(F7.1,X))') (V)

  ! Вычисление обратной матрицы

  call InvSVD(AInv, U, D, V)
  write(*,*) "A inversed"
  write(*,'(7(F7.1,X))') transpose(AInv)


  UDV1=matmul(D,U);!V;
  V = transpose(V);
  UDV=matmul(UDV1,V)
  write(*,*) "UDV inversed"
  write(*,'(7(F7.1,X))') transpose(UDV)

end program main

! Единичная матрицы
subroutine ident(A)
    implicit none
    double precision, dimension(:, :) :: A

    integer i, j, N
    N = size(A,1)
    do i=1,N
        do j=1,N
            if (i == j) then
                A(i,j)=1.0
            else
                A(i,j)=0.0
            end if
        end do
    end do
end subroutine ident

! Квадрат евклидовой нормы вектора
function n2(v) result (n)
    implicit none
    double precision, dimension(:) :: v
    double precision :: n
    integer i

    do i=1, size(v)
        n = n + v(i) * v(i)
    end do
end function n2

! Сингулярное разложение
subroutine SVD(A0, U0, D0, V0)
    implicit none

    double precision, dimension(:,:) :: A0, U0, D0, V0
    integer N, M
    double precision, dimension(:,:,:), allocatable :: A, B, H, G
    double precision, dimension(:,:), allocatable :: v, vt, w, wt;
    double precision, dimension(:,:), allocatable :: Tv, Tw
    double precision s, norm
    integer i, j, k, r

    interface
      subroutine ident(a)
        double precision, dimension(:,:) :: a
      end subroutine ident

      double precision function n2(v)
        double precision, dimension(:) :: v
      end function n2
    end interface

    N = size(A0, 1);    M = size(A0, 2);    ! Размерность исходной матрицы
    allocate(A(0:M-1,N,M))    ! Матрицы A
    allocate(B(M-1,N,M))    ! Матрицы B
    allocate(H(M-1,N,N))    ! Матрицы H
    allocate(G(M-1,M,M))    ! Матрицы G

    allocate(v(N,1))    ! Вектор v
    allocate(vt(1,N))    ! Транспонированный вектор v
    allocate(Tv(N,N))    ! Матрица для умножения v*vt
    allocate(w(M,1))    ! Вектор w
    allocate(wt(1,M))    ! Транспонированный вектор w
    allocate(Tw(M,M))    ! Матрица для умножения w*wt


    A(0,:,:) = A0

    do k=1,M-1
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Построение вектора v
        do r=1,N
            if ( r < k ) then
                v(r,1) = 0.0
            elseif ( r == k ) then
                s = 0.0
                do i=k,N
                    s = s + A(k-1,i,k)**2
                enddo
                if ( A(k-1, k, k) > 0 ) then
                    v(r,1) = A(k-1,k,k) + sqrt(s)
                else
                    v(r,1) = A(k-1,k,k) - sqrt(s)
                end if
            else
                v(r,1) = A(k-1,r,k)
            end if
        end do

        norm = n2(v(:,1))   ! Вычисляем норму
        if ( norm == 0.0 ) then
            exit
        end if

        vt = transpose(v)   ! Транспонированный v
        Tv = matmul(v, vt)  ! Умножаем v*vt

        call ident(H(k,:,:))    ! Единичная матрица
        H(k,:,:) = H(k,:,:)-(2.0/norm)*Tv   ! Матрица H(k)

        B(k,:,:) = matmul(H(k,:,:), A(k-1,:,:))

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Построение вектора w
        do r=1,M
            if ( r <= k ) then
                w(r,1) = 0.0
            elseif ( r == (k+1) ) then
                s = 0.0
                do j=k+1,M
                    s = s + B(k,k,j)**2
                enddo
                if ( B(k, k, k+1) > 0 ) then
                    w(r,1) = B(k,k,k+1) + sqrt(s)
                else
                    w(r,1) = B(k,k,k+1) - sqrt(s)
                end if
            else
                w(r,1) = B(k,k,r)
            end if
        end do

        norm = n2(w(:,1))   ! Вычисляем норму
        if ( norm == 0.0 ) then
            exit
        end if

        wt = transpose(w)   ! Транспонированный w
        Tw = matmul(w, wt)  ! Умножаем w*wt

        call ident(G(k,:,:))    ! Единичная матрица
        G(k,:,:) = G(k,:,:)-(2.0/norm)*Tw   ! Матрица G(k)

        A(k,:,:) = matmul(B(k,:,:), G(k,:,:))
    end do

    D0 = A(k-1,:,:)

    call Ident(U0)
    call Ident(V0)

    do k=1,M-2
        U0 = matmul(U0, H(k,:,:))   ! Накопление матрицы U
        V0 = matmul(V0, G(k,:,:))   ! Накопление матрицы Vt
    end do
    U0 = matmul(U0, H(M-1,:,:))

    ! Очистка памяти
    deallocate (A);    deallocate (B);    deallocate (G);    deallocate (H)
    deallocate (v);    deallocate (vt);   deallocate (Tv);
    deallocate (w);    deallocate (wt);   deallocate (Tw);
end subroutine SVD

subroutine InvSVD(AInv, U, D, VT)
    implicit none
    double precision, dimension(:,:) :: AInv, U, D, VT
    double precision, dimension(size(D,1),size(D,2)) :: UT, V, DInv

    UT = transpose(U)
    V = transpose(VT)
    call Inverse(D, DInv, size(D, 1))

    AInv = matmul(matmul(V, DInv), UT)
end subroutine InvSVD



! Обращение матрицы
subroutine Inverse(a,c,n)
    implicit none
    integer n
    double precision a(n,n), c(n,n)
    double precision L(n,n), U(n,n), b(n), d(n), x(n)
    double precision coeff
    integer i, j, k

    L = 0.0
    U = 0.0
    b = 0.0

    do k = 1,n-1
        do i = k+1,n
            coeff = a(i,k)/a(k,k)
            L(i,k) = coeff
            do j = k+1,n
                a(i,j) = a(i,j)-coeff*a(k,j)
            end do
        end do
    end do

    do i = 1,n
        L(i,i) = 1.0
    end do

    do j=1,n
        do i=1,j
            U(i,j) = a(i,j)
        end do
    end do

    do k = 1,n
        b(k) = 1.0
        d(1) = b(1)

        do i = 2,n
            d(i)=b(i)
            do j=1,i-1
                d(i) = d(i) - L(i,j)*d(j)
            end do
        end do

        x(n)=d(n)/U(n,n)

        do i = n-1,1,-1
            x(i) = d(i)
            do j=n,i+1,-1
                x(i)=x(i)-U(i,j)*x(j)
            end do

            x(i) = x(i)/u(i,i)
        end do
        do i=1,n
            c(i,k) = x(i)
        end do
        b(k)=0.0
    end do
end subroutine Inverse


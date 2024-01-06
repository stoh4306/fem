subroutine sparse_conj_gradient_sp(N, M, A_Ind, A, b, eps, x, err, iter)
    implicit none

    integer, intent(in) :: N, M
    integer, intent(in) :: A_Ind(N, M)
    real, intent(in)    :: A(N,M)
    real, intent(in)    :: b(N)
    real, intent(in)    :: eps
    real, intent(inout) :: x(N)
    real, intent(out)   :: err
    integer, intent(out) :: iter

    real        :: r(N), p(N), Ap(N)
    real        :: tau, mu, rho(2)
    integer     :: k = 0
    real        :: tv(N)

    call sparse_matvecmul(N, M, A_Ind, A, x, tv)

    r = b-tv
    rho(1) = dot_product(r, r)

    do 
        if ( rho(1) .lt. eps ) then
            exit
        end if

        k = k + 1

        if ( k .eq. 1 ) then
            p = r
        else
            tau = (rho(1)/rho(2))
            p = r + tau*p
        end if

        call sparse_matvecmul(N, M, A_Ind, A, p, Ap)
        mu = rho(1) / dot_product(p, Ap)
        x = x + mu * p
        r = r - mu * Ap
        rho(2) = rho(1)
        rho(1) = dot_product(r, r)

        if ( k .gt. 2*N ) then
            exit
        end if
    
    end do

    iter = k
    err = rho(1)

end subroutine

subroutine sparse_matvecmul(N, M, A_Ind, A, x, Ax)
    implicit none

    integer, intent(in) :: N, M
    integer, intent(in) :: A_Ind(N, M)
    real, intent(in) :: A(N, M)
    real, intent(in) :: x(N)
    real, intent(out) :: Ax(N)

    integer :: i, j

    do i = 1, N
        Ax(i) = 0.0
        do j = 2, A_Ind(i,1)+1
            Ax(i) = Ax(i) + A(i, j)*x(A_Ind(i,j))
        end do 
    end do
end subroutine
!---------------------------------------------------
! Compute K and b for Kx=b of Laplace equation
! Assume :
!   . Linear basis function : {1-(x+y), x, y}
!   . Standard triangle : {(0,0), (1,0), (0,1)}
!---------------------------------------------------
subroutine compute_fem_matrix_equation(domin, nsubd, grid_spacing, &
                nnode, node, ntriangle, triangle, f, g, &
                nmaxnonzeros, K_sparse_index, K_sparse, b)    
    implicit none
    
    real, intent(in) :: domin(2)
    integer, intent(in) :: nsubd(2)
    real, intent(in) :: grid_spacing(2)
    integer, intent(in) :: nnode, ntriangle
    real, intent(in) :: node(nnode, 2)
    integer, intent(in) :: triangle(ntriangle, 3)
    real, intent(in) :: f(nnode), g(nnode)
    integer, intent(in) :: nmaxnonzeros
    
    integer, intent(out) :: K_sparse_index(nnode, nmaxnonzeros)
    real, intent(out) :: K_sparse(nnode, nmaxnonzeros)
    real, intent(out) :: b(nnode)
    
    integer :: i
    real    :: v(3,2), tv(3, 2)
    real    :: grad_phi(3,2)
    
    ! Initialize K 
    K_sparse_index(:, 1) = 0 ! number of nonzero elements
    K_sparse = 0.0
    
    ! Iterate all triangles to compute the matrix K
    
    ! Set gradients of standard basis functions defined above
    grad_phi(1,:) = [-1.0, -1.0]
    grad_phi(2,:) = [ 1.0,  0.0]
    grad_phi(3,:) = [ 0.0,  1.0]
    
    do i = 1, ntriangle
        v(1, :) = node(triangle(i, 1), :)
        v(2, :) = node(triangle(i, 2), :)
        v(3, :) = node(triangle(i, 3), :)
        
        tv = v
        tv(1,:) = tv(1,:) - v(1,:)
        tv(2,:) = tv(2,:) - v(1,:)
        tv(3,:) = tv(3,:) - v(1,:)
        
        print *, tv(1,:)
        print *, tv(2,:)
        print *, tv(3,:)
        
    end do
    
    
end subroutine 
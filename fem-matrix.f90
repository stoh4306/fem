subroutine compute_fem_matrix_equation(domin, nsubd, grid_spacing, &
                nnode, node, ntriangle, triangle, f, g, &
                nmaxnonzeros, K_sparse_index, K_sparse, b)
    implicit none
    
    real, intent(in) :: domin(2)
    integer, intent(in) :: nsubd(2)
    real, intent(in) :: grid_spacing(2)
    integer, intent(in) :: nnode, ntriangle
    integer, intent(in) :: node(2, nnode), triangle(3, ntriangle)
    real, intent(in) :: f(nnode), g(nnode)
    integer, intent(in) :: nmaxnonzeros
    
    integer, intent(out) :: K_sparse_index(nnode, nmaxnonzeros)
    real, intent(out) :: K_sparse(nnode, nmaxnonzeros)
    real, intent(out) :: b(nnode)
    
    
end subroutine 
!  fem.f90 
!
!  FUNCTIONS:
!  fem - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: fem
!
!  PURPOSE:  Entry point for the console application.
!            -div(grad U) = f, U = g on boundary
!
!****************************************************************************

program fem
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    
    implicit none

    ! Variables
    real        :: domin(2), domax(2)
    integer     :: nsubd(2), grid_size(2)
    real        :: grid_spacing(2)
    
    integer, parameter      :: GRIDNX = 3, GRIDNY = 3
    integer, parameter      :: NNODE = (GRIDNX + 1) * (GRIDNY + 1)
    integer, parameter      :: NTRIANGLE = GRIDNX * GRIDNY * 2
    
    real        :: node(NNODE, 2)
    integer     :: triangle(NTRIANGLE, 3)
    
    real        :: f(NNODE), g(NNODE) ! -div(grad u) = f, u = g on boundary
    
    integer, parameter:: NMAXNONZEROS = 10
    integer     :: K_sparse_index(NNODE, NMAXNONZEROS)
    real        :: K_sparse(NNODE, NMAXNONZEROS), b(NNODE)
    
    ! Set domain
    domin = [0.0, 0.0]
    domax = [1.0, 1.0]
    nsubd = [GRIDNX, GRIDNY]
    
    grid_spacing = (domax-domin)/nsubd
    grid_size = nsubd + 1
    
    ! Set Mesh
    call set_mesh_on_rectangular_domain(domin, nsubd, grid_spacing, &
                nnode, node, ntriangle, triangle)
    
    call set_source_bdry_func(domin, nsubd, grid_spacing, nnode, node, f, g)
    
    ! Compute K-matrix and b
    call compute_fem_matrix_equation(domin, nsubd, grid_spacing, nnode, node, ntriangle, triangle, f, g, &
                NMAXNONZEROS, K_sparse_index, K_sparse, b)
    
    ! Body of fem
    print *, 'Domain : ', grid_spacing, grid_size, nnode, ntriangle
    print *, 'Nodes : ', node(2, :)
    print *, 'Triangles : ', triangle
    print *, 'f : ', f(1:nnode)
    print *, 'g : ', g(1:nnode)
    
end program fem
    
subroutine set_source_bdry_func(domin, nsubd, grid_spacing, nnode, node, f, g)
    implicit none
    real, intent(in) :: domin(2)
    integer, intent(in) :: nsubd(2)
    real, intent(in)    :: grid_spacing(2)
        
    integer, intent(in) :: nnode
    real, intent(in)    :: node(nnode, 2)
        
    real :: source_func, bdry_func
        
    real, intent(out)   :: f(nnode), g(nnode)
        
    integer     :: i
        
    do i = 1, nnode
        f(i) = source_func(node(i, 1), node(i, 2))
        g(i) = bdry_func(node(i, 1), node(i, 2))
    end do
end subroutine
    
!----------------------------------
! Source function
!----------------------------------
real function source_func(x, y)
    implicit none
    real :: x, y
    source_func = 2.0*x*(1.0-x)*y*(1.0-y)
end function source_func
    
!----------------------------------
! Boundary value function
!----------------------------------
real function bdry_func(x, y)
    implicit none
    real :: x, y
    bdry_func = 0.0
end function bdry_func
    
subroutine set_mesh_on_rectangular_domain(domin, nsubd, grid_spacing, &
                                        nnode, node, ntriangle, triangle)
    implicit none
    real, intent(in) :: domin(2)
    integer, intent(in) :: nsubd(2)
    real, intent(in)    :: grid_spacing(2)
    integer, intent(in) :: nnode, ntriangle
    real, intent(out) :: node(nnode, 2)
    integer, intent(out) :: triangle(ntriangle, 3)
        
    integer :: i, j, nindex, tindex
    real    :: x, y
        
    ! Set nodes
    x = domin(1)
    nindex = 0
    do i = 1, nsubd(1)+1
        y = domin(2) 
        do j = 1, nsubd(2)+1
            nindex = nindex + 1   
                
            node(nindex, 1) = x
            node(nindex, 2) = y
                
            y = y + grid_spacing(2)
        end do
        x = x + grid_spacing(1)
    end do
            
    if (nindex .ne. nnode) then
        print*, "ERROR! Node index computation is incorrect."
    endif
        
    ! Set triangles
    tindex = 0
    do j = 1, nsubd(2)
        do i = 1, nsubd(1)
            tindex = tindex + 1
            triangle(tindex, 1) =  (j-1)*(nsubd(1)+1) + i   ! (i,j)
            triangle(tindex, 2) =  (j-1)*(nsubd(1)+1) + i+1 ! (i+1,j)
            triangle(tindex, 3) =      j*(nsubd(1)+1) + i+1 ! (i+1,j+1)
                
            tindex = tindex + 1
            triangle(tindex, 1) =  (j-1)*(nsubd(1)+1) + i   ! (i,j)
            triangle(tindex, 2) =      j*(nsubd(1)+1) + i+1 ! (i+1,j+1)
            triangle(tindex, 3) =      j*(nsubd(1)+1) + i   ! (i,j+1)
        end do
    end do
            
end subroutine


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
    
    integer, parameter      :: GRIDNX = 2, GRIDNY = 2
    integer, parameter      :: NNODE = (GRIDNX + 1) * (GRIDNY + 1)
    integer, parameter      :: NTRIANGLE = GRIDNX * GRIDNY * 2
    integer, parameter      :: NMAXNONZEROS = 10
    
    real        :: node(NNODE, 2)
    integer     :: triangle(NTRIANGLE, 3), node_neighbor(NNODE, NMAXNONZEROS)
    logical     :: isbdry(NNODE)
    
    real        :: f(NNODE), g(NNODE) ! -div(grad u) = f, u = g on boundary
    
    integer     :: K_sparse_index(NNODE, NMAXNONZEROS)
    real        :: K_sparse(NNODE, NMAXNONZEROS), b(NNODE)
    
    integer :: i
    real        :: tempValue
    
    ! Set domain
    domin = [0.0, 0.0]
    domax = [1.0, 1.0]
    nsubd = [GRIDNX, GRIDNY]
    
    grid_spacing = (domax-domin)/nsubd
    grid_size = nsubd + 1
    
    ! Set Mesh
    call set_mesh_on_rectangular_domain(domin, nsubd, grid_spacing, &
                nnode, node, ntriangle, triangle, &
                NMAXNONZEROS, node_neighbor, isbdry)
    
    call set_source_bdry_func(domin, nsubd, grid_spacing, nnode, node, f, g)
    
    ! Compute K-matrix and b
    call compute_fem_matrix_equation(domin, nsubd, grid_spacing, &
                nnode, node, ntriangle, triangle, node_neighbor, isbdry, f, g, &
                NMAXNONZEROS, K_sparse_index, K_sparse, b)
    
    ! Body of fem
    !print *, 'Domain : ', grid_spacing, grid_size, nnode, ntriangle
    
    !print *, 'Nodes : '
    !do i = 1, nnode
    !    print*, node(i,:)
    !end do
    
    !print *, 'Triangles : '
    !do i = 1, ntriangle
    !    print *, triangle(i,:)
    !end do
    
    print *, 'Neighbors : '
    do i = 1, nnode
        if (node_neighbor(i, 1) > 0) then
            print*, 'i=', i, " : ", node_neighbor(i, 1:node_neighbor(i,1)+1)
        else
            print*, 'i=', i, " : ", node_neighbor(i, 1)
        end if
    end do
    
    print *, 'f : ', f(1:nnode)
    print *, 'g : ', g(1:nnode)
    print *, 'bdry', isbdry
    
    print *, 'K :'
    do i = 1, nnode
        if (K_sparse_index(i, 1) > 0) then
            print*, 'i=', i, " : ", K_sparse(i, 1:K_sparse_index(i,1)+1)
        else
            print*, 'i=', i, " : "
        end if
    end do
    
    print*, 'b :', b
    
    tempValue = 0.0
    do i = 3,  K_sparse_index(5,1)+1
        tempValue = tempValue + K_sparse(5,i) * b(K_sparse_index(5,i))
    end do
    print*, tempValue
   
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
    source_func = -4.0 !2.0*(x*(1.0-x)+y*(1.0-y))
end function source_func
    
!----------------------------------
! Boundary value function
!----------------------------------
real function bdry_func(x, y)
    implicit none
    real :: x, y
    bdry_func = x*x+y*y
end function bdry_func
    
subroutine set_mesh_on_rectangular_domain(domin, nsubd, grid_spacing, &
                                        nnode, node, ntriangle, triangle, &
                                        nmaxnonzeros, node_neighbor, isbdry)
    implicit none
    real, intent(in) :: domin(2)
    integer, intent(in) :: nsubd(2)
    real, intent(in)    :: grid_spacing(2)
    integer, intent(in) :: nnode, ntriangle, nmaxnonzeros
    real, intent(out) :: node(nnode, 2)
    integer, intent(out) :: triangle(ntriangle, 3)
    integer, intent(out) :: node_neighbor(nnode, nmaxnonzeros)
    logical, intent(out) :: isbdry(nnode)
        
    integer :: i, j, k, nindex, tindex, ci(3)
    real    :: x, y
    logical :: exist = .false.
    integer :: tnnei(nnode, 21) ! nmaxnonzeros = 10, 2*nmaxnonzeros+1
    integer :: temp, count
    integer :: tmp_neighbor(nmaxnonzeros), numit
    integer :: bSort
    
    ! Set nodes
    y = domin(2)
    nindex = 0
    do j = 1, nsubd(2)+1
        x = domin(1) 
        do i = 1, nsubd(1)+1
            nindex = nindex + 1   
                
            node(nindex, 1) = x
            node(nindex, 2) = y
                
            x = x + grid_spacing(1)
        end do
        y = y + grid_spacing(2)
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
    
    ! Find node neighbors
    tnnei(:, 1) = 0
    
    do i = 1, ntriangle
        ci = triangle(i,:)
        
        ! For local node(1)
        temp = tnnei(ci(1), 1)
        tnnei(ci(1), 1) = tnnei(ci(1),1) + 2
        tnnei(ci(1), temp+2) = ci(2)
        tnnei(ci(1), temp+3) = ci(3)
        
        ! For local node(2)
        temp = tnnei(ci(2), 1)
        tnnei(ci(2), 1) = tnnei(ci(2),1) + 2
        tnnei(ci(2), temp+2) = ci(1)
        tnnei(ci(2), temp+3) = ci(3)
        
        ! For local node(3)
        temp = tnnei(ci(3), 1)
        tnnei(ci(3), 1) = tnnei(ci(3),1) + 2
        tnnei(ci(3), temp+2) = ci(1)
        tnnei(ci(3), temp+3) = ci(2)
    end do            
    
    !do i = 1, nnode
    !    if (tnnei(i, 1) > 0) then
    !        print*, 'i=', i, " : ", tnnei(i, 1:tnnei(i,1)+1)
    !    else
    !        print*, 'i=', i, " : ", tnnei(i, 1)
    !    end if
    !end do
    
    do i = 1, nnode
        count = 0
        if (tnnei(i,1) > 0) then
            do j = 2, tnnei(i,1) + 1
                temp = tnnei(i, j)
                exist = .false.
                do k = 2, count + 1
                    if (temp == node_neighbor(i, k)) then
                        exist = .true.
                        exit
                    end if
                end do
                
                if (exist == .false.) then
                    count = count + 1
                    node_neighbor(i,count+1) = temp
                end if
            end do
            node_neighbor(i,1) = count
            
            tmp_neighbor(1:count) = node_neighbor(i, 2:count+1)
            numit = bSort(tmp_neighbor, count)
            node_neighbor(i, 2:count+1) = tmp_neighbor(1:count)
        end if
    end do
    
    ! Set boundary nodes
    isbdry = .false.
    do j = 1, nsubd(2)+1
        do i = 1, nsubd(1)+1
            if ( i == 1 .or. i==nsubd(1)+1 .or. j==1 .or. j==nsubd(2)+1) then
                isbdry((j-1)*(nsubd(1)+1)+i) = .true.
            end if
        end do
    end do
    
end subroutine
    
INTEGER FUNCTION bSort(data, size)
IMPLICIT NONE

INTEGER, INTENT(IN) :: size
INTEGER, INTENT(OUT) :: data(size)

        INTEGER :: temp
        INTEGER :: i
        INTEGER :: j
        INTEGER :: num_it = 0
        INTEGER :: ns = 0

        DO i = SIZE, 1, -1
                ns = 0
                DO j = 1, i-1, 1
                        num_it = num_it + 1
                        IF (data(j) > data(j+1)) THEN
                        ns = 1
                        temp = data(j)
                        data(j) = data(j+1)
                        data(j+1) = temp
                        END IF
                END DO
                IF (ns == 0) EXIT
        END DO

bSort = num_it

END FUNCTION bSort

    


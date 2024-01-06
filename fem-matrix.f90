!---------------------------------------------------
! Compute K and b for Kx=b of Laplace equation
! Assume :
!   . Linear basis function : {1-(x+y), x, y}
!   . Standard triangle : {(0,0), (1,0), (0,1)}
!---------------------------------------------------
subroutine compute_fem_matrix_equation(domin, nsubd, grid_spacing, &
                nnode, node, ntriangle, triangle, node_neighbor, isbdry, f, bg, &
                nmaxnonzeros, K_sparse_index, K_sparse, b)    
    implicit none
    
    real, intent(in)    :: domin(2)
    integer, intent(in) :: nsubd(2)
    real, intent(in)    :: grid_spacing(2)
    integer, intent(in) :: nnode, ntriangle
    real, intent(in)    :: node(nnode, 2)
    integer, intent(in) :: triangle(ntriangle, 3)
    logical, intent(in) :: isbdry(nnode)
    real, intent(in)    :: f(nnode), bg(nnode)
    integer, intent(in) :: nmaxnonzeros
    integer, intent(in) :: node_neighbor(nnode, nmaxnonzeros)
    
    integer, intent(out) :: K_sparse_index(nnode, nmaxnonzeros)
    real, intent(out) :: K_sparse(nnode, nmaxnonzeros)
    real, intent(out) :: b(nnode)
    
    integer :: i, j, k, l
    integer :: nidx(3), tidx
    real    :: v(3,2), tv(3, 2)
    real    :: G(3,3), dphi(3,2)
    real    :: T(2,2), invT(2,2), W(2,2), det
    real    :: dphiw(3,2)
    real    :: ta(3,3), tG
    
    ! Initialize K and b
    K_sparse_index(:, 1) = node_neighbor(:,1) + 1 ! number of neighbor nodes + 1 (itself)
    do i = 1, nnode
        K_sparse_index(i, 2) = i
        K_sparse_index(i, 3:node_neighbor(i,1)+2) = &
            node_neighbor(i,2:node_neighbor(i,1)+1)
    end do
    
    K_sparse = 0.0
    K_sparse(:,1) = K_sparse_index(:,1)
    ! print*, 'K(:,1)=', K_sparse(:,1)
    
    b = 0.0
    
    ! Iterate all triangles to compute the matrix K
    
    ! Set integral of (phi_i * phi_j) over standard triangle
    G(1,1) = 1.0/12.0
    G(1,2) = 1.0/24.0
    G(1,3) = G(1,2)
    G(2,2) = G(1,1)
    G(2,3) = G(1,2)
    G(3,3) = G(1,1)
    G(2,1) = G(1,2)
    G(3,1) = G(1,3)
    G(3,2) = G(2,3)
    
    ! Set gradients of standard basis functions defined above
    dphi(1,:) = [-1.0, -1.0]
    dphi(2,:) = [ 1.0,  0.0]
    dphi(3,:) = [ 0.0,  1.0]
    
    ! Compute the element matrix and add it to K
    do i = 1, ntriangle
        nidx = triangle(i,:)
        
        v(1, :) = node(nidx(1), :)
        v(2, :) = node(nidx(2), :)
        v(3, :) = node(nidx(3), :)
        
        tv = v
        tv(1,:) = tv(1,:) - v(1,:)
        tv(2,:) = tv(2,:) - v(1,:)
        tv(3,:) = tv(3,:) - v(1,:)
        
        !print *, tv(1,:)
        !print *, tv(2,:)
        !print *, tv(3,:)
        
        T(:,1) = tv(2,:)
        T(:,2) = tv(3,:)
        
        invT(1,1) = T(2,2)
        invT(1,2) = -T(1,2)
        invT(2,1) = -T(2,1)
        invT(2,2) = T(1,1)
        
        det = abs(T(1,1)*T(2,2)-T(1,2)*T(2,1))
        invT = invT / det
        
        !if (i==1) then
            !print *, 'T='
            !print *, T(1, :)
            !print *, T(2, :)
            !print *, 'invT='
            !print *, invT(1, :)
            !print *, invT(2, :)
        !end if
                
        det = abs(T(1,1)*T(2,2)-T(1,2)*T(2,1))
        !print*, 'det=', det
        
        W = matmul(invT, transpose(invT))
        ! if (i==1) then
        !     print*, W(1, :)
        !     print*, W(2, :)
        ! end if
                
        !print*, T
        !print*, W
        
        !print*, matmul(dphi(1,:),W)
        !print*, dot_product(matmul(dphi(1,:),W),dphi(1,:))
        do j = 1, 3
            dphiw(j,:) = matmul(dphi(j,:),W)
            do k = j, 3
                ta(j,k) = dot_product(dphiw(j,:), dphi(k,:))*0.5*det
                ! if (i==1) then
                !     print*, j, k, ta(j,k)
                ! end if
                
                tG = G(j,k)*det
                ! if (i==1) then
                !     print*, j, k, tG, f(nidx(k))
                ! end if
                
                ! nidx(j)
                do l = 2, K_sparse_index(nidx(j),1) + 1
                    if (K_sparse_index(nidx(j),l) == nidx(k)) then
                        K_sparse(nidx(j),l) = K_sparse(nidx(j),l) + ta(j,k)
                        exit
                    end if    
                end do
                
                b(nidx(j)) = b(nidx(j)) + f(nidx(k)) * tG
                
                ! nidx(k)
                if ( j .ne. k ) then
                    do l = 2, K_sparse_index(nidx(k),1) + 1
                        if (K_sparse_index(nidx(k), l) == nidx(j)) then
                            K_sparse(nidx(k), l) = K_sparse(nidx(k),l) +ta(j,k)
                            exit
                        end if    
                    end do
                    
                    b(nidx(k)) = b(nidx(k)) + f(nidx(j)) * tG
                    
                end if
            end do
        end do        
    end do    
    
    ! print*, "b=", b
    ! print*, "K="
    ! do i = 1, nnode
    !     print*, "i=", K_sparse(i,:)
    ! end do

    ! For boundary condition
    do i = 1, nnode
        if (isbdry(i)) then
            K_sparse(i,2) = 1.0
            K_sparse(i, 3:K_sparse_index(i,1)+1) = 0.0
            b(i) = bg(i)
        else
            do j = 2, K_sparse_index(i,1)+1
                if (isBdry(K_sparse_index(i,j))) then
                    b(i) = b(i) - K_sparse(i,j) * bg(K_sparse_index(i,j))
                    K_sparse(i,j) = 0.0
                end if
            end do
        end if
    end do  
    
end subroutine 
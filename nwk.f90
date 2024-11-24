module nwk
    implicit none

    integer, parameter :: node_num = 1000
    real(8) :: avg_deg
    integer :: edge_num
    integer, dimension (:), allocatable :: n_deg
    integer, dimension (:,:), allocatable :: adj_mat, node_id

    real(8) :: rd1, rd2

    contains

    ! ---------------------------------------------------------------------------
    ! constructing the network

        subroutine const_nwk ()
            implicit none
            integer :: add_edge
            integer :: i, j

            edge_num = node_num * avg_deg / 2.0
            add_edge = 0

            allocate(n_deg(node_num))
            allocate(adj_mat(node_num, node_num))
            allocate(node_id(node_num, node_num))

            adj_mat = 0
            node_id = 0
            n_deg = 0

            do while (add_edge < edge_num)
                call random_number(rd1)
                i = ceiling(rd1 * node_num) + 1
                if (i > node_num) then
                    i = node_num
                end if

                call random_number(rd2)
                j = ceiling(rd2 * node_num) + 1
                if (j > node_num) then
                    j = node_num
                end if

                if (i /= j .and. adj_mat(i, j) == 0) then
                    adj_mat(i, j) = 1
                    adj_mat(j, i) = 1
                    add_edge = add_edge + 1
                end if
            end do

            do i = 1, node_num
                do j = 1, node_num
                    if (adj_mat(i,j) == 1) then
                        n_deg(i) = n_deg(i) + 1
                        node_id(i, n_deg(i)) = j
                    end if
                end do
            end do


        end subroutine

    ! ---------------------------------------------------------------------------
        subroutine clean_nwk () ! just deallocating
            deallocate(n_deg)
            deallocate(adj_mat)
            deallocate(node_id)
        end subroutine

end module nwk
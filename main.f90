program opinion_game_model

    use nwk
    implicit none

    ! 参数
    integer, parameter :: initials = 1000, net_init = 1, sim_time = 100000
    real(8), parameter :: coop_init = 0.5, sigma = 2.0

    ! 变量
    real(8) :: K, gamma, alpha, beta, epsilon, eta

    ! 中间量数组
    real(8), dimension(node_num): node_x_t_curr, node_x_t_last, node_x_tmp
    real(8), dimension(node_num, node_num): weight_mat
    real(8), dimension(node_num) :: mag

    ! 记录数组
    real(8), dimension(sim_time) :: coop_freq, benefit_coop, benefit_def
    real(8), dimension(sim_time) :: ord_para_global
    real(8), dimension(node_num) :: ord_para_local
    real(8), dimension(node_num) :: cost_local
    real(8), dimension(node_num) :: payoff
    real(8) :: ord_para_neighbor, gamma
    real(8) :: ord_para_0, x_mean, x_var

    integer, dimension(node_num) :: strategy

    character(len=128) :: arg
    character(200) :: filename1

    ! 临时变量
    real(8) :: rd, i, j, delta
    integer :: kk, k_start, k_end, net, kh
    integer :: coop_node, add_coop_node

    if (command_argument_count().ne.5) then
        write (*,*) 'usage: ./main [avg_deg] [K] [gamma] [alpha] [beta]'
        call exit(1)
    end if

    ! 参数
    call get_command_argument(1, arg)
    read (arg, *) avg_deg
    call get_command_argument(2, arg)
    read (arg, *) K
    call get_command_argument(3, arg)
    read (arg, *) gamma
    call get_command_argument(4, arg)
    read (arg, *) alpha
    call get_command_argument(5, arg)
    read (arg, *) beta


    ! 生成文件名
    write(filename1, "(A, F0.2, A, F0.2, A, F0.2, A, F0.2, A, F0.2, A, I0, A, I0, A)") "final_avg=", avg_deg, "_phi=", &
    phi, "_beta=", beta, "_b=", boy_init, "_c=", cat_init, "_as=", alpha_start, "_ae=", alpha_end, ".txt"
    open(unit=13, file = filename1)


    do kk = k_start, k_end

        print *, kk

        K = kk * 1.0 / 100.0

        ! initials
        node_x_tmp = 0.0
        node_x_t_curr = 0.0
        node_x_t_last = 0.0
        weight_mat = 0.0
        strategy = 0

        cc = 0
        do net = 1, net_init

            call const_nwk ()        
            
            ! 初始化net_init次节点状态
            do kh = 1, initials

                ! 初始化coop节点数量, 初始化50%
                coop_node = floor(coop_init * node_num)
                add_coop_node = 0
                do while (add_coop_node < coop_node)
                    call random_number(rd)
                    copy_node_id = floor(rd * node_num) + 1
                    if (strategy(copy_node_id) == 0) then
                        strategy(copy_node_id) = 1
                        add_coop_node = add_coop_node + 1
                    end if
                end do 

                ! 初始化节点意见值
                do i = 1, node_num
                    call random_number(rd)
                    node_x_t_curr(i) = rd * 2.0 - 1.0
                end do

                ! 初始的序参量
                ord_para_0 = 0.0
                x_mean = 0.0
                x_var = 0.0
                do i = 1, node_num
                    x_mean = x_mean + node_x_t_curr(i)
                end do
                x_mean = x_mean * 1.0 / node_num
                do i = 1, node_num
                    x_var = x_var + (node_x_t_curr(i) - x_mean) ** 2
                end do
                x_var = (x_var * 1.0 / node_num) ** (0.5)
                ord_para_0 = 1 - x_var

                ! 初始化权重
                mag = 0.0
                do i = 1, node_num
                    mag = abs(node_x_t_curr(i) - node_x_t_curr + epsilon)
                    weight_mat(i, :) = mag(i) ** (-1.0 * beta)
                    weight_mat(i, :) = weight_mat(i, :) / sum(weight_mat(i, :))
                end do
                !-----------------------------------------------------

                do tt = 1, sim_time
                    ! 假设不记录opinion
                    node_x_t_last = node_x_t_curr
                    node_x_tmp = 0.0

                    ! 意见更新
                    do i = 1, node_num
                        ! 策略0
                        if (strategy(i) == 0) then
                            node_x_tmp(i) = gamma * node_x_t_curr(i)
                        ! 策略1
                        elseif (strategy(i) == 1) then
                            node_x_tmp(i) = gamma * node_x_t_curr(i)
                            do j = 1, node_num
                                node_x_tmp(i) = node_x_tmp(i) + K * (adj_mat(i, j) * weight_mat(i, j) * tanh(alpha * node_x_t_curr(j)))
                            end do
                        endif
                    end do

                    node_x_t_curr = node_x_tmp

                    ! 权重更新
                    mag = 0.0
                    weight_mat = 0.0
                    do i = 1, node_num
                        mag = abs(node_x_t_curr(i) - node_x_t_curr + epsilon)
                        weight_mat(i, :) = mag(i) ** (-1.0 * beta)
                        weight_mat(i, :) = weight_mat(i, :) / sum(weight_mat(i, :))
                    end do

                    ! 计算每个节点的payoff
                    ord_para_local = 0.0
                    cost_local = 0.0
                    do i = 1, node_num
                        if (n_deg(i) > 0) then
                            do j = 1, n_deg(i)
                                delta = node_x_t_curr(i) - node_x_t_curr(node_id(i, j))
                                call heaviside(delta, sigma, gamma)
                                ord_para_neighbor = (1 - 0.5 * abs(delta)) * gamma
                                ord_para_local(i) = ord_para_local(i) + ord_para_neighbor
                            end do
                            ord_para_local(i) = ord_para_local(i) * 1.0 / n_deg(i)
                        end if
                        cost_local(i) = abs(node_x_t_curr(i) - node_x_t_last(i))
                        payoff(i) = ord_para_local(i) - eta * cost_local(i)
                    end do

                    ! 策略更新
                    do i = 1, node_num
                        
                        call random_number(rd)
                    end do
                    

                    
                end do

               

        ! 求均值
        do m = 1, net_init * initials
            sus_avg = sus_avg + infect_final(m)
        end do

        sus_avg = sus_avg / real(net_init * initials) 

        ! 求方差
        do m = 1, net_init * initials
            sus_var = sus_var + (infect_final(m)-sus_avg)**2
        end do
        sus_var = sus_var / real(net_init * initials)

        ! 计算磁化率
        susceptibility = real(node_num) * sus_var/sus_avg

        print *, rr, sus_avg, sus_var, (sus_var**2)/(sus_avg**2), susceptibility

        do m = 1, sim_time
            write(13, *) m, avg_deg, alpha, beta, boy_init, avg_infect_f(m), avg_infect_f_normal(m), avg_infect_f_im(m)
        end do

        write(14, *) avg_deg, alpha, beta, boy_init, susceptibility, (sus_var**2)/(sus_avg**2)

    end do

    close(unit=13)
    close(unit=14)


contains

    ! 阶跃函数定义
    subroutine heaviside(x, y, result)
        real(8), intent(in) :: x, y
        real(8), intent(out) :: result
        real(8) :: h1, h2

        if (x + y >= 0) then
            h1 = 1.0
        elseif
            h1 = 0.0
        end if

        if (x - y >= 0) then
            h2 = 1
        elseif
            h2 = 0
        end if

        result = h1 - h2

    end subroutine heaviside

end program opinion_game_model
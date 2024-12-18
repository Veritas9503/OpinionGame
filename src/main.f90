program opinion_game_model

    use nwk
    implicit none

    ! 参数
    integer, parameter :: initials = 1000, net_init = 1, sim_time = 100000
    real(8), parameter :: coop_init = 0.5, sigma = 2.0

    ! 变量
    real(8) :: K, gamma, alpha, beta, epsilon, eta, kappa

    ! 中间量数组
    real(8), dimension(node_num): node_x_t_curr, node_x_t_last, node_x_tmp
    real(8), dimension(node_num, node_num): weight_mat
    real(8), dimension(node_num) :: mag

    ! 临时数组
    integer, dimension(sim_time) :: transition_frac_dd, transition_frac_cc, transition_frac_cd, transition_frac_dc
    real(8), dimension(sim_time) :: coop_freq, benefit_avg
    real(8), dimension(sim_time) :: ord_para_global
    real(8), dimension(node_num) :: ord_para_local
    real(8), dimension(node_num) :: cost_local
    real(8), dimension(node_num) :: payoff
    real(8) :: ord_para_neighbor, gamma, imitate_prob
    real(8) :: ord_para_0, x_mean, x_var
    integer, dimension(node_num) :: strategy

    ! 记录稳态值和动态值
    integer, parameter :: record_time = 1000
    real(8) :: coop_freq_record, benefit_avg_record, ord_para_global_record
    real(8), dimension(sim_time) :: coop_freq_dyn, benefit_avg_dyn, ord_para_global_dyn
    real(8), dimension(sim_time) :: transition_frac_dd_dyn, transition_frac_cc_dyn, transition_frac_cd_dyn, transition_frac_dc_dyn

    character(len=128) :: arg
    character(200) :: filename1

    ! 临时变量
    real(8) :: rd, i, j, delta, p
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

        coop_freq_record = 0.0
        benefit_avg_record = 0.0
        ord_para_global_record = 0.0

        coop_freq_dyn = 0.0
        benefit_avg_dyn = 0.0
        ord_para_global_dyn = 0.0

        transition_frac_dd_dyn = 0.0
        transition_frac_cc_dyn = 0.0
        transition_frac_cd_dyn = 0.0
        transition_frac_dc_dyn = 0.0

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
                x_mean = sum(node_x_t_curr) * 1.0 / node_num
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

                ! coop_freq, benefit_avg
                coop_freq = 0.0
                benefit_avg = 0.0

                transition_frac_dd = 0.0
                transition_frac_cc = 0.0
                transition_frac_cd = 0.0
                transition_frac_dc = 0.0

                do tt = 1, sim_time
                    ! 假设不记录opinion
                    node_x_t_last = node_x_t_curr
                    node_x_tmp = 0.0
                    transition_frac = 0

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
                    ! DD : 0, CC : 1, CD : 2, DC : 3
                    do i = 1, node_num
                        call random_number(rd)
                        j = floor(rd * n_deg(i)) + 1
                        imitate_prob = 1.0 / (1 + exp(-kappa * (payoff(j) - payoff(i))))

                        call random_number(p)
                        ! 如果模仿
                        if (p <= imitate_prob) then
                            if (strategy(i) == 0 .and. strategy(j) == 0) then
                                transition_frac_dd(tt) = transition_frac_dd(tt) + 1.0 / node_num
                            elseif (strategy(i) == 1 .and. strategy(j) == 1) then
                                transition_frac_cc(tt) = transition_frac_cc(tt) + 1.0 / node_num
                            elseif (strategy(i) == 1 .and. strategy(j) == 0) then
                                transition_frac_cd(tt) = transition_frac_cd(tt) + 1.0 / node_num
                            elseif (strategy(i) == 0 .and. strategy(j) == 1) then
                                transition_frac_dc(tt) = transition_frac_dc(tt) + 1.0 / node_num
                            end if
                            strategy(i) = strategy(j)
                        ! 如果不模仿
                        else
                            if (strategy(i) == 0) then 
                                transition_frac_dd(tt) = transition_frac_dd(tt) + 1.0 / node_num
                            elseif (strategy(i) == 1) then
                                transition_frac_cc(tt) = transition_frac_cc(tt) + 1.0 / node_num
                            end if
                        end if
                    end do
                    
                    ! 统计合作比例，平均社会效益
                    coop_freq(tt) = sum(strategy) * 1.0 / node_num
                    benefit_avg(tt) = sum(payoff) * 1.0 / node_num

                    ! 计算全局序参量
                    do i = 1, node_num
                        ord_para_global(tt) = ord_para_global(tt) + &
                        (node_x_t_curr(i) - sum(node_x_t_curr) * 1.0 / node_num) ** 2
                    end do
                    ord_para_global(tt) = 1 - (ord_para_global(tt) * 1.0 / node_num) ** 0.5

                    ! 记录动态数据
                    coop_freq_dyn(tt) = coop_freq_dyn(tt) + coop_freq(tt) / (net_init * initials)
                    benefit_avg_dyn(tt) = benefit_avg_dyn(tt) + benefit_avg(tt) / (net_init * initials)
                    ord_para_global_dyn(tt) = ord_para_global_dyn(tt) + ord_para_global(tt) / (net_init * initials)

                    transition_frac_dd_dyn(tt) = transition_frac_dd_dyn(tt) + transition_frac_dd(tt) / (net_init * initials)
                    transition_frac_cc_dyn(tt) = transition_frac_cc_dyn(tt) + transition_frac_cc(tt) / (net_init * initials)
                    transition_frac_cd_dyn(tt) = transition_frac_cd_dyn(tt) + transition_frac_cd(tt) / (net_init * initials)
                    transition_frac_dc_dyn(tt) = transition_frac_dc_dyn(tt) + transition_frac_dc(tt) / (net_init * initials)
                end do

                ! 记录稳态数据
                coop_freq_record = coop_freq_record + sum(coop_freq(sim_time-record_time:sim_time)) & 
                * 1.0 / record_time
                benefit_avg_record = benefit_avg_record + sum(benefit_avg(sim_time-record_time:sim_time)) &
                * 1.0 / record_time
                ord_para_global_record = ord_para_global_record + &
                sum(ord_para_global(sim_time-record_time:sim_time)) * 1.0 / record_time
            end do

        end do

        ! 稳态数据
        coop_freq_record = coop_freq_record * 1.0 / (net_init * initials)    
        benefit_avg_record = benefit_avg_record * 1.0 / (net_init * initials)
        ord_para_global_record  = ord_para_global_record * 1.0 / (net_init * initials)
        
        do m = 1, sim_time
            write(13, *) m, coop_freq_dyn(m), benefit_avg_dyn(m), ord_para_global_dyn(m), &
            transition_frac_dd_dyn(m), transition_frac_cc_dyn(m), transition_frac_cd_dyn(m), transition_frac_dc_dyn(m)
        end do

        write(14, *) coop_freq_record, benefit_avg_record, ord_para_global_record

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
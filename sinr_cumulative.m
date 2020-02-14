%This file is for DOA estimation in LS-MIMO angular spread scenarios

clear;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nt = 500;%number of snapshots
Ite_num = 1e3;%number of iterations
u = 2 * pi * 0.5;%d/lambda = 0.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K = 10;%user number
Nk = 50;%number of paths
snr = 10;%dB, received SNR
M = 100;
Iteg_step = 0.1;
Mx = sqrt(M);
My = Mx;
theta_d = 1;%angular deviation
phi_d = 0.5;%angular deviation
a = zeros(M, 1);%steering vector
N = zeros(M, Nt);%received noise at the BS
L = 3;%cell number
N_total = Ite_num * Nt * L;
SINR_store = zeros(N_total, 7);
r = 300;%center to edge distance(m)
rc = r * 0.8;
rh = 30;%minimum terminal radius of the cell(m)
ra = rc / rh - 1;
BS_height = 30;%m
gamma = 3.8;%decay exponent
sigma = 10^0.8;
pos = zeros(K, L);
base(1:L,1) = [0;-1i * 2 * r;(sqrt(3) * r - 1i * r);];
r_inside = zeros(K, L);
dis_inside = zeros(K, L);
for ii = 1 : Ite_num
    dis(1:K,1:L) = (rem(rand(K,L) * ra, ra) + 1) * rh;
    ang(1:K,1) = rand(K, 1) * 2 * pi / 3 + 4 * pi / 3;
    ang(1:K,2) = rand(K, 1) * 2 * pi / 3;
    ang(1:K,3) = rand(K, 1) * 2 * pi / 3 + 2 * pi / 3;
    pos(1:K,1:L) = dis .* (exp(1i * ang));
    for ll = 1 : L
        pos(:,ll) = pos(:,ll) + base(ll,1);
    end
    amp_k = sqrt(10 ^ (snr * 0.1) / Nk);
    Sk = (sign(randn(L, K*Nt))+1i*sign(randn(L, K*Nt))) / sqrt(2) * amp_k;%QPSK modulation
    for j = 1 : L
        r_inside(:, j) = abs(pos(:, j) - base(j, 1));
        dis_inside(:, j) = sqrt(r_inside(:, j).^2 + BS_height^2);%the distance inside the cell
    end
    shadow_amp = 10.^(randn(K*L, L) * sigma * 0.1);%the shadow fading from all the KL users to all the L cells
    for j = 1 : L
        Y1 = zeros(M, Nt);%received signal
        YS1 = zeros(M, K*Nt);%received useful signal
        N = (randn(M, Nt) + 1i * randn(M, Nt)) / sqrt(2);
        Y1 = Y1 + N;
        theta_store = zeros(K, L);
        phi_store= zeros(K, L);
        store_count = 1;
        for l = 1 : L
            %generate the received signal
            theta = zeros(K, 1);
            theta(:, 1) = atan2(imag(pos(:, l) - base(j, 1)), real(pos(:, l) - base(j, 1)));%-pi~pi
            switch j
                case 1
                    theta(:, 1) = theta(:, 1) + 2*pi/3;
                case 3
                    for k = 1 : K
                        if theta(k, 1)>0
                            theta(k, 1) = theta(k, 1) - 2*pi/3;
                        else
                            theta(k, 1) = theta(k, 1) + 2*pi - 2*pi/3;%0~2*pi/3
                        end
                    end
                otherwise
            end
            theta = theta + pi / 6;%pi/6~5*pi/6
            r_jl = abs(pos(:, l) - base(j, 1));
            dis_jl = sqrt(r_jl.^2 + BS_height^2);
            phi = atan(r_jl / BS_height);%the angle between the propogation path and the vertical line
            theta_store(:, l) = theta;%the angles from terminals in the lth cell to the BS in the jth cell
            phi_store(:, l) = phi;
            if l == j
                for k = 1 : K
                    for jj = 1 : Nt
                        alpha_k = (randn(Nk, 1) + 1i * randn(Nk, 1)) / sqrt(2);%small scale fading
                        theta_k = theta(k, 1) + randn(Nk, 1) * theta_d / 180 * pi;%Gaussian distribution
                        phi_k = phi(k, 1) + randn(Nk, 1) * phi_d / 180 * pi;%Gaussian distribution
                        for ll = 1 : Nk
                            for mx = 1 : Mx
                                for my = 1 : My
                                    m = mx + (my - 1) * Mx;
                                    a(m, 1) = exp(1i * u * sin(phi_k(ll,1)) * ((mx - 1) * cos(theta_k(ll,1)) + (my - 1) * sin(theta_k(ll,1))));
                                end
                            end
                            YS1(:, (k-1)*Nt+jj) = YS1(:, (k-1)*Nt+jj) + alpha_k(ll, 1) * a * Sk(l, (k-1)*Nt+jj);
                        end
                    end
                    YS1(:, (k-1)*Nt+1 : k*Nt) = YS1(:, (k-1)*Nt+1 : k*Nt);%power control
                    Y1 = Y1 + YS1(:, (k-1)*Nt+1 : k*Nt);
                end
            else
                for k = 1 : K
                    temp = zeros(M, Nt);
                    for mx = 1 : Mx
                        for my = 1 : My
                            m = mx + (my - 1) * Mx;
                            a(m, 1) = exp(1i * u * sin(phi(k, 1)) * ((mx - 1) * cos(theta(k, 1)) + (my - 1) * sin(theta(k, 1))));
                        end
                    end
                    for jj = 1 : Nt
                        alpha_k = (randn(1, 1) + 1i * randn(1, 1)) / sqrt(2);%small scale fading
                        temp(:, jj) =  alpha_k * a * sqrt(Nk) * Sk(l, (k-1)*Nt+jj);
                    end
                    Y1 = Y1 + temp * sqrt(dis_jl(k, 1)^(-gamma)) * shadow_amp(K*(l-1)+k, j) / sqrt(dis_inside(k, l)^(-gamma)) / shadow_amp(K*(l-1)+k, l);%power control
                end
            end
        end
        %%%%%%%%%%
        %%beamforming
        Q_I = zeros(M, M);
        a_temp = zeros(M, 1);
        N_I_a = ceil(theta_d / Iteg_step) * 2;
        N_I_e = ceil(phi_d / Iteg_step) * 2;
        for k1 = 2 : K
            for n_a = 1 : N_I_a
                theta_temp = theta_store(2, j) / pi * 180 + (n_a - 1 - N_I_a * 0.5) * Iteg_step;%the second user is the interfering user
                for n_e = 1 : N_I_e
                    phi_temp = phi_store(2, j) / pi * 180 + (n_e - 1 - N_I_e * 0.5) * Iteg_step;
                    for mx = 1 : Mx
                        for my = 1 : My
                            m = mx + (my - 1) * Mx;
                            a_temp(m, 1) = exp(1i * u * sin(phi_temp / 180 * pi) * ((mx - 1) * cos(theta_temp / 180 * pi) + (my - 1) * sin(theta_temp / 180 * pi)));
                        end
                    end
                    Q_I = Q_I + a_temp * a_temp' * Iteg_step^2;
                end
            end
        end
        [Q_E, Q_D] = eig(Q_I);
        Q_D_d = diag(Q_D);
        [Q_D_d_s, Ind] = sort(Q_D_d, 'ascend');
        %choose N_tau according to the eigenvalues
        sum_Q_d = 0;
        for N_tau = 1 : M
            sum_Q_d = sum_Q_d + Q_D_d_s(M-N_tau+1);
            if sum_Q_d > 0.9 * sum(Q_D_d_s)
                break;
            else
            end
        end
        Q_s1 = Q_E(:, Ind(M-N_tau+1:M));
        Q_s2 = zeros(M, K*(L-1));
        q_count = 0;
        for l = 1 : L
            if j == l
                continue;
            else
                for k = 1 : K
                    q_count = q_count + 1;
                    for mx = 1 : Mx
                        for my = 1 : My
                            m = mx + (my - 1) * Mx;
                            Q_s2(m, q_count) = exp(1i * u * sin(phi_store(k, l))  ...
                                * ((mx - 1) * cos(theta_store(k, l)) + (my - 1) * sin(theta_store(k, l))));
                        end
                    end
                end
            end
        end
        Q_s = [Q_s1, Q_s2];
        [Q_se, Q_sd] = eig(Q_s * Q_s');
        Q_sd_d = diag(Q_sd);
        [Q_sd_d_s, Ind] = sort(Q_sd_d, 'ascend');
        for Ind_count = 1 : M
            if Q_sd_d_s(Ind_count, 1) > 1e-5
                break;
            else
            end
        end
        Q_Null = Q_se(:, Ind(1:Ind_count-1));
        S_d = zeros(M, 1);
        for mx = 1 : Mx
            for my = 1 : My
                m = mx + (my - 1) * Mx;
                S_d(m, 1) = exp(1i * u * sin(phi_store(1, j)) * ((mx - 1) * cos(theta_store(1, j)) + (my - 1) * sin(theta_store(1, j))));
            end
        end
        w = Q_Null * Q_Null' * S_d;
        w = w / norm(w);
        Y1_b = w' * Y1;
        S1_b = w' * YS1(:, 1:Nt);
        SINR_store((ii-1)*Nt*L + (j-1)*Nt + 1: (ii-1)*Nt*L + j*Nt, 1) = (abs(S1_b).^2 ./ (abs(Y1_b - S1_b).^2))';
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Q_s = zeros(M, M);
        a_temp = zeros(M, 1);
        N_s_a = ceil(theta_d / Iteg_step) * 2;
        N_s_e = ceil(phi_d / Iteg_step) * 2;
        for n_a = 1 : N_s_a
            theta_temp = theta_store(1, j) / pi * 180 + (n_a - 1 - N_s_a * 0.5) * Iteg_step;%the second user is the interfering user
            for n_e = 1 : N_s_e
                phi_temp = phi_store(1, j) / pi * 180 + (n_e - 1 - N_s_e * 0.5) * Iteg_step;
                for mx = 1 : Mx
                    for my = 1 : My
                        m = mx + (my - 1) * Mx;
                        a_temp(m, 1) = exp(1i * u * sin(phi_temp / 180 * pi) * ((mx - 1) * cos(theta_temp / 180 * pi) + (my - 1) * sin(theta_temp / 180 * pi)));
                    end
                end
                Q_s = Q_s + a_temp * a_temp' * Iteg_step^2;
            end
        end
        [Q_E, Q_D] = eig(Q_s);
        Q_D_d = diag(Q_D);
        [Q_s_d, Ind] = sort(Q_D_d, 'descend');
        Q_sa = Q_E(:, Ind(1));
        Q_s3 = zeros(M, 1);
        for l = 1 : L
            if j == l
                continue;
            else
                for k = 1 : K
                    for mx = 1 : Mx
                        for my = 1 : My
                            m = mx + (my - 1) * Mx;
                            Q_s3(m, 1) = exp(1i * u * sin(phi_store(k, l))  ...
                                * ((mx - 1) * cos(theta_store(k, l)) + (my - 1) * sin(theta_store(k, l))));
                        end
                    end
                    pcjl = sqrt(dis_jl(k, 1)^(-gamma)) * shadow_amp(K*(l-1)+k, j) / sqrt(dis_inside(k, l)^(-gamma)) / shadow_amp(K*(l-1)+k, l);
                    Q_I = Q_I + Q_s3 * Q_s3' * pcjl^2;
                end
            end
        end
        [Q_E, Q_D] = eig(Q_I);
        Q_D_d = diag(Q_D);
        [Q_D_d_s, Ind] = sort(Q_D_d, 'ascend');
        %choose N_tau according to the eigenvalues
        sum_Q_d = 0;
        for N_tau = 1 : M
            if Q_D_d_s(M-N_tau+1) < 1e-2
                break;
            else
            end
        end
        Q_Null = Q_E(:, Ind(1:M-N_tau+1));
        w = Q_Null * Q_Null' * Q_sa;
        w = w / norm(w);
        Y1_b = w' * Y1;
        S1_b = w' * YS1(:, 1:Nt);
        SINR_store((ii-1)*Nt*L + (j-1)*Nt + 1: (ii-1)*Nt*L + j*Nt, 2) = (abs(S1_b).^2 ./ (abs(Y1_b - S1_b).^2))';
    end
    display(ii);
end
SINR = mean(SINR_store);
SINR_store = sort(SINR_store);
SINR_store = log10(SINR_store) * 10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = figure;
set(h,'PaperType','A4');
xx = axes('FontSize',16);
plot(SINR_store(:,1),(1:N_total)/N_total,'--','LineWidth',2,'MarkerSize',14)
hold on
plot(SINR_store(:,2),(1:N_total)/N_total, '-k','LineWidth',2,'MarkerSize',14)
plot(SINR_store(0.05*N_total,1),0.05,'o','MarkerSize',10)
plot(SINR_store(0.05*N_total,2),0.05,'ko','MarkerSize',10)
grid on
le = legend('the approach in [15]','the proposed approach', 'Location', 'northwest');
set(le,'Fontsize',16,'Fontname','Times')
xlabel('The SINR with beamforming (dB)','Fontsize',20,'Fontname','Times')
ylabel('The cumulative distribution','Fontsize',20,'Fontname','Times')
%print(h,'-dpdf','sinr_cumulative')
%%%%%%%%%%%%%%%%%%%




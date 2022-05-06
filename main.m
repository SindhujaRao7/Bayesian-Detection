clc;
%As per requirement
sigma2_s = 1;
N_d = 32;
N_c = 8;
K = 50;
P_FA = 0.05;
snr_start = -20;
snr_step = 2;
snr_end = 4;
snr = 1;
N = (K + 1)*(N_c + N_d);
P_H1 = 0.2;
STAT_COUNT = 1000;
snr_list = transpose(snr_start:snr_step:snr_end);
sigma2_w = snr_list(5,1); %for a fixed SNR to generate T(y) for Implementation 1
%required functions
[T_h0, T_h1] = generate_stat(STAT_COUNT,sigma2_s,N,sigma2_w); %generate 1000 test stats for H1 hypothesis Implementation 1.
plot_stat_threshold(STAT_COUNT,sigma2_s,N,P_FA,snr_list);%Plot and Calculate Pfa and PD vs SNR Implementation 1.(a)
plot_stat_threshold_tilda(STAT_COUNT,sigma2_s,N,P_FA,snr_list); %Plot and Calculate Pfa and PD vs SNR with sigma2_w_tilda Implementation 1.(b)
plot_snr_threshold(snr_list,N, sigma2_s,P_H1,P_FA); %Implementation 1.(c)

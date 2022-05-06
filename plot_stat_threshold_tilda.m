function plot_stat_threshold_tilda(STAT_COUNT,sigma2_s,N,P_FA,snr_list)

P_fa_theoritical = zeros(size(snr_list));
P_d_theoritical = zeros(size(snr_list));
P_fa_numerical = zeros(size(snr_list));
P_d_numerical = zeros(size(snr_list));

% Evaluate the performance on H0 and H1

FA_COUNT = 0;
DET_COUNT = 0;
for i = 1: size(snr_list,1)
    %Generate the inaccurate noise variance
    sigma2_w = sigma2_s/(10 ^(snr_list(i,1) / 10));
    sigma_w_tilda = sigma2_w/(10 ^(snr_list(i,1) / 10));
    gamma = sigma_w_tilda * chi2inv(1-P_FA,N);

    [T_h0,T_h1] = generate_stat(STAT_COUNT,sigma2_s,N,sigma_w_tilda);

    for j=1: size(T_h0,1)
        if T_h0(j,1) >= gamma
            FA_COUNT = FA_COUNT + 1;
        end
    end
    for j=1: size(T_h1,1)
        if T_h1(j,1) >= gamma
            DET_COUNT = DET_COUNT + 1;
        end
    end

    %FA_COUNT
    P_fa_numerical(i) = FA_COUNT / STAT_COUNT;
    P_d_numerical(i) = DET_COUNT / STAT_COUNT;

    P_fa_theoritical(i) = P_FA;
    %P_d_theoritical(i) = chi2cdf((gamma/(sigma2_s + sigma_w_tilda)), nu,'upper');
    P_d_theoritical(i) = chi2cdf((gamma/(sigma2_s + sigma2_w)),N, 'upper');
    %     for k = 1: nu/2
    %     P_d_theoritical(i) = (exp(-1/2)*(gamma/(sigma2_s + sigma2_w))) * sum((gamma/(sigma2_s + sigma2_w))/2)^k /factorial(k);
    %     end
end

% display(P_fa_numerical);
% display(P_d_numerical);
% display(P_fa_theoritical);
% display(P_d_theoritical);


%Plot SNR vs P_fa and P_d with inaccurate variance
figure
plot(snr_list, P_fa_numerical);
hold on;
plot(snr_list, P_d_numerical);
hold on;
plot(snr_list, P_fa_theoritical);
hold on;
plot(snr_list, P_d_theoritical);
xlabel('SNR');
xlim([-20 5])
ylabel('Probability');
ylim([0 1])
title('Detector with inaccurate noise variance')
legend('Numerical P_{FA}','Numerical P_{D}','Theoretical P_{FA}','Theoretical P_{D}');


end
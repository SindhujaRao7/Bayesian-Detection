function plot_snr_threshold(snr_list,N,sigma2_s,P_H1,P_FA)

%Initialize thresholds for both detectors
gamma_NP = zeros(size(snr_list));
gamma_BD = zeros(size(snr_list));

for i = 1: size(snr_list)
    sigma2_w = sigma2_s/(10 ^(snr_list(i,1) / 10));
    np_temp = chi2inv(1-P_FA,N) * (sigma2_w);
    bd_temp = 2 * (sigma2_s + sigma2_w) * (sigma2_w /sigma2_s)*(N/2*log(1 + sigma2_s / sigma2_w)+ log(1 - P_H1) - log(P_H1));
    gamma_NP(i) = np_temp;
    gamma_BD(i) = bd_temp;
end

%Compare both detectors
figure
plot(snr_list, gamma_NP,'b');
hold on;
plot(snr_list, gamma_BD,'r');
xlabel('SNR');
ylabel('Threshold');
title('Detectors Comparison')
legend('NP Threshold','Bayes Threshold');

end
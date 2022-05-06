function [T_h0, T_h1] = generate_stat(STAT_COUNT,sigma2_s,N,sigma2_w)

%     Generate H0 test statistics
T_h0 = zeros(STAT_COUNT,1);

for i = 1:STAT_COUNT
    w = sqrt(sigma2_w/2)*(randn(N,1)+1j*randn(N,1));
    y = w;
    %Calculate test statistic
    %using abs when using complex signal
    T_h0(i,1) = sum(abs(y).^2);
end
%disp(T_h0);

%     Generate H1 test statistics
T_h1 = zeros(STAT_COUNT,1);

for i = 1:STAT_COUNT
    w = sqrt(sigma2_w/2)*(randn(N,1)+1j*randn(N,1));
    x = sqrt(sigma2_s)*randn(N, 1);
    y = x + w;
    %Calculate test statistic
    T_h1(i,1) = sum(abs(y).^2);
end

end
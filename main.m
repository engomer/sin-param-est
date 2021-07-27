clear;clc;
%% 1)
A = 1;
fs = 500;
sigma_sq = 1;
phi = pi/4;
fc = 100;
f0 = fc/fs;
ts = 1/f0;
N = 4*ts;
x = zeros();
s = zeros();
for n = 1:N
    s(n,1) = A*cos(2*pi*f0*(n-1)+ phi);
    x(n,1) = s(n,1) + sigma_sq*randn(1,1);
end

figure
plot(x);
hold on;
plot(s);
title("Noisy sinusoidal & Exact signal");
legend(["Noisy Signal","Exact Signal"]);
grid on;
fprintf("N is %d which is 4 times of the period of the signal and the period is %d\n",N,ts);
%% 2)
L = N;
Y = fft(x);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = (0:(L/2))/L;
figure;
plot(f,P1);
title("FFT of x[n]");
grid on;
xlabel("Frequency (Hz)");
ylabel("Amplitude");
f0_est = f(P1 == max(P1));
%% 3)
alpha = [A*cos(phi) -A*sin(phi)]';
H = zeros();
for i=1:N
    H(i,1) = cos(2*pi*f0_est*(i-1));
    H(i,2) = sin(2*pi*f0_est*(i-1));
end

alpha_est = (H'*H)^(-1)*H'*x;
A_est = sqrt(alpha_est(1,1)^2 + alpha_est(2,1)^2);
phi_est = atan(-alpha_est(2,1)/alpha_est(1,1));

%% 4)
x_est = zeros();
for n=1:N
    x_est(n,1) = A_est*cos(2*pi*f0_est*(n-1)+ phi_est);
end

figure
plot(x);
hold on;
plot(x_est);
grid on;
plot(s);
title("Estimated Signal");
ylabel("Amplitude");
xlabel("Samples");
legend(["Noisy Signal","Reconstructed Signal","Exact Signal"]);

%% 5)
A_est_arr = zeros();
phi_est_arr = zeros();
f0_est_arr = zeros();

for ii = 1: 100
    for n = 1:N
        x(n,1) = A*cos(2*pi*f0*(n-1) + phi) + sigma_sq*randn(1,1);
    end
    
    L = N;
    Y = fft(x);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = (0:(L/2))/L;
    f0_est_arr(ii) = f(P1 == max(P1));
    
    alpha = [A*cos(phi) -A*sin(phi)]';
    H = zeros();
    for i=1:N
        H(i,1) = cos(2*pi*f0_est_arr(ii)*(i-1));
        H(i,2) = sin(2*pi*f0_est_arr(ii)*(i-1));
    end
    
    alpha_est = (H'*H)^(-1)*H'*x;
    A_est_arr(ii) = sqrt(alpha_est(1,1)^2 + alpha_est(2,1)^2);
    phi_est_arr(ii) = atan(-alpha_est(2,1)/alpha_est(1,1));
end

mse_A = 0;
mse_phi = 0;
mse_f0 = 0;
for ii=1:100
    mse_A = mse_A + (A - A_est_arr(ii))^2;
    mse_phi = mse_phi + (phi - phi_est_arr(ii))^2;
    mse_f0 = mse_f0 + (f0 - f0_est_arr(ii))^2;
end

mse_A = 1/100 .* mse_A;
mse_f0 = 1/100 .* mse_f0;
mse_phi = 1/100 .* mse_phi;

%% 6)
ii = 1;
mse_arr = zeros();
ssq = 0.01:0.1:1;
for iii = ssq
    ests = zeros();
    for jj = 1: 100
        ests(jj,1:3) = sinusoidal_est(A,fc,fs,phi,N,iii);
    end
    mse_arr(ii,1:3) = mse_calc(ests(:,1)',ests(:,2)',ests(:,3)',A,phi,f0);
    ii = ii + 1;
end
figure
plot(1./ssq,mse_arr(:,1));
hold on;
plot(1./ssq,mse_arr(:,2));
plot(1./ssq,mse_arr(:,3));
grid on;
xlabel("1/\sigma^2");
ylabel("MSE");
title("MSE with varying \sigma^2");
legend(["MSE F_0","MSE A","MSE \phi"]);

%% 7)
ii = 1;
mse_arr = zeros();
NN = 50:10:1000;
for iii = NN
    ests = zeros();
    for jj = 1: 100
        ests(jj,1:3) = sinusoidal_est(A,fc,fs,phi,iii,1);
    end
    mse_arr(ii,1:3) = mse_calc(ests(:,1)',ests(:,2)',ests(:,3)',A,phi,f0);
    ii = ii + 1;
end
figure
plot(NN,mse_arr(:,1));
hold on;
plot(NN,mse_arr(:,2));
plot(NN,mse_arr(:,3));
grid on;
title("MSE with varying N");
xlabel("N");
ylabel("MSE");
legend(["MSE F_0","MSE A","MSE \phi"]);
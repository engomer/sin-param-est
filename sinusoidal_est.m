function res = sinusoidal_est(A,fc,fs,phi,N,sigma_sq)

x = zeros();
for n = 1:N
    x(n,1) = A*cos(2*pi*fc/fs*(n-1) + phi) + sigma_sq*randn(1,1);
end

L = N;
Y = fft(x);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = (0:(L/2))/L;
f0_est = f(P1 == max(P1));

alpha = [A*cos(phi) -A*sin(phi)]';
H = zeros();
for i=1:N
    H(i,1) = cos(2*pi*f0_est*(i-1));
    H(i,2) = sin(2*pi*f0_est*(i-1));
end

alpha_est = (H'*H)^(-1)*H'*x;
A_est = sqrt(alpha_est(1,1)^2+ alpha_est(2,1)^2);
phi_est = atan(-alpha_est(2,1)/alpha_est(1,1));

res = [f0_est A_est phi_est];
end
function res = mse_calc(f0_arr,A_arr,phi_arr,A,phi,f0)
[~,s] = size(A_arr);

mse_A = 0;
mse_phi = 0;
mse_f0 = 0;
for ii=1:s
    mse_A = mse_A + (A - A_arr(ii))^2;
    mse_phi = mse_phi + (phi - phi_arr(ii))^2;
    mse_f0 = mse_f0 + (f0 - f0_arr(ii))^2;
end

res = [1/s.*mse_f0 1/s.*mse_A 1/s.*mse_phi];

end
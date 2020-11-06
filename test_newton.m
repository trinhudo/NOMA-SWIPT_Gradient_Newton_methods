% This is a test code of Newton's method
% but I forgot for what problem
% Tri Nhu

format long
k = 0;
x_k = 0.001;
d_k = 100;
s_th = 10^-6;
a_k = 1;
while d_k > s_th
    G_k = 7 - 1/x_k;
    H_k = 1/(x_k^2);
    d_k = - G_k/H_k; % Gradient method does not work for this problem
    k = k+1;
    x_k = x_k + a_k*d_k;
end

disp(['k = ' num2str(k)])
disp(['x_k = ' num2str(x_k)])
    
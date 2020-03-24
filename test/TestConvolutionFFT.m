%
% Script to test FFT implementation.
%

clear all

% test vector of length 2^n
vd_x = [1:4 1:4];
vd_y = [5:8 5:8];
i_n  = numel(vd_x);

% cyclic convolution through FFT
vd_trans1 = 1/i_n * FFT(i_n^2 * (FFT(vd_x, true) .* FFT(vd_y, true)), false);
display(['FFT implementation: ' num2str(vd_trans1)])

% naive summation
vd_indices = 0:(i_n-1);
vd_trans2  = zeros(1, i_n);
for l = 1 : numel(vd_indices)
    tmp = 0;
    for k = 1 : i_n
       tmp = tmp + vd_x(k) * vd_y(mod((l-k), i_n)+1); 
    end
    vd_trans2(l) = tmp;
end
display(['naive summation: ' num2str(vd_trans2)])

% plot results
figure
plot(real(vd_trans1))
hold on
plot(real(vd_trans2))
title('convolution through FFT test: plot of real parts convolutions')
legend({'FFT implementation', 'naive summation'})

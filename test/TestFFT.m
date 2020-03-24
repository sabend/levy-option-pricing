%
% Script to test FFT implementation.
%

clear all

% test vector of length 2^n
vd_test = 1:8;
i_n     = numel(vd_test);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FFT

% FFT implementation
vd_trans1 = FFT(vd_test, false);
display(['FFT implementation: ' num2str(vd_trans1)])

% naive summation
vd_indices = 0:(i_n-1);
vd_trans2 = arrayfun(@(d_l) sum(vd_test .* exp(2*pi*1i/i_n*vd_indices*d_l)), vd_indices);
display(['naive summation: ' num2str(vd_trans2)])

% plot results
figure
plot(real(vd_trans1))
hold on
plot(real(vd_trans2))
title('FFT implementation test: plot of real parts of transforms')
legend({'FFT implementation', 'naive summation'})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inverse FFT

% iFFT implementation
vd_trans3 = FFT(vd_test, true);
display(['iFFT implementation: ' num2str(vd_trans3)])

% naive summation
vd_indices = 0:(numel(vd_test)-1);
vd_trans4 = 1/i_n * arrayfun(@(d_l) sum(vd_test .* exp(-2*pi*1i/i_n*vd_indices*d_l)), vd_indices);
display(['naive summation: ' num2str(vd_trans4)])

% plot results
figure
plot(real(vd_trans3))
hold on
plot(real(vd_trans4))
title('inverse FFT implementation test: plot of real parts of transforms')
legend({'iFFT implementation', 'naive summation'})

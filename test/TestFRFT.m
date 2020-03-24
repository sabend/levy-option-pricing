%
% Script to test FRFT implementation.
%

clear all

d_beta = -0.45;

% test vector of length 2^n
vd_test = 1:8;
i_n     = numel(vd_test);

% FRFT implementation
vd_trans1 = FRFT(vd_test, d_beta);
display(['FRFT implementation: ' num2str(vd_trans1)])

% naive summation
vd_indices = 0:(i_n-1);
vd_trans2 = arrayfun(@(d_l) sum(vd_test .* exp(2*pi*1i*vd_indices*d_l*d_beta)), vd_indices);
display(['naive summation: ' num2str(vd_trans2)])

% plot results
figure
plot(real(vd_trans1))
hold on
plot(real(vd_trans2))
title('FRFT implementation test: plot of real parts of transforms')
legend({'FRFT implementation', 'naive summation'})

%
% Script to test object model for characteristic functions.
%

clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Black-Scholes

o_bs = CF_BlackScholes(0.01, 2, 5, struct('sigma', 0.2));

v_x = linspace(-5, 5, 100);
v_y = o_bs.eval(v_x);

figure
hold on
plot(v_x, real(v_y))
plot(v_x, imag(v_y))
title('Characteristic function of Black-Scholes model')
legend({'real part' 'imaginary part'})
hold off

o_bs.calcCumulants();
display(['cumulant 1: ' num2str(o_bs.cumulant1())])
display(['cumulant 2: ' num2str(o_bs.cumulant2())])
display(['cumulant 4: ' num2str(o_bs.cumulant4())])


o_bs.reset(struct('sigma', 0.15'));
o_bs.calcCumulants();
display(['cumulant 1: ' num2str(o_bs.cumulant1())])
display(['cumulant 2: ' num2str(o_bs.cumulant2())])
display(['cumulant 4: ' num2str(o_bs.cumulant4())])











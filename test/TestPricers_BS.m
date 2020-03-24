%
% Script to test object model for pricers.
%

clear all

% create sample strike vector
vd_strikes = linspace(0.5, 4, 10);

% create sample characteristic function
o_bs = CF_BlackScholes(0.01, 2, 5, struct('sigma', 0.2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COS pricer

str_cos_params = struct('N', 64, ...
                        'L', 10);

o_cos_pricer = Pricer_COS(o_bs, str_cos_params);

display(['COS prices: ' num2str(o_cos_pricer.run(vd_strikes))]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analytic Black-Scholes pricer

o_bs_pricer = Pricer_BS(o_bs, struct());

display(['Black-Scholes prices: ' num2str(o_bs_pricer.run(vd_strikes))]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fourier transform pricer

str_ft_params = struct('N', 64,     ...
                       'eta', 0.25, ...
                       'alpha', -4);

o_ft_pricer = Pricer_FT(o_bs, str_ft_params);

display(['Fourier transform price: ' num2str([o_ft_pricer.run(vd_strikes(1)) o_ft_pricer.run(vd_strikes(end))])]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fast Fourier transform pricer

str_fft_params = struct('N', 128,    ...
                        'eta', 0.25, ...
                        'alpha', -4);

o_fft_pricer = Pricer_FFT(o_bs, str_fft_params);

display(['FFT prices: ' num2str(o_fft_pricer.run(vd_strikes))]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fractional Fast Fourier transform pricer

str_frft_params = struct('N', 64,     ...
                         'eta', 0.25, ...
                         'alpha', -4);

o_frft_pricer = Pricer_FRFT(o_bs, str_frft_params);

display(['FRFT prices: ' num2str(o_frft_pricer.run(vd_strikes))]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% quadgk based pricer

str_qgk_params = struct('N', 20000, 'alpha', -4);

o_qgk_pricer = Pricer_QGK(o_bs, str_qgk_params);

display(['QGK prices: '  num2str([o_qgk_pricer.run(vd_strikes(1)) o_qgk_pricer.run(vd_strikes(end))])]);






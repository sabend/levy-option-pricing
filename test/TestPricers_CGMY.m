%
% Script to test object model for pricers.
%

clear all

vd_strikes = [2 2.5];

% create sample characteristic function
o_cgmy = CF_CGMY(0.01, 2, 5, struct('C', 2, 'G', 20, 'M', 20, 'Y', 1.2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COS pricer

str_cos_params = struct('N', 64, ...
                        'L', 10);

o_cos_pricer = Pricer_COS(o_cgmy, str_cos_params);

display(['COS prices: ' num2str(o_cos_pricer.run(vd_strikes))]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fourier transform pricer

str_ft_params = struct('N', 64,     ...
                       'eta', 0.25, ...
                       'alpha', -4);

o_ft_pricer = Pricer_FT(o_cgmy, str_ft_params);

display(['Fourier transform price: ' num2str(o_ft_pricer.run(vd_strikes))]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fast Fourier transform pricer

str_fft_params = struct('N', 512,   ...
                       'eta', 0.25, ...
                       'alpha', -4);

o_fft_pricer = Pricer_FFT(o_cgmy, str_fft_params);

display(['FFT prices: ' num2str(o_fft_pricer.run(vd_strikes))]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fractional Fast Fourier transform pricer

str_frft_params = struct('N', 128,    ...
                         'eta', 0.25, ...
                         'alpha', -4);

o_frft_pricer = Pricer_FRFT(o_cgmy, str_frft_params);

display(['FRFT prices: ' num2str(o_frft_pricer.run(vd_strikes))]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% quadgk based pricer

str_qgk_params = struct('N', 20000, 'alpha', -4);

o_qgk_pricer = Pricer_QGK(o_cgmy, str_qgk_params);

display(['QGK prices: ' num2str([o_qgk_pricer.run(vd_strikes(1)) o_qgk_pricer.run(vd_strikes(end))])]);






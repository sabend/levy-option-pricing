%
% Script to compare performance of several pricing methods.
%

clear all

% script parameters
vi_N   = 2.^(2:15);
i_iter = 10;

% sample strikes
vd_ref_strike = linspace(1.8, 2.2, 100);

% create characteristic function
o_bs = CF_BlackScholes(0.01, 2, 3, struct('sigma', 0.2));

% create analytics Black-Scholes reference pricer
o_bs_pricer = Pricer_BS(o_bs, struct());

% calculate Black-Scholes reference price
vd_bs_price = o_bs_pricer.run(vd_ref_strike);

% initialise result cubes
m3_cos_results  = zeros(2, numel(vi_N), i_iter);
m3_ft_results   = zeros(2, numel(vi_N), i_iter);
m3_fft_results  = zeros(2, numel(vi_N), i_iter);
m3_frft_results = zeros(2, numel(vi_N), i_iter);
for j = 1 : i_iter
    for i = numel(vi_N) : -1 : 1

        % COS pricer
        o_cos_pricer = Pricer_COS(o_bs, struct('N', vi_N(i), 'L', 10));
        tic
        vd_cos_price = o_cos_pricer.run(vd_ref_strike);
        m3_cos_results(1, i, j) = toc;
        m3_cos_results(2, i, j) = sum(abs(vd_cos_price - vd_bs_price));

        % FT pricer
        o_ft_pricer = Pricer_FT(o_bs, struct('N', vi_N(i), 'eta', 0.25, 'alpha', -4));
        tic
        vd_ft_price  = o_ft_pricer.run(vd_ref_strike);
        m3_ft_results(1, i, j) = toc;
        m3_ft_results(2, i, j) = sum(abs(vd_ft_price - vd_bs_price));

        % FFT pricer
        o_fft_pricer = Pricer_FFT(o_bs, struct('N', vi_N(i), 'eta', 0.25, 'alpha', -4));
        tic
        vd_fft_price  = o_fft_pricer.run(vd_ref_strike);
        m3_fft_results(1, i, j) = toc;
        m3_fft_results(2, i, j) = sum(abs(vd_fft_price - vd_bs_price));

        % FRFT pricer
        o_frft_pricer = Pricer_FRFT(o_bs, struct('N', vi_N(i), 'eta', 0.25, 'alpha', -4));
        tic
        vd_frft_price  = o_frft_pricer.run(vd_ref_strike);
        m3_frft_results(1, i, j) = toc;
        m3_frft_results(2, i, j) = sum(abs(vd_frft_price - vd_bs_price));

    end
end

% average runs
m_cos_results  = mean(m3_cos_results(:, :, :), 3);
m_ft_results   = mean(m3_ft_results(:, :, :), 3);
m_fft_results  = mean(m3_fft_results(:, :, :), 3);
m_frft_results = mean(m3_frft_results(:, :, :), 3);

% plot cpu time vs N
figure
loglog(vi_N, m_ft_results(1, :), 'color', 'blue', 'linestyle', '-', 'marker', '+')
hold on
loglog(vi_N, m_fft_results(1, :), 'color', 'red', 'linestyle', '-', 'marker', 'o')
loglog(vi_N, m_frft_results(1, :), 'color', 'black', 'linestyle', '-', 'marker', 'square')
loglog(vi_N, m_cos_results(1, :), 'color', 'green', 'linestyle', '-', 'marker', 'x')
hold off
legend({'pFT', 'pFFT', 'pFRFT', 'pCOS'})
xlabel('N')
ylabel('cpu time')

% plot pricing error vs cpu time
[vd_ft_sorted, vd_ft_ind]     = sort(m_ft_results(1, :));
[vd_fft_sorted, vd_fft_ind]   = sort(m_fft_results(1, :));
[vd_frft_sorted, vd_frft_ind] = sort(m_frft_results(1, :));
[vd_cos_sorted, vd_cos_ind]   = sort(m_cos_results(1, :));
figure
loglog(vd_ft_sorted, m_ft_results(2, vd_ft_ind), 'color', 'blue', 'linestyle', '-', 'marker', '+')
hold on
loglog(vd_fft_sorted, m_fft_results(2, vd_fft_ind), 'color', 'red', 'linestyle', '-', 'marker', 'o')
loglog(vd_frft_sorted, m_frft_results(2, vd_frft_ind), 'color', 'black', 'linestyle', '-', 'marker', 'square')
loglog(vd_cos_sorted, m_cos_results(2, vd_cos_ind), 'color', 'green', 'linestyle', '-', 'marker', 'x')
hold off
legend({'pFT', 'pFFT', 'pFRFT', 'pCOS'})
xlabel('cpu time')
ylabel('pricing error')

% plot pricing error vs N
figure
loglog(vi_N, m_ft_results(2, :), 'color', 'blue', 'linestyle', '-', 'marker', '+')
hold on
loglog(vi_N, m_fft_results(2, :), 'color', 'red', 'linestyle', '-', 'marker', 'o')
loglog(vi_N, m_frft_results(2, :), 'color', 'black', 'linestyle', '-', 'marker', 'square')
loglog(vi_N, m_cos_results(2, :), 'color', 'green', 'linestyle', '-', 'marker', 'x')
hold off
legend({'pFT', 'pFFT', 'pFRFT', 'pCOS'})
xlabel('N')
ylabel('pricing error')














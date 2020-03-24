%
% Analysis of dependence of COS method's performance on the truncation
% bounds
% Models: Black-Scholes and CGMY
%

clear all
rng(309)
b_save = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% > PREPARATION

i_N_param = 100;

% market parameters
d_r     = 0.01;
d_stock = 5;

% COS parameters
d_L      = 10;
vi_N_cos = [10:10:100 150:50:950 1000:500:10000 20000];

% QGK parameters
str_qgk_params = struct('N', 1000000, 'alpha', -1.1);

% create random parameter samples
vd_sigma  = 0.05 + 0.9*rand(1, i_N_param);
vd_C      = 0.01 + 1.99*rand(1, i_N_param);
vd_G      = 1 + 4*rand(1, i_N_param);
vd_M      = 1 + 4*rand(1, i_N_param);
vd_Y      = 0.1 + 0.8*rand(1, i_N_param);
vd_T      = 1 + 2*rand(1, i_N_param);
vd_strike = 4 + 2*rand(1, i_N_param);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% > CALCULATION

% 1. calculate truncation bounds for each parameter set for both
%    Black-Scholes and CGMY
% 2. evaluate analytical Black-Scholes formula
% 3. run QGK pricer for CGMY with high precision
vd_zeros  = zeros(1, i_N_param);
vd_a_bs   = vd_zeros;
vd_b_bs   = vd_zeros;
vd_a_cgmy = vd_zeros;
vd_b_cgmy = vd_zeros;
vd_bs_price       = vd_zeros;
vd_cgmy_qgk_price = vd_zeros;
for i = 1 : i_N_param
    
    d_T      = vd_T(i);
    d_strike = vd_strike(i);
    d_tmp    = log(d_stock/d_strike);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1.
    
    % Black-Scholes
    d_sigma = vd_sigma(i);
    
    o_bs = CF_BlackScholes(d_r, d_stock, d_T, struct('sigma', d_sigma));
    o_bs.calcCumulants();
    
    vd_cumulants_bs = o_bs.cumulants();
    d_c1  = vd_cumulants_bs(1);
    d_c2  = vd_cumulants_bs(2);
    d_c4  = vd_cumulants_bs(4);
    d_c6  = vd_cumulants_bs(6);
    d_c8  = vd_cumulants_bs(8);
    d_c10 = vd_cumulants_bs(10);
    d_c12 = vd_cumulants_bs(12);
    d_c14 = vd_cumulants_bs(14);
    d_c16 = vd_cumulants_bs(16);
    
    d_tmp2 = d_L*sqrt(d_c2+sqrt(d_c4+sqrt(d_c6+sqrt(d_c8+sqrt(d_c10+ ...
             sqrt(d_c12+sqrt(d_c14+sqrt(d_c16))))))));
    vd_a_bs(i) = d_c1 + d_tmp - d_tmp2;
    vd_b_bs(i) = d_c1 + d_tmp + d_tmp2;
    
    % CGMY
    d_C = vd_C(i);
    d_G = vd_G(i);
    d_M = vd_M(i);
    d_Y = vd_Y(i);
    
    o_cgmy = CF_CGMY(d_r, d_stock, d_T, struct('C', d_C, 'G', d_G, 'M', d_M, 'Y', d_Y));
    o_cgmy.calcCumulants();
    
    vd_cumulants_cgmy = o_cgmy.cumulants();
    d_c1  = vd_cumulants_cgmy(1);
    d_c2  = vd_cumulants_cgmy(2);
    d_c4  = vd_cumulants_cgmy(4);
    d_c6  = vd_cumulants_cgmy(6);
    d_c8  = vd_cumulants_cgmy(8);
    d_c10 = vd_cumulants_cgmy(10);
    d_c12 = vd_cumulants_cgmy(12);
    d_c14 = vd_cumulants_cgmy(14);
    d_c16 = vd_cumulants_cgmy(16);

    d_tmp2 = d_L*sqrt(d_c2+sqrt(d_c4+sqrt(d_c6+sqrt(d_c8+sqrt(d_c10+ ...
             sqrt(d_c12+sqrt(d_c14+sqrt(d_c16))))))));
    vd_a_cgmy(i) = d_c1 + d_tmp - d_tmp2;
    vd_b_cgmy(i) = d_c1 + d_tmp + d_tmp2;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 2.
    
    % Black-Scholes analytical value
    o_bs_pricer = Pricer_BS(o_bs, struct());
    vd_bs_price(i) = o_bs_pricer.run(d_strike);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 3.
    
    % QGK pricing for CGMY
    o_qgk_pricer = Pricer_QGK(o_cgmy, str_qgk_params);
    vd_cgmy_qgk_price(i) = o_qgk_pricer.run(d_strike);   

end

% save FT reference prices for CGMY to file
if b_save
    save('Save\FT_CGMY_Prices.mat', 'vd_cgmy_qgk_price')
end
    
% 4. run COS pricer for Black-Scholes and CGMY with varying N
m_bs_cos_price   = zeros(numel(vi_N_cos), i_N_param);
m_cgmy_cos_price = zeros(numel(vi_N_cos), i_N_param);

for j = 1 : numel(vi_N_cos)
    
    i_N_cos = vi_N_cos(j);
    
    for i = 1 : i_N_param

        d_T      = vd_T(i);
        d_strike = vd_strike(i);

        % Black-Scholes COS pricing
        o_bs = CF_BlackScholes(d_r, d_stock, d_T, struct('sigma', vd_sigma(i)));
        o_cos_pricer_bs = Pricer_COS(o_bs, struct('N', i_N_cos, 'L', d_L));
        m_bs_cos_price(j, i) = o_cos_pricer_bs.run(d_strike);

        % CGMY COS pricing
        o_cgmy = CF_CGMY(d_r, d_stock, d_T, struct('C', vd_C(i), 'G', vd_G(i), 'M', vd_M(i), 'Y', vd_Y(i)));
        o_cos_pricer_cgmy = Pricer_COS(o_cgmy, struct('N', i_N_cos, 'L', d_L));
        m_cgmy_cos_price(j, i) = o_cos_pricer_cgmy.run(d_strike);

    end

end

% 5. calculate maximal, average and minimal pricing errors
vd_bs_error_max   = transpose(max(abs(bsxfun(@minus, m_bs_cos_price, vd_bs_price)), [], 2));
vd_cgmy_error_max = transpose(max(abs(bsxfun(@minus, m_cgmy_cos_price, vd_cgmy_qgk_price)), [], 2));

vd_bs_error_min   = transpose(min(abs(bsxfun(@minus, m_bs_cos_price, vd_bs_price)), [], 2));
vd_cgmy_error_min = transpose(min(abs(bsxfun(@minus, m_cgmy_cos_price, vd_cgmy_qgk_price)), [], 2));

vd_bs_error_avg   = transpose(mean(abs(bsxfun(@minus, m_bs_cos_price, vd_bs_price)), 2));
vd_cgmy_error_avg = transpose(mean(abs(bsxfun(@minus, m_cgmy_cos_price, vd_cgmy_qgk_price)), 2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% > PLOTTING

% BS and CGMY maximal pricing errors vs number of COS summands
figure
loglog(vi_N_cos, vd_bs_error_max, 'color', 'blue', 'linestyle', '-')
hold on
loglog(vi_N_cos, vd_bs_error_avg, 'color', 'blue', 'linestyle', '--')
loglog(vi_N_cos, vd_bs_error_min, 'color', 'blue', 'linestyle', '-.')
loglog(vi_N_cos, vd_cgmy_error_max, 'color', 'magenta', 'linestyle', '-')
loglog(vi_N_cos, vd_cgmy_error_avg, 'color', 'magenta', 'linestyle', '--')
loglog(vi_N_cos, vd_cgmy_error_min, 'color', 'magenta', 'linestyle', '-.')
hold off
xlabel('N')
ylabel('pricing error')
legend({'BS: max error', 'BS: avg error', 'BS: min error', ...
        'CGMY: max error', 'CGMY: avg error', 'CGMY: min error'})

% pricing error accross random parameters sets for highest COS summand
figure
semilogy(abs(m_cgmy_cos_price(end, :)-vd_cgmy_qgk_price))
xlabel('parameter set')
ylabel('pricing error')

% pricing error for highest COS summand vs b-a
vd_bs_width = vd_b_bs - vd_a_bs;
[vd_bs_width, vd_index_bs] = sort(vd_bs_width);
vd_bs_last  = abs(m_bs_cos_price(end, :)-vd_bs_price);
vd_cgmy_width = vd_b_cgmy - vd_a_cgmy;
[vd_cgmy_width, vd_index_cgmy] = sort(vd_cgmy_width);
vd_cgmy_last  = abs(m_cgmy_cos_price(end, :)-vd_cgmy_qgk_price);
figure
semilogy(vd_bs_width, vd_bs_last(vd_index_bs))
hold on
semilogy(vd_cgmy_width, vd_cgmy_last(vd_index_cgmy))
hold off
xlabel('b - a')
ylabel('pricing error')
legend({'Black-Scholes', 'CGMY'})


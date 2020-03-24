%
% Script to plot the calibration target function for
% the Black-Scholes model.
%

clear all

% basic parameters
d_sigma       = 0.2;
d_r           = 0.01;
d_stock       = 3000;
vd_maturities = [1.7 2.4 3.4 4.5 6];
vd_strikes    = 2500:50:3500;

% create characteristic function
o_bs = CF_BlackScholes(d_r, d_stock, 0, struct('sigma', d_sigma));

% create pricer
o_cos_pricer = Pricer_COS(o_bs, struct('N', 128, 'L', 10));

% simulate option prices
i_count       = numel(vd_maturities);
c_option_data = cell(i_count, 1);
for i = 1 : i_count
        
    d_maturity = vd_maturities(i);
    o_cos_pricer.setMaturityCF(d_maturity);
    vd_tmp = ones(numel(vd_strikes), 1) * d_maturity;
    vd_prices = transpose(o_cos_pricer.run(vd_strikes));
    
    c_option_data{i, 1} = [vd_tmp transpose(vd_strikes) vd_prices];
    
end

% create loss function object
o_loss_function = Loss_Quadratic(c_option_data, 3);

% determine pricing error for varying sigma
vd_sigma = linspace(0, 1, 100);
vd_error = zeros(1, numel(vd_sigma));
for i = 1 : numel(vd_sigma)
    
    o_bs = CF_BlackScholes(d_r, d_stock, 0, struct('sigma', vd_sigma(i)));
    o_cos_pricer = Pricer_COS(o_bs, struct('N', 128, 'L', 10));
    
    c_tmp_prices = cell(i_count, 1);
    for j = 1 : i_count
        
        d_maturity = vd_maturities(j);
        o_cos_pricer.setMaturityCF(d_maturity);
        vd_tmp = ones(numel(vd_strikes), 1) * d_maturity;
        
        vd_prices = transpose(o_cos_pricer.run(vd_strikes));

        c_tmp_prices{j, 1} = [vd_tmp transpose(vd_strikes) vd_prices];
    
    end
    
    vd_error(i) = o_loss_function.calculate(c_tmp_prices);
    
end

% plot calibration target function
figure
plot(vd_sigma, vd_error)
xlabel('\sigma')
ylabel('pricing error')

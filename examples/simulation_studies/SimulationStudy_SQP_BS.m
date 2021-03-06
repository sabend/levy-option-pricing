%
% Script for simulation study on calibration implementation.
% Option prices are simulated from the Black-Scholes model.
%

clear all

b_plot_iters = true;

% basic parameters
d_sigma       = 0.1234;
d_r           = 0.01;
d_stock       = 3000;
vd_maturities = 1.7; %[1.7 2.4 3.4 4.5 6];
vd_strikes    = 2500:50:3500;
d_guess       = 0.4;

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

% create calibration engine
o_sqp = Calibration_SQP(o_cos_pricer, o_loss_function, ...
                        struct(), c_option_data, d_guess, true);

% run calibration
o_sqp.run();

% get and print calibration results
str_results = o_sqp.result();
display(str_results)

% plot market data and calibrated option prices
plotCalibrationResult_Single(c_option_data, o_sqp.surface(), b_plot_iters);

% save iteration results
t_bs_log = o_sqp.log();
writetable(t_bs_log);

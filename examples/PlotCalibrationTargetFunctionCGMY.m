%
% Script to plot the calibration target function for
% the Black-Scholes model.
%

clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% basic parameters
d_C           = 0.8;
d_G           = 5;
d_M           = 7;
d_Y           = 0.6;
d_r           = 0.01;
d_stock       = 3000;
vd_maturities = [1.7 2.4 3.4 4.5 6];
vd_strikes    = 2500:50:3500;

% create characteristic function
str_params = struct('C', d_C, 'G', d_G, 'M', d_M, 'Y', d_Y);
o_cgmy = CF_CGMY(d_r, d_stock, 0, str_params);

% create pricer
o_cos_pricer = Pricer_COS(o_cgmy, struct('N', 512, 'L', 10));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulate option prices
i_count       = numel(vd_maturities);
c_option_data = cell(i_count, 1);
for i = 1 : i_count
        
    d_maturity = vd_maturities(i);
    o_cos_pricer.setMaturityCF(d_maturity);
    vd_tmp = ones(numel(vd_strikes), 1) * d_maturity;
    vd_prices     = transpose(o_cos_pricer.run(vd_strikes));
    
    c_option_data{i, 1} = [vd_tmp transpose(vd_strikes) vd_prices];
    
end

% create loss function object
o_loss_function = Loss_Quadratic(c_option_data, 3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% determine pricing error for varying C
vd_C     = linspace(0.01, 2, 100);
vd_error = zeros(1, numel(vd_C));
for i = 1 : numel(vd_C)
    
    o_cgmy = CF_CGMY(d_r, d_stock, 0, struct('C', vd_C(i), 'G', d_G, 'M', d_M, 'Y', d_Y));
    o_cos_pricer = Pricer_COS(o_cgmy, struct('N', 512, 'L', 10));
    
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
semilogy(vd_C, vd_error)
xlabel('C')
ylabel('pricing error')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% determine pricing error for varying G
vd_G     = linspace(0.01, 20, 100);
vd_error = zeros(1, numel(vd_G));
for i = 1 : numel(vd_G)
    
    o_cgmy = CF_CGMY(d_r, d_stock, 0, struct('C', d_C, 'G', vd_G(i), 'M', d_M, 'Y', d_Y));
    o_cos_pricer = Pricer_COS(o_cgmy, struct('N', 512, 'L', 10));
    
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
semilogy(vd_G, vd_error)
xlabel('G')
ylabel('pricing error')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% determine pricing error for varying Y
vd_Y     = [linspace(0.01, 0.99, 50) linspace(1.01, 1.99, 50)];
vd_error = zeros(1, numel(vd_Y));
for i = 1 : numel(vd_Y)
    
    o_cgmy = CF_CGMY(d_r, d_stock, 0, struct('C', d_C, 'G', d_G, 'M', d_M, 'Y', vd_Y(i)));
    o_cos_pricer = Pricer_COS(o_cgmy, struct('N', 512, 'L', 10));
    
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
semilogy(vd_Y, vd_error)
xlabel('Y')
ylabel('pricing error')

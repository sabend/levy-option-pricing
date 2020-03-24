%
% Script to test object model for calibration engine.
%

% load option data
s_asof = '20160324';
s_folder = [pwd '\Data\EuroSTOXX\' s_asof];
vs_maturities = ['201712'];% '201812'; '201912'; '202012'; '202112'];
c_maturities = cellstr(vs_maturities);
c_option_surface = cell(numel(c_maturities), 1);
for i = 1:numel(c_maturities)
  s_maturity   = c_maturities{i, 1};
  c_option_surface{i, 1} = txt2optTable([s_folder '\EuroSTOXX_Call_' ...
                                         s_maturity '.txt'], s_asof, ...
                                         s_maturity, 2800, 4000);
end

% create loss function object
o_loss_function = Loss_Quadratic(c_option_surface, 3);

% create characteristic function
o_cgmy = CF_CGMY(0.01, 2986.73, 1, struct('C', 1, 'G', 20, 'M', 20, 'y', 0.5));

% create pricer
o_cos_pricer = Pricer_COS(o_cgmy, struct('N', 128, 'L', 10));

% create calibration engine
o_sqp = Calibration_SQP(o_cos_pricer, o_loss_function, ...
                        struct('a', 1, 'b', 2), ...
                        c_option_surface, [3;40;30;0.4], true);

% run calibration
o_sqp.run();

% get and print calibration results
str_results = o_sqp.result();
display(str_results)

% plot market data and calibrated option prices
plotCalibrationResult_Single(c_option_surface, o_sqp.surface());

% ToDo: create object for option surface -> cell array with c{i, 1} = [maturity strike price]












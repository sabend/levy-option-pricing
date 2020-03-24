classdef Calibration_SQP < CalibrationEngine
    %Calibration_SQP: Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        % inherited: m_b_logging
        % inherited: m_c_option_data % {i, 1} = [maturity strike price]
        % inherited: m_o_pricer
        % inherited: m_o_loss_function
        % inherited: m_str_params
        % inherited: m_str_result
        
        m_c_surface % ToDo -> better: m_o_calibrated_surface
        m_d_initial_guess
        m_t_log
        
    end
    
    methods (Access = public)
        
        % constructor
        function obj = Calibration_SQP(o_pricer, o_loss_function, ...
                                       str_params, c_option_data, ...
                                       d_initial_guess, b_logging)
                                    
            obj@CalibrationEngine(o_pricer, o_loss_function, ...
                                  str_params, c_option_data, ...
                                  b_logging);
                              
            obj.m_d_initial_guess = d_initial_guess;
                                     
        end
        
        % run calibration
        function obj = run(obj)
            
            % get characteristic function from pricer
            o_cf = obj.m_o_pricer.cf();
            vd_guess = obj.m_d_initial_guess;
            [vd_g, ~, vd_h, ~] = o_cf.constraints(vd_guess);
            vd_initial_lambda  = [];
            vd_initial_mu      = [];
            if ~isempty(vd_g)
                vd_initial_lambda = zeros(numel(vd_g), 1);
            end
            if ~isempty(vd_h)
                vd_initial_mu = zeros(numel(vd_h), 1);
            end
            
            % optimization parameters
            str_params = struct('initial_x',      vd_guess,          ...
                                'initial_lambda', vd_initial_lambda, ...
                                'initial_mu',     vd_initial_mu,     ...
                                'alpha',          3,                 ...
                                'max_iter',       100,               ...
                                'beta',           0.5,               ...
                                'gamma',          0.1,               ...
                                'epsilon',        10^(-8));
                    
            % instantiate optimiser
            o_sqp = Optimiser_SQP(obj, str_params, obj.m_b_logging);
        
            % run optimiser
            o_sqp.run();

            % store results
            obj.m_str_result = o_sqp.results();
            obj.m_t_log      = o_sqp.log();

        end
        
        % optimisation constraints
        function [out_g, out_gF, out_h, out_hF] = constraints(obj, vd_x)
           
            [out_g, out_gF, out_h, out_hF] = obj.m_o_pricer.cf().constraints(vd_x);
            
        end
        
        % optimisation objective
        function out_error = objective(obj, vd_x, b_store_surface)
          
            o_pricer = obj.m_o_pricer;
            
            % initialise price cell array
            i_count      = numel(obj.m_c_option_data);
            c_cur_prices = cell(i_count, 1);
            
            % update characteristic function parameters
            str_params = o_pricer.cf().fillParameters(vd_x);
            o_pricer.setParamsCF(str_params);
            
            % fill price cell array
            for i = 1 : i_count
    
                d_maturity = obj.m_c_option_data{i, 1}(1, 1);
                o_pricer.setMaturityCF(d_maturity);
                
                vd_strikes    = obj.m_c_option_data{i, 1}(:, 2);
                vd_prices     = transpose(o_pricer.run(transpose(vd_strikes)));
                vd_maturities = ones(numel(vd_strikes), 1) * d_maturity;
                
                c_cur_prices{i, 1} = [vd_maturities vd_strikes vd_prices];

            end
            
            % return pricing error of current iteration
            out_error = obj.m_o_loss_function.calculate(c_cur_prices);
            
            % store current prices
            if b_store_surface
                obj.m_c_surface{numel(obj.m_c_surface)+1} = c_cur_prices;
            end
 
        end
        
        % get calibrated surface
        function s = surface(obj)
            
            s = obj.m_c_surface;
            
        end

        % get log table
        function out_log = log(obj)

            out_log = obj.m_t_log;

        end
        
    end
    
end


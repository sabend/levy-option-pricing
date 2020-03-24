classdef Pricer_BS < CF_Pricer
    %Pricer_BS: Implements the analytic Black-Scholes formula for European calls.
    
    properties
        
        % inherited: m_str_params (struct: empty)
        % inherited: m_o_cf
        
    end
    
    methods (Access = public)
        
        % constructor
        function obj = Pricer_BS(o_cf, str_params)

            obj@CF_Pricer(o_cf, str_params);

        end
        
        % run pricing calculation
        function out_value = run(obj, vd_strikes)
            
            d_stock    = obj.m_o_cf.stock();
            d_r        = obj.m_o_cf.r();
            d_sigma    = obj.m_o_cf.params().sigma;
            d_T        = obj.m_o_cf.T();
        
            d_d1 = (log(d_stock ./ vd_strikes) + (d_r + 0.5*d_sigma^2)*d_T) ...
                   / (d_sigma * sqrt(d_T));
            d_d2 = d_d1 - d_sigma * sqrt(d_T);
            
            out_value = d_stock * normcdf(d_d1) ...
                        - vd_strikes * exp(-d_r*d_T) .* normcdf(d_d2);
        
        end
        
    end
    
end


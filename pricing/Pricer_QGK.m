classdef Pricer_QGK < CF_Pricer
    %Pricer_QGK: Fourier pricing with integration by 'quadgk'.
    
    properties
        
        % inherited: m_str_params (struct: empty)
        % inherited: m_o_cf
        
    end
    
    methods (Access = public)
        
        % constructor
        function obj = Pricer_QGK(o_cf, str_params)

            obj@CF_Pricer(o_cf, str_params);

        end
        
        % run pricing calculation
        function out_value = run(obj, vd_strikes)
            
            d_r        = obj.m_o_cf.r();
            d_T        = obj.m_o_cf.T();
            d_stock    = obj.m_o_cf.stock();
            d_alpha    = obj.m_str_params.alpha;
            i_N        = obj.m_str_params.N;
            o_cf       = obj.m_o_cf;
            
            o_cf.setStock(1);
            
            f_eval = @(x) real(exp(1i*x*log(d_stock)) ...
                          .* o_cf.eval(x+d_alpha*1i) ...
                          .* vd_strikes.^(d_alpha+1-1i*x) ...
                          ./ (d_alpha+1-1i*x) ./ (d_alpha-1i*x));

            out_value = exp(-d_alpha*log(d_stock)-d_r*d_T)/(2*pi) ...
                        .* quadgk(f_eval, -Inf, Inf, 'AbsTol', 10^(-12), ...
                                  'RelTol', 10^(-12), 'MaxIntervalCount', i_N);

            o_cf.setStock(d_stock);
            
        end
        
    end
    
end


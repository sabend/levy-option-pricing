classdef Pricer_FT < CF_Pricer
    %Pricer_FT: Implements Fourier pricing for European calls.
    %   Uses trapezoidal rule.
    %   References:
    %   - P. Carr & D. B. Madan 1999.
    %   - E. Eberlein, K. Glau & A. Papapantoleon 2010.
    
    properties
        
        % inherited: m_str_params (struct: N, eta, alpha)
        % inherited: m_o_cf
        
    end
    
    methods (Access = public)
        
        % constructor
        function obj = Pricer_FT(o_cf, str_params)

            obj@CF_Pricer(o_cf, str_params);

        end
        
        % run pricing calculation
        function out_value = run(obj, vd_strikes)
            
            % initialise
            vd_sum     = 0;
            vd_k       = log(vd_strikes);
            d_r        = obj.m_o_cf.r();
            d_stock    = obj.m_o_cf.stock();
            d_T        = obj.m_o_cf.T();

            % set parameter stock to 1
            obj.m_o_cf.setStock(1);

            % Sum over grid. Multiply first and last summand by 0.5 to 
            % account for the trapezoidal approximation.
            for i = 0 : (obj.m_str_params.N-1)

                d_v = obj.m_str_params.eta * i;
                
                d_psi       = obj.m_o_cf.eval(d_v+1i*obj.m_str_params.alpha);
                d_psi       = d_psi ./ (1+obj.m_str_params.alpha-1i*d_v) ...
                              ./ (obj.m_str_params.alpha-1i*d_v);
                          
                vd_integrand = exp(1i*d_v*(log(d_stock)-vd_k))*d_psi;
                
                if (i == 0 || i == (obj.m_str_params.N-1)) 
                    vd_integrand = 0.5 * vd_integrand;
                end
                vd_sum = vd_sum + vd_integrand;

            end

            % reset parameter stock
            obj.m_o_cf.setStock(d_stock);

            % calculate option value
            out_value = obj.m_str_params.eta ./ pi ...
                      .* exp((1+obj.m_str_params.alpha).*vd_k ...
                      - obj.m_str_params.alpha*log(d_stock) ...
                      - d_r*d_T) .* real(vd_sum);

        end
        
    end
    
end


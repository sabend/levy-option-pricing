classdef Pricer_FRFT < CF_Pricer
    %Pricer_FRFT: Implements the FFT pricing methodology for European calls.
    %   Reference:  1999
    
    properties
        
        % inherited: m_str_params (struct: N, eta, alpha)
        % inherited: m_o_cf
        
    end
    
    methods (Access = public)
        
        % constructor
        function obj = Pricer_FRFT(o_cf, str_params)

            obj@CF_Pricer(o_cf, str_params);

        end
        
        % run pricing calculation
        function [out_values, out_strikes] = run(obj, vd_strikes)
            
            % get relevant parameters of characteristic function
            d_r     = obj.m_o_cf.r();
            d_stock = obj.m_o_cf.stock();
            d_T     = obj.m_o_cf.T();
            
            % set parameter stock to 1
            obj.m_o_cf.setStock(1);
            
            % set up strike vector
            i_N            = obj.m_str_params.N;
            d_log_lbound   = log(min(vd_strikes) * 0.98);
            d_log_ubound   = log(max(vd_strikes) * 1.02);
            d_lambda       = (d_log_ubound - d_log_lbound) / (i_N - 1);
            vi_indices     = 0:(i_N-1);
            vd_log_strikes = d_log_lbound + d_lambda * vi_indices;
            out_strikes    = exp(vd_log_strikes);
            
            % determine beta for Fractional FFT algorithm.
            d_eta  = obj.m_str_params.eta;
            d_beta = d_lambda * d_eta / 2 / pi;
            
            % evaluate dampened Fourier transform on grid
            % multiply first and last value by 0.5 due to trapezoidal approximation
            d_alpha   = obj.m_str_params.alpha;
            vd_grid   = d_eta * vi_indices;
            vd_phi    = obj.m_o_cf.eval(vd_grid+1i*d_alpha);
            vd_series = exp(1i*vd_grid*(log(d_stock) - d_log_lbound));
            vd_series = vd_series ./ (1+d_alpha-1i*vd_grid) ./ ...
                        (d_alpha-1i*vd_grid);
            vd_series = vd_series .* vd_phi;
            vd_series(1)   = 0.5 * vd_series(1);
            vd_series(i_N) = 0.5 * vd_series(i_N);

            % perform inverse FRFT.
            vd_transform = real(FRFT(vd_series, -d_beta));

            % Calculate option prices. Multiply by N due to the inverse FRFT.
            vd_values = d_eta / pi * exp((1+d_alpha)*vd_log_strikes-d_alpha ...
                        * log(d_stock)-d_r*d_T) .* vd_transform;
            
            % reset parameter stock
            obj.m_o_cf.setStock(d_stock);

            % interpolate option prices of interest
            out_values = real(extractInterpolatedElements(vd_strikes, out_strikes, vd_values, 'loglinear'));
            %out_values = interp1(out_strikes, vd_values, vd_strikes, 'cubic');
            
        end
        
    end
    
end


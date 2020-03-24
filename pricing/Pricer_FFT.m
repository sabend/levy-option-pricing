classdef Pricer_FFT < CF_Pricer
    %Pricer_FFT: Implements the FFT pricing methodology for European calls.
    %   Reference: P. Carr & D. B. Madan 1999
    
    properties
        
        % inherited: m_str_params (struct: N, eta, alpha)
        % inherited: m_o_cf
        
    end
    
    methods (Access = public)
        
        % constructor
        function obj = Pricer_FFT(o_cf, str_params)

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
            d_lambda       = 2*pi / (obj.m_str_params.eta*i_N);
            d_lbound       = log(d_stock) - 0.5*d_lambda*i_N;
            vi_indices     = 0:(i_N-1);
            vd_log_strikes = d_lbound + d_lambda * vi_indices;
            out_strikes    = exp(vd_log_strikes);

            % evaluate dampened Fourier transform on grid
            % multiply first and last value by 0.5 due to trapezoidal approximation
            d_alpha   = obj.m_str_params.alpha;
            vd_grid   = obj.m_str_params.eta * vi_indices;
            vd_phi    = obj.m_o_cf.eval(vd_grid+1i*d_alpha);
            vd_series = exp(1i*vd_grid*(log(d_stock) - d_lbound));
            vd_series = vd_series ./ (1+d_alpha-1i*vd_grid) ./ ...
                        (d_alpha-1i*vd_grid);
            vd_series = vd_series .* vd_phi;
            vd_series(1)   = 0.5 * vd_series(1);
            vd_series(i_N) = 0.5 * vd_series(i_N);

            % perform inverse FFT
            vd_transform = real(FFT(vd_series, true));

            % Calculate option prices. Multiply by N due to the inverse FFT.
            vd_values = obj.m_str_params.eta / pi * i_N ...
                        * exp((1+d_alpha)*vd_log_strikes-d_alpha ...
                        * log(d_stock)-d_r*d_T) .* vd_transform;

            % reset parameter stock
            obj.m_o_cf.setStock(d_stock);
                    
            % Interpolate option prices of interest.
            out_values = real(extractInterpolatedElements(vd_strikes, out_strikes, vd_values, 'loglinear'));
            %out_values = interp1(out_strikes, vd_values, vd_strikes, 'cubic');
            
        end
        
    end
    
end


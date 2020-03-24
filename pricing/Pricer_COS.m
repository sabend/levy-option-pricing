classdef Pricer_COS < CF_Pricer
    %Pricer_COS: Implements the COS pricing algorithm for European calls.
    %   Reference: F. Fang and C. W. Oosterlee.
    %   A novel pricing method for european options based on
    %   fourier-cosine series expansions.
    %   SIAM Journal on Scientific Computing, 31(2):826--848, 2008.
    
    properties
    
        % inherited: m_str_params (struct: N, L)
        % inherited: m_o_cf
        
        m_vd_lbound
        m_vd_ubound
        
        m_b_user_bounds
        
    end
    
    methods (Access = public)
    
        % constructor
        function obj = Pricer_COS(o_cf, str_params)

            obj@CF_Pricer(o_cf, str_params);
            
            obj.m_b_user_bounds = false;
            if isfield(str_params, 'user_bounds')
               obj.m_b_user_bounds = str_params.user_bounds;
               obj.m_vd_lbound      = str_params.lbound;
               obj.m_vd_ubound      = str_params.ubound;
            end

        end
        
        % run pricing calculation
        function out_values = run(obj, vd_strikes)
            
            % read relevant parameters from characteristic function
            d_r     = obj.m_o_cf.r();
            d_stock = obj.m_o_cf.stock();
            d_T     = obj.m_o_cf.T();
            
            % pre-processing: calculate truncation bounds
            if ~obj.m_b_user_bounds
                obj.m_o_cf.calcCumulants();
                vd_cumulants = obj.m_o_cf.cumulants();
                d_c1  = vd_cumulants(1);
                d_c2  = vd_cumulants(2);
                d_c4  = vd_cumulants(4);
                d_c6  = vd_cumulants(6);
                d_c8  = vd_cumulants(8);
                d_c10 = vd_cumulants(10);
                d_c12 = vd_cumulants(12);
                d_c14 = vd_cumulants(14);
                d_c16 = vd_cumulants(16);
                
                vd_adj = log(d_stock ./ vd_strikes);
                
                d_tmp = obj.m_str_params.L*sqrt(d_c2+sqrt(d_c4+sqrt(d_c6+...
                        sqrt(d_c8+sqrt(d_c10+sqrt(d_c12+sqrt(d_c14+sqrt(d_c16))))))));
                obj.m_vd_lbound = d_c1 + vd_adj - d_tmp;
                obj.m_vd_ubound = d_c1 + vd_adj + d_tmp;
            end
            
            % create grid for evaluation of characteristic function
            vi_k    = 0:(obj.m_str_params.N-1);
            vd_diff = obj.m_vd_ubound - obj.m_vd_lbound;
            m_grid  = (pi ./ transpose(vd_diff)) * vi_k;
            
            % set parameter stock to 1
            obj.m_o_cf.setStock(1);

            % evaluate characteristic function on grid
            m_cf_values = obj.m_o_cf.eval(m_grid);

            % prepare special cosine series coefficients
            vd_arg1 = vi_k*pi;
            m_arg2  = -pi .* transpose(obj.m_vd_lbound./vd_diff) * vi_k;
            m_tmp   = transpose(1./vd_diff) * (vi_k.*pi);
            m_chi   = 1 ./ (1 + m_tmp.^2) ...
                      .* (transpose(exp(obj.m_vd_ubound))*cos(vd_arg1) ...
                      - cos(m_arg2) + m_tmp.*(transpose(exp(obj.m_vd_ubound))*sin(vd_arg1) ...
                      - sin(m_arg2)));
            m_psi   =  (repmat(sin(vd_arg1), numel(vd_strikes), 1) - sin(m_arg2)) ./ m_tmp;
            m_psi(:, 1) = obj.m_vd_ubound;
            m_coeff = m_chi - m_psi;
            m_coeff(:, 1) = 0.5 * m_coeff(:, 1);

            % create summands
            vd_x    = log(d_stock./vd_strikes);
            m_sums  = real(m_cf_values .* exp(pi.*1i.*transpose(((vd_x-obj.m_vd_lbound)./vd_diff)) * vi_k));
            m_sums  = m_sums .* m_coeff;
            vd_sums = sum(m_sums, 2);

            % reset parameter stock
            obj.m_o_cf.setStock(d_stock);

            % calculate and return option values
            out_values = exp(-d_r*d_T).*2./vd_diff.*vd_strikes.*transpose(vd_sums);
                    
        end
        
    end
    
end


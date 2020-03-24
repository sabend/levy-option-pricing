classdef CF_CGMY < CharacteristicFunction
    %CF_BlackScholes: Characteristic function for CGMY model.
    %   Reference: Peter Carr, Helyette Geman, Dilip B. Madan, and Marc Yor.
    %              The fine structure of asset returns: 
    %              An empirical investigation.
    %              Journal of Business
    %              75(2):305-332, 2002
    
    properties
        
        % inherited: m_d_r
        % inherited: m_d_stock
        % inherited: m_d_T
        % inherited: m_str_params (struct: C, G, M, Y)
        % inherited: m_vd_cumulants
        
        m_d_epsilon = 10^(-4);
        
    end
    
    methods (Access = public)
        
        % constructor
        function obj = CF_CGMY(d_r, d_stock, d_T, str_params)
        
            obj = obj@CharacteristicFunction(d_r, d_stock, d_T, str_params);
            
        end
        
        % evaluate characteristic function
        function out_value = eval(obj, vd_u)
            
            d_r     = obj.m_d_r;
            d_stock = obj.m_d_stock;
            d_T     = obj.m_d_T;
            d_C     = obj.m_str_params.C;
            d_G     = obj.m_str_params.G;
            d_M     = obj.m_str_params.M;
            d_Y     = obj.m_str_params.Y;
            
            d_omega = -1/d_T*log(obj.CGMY((-1i), d_T, d_C, d_G, d_M, d_Y));

            out_value = exp(1i*vd_u*(log(d_stock) + (d_r+d_omega)*d_T)) ...
                        .* obj.CGMY(vd_u, d_T, d_C, d_G, d_M, d_Y);
           
        end
        
        % density function
        % this is just the density of the Levy measure
        function out_values = density(obj, vd_x)
            
            d_C     = obj.m_str_params.C;
            d_G     = obj.m_str_params.G;
            d_M     = obj.m_str_params.M;
            d_Y     = obj.m_str_params.Y;
            
            vi_negative = find(vd_x<0);
            vd_negative = vd_x(vi_negative);
            vi_positive = find(vd_x>0);
            vd_positive = vd_x(vi_positive);
            
            out_values = NaN(size(vd_x, 1), size(vd_x, 2));
            
            out_values(vi_negative) = d_C * exp(-d_G*abs(vd_negative)) ./ ...
                                      abs(vd_negative).^(1+d_Y);
            out_values(vi_positive) = d_C * exp(-d_M*vd_positive) ./ ...
                                      vd_positive.^(1+d_Y);
            
        end
        
        % calculate cumulants
        function obj = calcCumulants(obj)
         
            d_C     = obj.params().C;
            d_G     = obj.params().G;
            d_M     = obj.params().M;
            d_Y     = obj.params().Y;
            d_r     = obj.m_d_r;
            d_T     = obj.m_d_T;
            
            % the below values are absolute values of the actual cumulants
            d_c1  = d_T*(d_r-d_C*gamma(-d_Y)*((d_M-1)^d_Y ...
                    - d_M^d_Y+(d_G+1)^d_Y-d_G^d_Y)) ...
                    + d_C*d_T*gamma(1-d_Y)*(d_M^(d_Y-1) ...
                    - d_G^(d_Y-1));
            d_c2  = d_C*d_T*gamma(2-d_Y)*(d_M^(d_Y-2)+d_G^(d_Y-2));
            d_c3  = -d_C*d_T*gamma(3-d_Y)*(d_M^(d_Y-3)-d_G^(d_Y-3));
            d_c4  = d_C*d_T*gamma(4-d_Y)*(d_M^(d_Y-4)+d_G^(d_Y-4));
            d_c5  = -d_C*d_T*gamma(5-d_Y)*(d_M^(d_Y-5)-d_G^(d_Y-5));
            d_c6  = d_C*d_T*gamma(6-d_Y)*(d_M^(d_Y-6)+d_G^(d_Y-6));
            d_c7  = -d_C*d_T*gamma(7-d_Y)*(d_M^(d_Y-7)-d_G^(d_Y-7));
            d_c8  = d_C*d_T*gamma(8-d_Y)*(d_M^(d_Y-8)+d_G^(d_Y-8));
            d_c9  = -d_C*d_T*gamma(9-d_Y)*(d_M^(d_Y-9)-d_G^(d_Y-9));
            d_c10 = d_C*d_T*gamma(10-d_Y)*(d_M^(d_Y-10)+d_G^(d_Y-10));
            d_c11 = -d_C*d_T*gamma(11-d_Y)*(d_M^(d_Y-11)-d_G^(d_Y-11));
            d_c12 = d_C*d_T*gamma(12-d_Y)*(d_M^(d_Y-12)+d_G^(d_Y-12));
            d_c13 = -d_C*d_T*gamma(13-d_Y)*(d_M^(d_Y-13)-d_G^(d_Y-13));
            d_c14 = d_C*d_T*gamma(14-d_Y)*(d_M^(d_Y-14)+d_G^(d_Y-14));
            d_c15 = -d_C*d_T*gamma(15-d_Y)*(d_M^(d_Y-15)-d_G^(d_Y-15));
            d_c16 = d_C*d_T*gamma(16-d_Y)*(d_M^(d_Y-16)+d_G^(d_Y-16));
            
            obj.m_vd_cumulants = [d_c1, d_c2, d_c3, d_c4, d_c5, ...
                                  d_c6, d_c7, d_c8, d_c9, d_c10, ...
                                  d_c11, d_c12, d_c13, d_c14, d_c15, ...
                                  d_c16];
            
        end
        
        % return value and Jacobi matrix of parameter constraint function
        % 'g' represents the inequality constraints, 'h' the equality
        % constraints
        function [out_g, out_gF, out_h, out_hF] = constraints(obj, vd_x)
            
            d_epsilon = obj.m_d_epsilon;
            
            m_tmp = [eye(4); -eye(4)];
            
            vd_bound = [3; 20; 20; (2 - d_epsilon); ...
                        -d_epsilon; 0; 0; -d_epsilon];
            
            out_g  = m_tmp * vd_x - vd_bound;
            out_gF = m_tmp;
            
            out_h  = [];
            out_hF = [];
            
        end
        
        % fill parameter struct from vector
        function out_struct = fillParameters(obj, vd_params)
           
            out_struct = struct('C', vd_params(1), ...
                                'G', vd_params(2), ...
                                'M', vd_params(3), ...
                                'Y', vd_params(4));
            
        end
        
        % lower parameter bound
        function out_bound = lowerBound(obj)
            
            out_bound = [obj.m_d_epsilon; 0; 0; 0];
            
        end
        
        % upper parameter bound
        function out_bound = upperBound(obj)
           
            out_bound = [20; 100; 100; (2-obj.m_d_epsilon)];
            
        end
        
        % evaluate constraint function signaling whether two instances 
        % of models with this characteristic function are ordered with
        % respect to integral stochastic ordering
        % important: first half of 'vd_x' represents ask price parameters
        %            while second half represents bid price parameters
        function [out_g, out_gF, out_h, out_hF] = orderingConstraints(obj, vd_x)

            if (numel(vd_x) ~= 8)
                error('[CF] needs 8 parameters to represent two CGMY models!')
            end
         
            if ((abs(vd_x(1)-vd_x(5))<obj.m_d_epsilon) && ...
                 abs(vd_x(4)-vd_x(8))<obj.m_d_epsilon)
                
                out_gF = [0 -1 0 0 0 1 0 0; ...
                          0 0 1 0 0 0 -1 0; ...
                          0 0 1 0 0 -1 0 0; ...
                          0 1 0 0 0 0 -1 0];
                     
                out_hF = [1 0 0 0 -1 0 0 0; ...
                          0 0 0 1 0 0 0 -1];
                     
                out_g = out_gF * vd_x;
                out_h = out_hF * vd_x;

            else
            
                display(['aC=' num2str(vd_x(1))])
                display(['bC=' num2str(vd_x(5))])
                display(['aY=' num2str(vd_x(4))])
                display(['bY=' num2str(vd_x(8))])
                error('[CF] C and Y have each to be identical for both CGMY processes!')
                
            end
            
        end
            
    end
    
    methods (Access = private)

        % basic CGMY characteristic function
        function out_value = CGMY(obj, vd_u, d_T, d_C, d_G, d_M, d_Y)

            out_value = exp(d_T*d_C*gamma(-d_Y)*((d_M-1i*vd_u).^d_Y ...
                        - d_M^d_Y + (d_G+1i*vd_u).^d_Y - d_G^d_Y));

        end

    end
    
end



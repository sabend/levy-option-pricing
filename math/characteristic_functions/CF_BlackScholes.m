classdef CF_BlackScholes < CharacteristicFunction
    %CF_BlackScholes: Characteristic function for Black-Scholes model.
    %   The Black-Scholes model uses the characteristic function of
    %   a normally distributed random variable.
    
    properties
        
        % inherited: m_d_r
        % inherited: m_d_stock
        % inherited: m_d_T
        % inherited: m_str_params (struct: sigma)
        % inherited: m_vd_cumulants
        
    end
    
    methods (Access = public)
        
        % constructor
        function obj = CF_BlackScholes(d_r, d_stock, d_T, str_params)
        
            obj = obj@CharacteristicFunction(d_r, d_stock, d_T, str_params);
            
        end
        
        % evaluate characteristic function
        function out_value = eval(obj, vd_u)
            
            d_r     = obj.m_d_r;
            d_stock = obj.m_d_stock;
            d_T     = obj.m_d_T;
            d_sigma = obj.m_str_params.sigma;
           
            out_value = exp(1i * vd_u * ((d_r - 0.5*d_sigma^2) * d_T ...
                        + log(d_stock)) - 0.5*d_sigma^2 * vd_u.^2 * d_T);
            
        end
        
        % calculate cumulants
        function obj = calcCumulants(obj)
           
            d_sigma = obj.params().sigma;
            d_r     = obj.m_d_r;
            d_T     = obj.m_d_T;
            
            d_c1 = (d_r - 0.5*d_sigma^2) * d_T;
            d_c2 = d_sigma^2 * d_T;
            
            obj.m_vd_cumulants = [d_c1, d_c2, zeros(1, 14)];
            
        end
        
        % return value and Jacobi matrix of parameter constraint function
        % 'g' represents the inequality constraints, 'h' the equality
        % constraints
        function [out_g, out_gF, out_h, out_hF] = constraints(obj, vd_x)
            
            out_g  = -vd_x;
            out_gF = -1;
            
            out_h  = [];
            out_hF = [];
            
        end
        
        % fill parameter struct from vector
        function out_struct = fillParameters(obj, vd_params)
           
            out_struct = struct('sigma', vd_params);
            
        end
        
        % lower parameter bound
        function out_bound = lowerBound(obj)
            
            out_bound = eps;
            
        end
        
        % upper parameter bound
        function out_bound = upperBound(obj)
            
            out_bound = 1;
            
        end
        
        % evaluate constraint function signaling whether two instances 
        % of models with this characteristic function are ordered with
        % respect to integral stochastic ordering
        function [out_g, out_gF, out_h, out_hF] = orderingConstraints(obj, vd_x)
        
            warning('[CF] stochastic ordering of Black-Scholes prices is not reasonable!')
            
            out_g  = [];
            out_gF = [];
            
            out_h  = [];
            out_hF = [];
            
        end
            
    end
    
end


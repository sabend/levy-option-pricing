classdef (Abstract) CharacteristicFunction < handle
    %CharacteristicFunction: Abstract characteristic function.
    %   Defines interface for functionality of concrete
    %   characteristic functions.
    
    properties (Access = protected)
    
        m_d_r
        m_d_stock
        m_d_T
        
        m_str_params
        
        m_vd_cumulants
        
    end
    
    methods (Access = protected)
        
        % constructor with parameters
        function obj = CharacteristicFunction(d_r, d_stock, d_T, str_params)

            obj.m_d_r        = d_r;
            obj.m_d_stock    = d_stock;
            obj.m_d_T        = d_T;
            obj.m_str_params = str_params;

        end
        
    end
    
    methods (Access = public, Abstract = false)
        
        % get risk-free rate
        function r = r(obj)
            
            r = obj.m_d_r;
            
        end
        
        % get stock price
        function s = stock(obj)
           
            s = obj.m_d_stock;
            
        end
        
        % get maturity
        function m = T(obj)
           
            m = obj.m_d_T;
            
        end
        
        % get cumulants
        function c = cumulants(obj)
       
            c = obj.m_vd_cumulants;
            
        end
        
        % get parameters struct
        function p = params(obj)
           
            p = obj.m_str_params;
            
        end
        
        % set risk-free rate
        function obj = setRate(obj, d_r)
            
            obj.m_d_r = d_r;
            
        end
        
        % set stock price
        function obj = setStock(obj, d_stock)
           
            obj.m_d_stock = d_stock;
            
        end
        
        % set maturity
        function obj = setMaturity(obj, d_T)
           
            obj.m_d_T = d_T;
            
        end
        
        % reset parameters
        function obj = reset(obj, str_params)
          
            obj.m_str_params = str_params;
            
        end
        
        % density function
        density(obj, vd_x)
         
    end
    
    methods (Access = public, Abstract = true)
        
        % calculate cumulants
        calcCumulants(obj)
        
        % return value and Jacobi matrix of parameter constraint function
        % 'g' represents the inequality constraints, 'h' the equality
        % constraints
        constraints(obj, vd_x)
        
        % evaluate characteristic function
        eval(obj, vd_u)
        
        % fill parameter struct from vector
        fillParameters(obj, vd_params)
        
        % lower parameter bound
        lowerBound(obj)
        
        % upper parameter bound
        upperBound(obj)
        
        % evaluate constraint function signaling whether two instances 
        % of models with this characteristic function are ordered with
        % respect to integral stochastic ordering
        orderingConstraints(obj, vd_x)
        
    end
    
end


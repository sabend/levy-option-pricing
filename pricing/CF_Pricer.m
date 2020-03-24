classdef (Abstract) CF_Pricer < handle
    %CF_Pricer: Prices call options with characteristic function of
    %           underlying model.
    
    properties (Access = protected)
        
        m_str_params
        m_o_cf
        
    end
    
    methods (Access = protected)
    
        % constructor
        function obj = CF_Pricer(o_cf, str_params)

            obj.m_str_params = str_params;
            obj.m_o_cf       = o_cf;

        end
        
    end
    
    methods (Access = public, Abstract = true)
        
        % run pricing calculation
        run(obj, vd_strikes)
        
    end
    
    methods (Access = public, Abstract = false)
        
        % get characteristic function member
        function out_cf = cf(obj)
           
            out_cf = obj.m_o_cf;
            
        end
        
        % set parameters of characteristic function
        function obj = setParamsCF(obj, str_params)
           
            obj.m_o_cf.reset(str_params);
            
        end
        
        % set maturity member of characteristic function
        function obj = setMaturityCF(obj, d_T)
           
            obj.m_o_cf.setMaturity(d_T);
            
        end
        
    end
    
end


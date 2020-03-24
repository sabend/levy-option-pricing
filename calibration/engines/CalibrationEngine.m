classdef (Abstract) CalibrationEngine < handle
    %CalibrationEngine: Calibration engine prototype.
    %   Derived classes implement the actual calibration logic.
    %   Every calibration instance has to be instantiated with a pricer
    %   and a loss function.
    
    properties (Access = protected)
        
        m_b_logging
        m_c_option_data % {i, 1} = [maturity; strike; price]
        m_o_pricer
        m_o_loss_function
        m_str_params
        m_str_result
        
    end
    
    methods (Access = protected)
        
        % constructor
        function obj = CalibrationEngine(o_pricer, o_loss_function, ...
                                         str_params, c_option_data, ...
                                         b_logging)
            
            obj.m_b_logging       = b_logging;
            obj.m_o_pricer        = o_pricer;
            obj.m_o_loss_function = o_loss_function;
            obj.m_str_params      = str_params;
            obj.m_c_option_data   = c_option_data;
            
        end
        
    end
    
    methods (Access = public)
        
        % run calibration
        run(obj)
        
        % return results
        function r = result(obj)
           
            r = obj.m_str_result;
            
        end
        
    end
    
end


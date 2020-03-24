classdef (Abstract) LossFunction < handle
    %LossFunction: Assign loss to pair of input and target value.
    %   Derived classes implement the actual measurement of loss.
    
    properties (Access = protected)
        
        m_c_target
        m_i_col
        
    end
    
    methods (Access = protected)
        
        % constructor
        function obj = LossFunction(c_target, i_col)
           
            obj.m_c_target = c_target;
            obj.m_i_col    = i_col;
            
        end
        
    end
    
    methods (Access = public)
        
        % calculate loss
        calculate(obj, c_input)
        
    end
    
end


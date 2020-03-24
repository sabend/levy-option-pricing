classdef Loss_Quadratic < LossFunction
    %Loss_Quadratic: Squared-loss function.
    %   Sum of squared differences between input and target value for each
    %   element of structs.
    
    properties
        
        % inherited: m_c_target
        % inherited: m_i_col
        
    end
    
    methods (Access = public)
        
        % constructor
        function obj = Loss_Quadratic(c_target, i_col)
            
            obj@LossFunction(c_target, i_col);
            
        end
        
        % calculate loss
        function out_loss = calculate(obj, c_input)
           
            i_length = numel(obj.m_c_target);
            
            if (numel(c_input) ~= i_length)
               
                error('[LOSS] cell array sizes do not match!')
                
            else
                
                d_loss = 0;
                
                for i = 1 : i_length
                    
                    d_diff = obj.m_c_target{i, 1}(:, obj.m_i_col)  ...
                             - c_input{i, 1}(:, obj.m_i_col);
                    d_loss = d_loss + sum(d_diff.^2);
                    
                end
                
                out_loss =  d_loss;
                
            end

        end
        
    end
    
end


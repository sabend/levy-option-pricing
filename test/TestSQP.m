%
% Script to test external SQP implementation.
%

% % feasible initial guess
% % x0 = [0; 4];
% x0 = 1.7;
% 
% % objective function
% objF = @objectiveFunction;
% 
% % constraints function
% constF = @constraintsFunction;
% 
% % lower bounds
% % lb = [0; 0];
% lb = [];
% 
% % start optimization
% [xMin,fMin,perf] = modSQP(objF, x0, constF, lb);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test object oriented modification of external SQP algo

o_optimisable = TestTarget;

[d_x, d_y, str_info] = modSQP(o_optimisable, 2, [], [], eye(1), 100, 0.0001);
















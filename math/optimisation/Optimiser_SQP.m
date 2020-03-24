classdef Optimiser_SQP < handle
    %Optimiser_SQP: Implementation of sequential quadratic programming.
    %   The constructor takes as argument an object implementing the
    %   following methods:
    %   1. objective(vd_x)
    %   2. constraints(vd_x)
    %   Note that the optimisation target has to be derived from 'handle'.
    %   If there exist no equality constraints the function 'constraints'
    %   is required to return [] for 'vd_h' and 'm_hF'.
    
    properties
        
        m_b_verbose
        
        m_d_alpha
        m_d_beta
        m_d_gamma
        m_d_epsilon
        m_d_objective
        m_d_objective_prev
        
        m_i_iteration
        m_i_max_iter
        
        m_m_H

        m_t_log
        
        m_o_optimisable
        
        m_str_params
        m_str_results
        
        m_vd_x
        m_vd_x_prev
        m_vd_lambda
        m_vd_lambda_prev
        m_vd_mu
        m_vd_mu_prev
        m_vd_gradient
        
    end
    
    methods (Access = public)
        
        % constructor
        function obj = Optimiser_SQP(o_optimisable, str_params, b_verbose)
        
            obj.m_o_optimisable = o_optimisable;
            
            obj.m_b_verbose  = b_verbose;
            
            obj.m_d_alpha    = str_params.alpha;
            obj.m_d_epsilon  = str_params.epsilon;
            obj.m_i_max_iter = str_params.max_iter;
            obj.m_d_beta     = str_params.beta;
            obj.m_d_gamma    = str_params.gamma;

            obj.m_t_log      = table;
                                                       
            obj.m_str_params = str_params;
            
        end
        
        % run optimisation routine
        function run(obj)
            
            % initialisation
            obj.m_i_iteration      = 0;
            obj.m_vd_x             = obj.m_str_params.initial_x;
            obj.m_vd_x_prev        = obj.m_str_params.initial_x;
            obj.m_vd_lambda        = obj.m_str_params.initial_lambda;
            obj.m_vd_lambda_prev   = obj.m_str_params.initial_lambda;
            obj.m_vd_mu            = obj.m_str_params.initial_mu;
            obj.m_vd_mu_prev       = obj.m_str_params.initial_mu;
            obj.m_m_H              = eye(numel(obj.m_vd_x));
            obj.m_d_objective      = obj.m_o_optimisable.objective(obj.m_vd_x, true);
            obj.m_d_objective_prev = obj.m_d_objective;
            obj.approxGradient();
        
            while (obj.m_i_iteration < obj.m_i_max_iter)
              
                % check whether current iterates are a KKT pair
                if obj.isKKT()
                    warning('SQP algorithm terminated: reached KKT point!')
                    break;
                end

                % stop if deviation of current iterate from previous one
                % does not differ by more than chosen precision
                vd_diff = norm(obj.m_vd_x - obj.m_vd_x_prev);
                if (obj.m_i_iteration > 1 && vd_diff < obj.m_d_epsilon)
                    warning('SQP algorithm terminated: change in x small enough!')
                    break;
                end

                % stop if deviation of current error from previous one
                % does not differ by more than chosen precision
                vd_diff = abs(obj.m_d_objective - obj.m_d_objective_prev);
                if (obj.m_i_iteration > 1 && vd_diff < obj.m_d_epsilon)
                    warning('SQP algorithm terminated: change in objective small enough!')
                    break;
                end
                
                % solve the inner quadratic problem
                [vd_s, vd_lambda_qp, vd_mu_qp, b_terminated] = obj.solveQP();
                if b_terminated
                    warning('[CAL] terminated SQP algorithm due to attained iteration limit!')
                    break
                end
                
                % determine Armijo step-size
                d_penalised_x = obj.penalised(obj.m_vd_x);
                d_sigma = 1;

                while true

                    if ((obj.penalised(obj.m_vd_x+d_sigma*vd_s) - d_penalised_x) ...
                        <= -obj.m_d_gamma*d_sigma*vd_s'*obj.m_m_H*vd_s) % 0) %
                        break
                    else
                        d_sigma = d_sigma * obj.m_d_beta;
                    end

                end

                % update the iterates
                obj.m_vd_x_prev      = obj.m_vd_x;
                obj.m_vd_x           = obj.m_vd_x + d_sigma * vd_s;
                obj.m_vd_lambda_prev = obj.m_vd_lambda;
                obj.m_vd_lambda      = vd_lambda_qp;
                obj.m_vd_mu_prev     = obj.m_vd_mu;
                obj.m_vd_mu          = vd_mu_qp;
                
                % approximate gradient
                obj.approxGradient();
                
                % calculate the Hessian
                obj.updateHessian();
                
                % update objective
                obj.m_d_objective_prev = obj.m_d_objective;
                obj.m_d_objective = obj.m_o_optimisable.objective(obj.m_vd_x, true);
                
                % update iteration
                obj.m_i_iteration = obj.m_i_iteration + 1;
                
                % possible write results of current iteration
                if obj.m_b_verbose
                    fprintf('iteration= %5.0f,  ', obj.m_i_iteration)
                    fprintf('x=')
                    fprintf('%15.5f,  ', obj.m_vd_x)
                    fprintf('s=')
                    fprintf('%15.5f,  ', vd_s)
                    fprintf('sigma=')
                    fprintf('%15.5f,  ', d_sigma)
                    fprintf('lambda=')
                    fprintf('%15.5f,  ', obj.m_vd_lambda)
                    fprintf('mu=')
                    fprintf('%15.5f,  ', obj.m_vd_mu)
                    fprintf('objective=%e\n', obj.m_d_objective)
                end

                % log results
                t_tmp = table;
                t_tmp.iteration = obj.m_i_iteration;
                t_tmp.x         = obj.m_vd_x';
                t_tmp.s         = vd_s';
                t_tmp.sigma     = d_sigma;
                t_tmp.objective = obj.m_d_objective;
                obj.m_t_log     = [obj.m_t_log; t_tmp];

            end
            
            % fill result struct
            obj.m_str_results.x          = obj.m_vd_x;
            obj.m_str_results.objective  = obj.m_d_objective;
            obj.m_str_results.iterations = obj.m_i_iteration;
            
        end
        
        % return results
        function out_results = results(obj)
            
            out_results = obj.m_str_results;
            
        end

        % get result of iterations
        function out_table = log(obj)

            out_table = obj.m_t_log;

        end
        
    end
    
    methods (Access = private)
        
        % approximate gradient
        function obj = approxGradient(obj)
        
            i_N = numel(obj.m_vd_x);
            vd_grad = zeros(i_N, 1);
            
            d_h = sqrt(eps) * max(norm(obj.m_vd_x, inf), 1);
            m_h = d_h * eye(i_N);
       
            vd_objective = obj.m_o_optimisable.objective(obj.m_vd_x, false);
            for i = 1 : i_N
               
                vd_grad(i) = (obj.m_o_optimisable.objective(obj.m_vd_x + m_h(:, i), false) - ...
                              vd_objective) / d_h;
    
            end
            
            obj.m_vd_gradient = vd_grad;         

        end
        
        % evaluate the Lagrange function on either the current or previous
        % iterates
        function out_lagrange = nablaLagrange(obj, s_mode)
            
            switch s_mode
                case 'current'
                    vd_x      = obj.m_vd_x;
                    vd_lambda = obj.m_vd_lambda;
                    vd_mu     = obj.m_vd_mu;
                   
                case 'previous'
                    vd_x      = obj.m_vd_x_prev;
                    vd_lambda = obj.m_vd_lambda_prev;
                    vd_mu     = obj.m_vd_mu_prev;
                    
                otherwise
                    error(['[OPT] mode ' s_mode ' not implemented for ' ...
                           'evaluation of Langrange function!'])
                
            end

            % evaluate constraint functions
            [~, m_gF, ~, m_hF] = obj.m_o_optimisable.constraints(vd_x);

            % evaluate nabla Lagrange function
            out_lagrange = obj.m_vd_gradient + m_gF'*vd_lambda;
            if ~isempty(m_hF)
                out_lagrange = out_lagrange + m_hF'*vd_mu;
            end
                
        end
        
        % checks whether current iterate is a KKT point
        function out_is_kkt = isKKT(obj)
        
            out_is_kkt = false;

            % check KKT conditions
            if (sum(abs(obj.nablaLagrange('current')) > ...
                   (obj.m_d_epsilon*max(norm(obj.m_vd_x, inf), 1)))==0)
             
                [vd_g, ~, vd_h, ~] = obj.m_o_optimisable.constraints(obj.m_vd_x);

                if (sum(abs(vd_h)) < obj.m_d_epsilon*numel(vd_h))

                    if (sum(vd_g > 0) == 0)

                        if (sum(obj.m_vd_lambda < 0) == 0)

                            if (abs(obj.m_vd_lambda'*vd_g) < obj.m_d_epsilon)

                                out_is_kkt = true;
                                
                            end
                            
                        end
                        
                    end
                    
                end
        
            end
            
        end
        
        % performs Powell's modified BFGS update on the hessian
        function updateHessian(obj)
            
            return
            
            vd_s = obj.m_vd_x - obj.m_vd_x_prev;
            vd_y = obj.nablaLagrange('current') - obj.nablaLagrange('previous');
            
            d_sy  = vd_s' * vd_y;
            d_sHs = vd_s' * obj.m_m_H * vd_s;
            
            if (d_sy < 0.2*d_sHs)
                d_theta = 0.8*d_sHs / (d_sHs-d_sy);
            else
                d_theta = 1;
            end
        
            vd_y_mod = d_theta*vd_y + (1-d_theta)*obj.m_m_H*vd_s;
            
            obj.m_m_H = obj.m_m_H + vd_y_mod*vd_y_mod' / (vd_s'*vd_y_mod) ...
                        - obj.m_m_H*(vd_s*vd_s')*obj.m_m_H/d_sHs;
       
        end
        
        % penalisation function ensuring global convergence of SQP
        function out_val = penalised(obj, vd_x)
            
            [vd_g, ~, vd_h, ~] = obj.m_o_optimisable.constraints(vd_x);
            
            out_val = obj.m_o_optimisable.objective(vd_x, false) + obj.m_d_alpha * ...
                      (sum(max(vd_g, 0)) + sum(abs(vd_h)));
            
        end
        
        % determine solution triplet for KKT problem
        function [out_x2, out_lambda2, out_mu, out_terminated] = calcKKT(obj, vi_active)
         
            vi_active = sort(vi_active);
            out_terminated = false;
           
            % prepare equation system
            [~, m_gF, ~, m_hF] = obj.m_o_optimisable.constraints(obj.m_vd_x);
            
            m_gA = m_gF(vi_active, :);
            i_g1 = 0;
            i_h1 = 0;
            if ~isempty(m_gA)
                i_g1 = size(m_gA, 1);
            else
                m_gA = [];
            end
            if ~isempty(m_hF)
                i_h1 = size(m_hF, 1);
            end
            m_zeros1 = [];
            m_zeros2 = [];
            if ~isempty(m_gA)
                m_zeros1 = zeros(i_g1, i_g1+i_h1);
            end
            if ~isempty(m_hF)
                m_zeros2 = zeros(i_h1, i_g1+i_h1);
            end
            m_left = [obj.m_m_H m_gA' m_hF'; ...
                      m_gA m_zeros1;         ...
                      m_hF m_zeros2];
            vd_right = [-obj.m_vd_gradient; zeros((i_g1+i_h1), 1)];

            % solve equation system
            if (sum(eig(m_left)<=0) == 0)               
                m_R    = chol(m_left);
                vd_sol = m_R \ (m_R \ vd_right);
            else
                vd_sol = m_left \ vd_right;
            end

            % ToDo: fix
            if ((sum(isnan(vd_sol))+sum(isinf(vd_sol)))>0)
                %out_terminated = true;
                vd_sol([find(isnan(vd_sol)); find(isinf(vd_sol))]) = 0;
            end

            % return solution
            i_M1 = size(obj.m_m_H, 1);
            i_M2 = numel(vi_active);
            out_x2      = vd_sol(1:i_M1);
            out_lambda2 = [];
            out_mu      = zeros(numel(obj.m_vd_mu), 1);
            if ~isempty(m_gA)
                out_lambda2 = vd_sol((i_M1+1):(i_M1+i_M2));
            end
            if ~isempty(m_hF)
                out_mu = vd_sol((i_M1+i_M2+1):end);
            end

        end
        
        % solver for the inner quadratic problem of SQP
        function [out_s, out_lambda_qp, out_mu_qp, out_terminated] = solveQP(obj)
         
            [vd_g, m_gF, vd_h, m_hF] = obj.m_o_optimisable.constraints(obj.m_vd_x);
            vd_nf = obj.m_vd_gradient;
            
            % find initial feasible starting point
            vd_s      = zeros(numel(obj.m_vd_x), 1);
            vd_lambda = zeros(numel(obj.m_vd_lambda), 1);
            vd_mu     = zeros(numel(obj.m_vd_mu), 1);
            
            % initialise the set of active constraints
            vi_active = find(abs(m_gF*vd_s+vd_g) < obj.m_d_epsilon);

            % number of inequality constraints
            vi_all = 1:numel(vd_g);

            i_count = 0;
            out_terminated = false;
            while true

                i_count = i_count + 1;
                if (i_count == 100*obj.m_i_max_iter) 
                    out_terminated = true;
                    break
                end
                
                % ToDo: fix
                vi_active = unique(vi_active);
                
                % check for KKT
                vd_tmp1 = zeros(numel(obj.m_vd_x), 1);
                vd_tmp2 = zeros(numel(obj.m_vd_x), 1);
                if ~isempty(m_gF)
                    vd_tmp1 = m_gF'*vd_lambda;
                end
                if ~isempty(m_hF)
                    vd_tmp2 = m_hF'*vd_mu;
                end

                if (sum(abs(vd_nf+obj.m_m_H*vd_s+vd_tmp1+vd_tmp2)) ...
                        < obj.m_d_epsilon*numel(vd_nf))
                    if (sum(vd_lambda<0)==0)
                        if ~isempty(m_gF)
                            if ((sum(vd_g+m_gF*vd_s>0)==0) && (abs(vd_lambda'*vd_g)<obj.m_d_epsilon))
                                if ~isempty(m_hF)
                                    if (abs(vd_h+m_hF*vd_s)<obj.m_d_epsilon*numel(vd_h))
                                        break
                                    end
                                else
                                    break
                                end
                            end
                        else
                            break
                        end
                    end
                end

                % identify active and inactive indices
                vi_inactive = setdiff(vi_all, vi_active);
                vd_lambda(vi_inactive) = 0;
                [vd_ds, vd_lambda2, vd_mu, b_terminated] = obj.calcKKT(vi_active);
                if b_terminated
                    out_terminated = true;
                    break
                end
                vd_lambda(vi_active) = vd_lambda2;

                % handle cases as in Geiger
                if (sum(abs(vd_ds)>obj.m_d_epsilon)==0)               
                    if (sum(vd_lambda2<0) == 0)                       
                        break
                    end
                    d_min = min(vd_lambda2);
                    if (d_min < 0)
                        vi_negative = find(vd_lambda==d_min);
                        vi_active   = setdiff(vi_active, vi_negative(1));
                        continue
                    end
                else
                    b_h_cond = true;
                    if (~isempty(vd_h) && ~isempty(m_hF))
                        b_h_cond = (sum(abs(m_hF*(vd_s+vd_ds)+vd_h)>=obj.m_d_epsilon)==0);
                    end
                    if (sum((m_gF*(vd_s+vd_ds)+vd_g)>0)==0) && b_h_cond
                        vd_s = vd_s + vd_ds;
                        continue
                    else
                        vd_cond  = m_gF*vd_ds;
                        vi_which = find(vd_cond>0);
                        if isempty(vi_which)
                            error('[CAL] empty condition set in SQP algorithm!')
                        end
                        vd_tmp    = (-vd_g-m_gF*vd_s) ./ vd_cond;
                        vi_which  = intersect(vi_which, vi_inactive);
                        d_sigma   = min(vd_tmp(vi_which));

                        vd_s      = vd_s + d_sigma*vd_ds;
                        vi_new    = find(vd_tmp==d_sigma);
                        vi_active = [vi_active; vi_new(1)];
                    end
                end

            end

            out_s         = vd_s;
            out_lambda_qp = vd_lambda;
            out_mu_qp     = vd_mu;
   
        end
        
    end
    
end


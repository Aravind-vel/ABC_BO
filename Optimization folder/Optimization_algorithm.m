classdef Optimization_algorithm

    % This class have functions to get next point for optimization from
    % different algorithms

    properties
        
    end

    methods (Static)

        function next_point = bayesopt_algo(problem_details, varargin)

            % varargin options

            % log obj - by default obj values were not logged, set this  value to 'true' if needed
            % acq_func - acquisition function
            % epsilon - exploration factor
            % objective constrain

            % saving input arguments as structure - optimizaiton_details
            for i = 1:2:size(varargin{1,1}, 2)
                key = varargin{1,1}{1,i};
                value = varargin{1,1}{1,i+1};

                optimization_details.(key) = value;
            end

            %% var for bayesopt
            var_labels = [problem_details.cont_var_labels, problem_details.dis_num_var_labels, problem_details.cat_var_labels];

            vars = [];

            if ~isempty(problem_details.n_cont_var) && problem_details.n_cont_var
                for i = 1:problem_details.n_cont_var
                    name = problem_details.cont_var_labels(i);
                    name = char(name);
                    limit = [problem_details.cont_var_bounds(1,i), problem_details.cont_var_bounds(2,i)];
                    vars = [vars, optimizableVariable(name,limit,"Type","real")];
                end
            else
                problem_details.n_cont_var = 0;
            end

            if ~isempty(problem_details.n_dis_num_var) && problem_details.n_dis_num_var
                for i = 1:problem_details.n_dis_num_var
                    name = problem_details.dis_num_var_labels(i);
                    name = char(name);
                    num_points = problem_details.dis_num_var_lists{i,1};
                    levels = numel(num_points);
                    int_points = 1:1:levels; % for mapping
                    dis_num_var_map_to_val.(name) = dictionary(int_points,num_points);
                    dis_num_var_map_to_int.(name) = dictionary(num_points,int_points);
                    limit = [1,levels];
                    vars = [vars, optimizableVariable(name, limit, "Type","integer")];
                end
            else
                problem_details.n_dis_num_var = 0;
            end

            if ~isempty(problem_details.n_cat_var) && problem_details.n_cat_var
                for i = 1:problem_details.n_cat_var
                    name = problem_details.cat_var_labels(i);
                    name = char(name);
                    levels = problem_details.cat_var_lev_labels(i);
                    vars = [vars, optimizableVariable(name,levels{1,1}, "Type","categorical")];
                end
            else
                problem_details.n_cat_var = 0;
            end

            %% Initial variables and objectives

            X_initial_bayes = table();

            % add continuous variables
            if problem_details.n_cont_var >0
                for i = 1:problem_details.n_cont_var
                    field_name = problem_details.cont_var_labels(i);
                    field_name = char(field_name);
                    X_initial_bayes.(field_name) = problem_details.variables(:,i);
                end
            end

            % add discrete numeric variables
            % Map to corresponding integers
            if problem_details.n_dis_num_var > 0
                for i = 1:problem_details.n_dis_num_var
                    field_name = problem_details.dis_num_var_labels(i);
                    field_name = char(field_name);
                    % integer values
                    integer_values = dis_num_var_map_to_int.(field_name)(problem_details.variables(:,problem_details.n_cont_var+i));
                    X_initial_bayes.(field_name) = integer_values;
                end
            end

            % add categorical variables
            if problem_details.n_cat_var >0
                for i = 1:problem_details.n_cat_var
                    field_name = problem_details.cat_var_labels(i);
                    field_name = char(field_name);
                    values = problem_details.variables(:,problem_details.n_cont_var+problem_details.n_dis_num_var+i);
                    values_as_labels = cell(numel(values),1);
                    for j = 1:numel(values)
                        values_as_labels{j} = problem_details.cat_var_num_to_label_map.(field_name)(values(j));
                    end

                    X_initial_bayes.(field_name) = categorical(string(values_as_labels));
                end
            end


            % log of objectives
            if isfield(optimization_details,'log_obj')
                if optimization_details.log_obj
                    for i = 1:size(problem_details.objectives,1)
                        if problem_details.objectives(i) ==0
                            problem_details.objectives(i) = 1e-12;
                        end
                        problem_details.objectives(i) = log(problem_details.objectives(i));
                    end
                end
            end


            Y_initial_bayes = problem_details.objectives;

            if strcmp(problem_details.obj_criteria{1},'max')
                % bayesopt by default minimize objective
                Y_initial_bayes = -1.*Y_initial_bayes;
            elseif strcmp(problem_details.obj_criteria{1},'min')
                Y_initial_bayes = Y_initial_bayes;
            end

            if size(X_initial_bayes,1) ~= size(Y_initial_bayes,1)
                error('variables and objectives size mismatch')
            end

            fprintf('Number of data loaded %d\n', numel(Y_initial_bayes))

            %% optimization phase

            % dummy function to complete bayesopt loop
            obj_fun = @(x) 0;
            

            % budget
            n_point = numel(Y_initial_bayes) + 1;

            % acquisition function default - 'EIP'
            if ~isfield(optimization_details,'acq_func')
                optimization_details.acq_func = 'expected-improvement-plus';
            end

            % default epsilon value
            if ~isfield(optimization_details,'epsilon')
                optimization_details.epsilon = 0.5;
            end


            % Initialize bayesopt
            if ~isfield(optimization_details,'constraints')
                % normal BO
                bayes_obj = bayesopt(obj_fun, vars, 'MaxObjectiveEvaluations', n_point, 'InitialX', X_initial_bayes,...
                    'InitialObjective',Y_initial_bayes,'AcquisitionFunctionName',optimization_details.acq_func,...
                    'ExplorationRatio', optimization_details.epsilon,'Verbose',0,PlotFcn=[]);
            end

            %% objective constraints
            n_obj_constraints = 0;
            constraint_violations = [];
            if isfield(optimization_details,'objective_constraints')
                if isempty(optimization_details.objective_constraints)
                    error('objective constraints structure is empty.  Remove objective_constraints in function call if you dont want to apply constraints on objectives')
                end
                obj_constraints = fieldnames(optimization_details.objective_constraints);
                n_obj_constraints = numel(obj_constraints);
                
                        % calculating constraint violations
                        results_struct = table2struct(problem_details.results);
                        for i = 1:size(X_initial_bayes,1)
                            for j = 1:n_obj_constraints
                                text_name = char(obj_constraints{j});
                                constraint_violations(i,j) = optimization_details.objective_constraints.(text_name)(results_struct(i));
                            end
                        end
                        obj_fun = @(x) dummy_fun(x, n_obj_constraints);
                
            end

            %% variable constraints
            variable_constraint = [];
            % variable details
            variable_details.n_dis_num_var = problem_details.n_dis_num_var;
            variable_details.dis_num_var_labels = problem_details.dis_num_var_labels;
            variable_details.n_cont_var = problem_details.n_cont_var;
            variable_details.cont_var_labels = problem_details.cont_var_labels;
            % if  discrete numeric variable is there send the mapping
            % in variable details
            if problem_details.n_dis_num_var>0
                variable_details.dis_num_var_map_to_val = dis_num_var_map_to_val;
            end
            if isfield(optimization_details,'variable_constraints')
                if isempty(optimization_details.variable_constraints)
                    error('variable constraints structure is empty. Remove variable_constraints in function call if you dont want to apply constraints on variables');
                end


                
                variable_constraints = optimization_details.variable_constraints;
                variable_constraint = @(x) variable_constraint_func(x, variable_constraints, variable_details);
            end

            %% ABC BO
            if isfield(optimization_details,'ABC_BO')
               if isempty(optimization_details.ABC_BO)
                   error('ABC_BO structure is emtpy. Remove ABC_BO in function call if you dont want to use Adaptive Boundary constraint - BO');
               end
               % if isfield(optimization_details,'variable_constraints')
               %     error('provide variable constraints in ABC_BO constraint itself. ABC_BO also involves constraints in variables. HINT: Provide multiple constraints in a function')
               % end

               
                        if ~isfield(optimization_details.ABC_BO, 'constraint') || ...
                                ~isfield(optimization_details.ABC_BO, 'boundary_update') 
                            error(['One or more of the ABC_BO dependents are missing - provide constraints, ' ...
                                ' boundary update equations']);
                        end
                        
                        % add constraints to variable constraints if any

                        variable_constraints.ABC_BO = @(x) optimization_details.ABC_BO.constraint(x);
                        variable_constraint = @(x) variable_constraint_func(x, variable_constraints, variable_details);

                        % update boundary
                        for i = 1:numel(optimization_details.ABC_BO.boundary_update.variables)

                            var_to_update = optimization_details.ABC_BO.boundary_update.variables{i};
                            % index of var
                            index = find(strcmp({vars.Name}, var_to_update));

                            % continuous variables
                            if ismember(var_to_update,problem_details.cont_var_labels)
                                if strcmp(optimization_details.ABC_BO.boundary_update.which_limit{i},'upper')
                                    max_limit = optimization_details.ABC_BO.boundary_update.limit.(var_to_update);
                                    start_limit = vars(index).Range(1);
                                    end_limit = vars(index).Range(end);
                                    if max_limit < end_limit
                                        if max_limit < start_limit
                                            % if max_limit is less than
                                            % lower boundary; fixing the
                                            % variable to lower boundary                                          
                                            % fprintf('Fixing the variable %s to %d\n', var_to_update, start_limit);
                                            % fprintf('Maximum limit of %s is %d which is less than the specified lower boundary %d\n', var_to_update, max_limit, start_limit);
                                            % vars(index) = [];
                                            % % removing the existing results of variables
                                            % X_initial_bayes.(var_to_update) = [];
                                            % % adding the variable to fixed
                                            % % variable list
                                            % fixed_variables.(var_to_update) = start_limit;
                                            text = sprintf('Optimization not possible using ABC-BO; not possible to improve the objective from existing best;max limit %d of variable %s is less than the specified lower boundary %d ', max_limit, var_to_update, start_limit);
                                            error(text);
                                        else
                                            % update upper boundary if max limit is greater than lower boundary
                                            vars(index).Range(end) = max_limit;
                                            fprintf('Adjusted upper boundary of %s from %d to %d\n', var_to_update, end_limit, max_limit);
                                            fprintf('Updated limits for %s [%d, %d]\n', var_to_update, vars(index).Range(1), vars(index).Range(end));
                                        end
                                    else
                                        fprintf('No boundary update for %s; max limit %d is more than specified upper boundary %d\n', var_to_update, max_limit, vars(index).Range(end));
                                    end

                                elseif strcmp(optimization_details.ABC_BO.boundary_update.which_limit{i},'lower')
                                    min_limit = optimization_details.ABC_BO.boundary_update.limit.(var_to_update);
                                    start_limit = vars(index).Range(1);
                                    end_limit = vars(index).Range(2);
                                    if min_limit > start_limit
                                        if min_limit > end_limit
                                            % % if min_limit is more than upper boundary; fixing the variable to upper boundary                                          
                                            % fprintf('Fixing the variable %s to %d\n', var_to_update, vars(index).Range(end));
                                            % fprintf('Minimum limit of %s is %d which is more than the specified upper boundary %d\n', var_to_update, min_limit, vars(index).Range(end));
                                            % vars(index) = [];
                                            % % removing the existing results of variables
                                            % X_initial_bayes.(var_to_update) = [];
                                            % % adding the variable to fixed
                                            % % variable list
                                            % fixed_variables.(var_to_update) = end_limit;
                                            text = sprintf('Optimization not possible using ABC-BO; not possible to improve the objective from existing best;min limit %d of variable %s is more than the specified upper boundary %d ', min_limit, var_to_update, end_limit);
                                            error(text);
                                        else
                                            % update lower boundary if max limit is less than upper boundary
                                            vars(index).Range(1) = min_limit;
                                            fprintf('Adjusted lower boundary of %s from %d to %d\n', var_to_update, start_limit, min_limit);
                                            fprintf('Updated limits for %s [%d, %d]\n', var_to_update, vars(index).Range(1), vars(index).Range(end));
                                        end
                                    else
                                        fprintf('No boundary update for %s; min limit %d is less than specified lower boundary %d\n', var_to_update, min_limit, vars(index).Range(1));
                                    end
                                end
                            end

                            % discrete numeric lists update
                            if ismember(var_to_update,problem_details.dis_num_var_labels)
                                if strcmp(optimization_details.ABC_BO.boundary_update.which_limit{i},'upper')
                                    index_dn = find(strcmp(var_to_update,problem_details.dis_num_var_labels));
                                    end_limit = problem_details.dis_num_var_lists{index_dn}(end);
                                    start_limit = problem_details.dis_num_var_lists{index_dn}(1);
                                    max_limit = optimization_details.ABC_BO.boundary_update.limit.(var_to_update);
                                    if max_limit < end_limit
                                        list = problem_details.dis_num_var_lists{index_dn};
                                        if max_limit < start_limit
                                            % % if max_limit is less than lower boundary; fixing the upper boundary same as lower boundary
                                            % fprintf('Fixing the variable %s to %d\n', var_to_update, start_limit);
                                            % fprintf('Maximum limit of %s is %d which is less than the specified lower boundary %d\n', var_to_update, max_limit, start_limit);
                                            % vars(index) = [];
                                            % % removing the existing results of variables
                                            % X_initial_bayes.(var_to_update) = [];
                                            % % adding the variable to fixed
                                            % % variable list
                                            % fixed_variables.(var_to_update) = start_limit;
                                            text = sprintf('Optimization not possible using ABC-BO; not possible to improve the objective from existing best;max limit %d of variable %s is less than the specified lower boundary %d ', max_limit, var_to_update, start_limit);
                                            error(text);
                                            
                                        else
                                            % update upper boundary if max limit is greater than lower boundary
                                            new_list = list(list<=max_limit);

                                            if isscalar(new_list)
                                                % defining single list in discrete numeric is not possible
                                                % providing 2 list
                                                % anyway futile conditions can be avoided by constraints
                                                new_list = list(1:2);
                                            end
                                            vars(index).Range(end) = numel(new_list);
                                            formatted_list = sprintf('%g ', new_list);
                                            fprintf('Updated list for %s: %s\n', var_to_update, formatted_list);
                                        end

                                    else
                                        fprintf('No update for %s; max limit %d is more than any other number in the list\n', var_to_update, max_limit);
                                    end

                                elseif strcmp(optimization_details.ABC_BO.boundary_update.which_limit{i},'lower')
                                    index_dn = find(strcmp(var_to_update,problem_details.dis_num_var_labels));
                                    start_limit = problem_details.dis_num_var_lists{index_dn}(1);
                                    end_limit = problem_details.dis_num_var_lists{index_dn}(end);
                                    min_limit = optimization_details.ABC_BO.boundary_update.limit.(var_to_update);
                                    if min_limit > start_limit
                                        list = problem_details.dis_num_var_lists{index_dn};
                                        if min_limit > end_limit
                                            % % if min_limit is more than upper boundary; fixing the variable to upper boundary                                          
                                            % fprintf('Fixing the variable %s to %d\n', var_to_update, end_limit);
                                            % fprintf('Minimum limit of %s is %d which is more than the specified upper boundary %d\n', var_to_update, min_limit, end_limit);
                                            % vars(index) = [];
                                            % % removing the existing results of variables
                                            % X_initial_bayes.(var_to_update) = [];
                                            % % adding the variable to fixed
                                            % % variable list
                                            % fixed_variables.(var_to_update) = end_limit;
                                            text = sprintf('Optimization not possible using ABC-BO; not possible to improve the objective from existing best;min limit %d of variable %s is more than the specified upper boundary %d ', min_limit, var_to_update, end_limit);
                                            error(text);
                                        else
                                            % update upper boundary if max limit is greater than lower boundary
                                            new_list = list(list>=min_limit);

                                            if isscalar(new_list)
                                                % defining single list in discrete numeric is not possible
                                                % providing 2 list
                                                % anyway futile conditions can be avoided by constraints
                                                new_list = list(end-1:end);
                                            end
                                            start_new_list = new_list(1);
                                            int_new_start = dis_num_var_map_to_int.(var_to_update)(start_new_list);
                                            vars(index).Range(1) = int_new_start;
                                            formatted_list = sprintf('%g ', new_list);
                                            fprintf('Updated list for %s: %s\n', var_to_update, formatted_list);
                                        end

                                    else
                                        fprintf('No update for %s; min limit %d is less than any other number in the list\n', var_to_update, min_limit);
                                    end
                                end
                            end

                            % categorical 
                            % removing levels that are futile                            
                            if ismember(var_to_update,problem_details.cat_var_labels)
                                % Ensuring at least 2 levels should be there in
                                % categorical variable list
                                % Anyway constraint take care of futile
                                indices_to_remove = find(ismember(vars(index).Range, optimization_details.ABC_BO.boundary_update.limit.(var_to_update)));
                                if numel(vars(index).Range)-numel(indices_to_remove) <2
                                    while ~(numel(vars(index).Range)-numel(indices_to_remove)>=2)
                                        indices_to_remove(end) = [];
                                    end
                                end
                                if indices_to_remove                                   
                                    vars(index).Range(indices_to_remove) = [];
                                    rangeStr = strjoin(vars(index).Range, ', ');
                                    fprintf('Updated list for %s: %s\n', var_to_update, rangeStr);
                                end
                            end


                             

                        end
                
            end
%% initiate bayesopt
            bayes_obj = bayesopt(obj_fun, vars, 'MaxObjectiveEvaluations', n_point, 'InitialX', X_initial_bayes,...
                'InitialObjective',Y_initial_bayes,'AcquisitionFunctionName',optimization_details.acq_func, ...
                'NumCoupledConstraints', n_obj_constraints,'InitialConstraintViolations',constraint_violations,...
                'XConstraintFcn', variable_constraint,...
                'ExplorationRatio', optimization_details.epsilon,'Verbose',0,PlotFcn=[]);

            % save next point
            % continuous variables
            if problem_details.n_cont_var >0
                for i = 1:problem_details.n_cont_var
                    field_name = problem_details.cont_var_labels(i);
                    field_name = char(field_name);
                    next_point.point_1.(field_name) = bayes_obj.XTrace.(field_name)(end);
                end
            end

            % discrete numeric variables
            if problem_details.n_dis_num_var >0
                for i = 1:problem_details.n_dis_num_var
                    field_name = problem_details.dis_num_var_labels(i);
                    field_name = char(field_name);
                    int_point = bayes_obj.XTrace.(field_name)(end);
                    % map to real value
                    real_point = dis_num_var_map_to_val.(field_name)(int_point);
                    next_point.point_1.(field_name) = real_point;
                end
            end

            % categorical variables
            if problem_details.n_cat_var >0
                for i = 1:problem_details.n_cat_var
                    field_name = problem_details.cat_var_labels(i);
                    field_name = char(field_name);
                    next_point.point_1.(field_name) = char(bayes_obj.XTrace.(field_name)(end));
                end
            end


            
            function [obj, cons] = dummy_fun(x, n_constraints)

                % dummy value to get the next point
                obj = 0;

                % some negative value to indicate coupled constraints are
                % satisfied 
                if exist("n_constraints","var")
                    cons = -20 *ones(1,n_constraints);
                end
            end

            % variable constraints
            function tf = variable_constraint_func(x,constraints,variable_details)

                n_samples = size(x,1);

                % converting discrete numeric variables to real values
                % before checking constraints
                if ~isempty(variable_details.n_dis_num_var) && variable_details.n_dis_num_var
                    var = variable_details.dis_num_var_labels; % dis num var names                    
                    for i = 1:variable_details.n_dis_num_var                        
                        for j = 1:n_samples
                            x.(char(var(i)))(j) = variable_details.dis_num_var_map_to_val.(char(var(i)))(x.(char(var(i)))(j));
                        end                        
                    end
                end

                %  checking variable satisfaction
                n_constraints = numel(fieldnames(constraints));
                names = fieldnames(constraints);
                tf = true(n_samples,1);
                for i = 1:n_constraints
                    tf = tf & constraints.(names{i})(x);
                end
               
                p = 1;

            end






        end

        

       
    end
end


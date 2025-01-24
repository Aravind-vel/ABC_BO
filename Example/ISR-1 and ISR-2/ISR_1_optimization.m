clc
clear all

% add optimization folder to path
addpath('C:\Users\senthilvel-a\..................................\Optimization folder');


%% Define optimization problem

problem_details = struct();

% Ensure the variable labels are consistent with those used elsewhere, such as in constraints, class file that calculate objective.

% continuous variables
problem_details.n_cont_var = 4;
problem_details.cont_var_labels = {'temp','time','reag_eq','cat_load'};
problem_details.cont_var_bounds = [25,1,1,5;50,10,2,20]; %  [lower bounds;upper bounds]

% objective
problem_details.n_obj = 1;
problem_details.obj_labels = {'throughput'}; % g/h
problem_details.obj_criteria = {'max'};

%% reaction and optimization instance

react_obj = ISR_reaction_class; % calculation of objectives
opt_obj = Optimization(problem_details); % to perform optimization

%% sampling

% center point sampling - center point in each combination of categorical variables
sampling_points = opt_obj.Generate_samples('CP');
% LHS sampling
% sampling_points = opt_obj.Generate_samples('LHS','n_samples',5); % n_samples per level of combination of categorical variables

% total number of sampling points
n_sampling_points = numel(fieldnames(sampling_points));
fprintf('Total number of sampling points: %d\n', n_sampling_points);


% objective calculation and updating results
field_names = fieldnames(sampling_points);
yield = zeros(n_sampling_points,1);
for i = 1:numel(field_names)

    % ith sampling point
    field_name = field_names{i};
    temp = sampling_points.(field_name).temp;
    time = sampling_points.(field_name).time;
    reag_eq = sampling_points.(field_name).reag_eq;
    cat_load = sampling_points.(field_name).cat_load;
    yield(i) = ISR_1_insilico(temp, time, reag_eq, cat_load);

    % objective calculation
    objective = react_obj.Objective_value_calculation(sampling_points.(field_name),yield(i)); % g h-1 euro-1
    % update result
    opt_obj = opt_obj.Update_result(sampling_points.(field_name),objective);
end

%% standard bayesian optimization

% n_opt = 30;
%
% for i = 1:n_opt
%
%     % generate next point
%     next_point = opt_obj.Next_points_to_opt('bayesopt','log_obj', false,'acq_func','expected-improvement-plus','epsilon',0.5);
%
%     temp = next_point.point_1.temp;
%     time = next_point.point_1.time;
%     reag_eq = next_point.point_1.reag_eq;
%     cat_load = next_point.point_1.cat_load;
%
%     % predict yield - in silico
%     yield(end+1) = ISR_1_insilico(temp, time, reag_eq, cat_load);
%
%     % objective calculation
%     objective = react_obj.Objective_value_calculation(next_point.point_1, yield(end));
%     % update result
%     opt_obj = opt_obj.Update_result(next_point.point_1, objective);
% end

%% ABC-BO
n_opt = 5;

for i = 1:n_opt
    % best objective value
    [cbo, pos] = max(opt_obj.Objective_trend);

    % yield corresponding to best objective value
    best_yield = yield(pos);

    % experimental condition (variables) corresponding to best objective value
    for i = 1:opt_obj.n_var
        field_name = opt_obj.var_labels{i};
        best_condition.(field_name) = opt_obj.results.(field_name)(pos);
    end

    yield_noise = 0; % noise to avoid strict constraint

    % best objective value updated with noise
    cbo_t  = react_obj.Objective_value_calculation(best_condition, max(best_yield-yield_noise, 0));

    % ABC BO info
    ABC_BO = struct();

    % max. achievable objective > best objective (corrected with noise)
    ABC_BO.constraint = @(var) react_obj.Objective_value_calculation(var,100) > cbo_t;

    % variables eligible for update
    ABC_BO.boundary_update.variables = {'time','reag_eq','cat_load'};

    % 'upper' or 'lower' limit to update for continuous and discrete numeric variables
    % 'levels' if categorical variable levels can be removed
    ABC_BO.boundary_update.which_limit = {'upper','upper','upper'};

    % boundary update function provide the limit for continuous and discrete numeric variables
    % above which the condition is futile
    % for categorical variables it provide the levels that are futile
    futile_limit = react_obj.Boundary_update(cbo_t);
    ABC_BO.boundary_update.limit.time = futile_limit.time;
    ABC_BO.boundary_update.limit.cat_load = futile_limit.cat_load;
    ABC_BO.boundary_update.limit.reag_eq = futile_limit.reag_eq;
    
    % It is not must to provide the variables limit.  In case of
    % difficulties it can be avoided. Constraint itself take care of futile
    % experiments.  However, if total number of futile space is
    % comparatively large - optimization may fail to find the feasible
    % point

    % initiate ABC-BO
    next_point = opt_obj.Next_points_to_opt('bayesopt','ABC_BO',ABC_BO,'log_obj', false,'acq_func','expected-improvement-plus','epsilon',0.5);

    temp = next_point.point_1.temp;
    time = next_point.point_1.time;
    reag_eq = next_point.point_1.reag_eq;
    cat_load = next_point.point_1.cat_load;

    % predict yield - in silico
    yield(end+1) = ISR_1_insilico(temp, time, reag_eq, cat_load);

    % objective calculation
    objective = react_obj.Objective_value_calculation(next_point.point_1, yield(end));
    % update result
    opt_obj = opt_obj.Update_result(next_point.point_1, objective);

end

%% Tracking best objective
best_obj_track = opt_obj.Objective_trend;
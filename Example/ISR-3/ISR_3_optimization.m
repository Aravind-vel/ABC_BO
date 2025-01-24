clc
clear all

% add optimization folder to path
addpath('C:\Users\senthilvel-a\Documents\_UNC\phython practice\personal files\matlab_bo\project constraints\github repository\Optimization folder');


%% Define optimization problem

problem_details = struct();

% Ensure the variable labels are consistent with those used elsewhere, such as in constraints, class file that calculate objective.

% continuous variables
problem_details.n_cont_var = 3;
problem_details.cont_var_labels = {'temp','time','cat_load'};
problem_details.cont_var_bounds = [30,1,0.5;110,10,2.5]; %  [lower bounds;upper bounds]

%categorical variables
problem_details.n_cat_var = 1;
problem_details.cat_var_labels = {'catalyst'}; % categorical variables
problem_details.cat_var_lev_labels = {{'1','2','3','4','5','6','7','8'}}; % levels of the categorical variables

% objective
problem_details.n_obj = 1;
problem_details.obj_labels = {'TON'};
problem_details.obj_criteria = {'max'};

%% reaction and optimization instance

react_obj = ISR_3_reaction_class;
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
    cat_load = sampling_points.(field_name).cat_load;
    catalyst = sampling_points.(field_name).catalyst;

    objective = ISR_3_insilico(time, temp, cat_load, str2num(catalyst));

    % update result
    opt_obj = opt_obj.Update_result(sampling_points.(field_name),objective);
end

%% standard bayesian optimization

n_opt = 10;

% for i = 1:n_opt
% 
%     % generate next point
%     next_point = opt_obj.Next_points_to_opt('bayesopt','log_obj', false,'acq_func','expected-improvement-plus','epsilon',0.5);
% 
%     temp = next_point.point_1.temp;
%     time = next_point.point_1.time;
%     cat_load = next_point.point_1.cat_load;
%     catalyst = next_point.point_1.catalyst;
% 
%     % predict yield - in silico
%     objective = ISR_3_insilico(time, temp, cat_load, str2num(catalyst));
% 
% 
%     % update result
%     opt_obj = opt_obj.Update_result(next_point.point_1, objective);
% end

%% ABC-BO
n_opt = 20;

for i = 1:n_opt

    % best objective value
    [cbo, pos] = max(opt_obj.Objective_trend);


    % ABC BO info
    ABC_BO = struct();

    % max. achievable objective > best objective (corrected with noise)
    ABC_BO.constraint = @(var) react_obj.Objective_value_calculation(var,100) > cbo;

    % variables eligible for update
    ABC_BO.boundary_update.variables = {'cat_load'};

    % 'upper' or 'lower' limit to update for continuous and discrete numeric variables
    % 'levels' if categorical variable levels can be removed
    ABC_BO.boundary_update.which_limit = {'upper'};

    % boundary update function provide the limit for continuous and discrete numeric variables
    % above which the condition is futile
    % for categorical variables it provide the levels that are futile
    futile_limit = react_obj.Boundary_update(cbo);
    
    ABC_BO.boundary_update.limit.cat_load = futile_limit.cat_load;
    disp(futile_limit.cat_load)

    % initiate ABC-BO
    next_point = opt_obj.Next_points_to_opt('bayesopt','ABC_BO',ABC_BO,'log_obj', false,'acq_func','expected-improvement-plus','epsilon',0.5);

    temp = next_point.point_1.temp;
    time = next_point.point_1.time;
    cat_load = next_point.point_1.cat_load;
    catalyst = next_point.point_1.catalyst;

    % objective calculation
    objective = ISR_3_insilico(time, temp, cat_load, str2num(catalyst));
    % update result
    opt_obj = opt_obj.Update_result(next_point.point_1, objective);

end

%% Tracking best objective
best_obj_track = opt_obj.Objective_trend;
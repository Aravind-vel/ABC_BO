clc 
clear all

% add optimization folder to path
addpath('C:\Users\senthilvel-a\Documents\_UNC\phython practice\personal files\matlab_bo\project constraints\github repository\Optimization folder');


%% Define optimization problem

problem_details = struct();

% make sure the variable labels are same as in the reaction class file

% continuous variables
problem_details.n_cont_var = 1; % cat. loading
problem_details.cont_var_bounds = [0.5;5]; % mol%  [lb;up]
problem_details.cont_var_labels = {'catalyst_conc'};

% discrete numeric variables
problem_details.n_dis_num_var = 3;
problem_details.dis_num_var_lists = {[50:10:110];[0.5:0.5:8];[1:0.1:2]};
problem_details.dis_num_var_labels = {'temp','time','elec_stoic'};

% categorical variables
problem_details.n_cat_var = 3;
problem_details.cat_var_labels = {'base','ligand','electrophile'}; % categorical variables
problem_details.cat_var_lev_labels = {{'tBuKO','Cs2CO3'};{'Johnphos','Xantphos', 'PPh3'};{'PhBr','PhCl'}}; % levels of the categorical variables

% objective
problem_details.n_obj = 1;
problem_details.obj_labels = {'objective'}; % productivity/cost
problem_details.obj_criteria = {'max'};

%% Reaction and optimization instance

react_obj = Buchwald_Hardwig_coupling;  % for reaction related calculations like objective function, boundary update details

opt_obj = Optimization(problem_details); % to perform optimization

%% sampling 
sampling_points = opt_obj.Generate_samples('CP');
yield = [99.4;1;88.5;0;5.4;0;94;0;9.4;0;7.9;0];

% objective calculation and updating results
field_names = fieldnames(sampling_points);
for i = 1:numel(field_names)
    % ith sampling point 
    field_name = field_names{i};    
    % objective calculation
    objective = react_obj.Objective_value_calculation(sampling_points.(field_name),yield(i)); % g h-1 euro-1
    % update result 
    opt_obj = opt_obj.Update_result(sampling_points.(field_name),objective);
end

%% standard optimization

% next_point = opt_obj.Next_points_to_opt('bayesopt','log_obj', true,'acq_func','expected-improvement-plus','epsilon',0.5);


%% ABC-BO

% best objective value
[cbo, pos] = max(opt_obj.Objective_trend);

% yield corresponding to best objective value
best_yield = yield(pos);

% experimental condition (variables) corresponding to best objective value
for i = 1:opt_obj.n_var
field_name = opt_obj.var_labels{i};
    best_condition.(field_name) = opt_obj.results.(field_name)(pos);
end

yield_noise = 5; % noise to avoid strict constraint

% best objective value updated with noise
cbo_t  = react_obj.Objective_value_calculation(best_condition, max(best_yield-yield_noise, 0));

% ABC BO info
ABC_BO = struct();

% max. achievable objective > best objective (corrected with noise)
ABC_BO.constraint = @(var) react_obj.Objective_value_calculation(var,100) > cbo_t;

% variables eligible for update
ABC_BO.boundary_update.variables = {'catalyst_conc','time','elec_stoic','base','ligand','electrophile'};

% 'upper' or 'lower' limit to update for continuous and discrete numeric variables
% 'levels' if categorical variable levels can be removed
ABC_BO.boundary_update.which_limit = {'upper','upper','upper','levels','levels','levels'};

% boundary update function provide the limit for continuous and discrete numeric variables
% above which the condition is futile
% for categorical variables it provide the levels that are futile
futile_limit = react_obj.Boundary_update(cbo_t);

ABC_BO.boundary_update.limit.catalyst_conc = futile_limit.catalyst_conc;
ABC_BO.boundary_update.limit.time = futile_limit.time;
ABC_BO.boundary_update.limit.elec_stoic = futile_limit.elec_stoic;
ABC_BO.boundary_update.limit.base = futile_limit.base;
ABC_BO.boundary_update.limit.ligand = futile_limit.ligand;
ABC_BO.boundary_update.limit.electrophile = futile_limit.electrophile;

% initiate ABC-BO
next_point = opt_obj.Next_points_to_opt('bayesopt','ABC_BO',ABC_BO,'log_obj', true,'acq_func','expected-improvement-plus','epsilon',0.5);



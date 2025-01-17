classdef Optimization %Comment added
    %{
optimization - initialize object with all the data need to initiate the
optimization and executes the required optimization algorithm thus provide
the next best point to evaluate

object is created by using the information about the variable and objective defined.
All the information needed (properties) for the optimization that are need to
be used by the different optimization algorithm are initialized in this
object. 

different optimization algorithms are defined as methods.

# Single-Objective Optimization 
bayesopt (MATLAB inbuilt) - bayesopt_algo


To initiate the optimization

1.  Crete object using the data (in struct)
        o = optimization(variables, objectives);
2.  This object is executed/called using the algorithm that we need which
    returns the next best point to evaluate.  
        next_point = o.bayesopt_algo;
    %}
    %   Detailed explanation goes here

    properties (Access = public)
        variables  % variables of the existing results in array
        objectives % objective values of the existing results in array

        results % overall results in table

        n_cont_var % no. of continuous variables
        cont_var_bounds % bounds of the continuous variables; first column lower bound, second column upper bound

        cont_var_labels % labels of continuous variables


        n_dis_num_var % no. of discrete numeric variable
        dis_num_var_lists % value of the discrete numeric in cell format

        dis_num_var_labels % labels of discrete numeric variables


        n_cat_var % no. of categorical variables
        cat_var_levels % levels of categorical variables
        cat_var_lev_labels % labels of categorical variables as cell
        cat_var_label_to_num_map % object that convert categorical label to appropriate number
        cat_var_num_to_label_map % object that convert appropriate number to categorical label

        cat_var_labels % labels of categorical variables

        n_obj      % no. of objectives to optimize
        obj_labels % labels of the objectives
        obj_criteria  % objective criteria 'min' or 'max'; by default 'max'

        n_var % total number of variables
        var_labels % labels of the variables

        n_exp_done % No. of experiments finished/done

        constraints % constraint to be satisfied as text
    end

    methods
        function this = Optimization(problem_details)

            % problem details - structure - contains info regarding the
            % optimization problem from which properties of the class are
            % generated

            % in problem details create only field for the variables
            % involved in the problem. For example, provide n_cat_var if
            % you have categorical variable else ignore creating the fields

            this.var_labels = [];
            this.n_var = 0;
            var_types = string([]);

            % continuous variables
            if isfield(problem_details,'n_cont_var')
                this.n_cont_var = problem_details.n_cont_var;
                this.cont_var_bounds = problem_details.cont_var_bounds;
                this.cont_var_labels = problem_details.cont_var_labels;

                this.n_var = this.n_cont_var;
                this.var_labels = [this.var_labels,this.cont_var_labels];
                var_types = [var_types, repmat({"double"},1,this.n_cont_var)];
            end

            % discrete numeric variables
            if isfield(problem_details,'n_dis_num_var')
                this.n_dis_num_var = problem_details.n_dis_num_var;
                this.dis_num_var_lists = problem_details.dis_num_var_lists;
                this.dis_num_var_labels = problem_details.dis_num_var_labels;

                this.n_var = this.n_var + this.n_dis_num_var;
                this.var_labels = [this.var_labels,this.dis_num_var_labels];
                var_types = [var_types, repmat({"double"},1,this.n_dis_num_var)];
            end

            % categorical variables
            if isfield(problem_details, 'n_cat_var')
                this.n_cat_var = problem_details.n_cat_var;
                this.cat_var_lev_labels = problem_details.cat_var_lev_labels;

                % no. of levels in categorical variables
                this.cat_var_labels = zeros(1,this.n_cat_var);
                for i = 1:this.n_cat_var
                    this.cat_var_levels(1,i) = numel(this.cat_var_lev_labels{i,1});
                end

                this.n_var = this.n_var + this.n_cat_var;
                this.cat_var_labels = problem_details.cat_var_labels;

                % mapping object for number to text and text to number
                for i = 1:this.n_cat_var
                    label_as_number = 1:1:this.cat_var_levels(i);
                    field_name = char(this.cat_var_labels(i));
                    this.cat_var_label_to_num_map.(field_name) = dictionary(this.cat_var_lev_labels{i,1},label_as_number);
                    this.cat_var_num_to_label_map.(field_name) = dictionary(label_as_number,this.cat_var_lev_labels{i,1});
                end

                this.var_labels = [this.var_labels,this.cat_var_labels];
                var_types = [var_types, repmat({"string"},1,this.n_cat_var)];
            end



            % labels of variables
            this.var_labels = [this.cont_var_labels,this.dis_num_var_labels,this.cat_var_labels];

            % objective
            this.n_obj = problem_details.n_obj;
            this.obj_labels = problem_details.obj_labels;
            this.obj_criteria = problem_details.obj_criteria;

            var_types = [var_types, repmat({"double"},1,this.n_obj)];


            % %% results table
            total_var_count = this.n_var+this.n_obj;

            % var types
            this.results = table('Size', [0, total_var_count],'VariableTypes',var_types,'VariableNames',[this.var_labels,this.obj_labels]);
            % this.results = table('Size', [0, total_var_count],'VariableNames',[this.var_labels,this.obj_labels]);

            % constraint
            if isfield(problem_details,'constraints')
                this.constraints = problem_details.constraints;
            end

        end

        function [sampling_points_struct, sampling_points_table,sampling_points] = Generate_samples(this,sample_type,varargin)

            % varargin options
            % n_samples - per level of categorical variables
            % saving input arguments as structure - sampling details
            for i = 1:2:size(varargin, 2)
                key = varargin{1,i};
                value = varargin{1,i+1};

                sampling_details.(key) = value;
            end


            switch sample_type
                case 'LHS'
                    sampling_points = Sampling.LHS(sampling_details.n_samples,this.n_cont_var, this.cont_var_bounds,...
                        this.n_cat_var, this.cat_var_levels,...
                        this.n_dis_num_var, this.dis_num_var_lists);

                    % sampling points in table
                    sampling_points_table = Optimization.sampling_results_in_table(sampling_points,this);
                    
                case 'CP' 
                    % Center point on each level of the categorical
                    % variable
                    sampling_points = Sampling.CP(this.n_cont_var, this.cont_var_bounds,...
                        this.n_cat_var, this.cat_var_levels,...
                        this.n_dis_num_var, this.dis_num_var_lists);

                     % sampling points in table
                    sampling_points_table = Optimization.sampling_results_in_table(sampling_points,this);                   

            end

            % sampling results in structure           
            for i = 1:size(sampling_points,1)
                field_name = sprintf('point_%d',i);
                for j = 1:this.n_var
                    field_name_var = this.var_labels{j};
                    if ~isnumeric(sampling_points_table.(field_name_var)(i))
                    sampling_points_struct.(field_name).(field_name_var) = char(sampling_points_table.(field_name_var)(i));
                    else
                        sampling_points_struct.(field_name).(field_name_var) = sampling_points_table.(field_name_var)(i);
                    end
                end
            end

            

        end

        function next_points = Next_points_to_opt(this,opt_algorithm, varargin)

            switch opt_algorithm
                case 'bayesopt'
                    %% varargin

                    % log_obj - true or false
                    % acq_func - 'expected-improvement-per-second-plus' | 'expected-improvement' | 'expected-improvement-plus'(default)  | 'expected-improvement-per-second' | 'lower-confidence-bound' | 'probability-of-improvement'
                    % epsilon - scalar real value - exploration factor - only for expected improvement plus and expected improvement per second plus
                    % objective_constraints - obj_constraints.constraint_1 @(par) par.x1+par.x2 
                    %    par can contain both variable and objective labels
                    %    if the value is negative (<=0) then the constraint is satisfied
                    %    it can be multiple constraints constraint_2, constraint_3
                    % variable_constraints - var_constraints.constraint_1 @(par) par.x1 +par.x2 <=100
                    %     par should contain only one variable labels
                    %     full constraint equation to be satisfied should be mentioned
                    %     multiple constraints constraint_1, constraint_2 can be specified
                    next_points = Optimization_algorithm.bayesopt_algo(this, varargin);

            end


        end

        function this = Update_result(this, variable, objective)

            %% update results table
            data = variable;
            for i = 1:this.n_obj
                field_name = this.obj_labels{i};
                field_name = char(field_name);
                data.(field_name) = objective(i);
            end

            data = struct2table(data);
            this.results = [this.results;data];

            %% update variables and objectives as array

            this.objectives(end+1,:) = objective;

            this.variables(end+1,:) = zeros(1,this.n_var);
            n = 1;
            if ~isempty(this.n_cont_var) && this.n_cont_var
                for i = 1:this.n_cont_var
                    field_name = this.cont_var_labels{i};
                    field_name = char(field_name);
                    this.variables(end,n) = variable.(field_name);
                    n = n+1;
                end
            end

            if ~isempty(this.n_dis_num_var) && this.n_dis_num_var
                for i = 1:this.n_dis_num_var
                    field_name = this.dis_num_var_labels{i};
                    field_name = char(field_name);
                    this.variables(end,n) = variable.(field_name);
                    n = n+1;
                end
            end

            if ~isempty(this.n_cat_var) && this.n_cat_var
                for i = 1:this.n_cat_var
                    field_name = this.cat_var_labels{i};
                    field_name = char(field_name);
                    value_as_label = {variable.(field_name)};
                    value_as_num = this.cat_var_label_to_num_map.(field_name)(value_as_label);
                    this.variables(end,n) = value_as_num;
                    n = n+1;
                end
            end

        end

        function obj_trend = Objective_trend(this)
            obj_trend = Results_analysis.Objective_trend(this.objectives, this.obj_criteria);
        end


        function [pareto_points, pareto_variables] = Pareto_calculation(this)

            [pareto_points, pareto_variables] = Results_analysis.find_pareto(this.objectives, this.variables, this.obj_criteria);

        end

        function HV_track = Hypervolume_track(this, Ref_ideal_pt, Ref_anti_ideal_pt)



            HV_track = Results_analysis.Hypervolume_track(this.objectives, this.obj_criteria, Ref_ideal_pt, Ref_anti_ideal_pt);
        end



        function  this = update_variables(this, data)

            data_array = zeros(size(data));

            % continuous and discrete numeric data in array
            data_array(:,1:this.n_var-this.n_cat_var) = cell2mat(data(:,1:this.n_var-this.n_cat_var));

            % categorical variables are in the last column
            % converting categorical variables to number - mapping
            if this.n_cat_var > 0
                col = this.n_var-this.n_cat_var+1; % first col. of cat. var
                cat_data = [];
                for i = 1:this.n_cat_var
                    map = this.map_label_to_num{1};
                    for j = 1:length(data) % all data in a cat. var
                        cat_data(j,i) = map(data{j,col});
                    end
                    col = col+1; % next cat. var
                end
                col = this.n_var-this.n_cat_var+1; % first col. of cat. var
                data_array(:,col:end) = cat_data;
            end

            % update variables in the object
            this.variables = [this.variables;data_array];

            % update tables in results
            % Nan values for the new row

            len_before_update = height(this.results);

            n_r = length(data);
            n_c = width(this.results);
            this.results(end+1:end+n_r,:) = array2table(repmat({NaN}, n_r, n_c), 'VariableNames', this.results.Properties.VariableNames);

            % update data
            % for i = 1:this.n_var
            this.results{len_before_update+1:end,1:this.n_var} = data;
            % end
            p =1;

        end


    end

    methods (Static)

        function table = sampling_results_in_table(array,prob_det)

            table = prob_det.results(:,1:prob_det.n_var); % removing columns for objectives



                    for i = 1:size(array,1)

                        n = 1;
                        if ~isempty(prob_det.n_cont_var) && prob_det.n_cont_var
                            for j = 1:prob_det.n_cont_var
                                field_name = prob_det.cont_var_labels{j};
                                field_name = char(field_name);
                                variable.(field_name) = array(i,n);
                                n = n+1;
                            end
                        end
                        if ~isempty(prob_det.n_dis_num_var) && prob_det.n_dis_num_var
                            for j = 1:prob_det.n_dis_num_var
                                field_name = prob_det.dis_num_var_labels{j};
                                field_name = char(field_name);
                                variable.(field_name) = array(i,n);
                                n = n+1;
                            end
                        end
                        if ~isempty(prob_det.n_cat_var) && prob_det.n_cat_var
                            for j = 1:prob_det.n_cat_var
                                field_name = prob_det.cat_var_labels{j};
                                field_name = char(field_name);
                                value_as_label = prob_det.cat_var_num_to_label_map.(field_name)(array(i,n));
                                variable.(field_name) = char(value_as_label);
                                n = n+1;
                            end

                        end

                        % add to table
                        table = [table; struct2table(variable)];
                    end



        end



    end
end


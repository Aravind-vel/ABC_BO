classdef Sampling


    % This Sampling class has static methods (functions) to get samples
    % from different sampling methods (LHS, random, and more)

    properties
        Property1
    end

    methods (Static)
        function sampling_points = LHS(n_samples,n_cont_var, cont_var_bounds, n_cat_var, cat_var_levels, n_dis_num_var, dis_num_var_lists)

            % n_samples - per level of categorical variables

            % to get lhs design for discrete numeric variable
            % it is treated as continuous variable
            n_var_lhs = 0;
            var_bounds_lhs = []; % [cont vars, dis num vars]

            if ~isempty(n_cont_var) && n_cont_var
                n_var_lhs = n_var_lhs+n_cont_var;
                for i = 1:n_cont_var
                    var_bounds_lhs(1,i) = cont_var_bounds(1,i);
                    var_bounds_lhs(2,i) = cont_var_bounds(2,i);

                end
            end

            % treating dis num as continuous variable
            if ~isempty(n_dis_num_var) && n_dis_num_var
                n_var_lhs = n_var_lhs+n_dis_num_var;


                for i = 1:n_dis_num_var
                    if ~isempty(n_cont_var) && n_cont_var, idx = n_cont_var+i; else, idx = i; end
                    var_bounds_lhs(1,idx) = min(dis_num_var_lists{i,1});
                    var_bounds_lhs(2,idx) = max(dis_num_var_lists{i,1});
                end
            end

            if ~isempty(n_cat_var)
                tot_n_samples = n_samples*prod(cat_var_levels);
            else
                tot_n_samples = n_samples;
            end

            % creating LHS design - for continuous space
            lhs_design = lhsdesign(n_samples, n_var_lhs,"criterion","maximin","iterations",2000);


            %lhs samples for the boundary of continuous space
            lhs_samples = ((var_bounds_lhs(2,:) - var_bounds_lhs(1,:)) .* lhs_design) + var_bounds_lhs(1,:);


            % real value of discrete numeric variable is moved to closest
            % value in dis num lists
            if ~isempty(n_dis_num_var) && n_dis_num_var
                for i  = 1:n_dis_num_var
                    if ~isempty(n_cont_var) && n_cont_var, idx = n_cont_var+i; else, idx = i; end
                    lhs_samples(:,idx) = interp1(dis_num_var_lists{i,1},dis_num_var_lists{i,1},lhs_samples(:,idx),'nearest');
                end
            end

            % repeat lhs samples for each level of categorical variables
            if ~isempty(n_cat_var) && n_cat_var

                %creating different discrete combinations
                % Generate a cell array of ranges for each element in arr
                ranges = arrayfun(@(n) 1:n, cat_var_levels, 'UniformOutput', false);
                % Generate grids for each range
                [grids{1:length(cat_var_levels)}] = ndgrid(ranges{:});
                % Concatenate the grids into a matrix of combinations
                % Each grid is converted to a column vector and concatenated horizontally
                combinations = cell2mat(cellfun(@(c) c(:), grids, 'UniformOutput', false));
                n_combinations = size(combinations,1);

                % Replicating combinations
                n = 0;
                combinations_rep = zeros(n_samples*n_combinations,n_cat_var);
                for i = 1:n_combinations
                    combinations_rep(n+1:n+n_samples,:) = repmat(combinations(i,:),n_samples,1);
                    n = n+n_samples;
                end
            else
                combinations = 1;
                n_combinations = 1;
            end

            % replicating lhs samples
            lhs_samples = repmat(lhs_samples, n_combinations,1);




            % concacting lhs samples with categorical variables
            if ~isempty(n_cat_var) && n_cat_var
                lhs_samples = [lhs_samples,combinations_rep];
            end
            sampling_points = lhs_samples;


        end

        function sampling_points = CP(n_cont_var, cont_var_bounds, n_cat_var, cat_var_levels, n_dis_num_var, dis_num_var_lists)

            sampling_points = [];

            % continuous variable center point
            if ~isempty(n_cont_var) && n_cont_var
                avg_points = (cont_var_bounds(1,:) + cont_var_bounds(2,:))/2;
                sampling_points(1,:) = avg_points;
            end

            % discrete numeric center point (median)
            if ~isempty(n_dis_num_var) && n_dis_num_var
                for i = 1:n_dis_num_var
                    median_point = Custom_median(dis_num_var_lists{i,1});
                    sampling_points(1,end+1) = median_point;
                end
            end

            % categorical variables - repeat the lists for all combo
            if ~isempty(n_cat_var) && n_cat_var
                %creating different discrete combinations
                % Generate a cell array of ranges for each element in arr
                ranges = arrayfun(@(n) 1:n, cat_var_levels, 'UniformOutput', false);
                % Generate grids for each range
                [grids{1:length(cat_var_levels)}] = ndgrid(ranges{:});
                % Concatenate the grids into a matrix of combinations
                % Each grid is converted to a column vector and concatenated horizontally
                combinations = cell2mat(cellfun(@(c) c(:), grids, 'UniformOutput', false));
                n_combinations = size(combinations,1);

                % Replicating combinations
                n = 0;
                combinations_rep = zeros(n_combinations,n_cat_var);
                for i = 1:n_combinations
                    combinations_rep(n+1:n+1,:) = repmat(combinations(i,:),1,1);
                    n = n+1;
                end

                sampling_points = repmat(sampling_points,n_combinations,1);

                % concact
                sampling_points = [sampling_points,combinations_rep];
            end


            function med = Custom_median(data)
                % Sort the data
                sortedData = sort(data);
                n = length(sortedData);

                % Check if the number of elements is even
                if mod(n, 2) == 0
                    % If even, return the first of the two middle values
                    med = sortedData(n/2);
                else
                    % If odd, return the middle value
                    med = sortedData(ceil(n/2));
                end
            end
        end
    end
end


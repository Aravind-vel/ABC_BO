classdef Results_analysis
    %RESULTS_ANALYSIS 

    % This class contain functions that are useful for processing results
    % after the optimization
    
    properties
        Property1
    end
    
    methods (Static)
        function obj_trend = Objective_trend(objectives, obj_criteria)

            obj_trend = zeros(size(objectives));

            for i = 1:size(objectives,2)
                if strcmp(obj_criteria{i},'min')
                    obj_trend(:,i) = cummin(objectives(:,i));
                end
                if strcmp(obj_criteria{i},'max')
                    obj_trend(:,i) = cummax(objectives(:,i));
                end
            end

        end
    
        function [paretopoints_obj,paretopoints_var] = find_pareto(Objectives,variables, obj_criteria)

            % following pareto calculation is made for maximizing objectives

            for i = 1:size(Objectives,2)
                if strcmp(obj_criteria{i},'min')
                    Objectives(:,i) = -1.*Objectives(:,i);
                end
            end

            % Initialize an empty logical array to store Pareto front information
            isPareto = true(size(Objectives, 1), 1);

            % Determine the Pareto front
            for i = 1:size(Objectives, 1)
                for j = 1:size(Objectives, 1)
                    % Check if solution i dominates solution j
                    if all(Objectives(i, :) <= Objectives(j, :)) && any(Objectives(i, :) < Objectives(j, :))
                        % Solution i does not dominate solution j
                        isPareto(i) = false;
                        break;  % No need to check further for this solution

                    end
                end
            end

            % restoring the sign
            for i = 1:size(Objectives,2)
                if strcmp(obj_criteria{i},'min')
                    Objectives(:,i) = -1.*Objectives(:,i);
                end
            end

            % Extract the Pareto-optimal solutions
            paretopoints_obj = Objectives(isPareto, :);
            paretopoints_var = variables(isPareto,:);
        end
    
        function HV_track = Hypervolume_track(objectives, obj_criteria, Ref_ideal_pt,Ref_anti_ideal_pt)
            
            dummy_var = zeros(size(objectives));

            % Hypervolume is coded for maximization
            % changing the sign of Reference point for maximization
            for i = 1:numel(obj_criteria)   
                if strcmp(obj_criteria{i},'min')
                    Ref_ideal_pt(i) = -1.*Ref_ideal_pt(i);
                    Ref_anti_ideal_pt(i) = -1.*Ref_anti_ideal_pt(i);
                end
            end


            % pareto and HV track
            for i = 1:size(objectives,1)

                pareto = Results_analysis.find_pareto(objectives(1:i,:),dummy_var,obj_criteria);

                for j = 1:numel(obj_criteria)
                    if strcmp(obj_criteria{j},'min')
                        pareto(:,j) = -1.*pareto(:,j);
                    end
                end



                HV_track(i,1) = Results_analysis.Hypervolume_calculation(pareto, Ref_ideal_pt, Ref_anti_ideal_pt);

            end



        end

        function hypervolume = Hypervolume_calculation(pareto,Ref_ideal_pt,Ref_anti_ideal_pt)


            [n_pareto,n_obj] = size(pareto);

            %% Editable parameters
            n_ran_samples = 100000*n_obj;

            %% Hypervolume calculation

            rand_samples = bsxfun(@plus,Ref_anti_ideal_pt,bsxfun(@times,(Ref_ideal_pt-Ref_anti_ideal_pt),rand(n_ran_samples,n_obj)));

            dominated = 0;
            for i = 1:n_pareto
                idx = sum(bsxfun(@ge,pareto(i,:),rand_samples),2)==n_obj;
                dominated = dominated+sum(idx);
                rand_samples(idx,:) = [];
            end

            hypervolume = (dominated/n_ran_samples)*100;



        end
    end
end


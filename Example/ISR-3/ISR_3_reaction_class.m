classdef ISR_3_reaction_class

    % Class file to calculate the objective and other parameters needed for
    % running ABC-BO. This is file is specific for each reaction
    % optimization problem. 

    properties

        cat_load_bound = [0.5, 2.5]; % mol%
        CA_0 = 0.167; % initial conc. of limiting reagent

    end

    methods
        function this = ISR_3_reaction_class

        end

        function [objective_value] = Objective_value_calculation(this, exp_condition, yield)

            % Return objective value TON

            cat_mol = exp_condition.cat_load;
            cat_eq = cat_mol*0.01;
            cat_conc = this.CA_0*cat_eq;

            TON = (this.CA_0*yield*0.01)./cat_conc;

            objective_value = TON;

        end

        function boundary_update = Boundary_update(this, objective_value)

            % Identifies the limit of variables that influence the
            % objective calculation to exclude futile region in the search
            % space

            % variable values favoring objective
            cat_load = min(this.cat_load_bound);

            %  catalyst loading boundary limit
            cat_load_limit = (this.CA_0*100*0.01)/(objective_value*this.CA_0*0.01);

            % boundary update

            boundary_update.cat_load = cat_load_limit;

        end

    end

end


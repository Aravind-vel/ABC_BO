classdef Buchwald_Hardwig_coupling

    % Class file to calculate the objective and other parameters needed for
    % running ABC-BO. This is file is specific for each reaction
    % optimization problem. 


    properties
        
        electrophile_options = {'PhBr','PhCl'};
        base_options = {'tBuKO','Cs2CO3'};
        ligand_options = {'Johnphos','Xantphos','PPh3'};

        time_bound = [0.5,8]; % h
        catalyst_bound = [0.5,5]; % mol%
        elec_stoic_bound = [1,2]; % eq.

        volume_solvent
        lim_reagent_ini
        base_ini
        ligand_conc_ini
        internal_std_conc

        molecular_weight
        density
        cost

    end

    methods
        function this = Buchwald_Hardwig_coupling
            % Reaction details



            %% chemistry details

            this.volume_solvent = 2*0.001; % l
            this.lim_reagent_ini = 500*0.001*0.001; %  mol
            this.base_ini = 3; % eq
            this.ligand_conc_ini = 4; % x cat. loading
            this.internal_std_conc = 0.25; % M ANISOLE

            this.molecular_weight = []; % g/mol
            this.molecular_weight.lim = 345.79;
            this.molecular_weight.PhBr = 157.01;
            this.molecular_weight.PhCl = 112.56;
            this.molecular_weight.Pd_catalyst = 915.73;
            this.molecular_weight.Johnphos = 298.41;
            this.molecular_weight.Xantphos = 578.62;
            this.molecular_weight.PPh3 = 262.292;
            this.molecular_weight.tBuKO = 112.21;
            this.molecular_weight.Cs2CO3 = 325.82;
            this.molecular_weight.internal_std = 108.14; %

            this.molecular_weight.product = 385.43;

            this.density = []; % g/l
            this.density.PhBr = 1.5*1000;
            this.density.PhCl = 1.11*1000;


            %% cost of reagents

            this.cost = [];
            this.cost.lim_reagent = 564/(50*0.001); % euro/g

            this.cost.PhBr = 134; % euro/l
            this.cost.PhCl = 469/25; % euro/l

            this.cost.Pd_catalyst = 1850/50; % euro/g

            this.cost.Johnphos = 287/25; % euro/g
            this.cost.Xantphos = 556/25; % euro/g
            this.cost.PPh3 = 1500/(25*1000); % euro/g

            this.cost.tBuKO = 462/(2.5*1000); % euro/g
            this.cost.Cs2CO3 = 268/500; % euro/g

            this.cost.toluene = 635/20; % euro/l




        end

        function [objective_value,productivity,cost_spent] = Objective_value_calculation(this, exp_condition, yield)

            % Return objective value [productivity/cost] g h-1 euro-1
            % Return productivity g h-1
            % Return cost_spent euro - [struct] - total cost, cost ligand,
            % cost electrophile, cost base, cost catalyst


            ligand = exp_condition.ligand;
            base = exp_condition.base;
            electrophile = exp_condition.electrophile;

            reaction_time = exp_condition.time; %hr

            %% amount of reagent used
            used = [];
            % amount of limiting reagent used
            used.lim_reagent = this.lim_reagent_ini*this.molecular_weight.lim; % g

            % amount of electrophile used
            used.electrophile = this.lim_reagent_ini*exp_condition.elec_stoic*this.molecular_weight.(electrophile); % g

            % amount of base used
            used.base = this.lim_reagent_ini*this.base_ini*this.molecular_weight.(base); %g

            % amount of catalyst used
            mol_cat_used = this.lim_reagent_ini*exp_condition.catalyst_conc*0.01;
            used.catalyst = mol_cat_used*this.molecular_weight.Pd_catalyst; %g

            % amount of ligand used
            mol_ligand_used = mol_cat_used*this.ligand_conc_ini;
            used.ligand = mol_ligand_used*this.molecular_weight.(ligand); %g

            %% cost of each reagents

            % cost of limiting reagent
            cost_spent.lim_reagent = used.lim_reagent * this.cost.lim_reagent; % euro

            % cost of electrophile
            cost_spent.electrophile = (used.electrophile * this.cost.(electrophile))/this.density.(electrophile); % euro

            % cost of base
            cost_spent.base = used.base*this.cost.(base); % euro

            %cost of catalyst
            cost_spent.catalyst = used.catalyst*this.cost.Pd_catalyst; % euro

            % cost of ligand
            cost_spent.ligand = used.ligand*this.cost.(ligand);  % euro

            % cost of solvent
            cost_spent.solvent = this.volume_solvent*this.cost.toluene;

            % cost_fixed = cost_spent.lim_reagent +cost_spent.solvent;
            cost_fixed = 0;

            cost_variable = cost_spent.electrophile + cost_spent.base +cost_spent.catalyst+...
                cost_spent.ligand;

            cost_spent.total = cost_fixed + cost_variable; % euro

            %% productivity

            product_formed = this.lim_reagent_ini *yield*0.01; % moles
            product_formed = product_formed * this.molecular_weight.product; %g

            productivity = product_formed/reaction_time; % g/h

            %% productivity/cost

            objective_value = productivity/cost_spent.total;

        end

        function boundary_update = Boundary_update(this, objective_value)

            % Identifies the limit of variables that influence the
            % objective calculation to exclude futile region in the search
            % space

            % exp condition associated with favoring calculation of objective
            exp_condition = struct();
            
            % variable values favoring objective 
            exp_condition.time = min(this.time_bound);
            exp_condition.catalyst_conc = min(this.catalyst_bound);
            exp_condition.elec_stoic = min(this.elec_stoic_bound);

            % cost associated with each combination of reagents
            % identifying reagents favoring objective
            n_combo = numel(this.ligand_options)*numel(this.electrophile_options)*numel(this.base_options);
            cost_combo = zeros(n_combo,1);
            obj_combo = zeros(n_combo,1);
            reag_combination = cell(n_combo,3);
            n = 1;
            for l = 1:numel(this.ligand_options)                
                for e = 1:numel(this.electrophile_options)
                    for b = 1:numel(this.base_options)
                        
                        % fixing reagent
                        exp_condition.ligand = char(this.ligand_options{l});
                        exp_condition.base = char(this.base_options{b});
                        exp_condition.electrophile = char(this.electrophile_options{e});

                        [reag_combination{n, 1}, reag_combination{n, 2}, reag_combination{n, 3}] = deal(exp_condition.ligand, exp_condition.base, exp_condition.electrophile);

                        [obj_combo(n),~,cost_all] = this.Objective_value_calculation(exp_condition,100);
                        cost_combo(n) = cost_all.total;
                        n = n+1;
                                             
                    end
                end
            end

            % minimum cost
            [min_cost,pos] = min(cost_combo,[],"all");

            % reagent associated with minimum cost
            exp_condition.ligand = reag_combination{pos,1};
            exp_condition.base = reag_combination{pos,2};
            exp_condition.electrophile = reag_combination{pos,3};

            % cost of individual reagent associated with minimum consumption
            [~,~,cost_all] = this.Objective_value_calculation(exp_condition, 100);

            % product formed (maximum) at 100% yield
            product_formed = this.lim_reagent_ini *100*0.01; % moles
            product_formed = product_formed * this.molecular_weight.product; %g

            %% time boundary limit
            time_limit = product_formed/((objective_value)*min_cost);
            %% catalyst boundary limit
            catalyst_conc_limit = ((product_formed/(objective_value*exp_condition.time))...
                - cost_all.base-cost_all.electrophile-cost_all.ligand)/(this.cost.Pd_catalyst*this.molecular_weight.Pd_catalyst*this.lim_reagent_ini*0.01);
            %% electrophile stoic limit
            elec_stoic_limit = ((product_formed/(objective_value*exp_condition.time))...
                -cost_all.base-cost_all.catalyst-cost_all.ligand)*this.density.(exp_condition.electrophile)/(this.cost.(exp_condition.electrophile)* this.lim_reagent_ini*this.molecular_weight.(exp_condition.electrophile));

            %% categorical variable levels that are futile

            % ligand 
            ligand_futile = cell(0,0);
            for l = 1:numel(this.ligand_options)
                % all combo of ligand
                row_numbers = find(strcmp(reag_combination(:,1),this.ligand_options{l}));
                % obj values corresponding to all combo of this ligand
                if all(obj_combo(row_numbers,:)<objective_value)
                    % mark ligand as futile
                    ligand_futile{end+1} = this.ligand_options{l};
                end
            end


            % base
            base_futile = cell(0,0);
            for b = 1:numel(this.base_options)
                % all combo of base
                row_numbers = find(strcmp(reag_combination(:,2),this.base_options{b}));
                % obj values corresponding to all combo of this base
                if all(obj_combo(row_numbers,:)<objective_value)
                    % mark base as futile
                    base_futile{end+1} = this.base_options{b};
                end
            end

            % electrophile
            electrophile_futile = cell(0,0);
            for e = 1:numel(this.electrophile_options)
                % all combo of base
                row_numbers = find(strcmp(reag_combination(:,3),this.electrophile_options{e}));
                % obj values corresponding to all combo of this electrophile
                if all(obj_combo(row_numbers,:)<objective_value)
                    % mark electrophile as futile
                    electrophile_futile{end+1} = this.electrophile_options{e};
                end
            end


           %% boundary update
           boundary_update.time = time_limit;
           boundary_update.catalyst_conc = catalyst_conc_limit;
           boundary_update.elec_stoic = elec_stoic_limit;
           boundary_update.electrophile = electrophile_futile;
           boundary_update.ligand = ligand_futile;
           boundary_update.base = base_futile;
            
        end
        


    end



end


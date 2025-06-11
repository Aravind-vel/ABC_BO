classdef ISR_reaction_class

    % Class file to calculate the objective and other parameters needed for
    % running ABC-BO. This is file is specific for each reaction
    % optimization problem. 

    properties

        time_bound = [1,10]; % min
        reag_eq_bound = [1,2]; % eq
        cat_load_bound = [5,20]; % mol %

        inj_volume = 214; %micro l
        lim_stock = 0.5; %M
        reag_stock = 0.5; %M excess reagent
        cat_stock = 0.05; %M
        reactor_vol = 5; %ml
        mol_weight_product = 246;

    end

    methods
        function this = ISR_reaction_class

        end

        function [objective_value] = Objective_value_calculation(this, exp_condition, yield)

            % Return objective value [Throughput] g h-1

            cat_mol = exp_condition.cat_load;
            cat_eq = cat_mol*0.01;
            reag_eq = exp_condition.reag_eq;
            time = exp_condition.time;

            % volume of limiting in a fixed total volume
            volume_lim = (this.reag_stock.*this.inj_volume)./(this.reag_stock+this.lim_stock.*(((cat_eq.*this.reag_stock)+(reag_eq.*this.cat_stock))./(this.cat_stock)));
            % final concentration of limiting in total volume
            conc_lim_aftermix = this.lim_stock .* volume_lim ./this.inj_volume;

            throughput = yield .* 0.01 .* conc_lim_aftermix .* this.reactor_vol./(time); %mol/min
            throughput = throughput .* this.mol_weight_product .* 0.001 / 0.01667; % g/h

            objective_value = throughput;

        end

        function boundary_update = Boundary_update(this, objective_value)
            
            % Class file to calculate the objective and other parameters needed for
            % running ABC-BO. This is file is specific for each reaction
            % optimization problem.


            % variable values favoring objective
            time = min(this.time_bound);
            cat_load = min(this.cat_load_bound);
            cat_eq = cat_load*0.01;
            reag_eq = min(this.reag_eq_bound);

            %% time boundary limit

            % volume of limiting in a fixed total volume
            volume_lim = (this.reag_stock*this.inj_volume)/(this.reag_stock+this.lim_stock*(((cat_eq*this.reag_stock)+(reag_eq*this.cat_stock))/(this.cat_stock)));
            % final concentration of limiting in total volume
            conc_lim_aftermix = this.lim_stock * volume_lim /this.inj_volume;

            throughput = (objective_value*0.01667)/(this.mol_weight_product*0.001); % mol/min
            time_limit = 100 * 0.01 * conc_lim_aftermix * this.reactor_vol/(throughput);
            %% reagent eq and catalyst loading boundary limit
            throughput = (objective_value*0.01667)/(this.mol_weight_product*0.001); % mol/min
            conc_lim_aftermix = (throughput*time)/(100*0.01*this.reactor_vol);
            volume_lim = (conc_lim_aftermix*this.inj_volume)/this.lim_stock;

            reag_eq_limit =  [(((this.reag_stock*this.inj_volume)/(volume_lim*this.lim_stock))-(this.reag_stock/this.lim_stock))]-[(cat_eq*this.reag_stock)/this.cat_stock];
            cat_eq_limit = [(((this.reag_stock*this.inj_volume)/(volume_lim*this.lim_stock))-(this.reag_stock/this.lim_stock))]*(this.cat_stock/this.reag_stock)-[(reag_eq*this.cat_stock)/this.reag_stock];
            cat_load_limit = cat_eq_limit*100;

            %% boundary update
            boundary_update.time = time_limit;
            boundary_update.reag_eq = reag_eq_limit;
            boundary_update.cat_load = cat_load_limit;

        end



    end



end


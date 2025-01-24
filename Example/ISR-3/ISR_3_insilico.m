function TON =ISR_3_insilico(Res_time,Temp,Cat_con,Cat)

% case_study = 1;  choose 1 for reizman 1, 2 for reizman 2, and so on. 
case_study = 1;

%%%% Editable parameters %%%%%%
CA_0 = 0.167; % M
CB_0 = 0.250; % M


%initial concentration
CA = CA_0;
CB = CB_0;
CR = 0;

[CA,CB,CR] = Concentration(Res_time,Temp,Cat_con,Cat,case_study,CA,CB,CR);

yield = (CR/CA_0) *100;

Cat_con = Cat_con*0.01*CA_0; % mol% of limiting reagent
TON = CR/(Cat_con);

end


function [CA,CB,CR] = Concentration(Res_time,Temp,Cat_con,Cat,case_study,CA,CB,CR)

%%%%% Editable parameters%%%%%
A_R = 3.1*10^7; % L^0.5 mol^-1.5 s^-1
EA_R = 55; % KJ mol^-1
R = 8.314 * 10^-3; %KJ/mol.K
dt = 0.001; % Define time step for Eucler integration


% default for all cases except case 3 and 4
A_s_1 = 0; % s^-1
EA_s_1 = 0; % KJ mol^-1
A_s_2 = 0; % L^0.5 mol^-1.5 s^-1
EA_s_2 = 0; % KJ mol^-1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Identifying parameters applicable for the current case study and the
% catalyst
switch case_study
    case 3
        A_s_1 = 1 * 10^12; % s^-1
        EA_s_1 = 100; %KJ mol^-1
    case 4
        A_s_2 = 3.1 * 10^5; % L^0.5 mol^-1.5 s^-1
        EA_s_2 = 50; %KJ mol^-1
end
EA_i_table = [0   0   0  0  -5.0
            0.3 0   0.3  0.3 0.7
            0.3 0.3 0.3  0.3 0.7
            0.7 0.7 0.7  0.7 0.7
            0.7 0.7 0.7  0.7 0.7
            2.2 2.2 2.2  2.2 2.2 
            3.8 3.8 3.8  3.8 3.8
            7.3 7.3 7.3  7.3 7.3];

EA_i = EA_i_table(Cat,case_study);

if Cat ==1 && case_study == 5 && Temp > 80
    EA_i = -5.0 + 0.3 * (Temp-80);
end

% converting temp to agree with other units
Temp = Temp+273 ; % K
%converting time to agree with othe units
Res_time = Res_time*60 ; % sec
% catalyst mol % to M
Cat_con = Cat_con*0.01*CA; % mol% of limiting reagent


    % Calculate kR, ks1, ks2 using the Arrhenius equation
    k_R = (Cat_con^(1/2)) * A_R * exp((-1 * (EA_R + EA_i)) / (R * Temp));
    k_s_1 = A_s_1 * exp((-1 * EA_s_1 )/ (R * Temp));
    k_s_2 = A_s_2 * exp((-1 * EA_s_2 )/ (R * Temp));



    for t = 0:dt:Res_time

        r_A = -k_R * CA *CB;
        r_B = -(k_R * CA *CB) - (k_s_1*CB) - (k_s_2*CR*CB);
        r_C = (k_R * CA *CB) - (k_s_2*CB*CR);

        CA = CA + r_A *dt;
        CB = CB + r_B *dt;
        CR = CR + r_C *dt;


    end

 
end

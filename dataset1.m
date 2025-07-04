% Dataset one 
clc;clear;close all;

% Constants
rho = 996.69;                                               % kg/m^3 (water at 80°F)
g = 9.81;                                                   % m/s^2 
Voltage = 460;                                              % Voltage (3-phase)
pf = 0.875;                                                 % Power factor
eta_motor = 0.90;                                           % Motor efficiency
eta_mech = 1;                                               % Mechanichal efficiency
% Elevation 1ft = 0.3048 m 
ft_to_m = 0.3048;
z_s = 1 * ft_to_m;                                          % Suction hight
z_d = 3 * ft_to_m;                                          % Discharge hight
d_pipe = 6 * 0.0254;                                        % Pipe diameter (inches to meters)
A_pipe = pi * (d_pipe/2)^2;                                 % Pipe cross-section area

% Data set 1
Q_gpm = [0 500 800 1000 1100 1200 1400 1500];               % Flow in GPM
P_s_psi = [0.65 0.25 -0.35 -0.92 -1.24 -1.62 -2.42 -2.89];  % suction
P_d_psi = [53.3 48.3 42.3 36.9 33.0 27.8 15.3 7.3];         % discharge
I = [18.0 26.2 31.0 33.9 35.2 36.3 38.0 39.0];              % motor current

% Convert pressure to Pa
psi_to_Pa = 6894.76;
P_s = P_s_psi * psi_to_Pa;
P_d = P_d_psi * psi_to_Pa;

% Convert GPM to m^3/s
Q = Q_gpm * 0.00006309;

% Velocity
V_s = Q ./ A_pipe;
V_d = Q ./ A_pipe;

% Calculating Head of the pump
h_pump = (P_d - P_s) ./ (rho * g) + (V_d.^2 - V_s.^2) ./ (2 * g) + (z_d - z_s);

% Hydraulic Power
P_hyd = rho * g .* Q .* h_pump;

% Motor input power
P_motor = sqrt(3) * Voltage * I * pf* eta_motor;  % in watts

% Mechanical Power
P_mech = eta_mech  * P_motor ;

% Overall Efficiency
eta_overall = P_hyd ./ P_motor;

% Hydrolic Efficiency
eta_hyd = P_hyd ./ P_mech ;

% Specific Speed (US Customary Units)
N = 1750;                                    % rpm
H_us = h_pump / 0.3048;                      % meters to feet
Ns_US = N * sqrt(Q_gpm) ./ (H_us.^(3/4));

% Specific Speed in SI (dimensionless)
omega = 2 * pi * N / 60;                     % rad/s
Ns_SI = N * sqrt(Q) ./ (h_pump.^(3/4));  % dimensionless


% Curve Fitting Section
% 4rth-degree polynomials
p_head = polyfit(Q_gpm, h_pump, 4);
p_power = polyfit(Q_gpm, P_hyd/1000, 4);                 % in kW
p_eta_overall = polyfit(Q_gpm, eta_overall * 100, 4);    % in %
p_eta_hyd = polyfit(Q_gpm, eta_hyd *100, 4);
p_Ns_US = polyfit(Q_gpm, Ns_US, 4);
p_Ns_SI = polyfit(Q_gpm,Ns_SI, 4);

% Create smooth Q range for plotting
Q_smooth = linspace(min(Q_gpm), max(Q_gpm), 200);

% Evaluate fitted curves
head_fit = polyval(p_head, Q_smooth);
power_fit = polyval(p_power, Q_smooth);
eta_overall_fit = polyval(p_eta_overall, Q_smooth);
eta_hyd_fit = polyval(p_eta_hyd, Q_smooth);
Ns_US_fit = polyval(p_Ns_US, Q_smooth);
Ns_SI_fit = polyval(p_Ns_SI, Q_smooth);

% Plotting the results
figure("Name","Rsults of Dataset 1",'Units','normalized','Position', [0.1 0.1 0.8 0.8])

subplot(2,3,1)
plot(Q_gpm, h_pump, '-o')
hold on
plot(Q_gpm, h_pump, 'ro', Q_smooth, head_fit, 'b-')
legend('Data','', 'Fitted Curve','Location','best')
xlabel('Flow Rate (GPM)')
ylabel('Pump Head (m)')
title('Pump Head vs Flow')
grid minor

subplot(2,3,4)
plot(Q_gpm, P_hyd/1000, '-o')
hold on 
plot(Q_gpm, P_hyd/1000, 'ro', Q_smooth, power_fit, 'b-')
legend('Data','','Fitted Curve','Location','best')
xlabel('Flow Rate (GPM)')
ylabel('Hydraulic Power (kW)')
title('Hydraulic Power vs Flow')
grid minor

subplot(2,3,2)
plot(Q_gpm, eta_overall*100, '-o')
hold on
plot(Q_gpm, eta_overall * 100, 'ro', Q_smooth, eta_overall_fit, 'b-')
legend('Data','','Fitted Curve','Location','best')
xlabel('Flow Rate (GPM)')
ylabel('Overall Efficiency (%)')
title('Overall Efficiency vs Flow')
hold on
grid minor

subplot(2,3,5)
plot(Q_gpm, Ns_US, '-o')
hold on
plot(Q_gpm, Ns_US, 'ro', Q_smooth, Ns_US_fit, 'b-')
legend('Data','','Fitted Curve','Location','best')
xlabel('Flow Rate (GPM)')
ylabel('Specific Speed (US)')
title('Specific Speed(US) vs Flow')
grid minor

subplot(2,3,3)
plot(Q_gpm, eta_hyd*100, '-o')
hold on
plot(Q_gpm, eta_hyd * 100, 'ro', Q_smooth, eta_hyd_fit, 'b-')
legend('Data','','Fitted Curve','Location','best')
xlabel('Flow Rate (GPM)')
ylabel('Hydrolic Efficiency (%)')
title('Hydrolic Efficiency vs Flow')
grid minor

subplot(2,3,6)
plot(Q_gpm, Ns_SI, '-o')
hold on
plot(Q_gpm, Ns_SI, 'ro', Q_smooth, Ns_SI_fit, 'b-')
legend('Data','','Fitted Curve','Location','best')
xlabel('Flow Rate (GPM)')
ylabel('Specific Speed (SI)')
title('Specific Speed(SI) vs Flow')
grid minor

% display coefficients and polynominal formulas
disp("Displaying the fitted polynomial equations")
display_poly(p_head, 'Q', 'Pump Head h');
display_poly(p_power, 'Q', 'Hydraulic Power P_{hyd}');
display_poly(p_eta_overall, 'Q', 'Overall Efficiency η_{overall}');
display_poly(p_eta_hyd, 'Q', 'Hydraulic Efficiency η_{hyd}');
display_poly(p_Ns_US, 'Q', 'Specific Speed (US) N_s^{US}');
display_poly(p_Ns_SI, 'Q', 'Specific Speed (SI) N_s^{SI}');

% Save to excel
T1 = table(Q_gpm',Q' , P_s_psi', P_d_psi',P_s' ,P_d', V_s' , V_d' , I', ...
    h_pump', P_hyd'/1000, P_motor'/1000, P_mech', eta_overall'*100, eta_hyd'*100,Ns_US' , Ns_SI', ...
    'VariableNames', {'Q_GPM','Q_m^3' , 'P_s_psi', 'P_d_psi','P_s_Pa' ,'P_d_Pa','V_s','V_d' , 'Motor_Current_A', ...
    'Pump_Head_m', 'P_hyd_kW', 'P_motor_kW','P_mech_kW', ...
    'Overrall_Efficiency_%','Hydrolic_Efficiency_%' ,'SpecificSpeed_US' ,'SpecificSpeed_SI'});
filename1 = fullfile(pwd , 'results_dataset1.xlsx');
writetable(T1, filename1);

function display_poly(coeffs, varname, label)
    fprintf('%s(Q) = ', label);
    n = length(coeffs);
    for i = 1:n
        c = coeffs(i);
        p = n - i;
        if abs(c) < 1e-10
            continue;  % Skip near-zero terms
        end
        if i > 1 && c >= 0
            fprintf('+ ');
        elseif i > 1
            fprintf('- ');
            c = -c;
        end
        if p == 0
            fprintf('%.4g ', c);
        elseif p == 1
            fprintf('%.4g*%s ', c, varname);
        else
            fprintf('%.4g*%s^%d ', c, varname, p);
        end
    end
    fprintf('\n\n');
end
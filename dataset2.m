% Dataset2
clc;clear;close all;
% Constants
rho = 997.0;                 % kg/m^3 (water at 25°C)
g = 9.81;                    % m/s^2    
N_rpm = 1100;                % RPM 
omega = 2 * pi * N_rpm / 60;
d_pipe = 0.15;               % Pipe diameter(m)
A_pipe = pi * (d_pipe/2)^2;  % Pipe cross-section area
% Elevation (ignoring the effect)
z_s = 0;           
z_d = 0;

% Data
Torque = [2.1 2 1.9 1.7 1.8 1.7 1.5 1.5 1.4 1.2 1.1 1.0 0.9];
Q_Lmin = [254 228 197 163 177 155 127 129 99 75 50 27 2];     
P_s_bar = [-0.08 -0.07 -0.05 -0.04 -0.05 -0.04 -0.03 -0.03 -0.02 -0.02 -0.02 -0.01 -0.01];
P_d_bar = [0.06 0.11 0.18 0.25 0.21 0.25 0.29 0.29 0.32 0.34 0.35 0.36 0.36];
P_motor = [0.52 0.5 0.48 0.46 0.47 0.45 0.42 0.42 0.4 0.38 0.36 0.33 0.31];

% Unit conversions
Q = Q_Lmin * (1e-3) / 60;               % L/min to m^3/s
P_s = P_s_bar * 1e5;                    % bar to Pa
P_d = P_d_bar * 1e5;                    % bar to Pa

% Velocities
V_s = Q ./ A_pipe;
V_d = Q ./ A_pipe;

% Head
h_pump = (P_d - P_s) ./ (rho * g) + (V_d.^2 - V_s.^2)./(2*g) + (z_d - z_s);

% Hydraulic Power
P_hyd = rho * g .* Q .* h_pump;      

% Mechanical Power 
P_mech = Torque .* omega;            

% Efficiencies
eta_hyd = P_hyd ./ P_mech;
eta_overall = P_hyd ./ (P_motor * 1000);

% Specific Speed (SI)
Ns_SI = N_rpm * sqrt(Q) ./ (h_pump.^(3/4));
Ns_US = N_rpm * sqrt(Q / 0.00006309) ./ ((h_pump / 0.3048).^(3/4));

% Curve Fitting Section
% 4rth-degree polynomials
p_head = polyfit(Q, h_pump, 4);
p_power = polyfit(Q, P_hyd/1000, 4);        
p_eta_overall = polyfit(Q, eta_overall * 100, 4);   
p_eta_hyd = polyfit(Q, eta_hyd *100, 4);
p_P_mech = polyfit(Q, P_mech, 4);
p_Ns_SI = polyfit(Q,Ns_SI, 4);

% Create smooth Q range for plotting
Q_smooth = linspace(min(Q), max(Q), 200);

% Evaluate fitted curves
head_fit = polyval(p_head, Q_smooth);
power_fit = polyval(p_power, Q_smooth);
eta_overall_fit = polyval(p_eta_overall, Q_smooth);
eta_hyd_fit = polyval(p_eta_hyd, Q_smooth);
P_mech_fit = polyval(p_P_mech, Q_smooth);
Ns_SI_fit = polyval(p_Ns_SI, Q_smooth);

% Plotting the results
figure("Name","Rsults of Dataset 2",'Units','normalized','Position', [0.1 0.1 0.8 0.8])

subplot(2,3,1)
plot(Q, h_pump, '-o')
hold on
plot(Q, h_pump, 'ro', Q_smooth, head_fit, 'b-')
legend('Data','', 'Fitted Curve','Location','best')
xlabel('Flow Rate (m^3/s)')
ylabel('Pump Head (m)')
title('Pump Head vs Flow')
grid minor

subplot(2,3,4)
plot(Q, P_hyd/1000, '-o')
hold on 
plot(Q, P_hyd/1000, 'ro', Q_smooth, power_fit, 'b-')
legend('Data','','Fitted Curve','Location','best')
xlabel('Flow Rate (m^3/s)')
ylabel('Hydraulic Power (kW)')
title('Hydraulic Power vs Flow')
grid minor

subplot(2,3,2)
plot(Q, eta_overall*100, '-o')
hold on
plot(Q, eta_overall * 100, 'ro', Q_smooth, eta_overall_fit, 'b-')
legend('Data','','Fitted Curve','Location','best')
xlabel('Flow Rate (m^3/s)')
ylabel('Overall Efficiency (%)')
title('Overall Efficiency vs Flow')
hold on
grid minor

subplot(2,3,5)
plot(Q, P_mech, '-o')
hold on
plot(Q, P_mech, 'ro', Q_smooth, P_mech_fit, 'b-')
legend('Data','','Fitted Curve','Location','best')
xlabel('Flow Rate (m^3/s)')
ylabel('Mechanical Power (W)')
title('Mechanical Power vs Flow')
grid minor

subplot(2,3,3)
plot(Q, eta_hyd*100, '-o')
hold on
plot(Q, eta_hyd * 100, 'ro', Q_smooth, eta_hyd_fit, 'b-')
legend('Data','','Fitted Curve','Location','best')
xlabel('Flow Rate (m^3/s)')
ylabel('Hydrolic Efficiency (%)')
title('Hydrolic Efficiency vs Flow')
grid minor

subplot(2,3,6)
plot(Q, Ns_SI, '-o')
hold on
plot(Q, Ns_SI, 'ro', Q_smooth, Ns_SI_fit, 'b-')
legend('Data','','Fitted Curve','Location','best')
xlabel('Flow Rate (m^3/s)')
ylabel('Specific Speed (SI)')
title('Specific Speed(SI) vs Flow')
grid minor

% display coefficients and polynominal formulas
disp("Displaying the fitted polynomial equations")
display_poly(p_head, 'Q', 'Pump Head h');
display_poly(p_power, 'Q', 'Hydraulic Power P_{hyd}');
display_poly(p_eta_overall, 'Q', 'Overall Efficiency η_{overall}');
display_poly(p_eta_hyd, 'Q', 'Hydraulic Efficiency η_{hyd}');
display_poly(p_P_mech, 'Q', 'Mechanical Power P_{US}');
display_poly(p_Ns_SI, 'Q', 'Specific Speed (SI) N_s^{SI}');


% Save to excel
T1 = table(Q_Lmin',Q' , P_s_bar', P_d_bar',P_s' ,P_d', V_s' , V_d' , Torque', ...
    h_pump', P_hyd'/1000, P_motor', P_mech'/1000, eta_overall'*100, eta_hyd'*100,Ns_US' , Ns_SI', ...
    'VariableNames', {'Q_Lmin','Q_m^3' , 'P_s_bar', 'P_d_bar','P_s_Pa' ,'P_d_Pa','V_s','V_d' , 'Torque_Nm', ...
    'Pump_Head_m', 'P_hyd_kW', 'P_motor_kW','P_mech_kW', ...
    'Overrall_Efficiency_%','Hydrolic_Efficiency_%' ,'SpecificSpeed_US' ,'SpecificSpeed_SI'});
filename1 = fullfile(pwd , 'results_dataset2.xlsx');
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
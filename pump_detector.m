clc;
clear;
close all;

%% 1. Importing digitlized pump performance data (Q-H curves) into matlab from excels

figure('Name', 'PERFORMANCE RANGE OF CENTRIFUGAL PUMPS','Units','normalized','Position', [0.1 0.1 0.8 0.8]);
hold on

% i and j for looping the pumps data and storing in cell
i_list = [32, 40, 50, 65, 80, 100, 125];
j_list = [125, 160, 200, 250];

% empty cell for data from exel
pump_data_cell = cell(length(i_list), length(j_list));

for i_index = 1:length(i_list)
    for j_index = 1:length(j_list)
        i_value = i_list(i_index);
        j_value = j_list(j_index);
        % files of data extracted from digitizer app
        filename = sprintf('%d-%d.xls', i_value, j_value);
        
        if isfile(filename)
            data = readmatrix(filename);
            
            Q = data(:, 2);  % Flow
            H = data(:, 3);  % Head
            
            % Store in cell array as a structure with fields of i and j and
            % Q and H date
            pump_data_cell{i_index, j_index} = struct('i', i_value, 'j', j_value, 'Q', Q, 'H', H);
            
            % Plot the performance curve of all the pumps
            plot(Q, H, 'DisplayName', sprintf('%d-%d', i_value, j_value));
        else
            continue
        end
    end
end

xlabel('Flow rate Q (m^3/h)');
ylabel('Head H (m)');
title('Pump Performance Curves');
legend show
grid minor
hold on

%% 2. Getting requirments from user and selecting the suitable pump

Q_user = input(' Please enter the required flow rate (m^3/h): ');
H_user = input(' Please enter the required head (m): ');
scatter(Q_user, H_user, 60, 'b', 'filled');  % Show user point
disp(" Starting the analysis ...")

% checking which pump can satisfy the requirments of user

for i_index = 1:length(i_list)
    for j_index = 1:length(j_list)

        pump = pump_data_cell{i_index, j_index};
        if isempty(pump)
            continue;  % Skip empty cells
        end

        % Extract data 
        Q = pump.Q;
        H = pump.H;

        % Check point inside polygon
        [in, on] = inpolygon(Q_user, H_user, Q, H);

        if in || on
            i_value = i_list(i_index);
            j_value = j_list(j_index);
            valid_pump = sprintf('%d-%d', i_value, j_value);
            fprintf(" Pump %s can achieve Q = %.1f m³/h and H = %.2f m\n", valid_pump, Q_user, H_user);
            break
        end
    end
end


%% 3. Determining the diameter of the impeller of the pump

figure('Name', 'PERFORMANCE-DIAMETER RANGE OF CENTRIFUGAL PUMPS','Units','normalized','Position', [0.1 0.1 0.8 0.8])
hold on

% getting the existing diameter data from excel files
filename = sprintf('d-%s.xlsx',valid_pump);
if isfile(filename)
    data_d = readmatrix(filename);
    d_list = data_d(:,1);
    non_nan_d = ~isnan(d_list);
    d_list = d_list(non_nan_d);
else
    error(" File '%s' not found.", filename);
end

% empty cell array for diameter datas
diameter_data_cell = cell(length(d_list));

for d_index = 1:length(d_list)
    d_value = d_list(d_index); 
    filename = sprintf('%d-%d-%d.xls', i_value, j_value,d_value);   
    if isfile(filename)
        data = readmatrix(filename);
        diameter = d_value ; % pump diameter
        Q = data(:, 2);  % Flow
        H = data(:, 3);  % Head
        diameter_data_cell{d_index} = struct('i', i_value, 'j', j_value, 'd', d_value, 'Q', Q, 'H', H);
        plot(Q, H, 'DisplayName', sprintf('%d-%d-%d', i_value, j_value, d_value));
    else
        continue
    end
end

xlabel('Flow rate Q (m^3/h)');
ylabel('Head H (m)');
title('Pump Performance-Diameter Curves');
legend show
grid minor
hold on
scatter(Q_user, H_user, 60, 'b', 'filled');  % Show user point


% determining the diameter using 1d interploation 
H_interp_list = [];
d_list_valid = [];

for d_index = 1:length(d_list)
    data = diameter_data_cell{d_index};
    if isempty(data)
        continue;
    end

    Q = data.Q;
    H = data.H;
    % clearing nan data
    non_nan = ~isnan(Q) & ~isnan(H);
    Q = Q(non_nan);
    H = H(non_nan);
    % sorting the data
    [Q_sorted, idx] = sort(Q);
    H_sorted = H(idx);

    % Interpolate head at Q_user for this diameter
    H_interp = interp1(Q_sorted, H_sorted, Q_user);
    
    % Store for interpolation
    H_interp_list(end+1) = H_interp;
    d_list_valid(end+1) = data.d;
end

% Remove any NaN or Inf values before interpolation
finite_mask = isfinite(H_interp_list) & isfinite(d_list_valid);
H_interp_list = H_interp_list(finite_mask);
d_list_valid = d_list_valid(finite_mask);

if length(H_interp_list) < 2
    error('Not enough valid data to interpolate impeller diameter.');
end


if length(H_interp_list) < 2
    disp(" Not enough diameter curves at this Q to interpolate.");
else
    % D = f(H) → interpolate to find diameter for H_user
    D_interp = interp1(H_interp_list, d_list_valid, H_user);

    fprintf(' Required diameter to achieve Q = %.2f and H = %.2f is: %.2f mm\n', Q_user, H_user, D_interp);
end
hold on

% Plotting the interpolated diameter curve

Q_interp = linspace(min(Q_user, 0), max(Q_user*1.5, 30), 200);
H_interp_curve = zeros(size(Q_interp));

for k = 1:length(Q_interp)
    q = Q_interp(k);

    H_at_q = [];
    D_at_q = [];

    for d_index = 1:length(d_list)
        data = diameter_data_cell{d_index};
        if isempty(data)
            continue;
        end

        Q = data.Q;
        H = data.H;

        valid = ~isnan(Q) & ~isnan(H);
        Q = Q(valid);
        H = H(valid);
        [Q_sorted, idx] = sort(Q);
        H_sorted = H(idx);

        if q < min(Q_sorted) || q > max(Q_sorted)
            continue;
        end

        hq = interp1(Q_sorted, H_sorted, q);
        H_at_q(end+1) = hq;
        D_at_q(end+1) = data.d;
    end

    % Interpolate in D direction to find H at D_interp
    if length(D_at_q) >= 2
        H_interp_curve(k) = interp1(D_at_q, H_at_q, D_interp, 'linear', 'extrap');
    else
        H_interp_curve(k) = NaN;
    end
end

% Plot the interpolated diameter curve
plot(Q_interp, H_interp_curve, 'k', 'LineWidth', 2.5, 'DisplayName', sprintf('Interp D = %.1f mm', D_interp));
scatter(Q_user, H_user, 60, 'b', 'filled');  % Show user point
plot([Q_user Q_user], [0 H_user], 'b--', 'LineWidth', 1.5);
plot([0 Q_user], [H_user H_user], 'b--', 'LineWidth', 1.5);


%% 4. Determining the Electrical Power needed for the motor of the pump

figure('Name', 'PERFORMANCE-POWER RANGE OF CENTRIFUGAL PUMPS','Units','normalized','Position', [0.1 0.1 0.8 0.8])
hold on

% empty cell array for Power datas
Power_data_cell = cell(length(d_list));
% getting the existing diameter-Power data from excel files
for d_index = 1:length(d_list)
    d_value = d_list(d_index); 
    filename = sprintf('%d-%d---%d.xls', i_value, j_value,d_value);   
    if isfile(filename)
        data = readmatrix(filename);
        diameter = d_value ; % pump diameter
        Q = data(:, 2);  % Flow
        P = data(:, 3);  % Power
        Power_data_cell{d_index} = struct('i', i_value, 'j', j_value, 'd', d_value, 'Q', Q, 'P', P);
        plot(Q, P, 'DisplayName', sprintf('%d-%d-%d', i_value, j_value, d_value));
    else
        continue
    end
end

xlabel('Flow rate Q (m^3/h)');
ylabel('Power  (kW)');
title('Pump Performance Curves');
legend show
grid minor
hold on


% determining the power using 1d interploation 
% Interpolate power at Q_user for each diameter
P_interp_list = [];
D_list_valid = [];

for d_index = 1:length(d_list)
    data = Power_data_cell{d_index};
    if isempty(data)
        continue;
    end

    Q = data.Q;
    P = data.P;

    non_nan = ~isnan(Q) & ~isnan(P);
    Q = Q(non_nan);
    P = P(non_nan);
    [Q_sorted, idx] = sort(Q);
    P_sorted = P(idx);

    if Q_user < min(Q_sorted) || Q_user > max(Q_sorted)
        continue;
    end

    % Interpolate power at Q_user
    P_interp = interp1(Q_sorted, P_sorted, Q_user);
    
    P_interp_list(end+1) = P_interp;
    D_list_valid(end+1) = data.d;
end

% Check if we can interpolate power across diameters
if length(P_interp_list) < 2
    disp(" Not enough valid diameter-power data at this flow rate.");
else
    % Interpolate across diameters to find power for desired H_user
    P_user = interp1(D_list_valid, P_interp_list, D_interp);  % D_interp is interpolated diameter
    fprintf(" Estimated shaft power at Q = %.2f, H = %.2f, D ≈ %.2f mm is: %.2f kW\n", ...
        Q_user, H_user, D_interp, P_user);
end

% Plotting the interpolated Power curve

Q_interp = linspace(min(Q_user, 0), max(Q_user*1.5, 30), 200);
P_interp_curve = zeros(size(Q_interp));

for k = 1:length(Q_interp)
    q = Q_interp(k);

    P_at_q = [];
    D_at_q = [];

    for d_index = 1:length(d_list)
        data = Power_data_cell{d_index}; 
        if isempty(data)
            continue;
        end

        Q = data.Q;
        P = data.P;

        valid = ~isnan(Q) & ~isnan(P);
        Q = Q(valid);
        P = P(valid);
        [Q_sorted, idx] = sort(Q);
        P_sorted = P(idx);

        if q < min(Q_sorted) || q > max(Q_sorted)
            continue;
        end

        pq = interp1(Q_sorted, P_sorted, q);
        P_at_q(end+1) = pq;
        D_at_q(end+1) = data.d;
    end

    if length(D_at_q) >= 2
        P_interp_curve(k) = interp1(D_at_q, P_at_q, D_interp, 'linear', 'extrap');
    else
        P_interp_curve(k) = NaN;
    end
end


% Plot the interpolated diameter curve
plot(Q_interp, P_interp_curve, 'k', 'LineWidth', 2.5, 'DisplayName', sprintf('Interp D = %.1f mm', D_interp));
scatter(Q_user, P_user, 60, 'b', 'filled');  % Show user point
plot([Q_user Q_user], [0 P_user], 'b--', 'LineWidth', 1.5);
plot([0 Q_user], [P_user P_user], 'b--', 'LineWidth', 1.5);


%% 5. Determining the pump efficiency at the point of work

figure('Name', 'PERFORMANCE-EFFICIENCY RANGE OF CENTRIFUGAL PUMPS','Units','normalized','Position', [0.1 0.1 0.8 0.8])
hold on

% getting the existing efficiency data from excel files
filename = sprintf('e-%s.xlsx',valid_pump);
if isfile(filename)
    data_e = readmatrix(filename);
    eta_list = data_e(:,1);
    non_nan_e = ~isnan(eta_list);
    eta_list = eta_list(non_nan_e);
else
    error(" File '%s' not found.", filename);
end


eta_data_cell = cell(length(eta_list));

for eta_index = 1:length(eta_list)
    eta_value = eta_list(eta_index); 
    filename = sprintf('%d-%d--%d.xls', i_value, j_value,eta_value);   
    if isfile(filename)
        data = readmatrix(filename);
        eta = eta_value ; % pump efficiency
        Q = data(:, 2);  % Flow
        H = data(:, 3);  % Head
        eta_data_cell{eta_index} = struct('i', i_value, 'j', j_value, 'eta', eta_value, 'Q', Q, 'H', H);
        plot(Q, H, 'DisplayName', sprintf('%d-%d--%d', i_value, j_value, eta_value));
    else
        continue
    end
end

xlabel('Flow rate Q (m^3/h)');
ylabel('Head H (m)');
title('Pump Performance-Efficiency Curves');
legend show
grid minor
hold on
scatter(Q_user,H_user,60,'b','filled')


% Gather all (Q, H, eta) triplets from curves
Q_all = [];
H_all = [];
eta_all = [];

for eta_index = 1:length(eta_list)
    data = eta_data_cell{eta_index};
    if isempty(data)
        continue;
    end

    Q = data.Q;
    H = data.H;
    eta = eta_list(eta_index);

    valid = ~isnan(Q) & ~isnan(H);
    Q_all = [Q_all; Q(valid)];
    H_all = [H_all; H(valid)];
    eta_all = [eta_all; repmat(eta, sum(valid), 1)];
end

% Use scatteredInterpolant
F = scatteredInterpolant(Q_all, H_all, eta_all, 'linear', 'none');

% Interpolate at user's point
eta_interp = F(Q_user, H_user);

if isnan(eta_interp)
    fprintf(" Q = %.2f, H = %.2f is outside the efficiency map.\n", Q_user, H_user);
else
    fprintf(" Estimated efficiency at Q = %.2f m³/h and H = %.2f m is: %1.1f %%\n", ...
        Q_user, H_user, eta_interp);
    disp(" Anyalysis successfull.")
end

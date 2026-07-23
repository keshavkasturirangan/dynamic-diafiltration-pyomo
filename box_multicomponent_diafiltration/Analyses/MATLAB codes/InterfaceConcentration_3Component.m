clear, clc, clf, close all

% From 0 to the boundary layer thickness.
xspan = linspace(0,1E-4,100); % [=] m

FilePath = 'C:\Users\lauri\OneDrive\Documents\Multicomponent Project\Diafiltration Raw Data';
ExcelFile = 'Rejection_Analysis.xlsx';
ExcelSheet = 'NF270_MC2 05.08.24_E1';
ExcelFile_Data = fullfile(FilePath, ExcelFile);

% Read in the data from the excel file.
Datapoints = readmatrix(ExcelFile_Data, 'Sheet', ExcelSheet, 'Range', 'A1:A1');
Cr_salt1 = readmatrix(ExcelFile_Data, 'Sheet', ExcelSheet, 'Range', 'B2:B11'); % mM
Cr_salt2 = readmatrix(ExcelFile_Data, 'Sheet', ExcelSheet, 'Range', 'C2:C11'); % mM
Ionic_strength = readmatrix(ExcelFile_Data, 'Sheet', ExcelSheet, 'Range', 'D2:D11'); % mM
salt_ratio = readmatrix(ExcelFile_Data, 'Sheet', ExcelSheet, 'Range', 'E2:E11');
rejection_obs_salt1 = readmatrix(ExcelFile_Data, 'Sheet', ExcelSheet, 'Range', 'F2:F11');
rejection_obs_salt2 = readmatrix(ExcelFile_Data, 'Sheet', ExcelSheet, 'Range', 'G2:G11');
Cp_salt1 = readmatrix(ExcelFile_Data, 'Sheet', ExcelSheet, 'Range', 'H2:H11'); % mM
Cp_salt2 = readmatrix(ExcelFile_Data, 'Sheet', ExcelSheet, 'Range', 'I2:I11'); % mM
Jw = readmatrix(ExcelFile_Data, 'Sheet', ExcelSheet, 'Range', 'J2:J11'); % um3/m2/s

% Initial conditions: c1(0), c2(0) of vial you are analyzing
c_r = [Cr_salt1(1); Cr_salt2(1)]; % [=] mM

% Parameters
params.DNa_s = 1.33E-9; % [=] m2/s
params.DCa_s = 0.79E-9; % [=] m2/s
params.DLa_s = 6.26E-10; % [=] m2/s
params.DCl_s = 2.03E-9; % [=] m2/s

params.Jw = Jw(1); % [=] m3/m2/s
params.cp1 = Cp_salt1(1);  % [=] mM
params.cp2 = Cp_salt2(1);  % [=] mM

% Solve ODE system for c1, c2
[x, c] = ode15s(@(x,c) ODEfunc(x, c, params), xspan, c_r);

% Recover c3 from electroneutrality
z1 = abs(1); z2 = abs(2); z3 = abs(-1);
c3 = (-z1.*c(:,1)-z2.*c(:,2))./z3;

% Plot all concentrations
hold on, axis square, grid on, box on, set(gcf,'Color','w'), set(gca,'LineWidth', 1.5), set(gca,'FontName','Arial'),set(gca,'FontSize',14), set(gca,"FontWeight","bold")
xlim([0 1E-4]), ylim([0 60])
xlabel('Boundary Layer [m]'), ylabel('Concentration [mM]')

Delta_Na = params.DNa_s/k(params.DNa_s);
Delta_Ca = params.DCa_s/k(params.DCa_s);
Delta_Cl = params.DCl_s/k(params.DCl_s);
Delta_avg = (Delta_Na+Delta_Ca+Delta_Cl)/3;

xline(Delta_avg, 'LineWidth', 1.5, 'LineStyle', '-.', Color='k');

p1 = plot(x, c(:,1), 'r', 'LineWidth', 2);
p2 = plot(x, c(:,2), 'b', 'LineWidth', 2);
p3 = plot(x, c3, 'k', 'LineWidth', 2);
legend([p1 p2 p3],'c_{Na}','c_{Ca}','c_{Cl}',location = 'northeast');
title('vial 1');

% Find and output where the boundary layer intersects the interfacial concentration profile.
c_int_salt1 = interp1(x, c(:,1), Delta_avg);
c_int_salt2 = interp1(x, c(:,2), Delta_avg);

ExportData(:,1) = c_int_salt1;
ExportData(:,2) = c_int_salt2;

function dcdx = ODEfunc(x, c, params)
    c1 = c(1); c2 = c(2);

    % Extract constants
    Jw = params.Jw;
    cp1 = params.cp1;
    cp2 = params.cp2;

    % Parameters
    chi = 0;
    D1 = 1.33E-9; % [=] m2/s solution phase diffusion coefficient for Na
    D2 = 0.793E-9; % [=] m2/s solution phase diffusion coefficient for Ca
    D3 = 2.03E-9; % [=] m2/s solution phase diffusion coefficient for Cl
    z1 = abs(1); z2 = abs(2); z3 = abs(-1); % [=] unitless

    % Denominators (with regularization)
    denom = c1*(z1^2*D1 + z1*z3*D3) + c2*(z2^2*D2 + z2*z3*D3); % flux in solution
    denom2 = denom - z3*D3*chi; % flux in the membrane

    % Compute Dij (reduced system)
    D11 = (c1*(-z1^2*D1*D3* - z1*z3*D1*D3) + c2*(-z2^2*D1*D2 - z2*z3*D1*D3)) / denom;
    D12 = (c1*(z1*z2*D1*D2 - z1*z2*D1*D3)) / denom;

    D21 = (c2*(z1*z2*D1*D2 - z1*z2*D2*D3)) / denom;
    D22 = (c1*(-z1^2*D1*D2 - z1*z3*D2*D3) + c2*(-z2*z3*D2*D3 - z2^2*D2*D3)) / denom;

    D = [D11, D12;
         D21, D22];

    % RHS vector for c1 and c2 only
    rhs = [Jw*(cp1 - c1);
           Jw*(cp2 - c2)];

    % Solve d[c1;c2]/dx = D⁻¹ * rhs
    dcdx = D \ rhs;
end

function MassTransferCoefficient = k(D_coeff) % units should be m/s
        b = 0.0254; % Units [=] m
        KinematicViscosity = 1.003*10^-6; % [=] m2/s
        AngularVelocity = 2*pi*350/60; % rotations per minute (i.e., RPM)
        v0 = AngularVelocity*b/2;
        MassTransferCoefficient = 0.23*(v0^0.57)*(D_coeff^0.67)/(KinematicViscosity^0.24)/(b^0.43);
end
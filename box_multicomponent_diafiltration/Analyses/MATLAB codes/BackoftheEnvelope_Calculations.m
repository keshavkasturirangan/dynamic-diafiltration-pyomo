clear, clc, clf, close all

% Define the experiment you are analyzing.
FilePath = 'C:\Users\lauri\OneDrive\Documents\Multicomponent Project\Diafiltration Raw Data';
ExcelFile = 'Rejection_Analysis.xlsx';
ExcelSheet = 'NF270_MC4 07.16.24_E16';
ExcelFile_Data = fullfile(FilePath, ExcelFile);

% Read in the data from the excel file.
Cr_salt1 = readmatrix(ExcelFile_Data, 'Sheet', ExcelSheet, 'Range', 'B2:B11'); % mM
Cr_salt2 = readmatrix(ExcelFile_Data, 'Sheet', ExcelSheet, 'Range', 'C2:C11'); % mM
Ionic_strength = readmatrix(ExcelFile_Data, 'Sheet', ExcelSheet, 'Range', 'D2:D11'); % mM
salt_ratio = readmatrix(ExcelFile_Data, 'Sheet', ExcelSheet, 'Range', 'E2:E11');
rejection_obs_salt1 = readmatrix(ExcelFile_Data, 'Sheet', ExcelSheet, 'Range', 'F2:F11');
rejection_obs_salt2 = readmatrix(ExcelFile_Data, 'Sheet', ExcelSheet, 'Range', 'G2:G11');
Cp_salt1 = readmatrix(ExcelFile_Data, 'Sheet', ExcelSheet, 'Range', 'H2:H11'); % mM
Cp_salt2 = readmatrix(ExcelFile_Data, 'Sheet', ExcelSheet, 'Range', 'I2:I11'); % mM
Jw = readmatrix(ExcelFile_Data, 'Sheet', ExcelSheet, 'Range', 'J2:J11'); % um3/m2/s

% Define known solution diffusion coefficients.
params.DNa_s = 1.33E-9; % [=] m2/s
params.DCa_s = 0.79E-9; % [=] m2/s
params.DLa_s = 6.26E-10; % [=] m2/s
params.DCl_s = 2.03E-9; % [=] m2/s

% Preallocate vector sizes for storing as you work through a for loop.
c1_int = zeros(length(Jw), 1);
c2_int = zeros(length(Jw), 1);

% Create a for loop, which will analyze c_int for each ion within each
% experimental vial collected (usually 10 vials).
for i = 1:length(Jw)
    c_bulk = [Cr_salt1(i); Cr_salt2(i)]; % define the bulk concentration values for ion 1 and ion 2.
    c_perm = [Cp_salt1(i); Cp_salt2(i)]; % define the permeate concentration values for ion 1 and ion 2.
    Current_Jw = Jw(i); % Define what the water flux is for the current vial that is being analyzed.

    % Calculate the Diffusion Coefficient Matrix, calling the function
    % described at the end of this code. The matrix is a function of both
    % the concentration of ion 1 and ion 2 of each vial.
    D_matrix = calculate_D(Cr_salt1(i), Cr_salt2(i));

    % Calculate the Boundary Layer (delta) for each ion using the function 
    % described at the end of this code. Then, take an average of each
    % delta_ion so a constant delta can be used for the following
    % interfacial conentration calculation.
    Delta_Na = params.DNa_s/k(params.DNa_s);
    Delta_Ca = params.DCa_s/k(params.DCa_s);
    Delta_La = params.DCa_s/k(params.DLa_s);
    Delta_Cl = params.DCl_s/k(params.DCl_s);
    Delta_avg = (Delta_Na+Delta_La+Delta_Cl)/3;

    % Find the sigma values, which are the eigenvalues of the transformed
    % diffusion coefficients.
    [t, sigma] = eig(D_matrix); % Sigma is the solution to a nonlinear equation, using the quadratic formula, stemming from the diffusion coefficients. Using "eig" tells MATLAB to solve this itseld, rather than us having to type in the individual, and gross-looking, equations.
    sigma_vector = diag(sigma); % Here, we put the sigma values we solved for above into a diagonal matrix. This is because sigma is the diagonal matrix of eigenvalues for the diffusion coefficient matrix.
    sigma_vector = sigma_vector(:); % We ensure that the matrix is in proper 2x1 form, that way the matrix multiplication works out.

    % Calculate the exponential term in the concentration equation after
    % performing the diffusion coefficient and concentration transform.
    exp_term = diag(exp(abs(Current_Jw*Delta_avg ./ real(sigma_vector))));

    % Calcualte the concentration differences, (deltaCa = c_int - c_perm; 
    % deltaC0 = c_bulk - c_perm), which are defined from the transformed
    % concentration back to the actual concentration,
    dc_0 = c_bulk - c_perm;
    dc_0 = dc_0(:); % creates a 2x1 matrix.
    dc_a = t*exp_term*(t\dc_0);
    c_int = dc_a + c_perm;

    % This takes the value from the first row of the c_int vector and stores
    % it into the i-th position of the c1_int results array. The same thing
    % is done for ion 1.
    c1_int(i) = c_int(1);
    c2_int(i) = c_int(2);

end

% Recover c3 from electroneutrality
z1 = abs(1); z2 = abs(3); z3 = abs(-1);
Cr_salt3(i)  = (-z1.*Cr_salt1(i)-z2.*Cr_salt2(i))./z3;

vial_index = 1:length(c1_int); % gives us the number of data points for us to plot. Mostly important for plotting interfacial concentrations for each ion per vial.

% Plots
figure(1);
scatter(Ionic_strength, c1_int);
hold on
scatter(Ionic_strength, c2_int);
hold off
xlabel('Ionic Strength (mM)');
ylabel('Interfacial Concentration (mM)');
legend('NaCl', 'LaCl3');

% Expore the data so it is easier to access.
ExportData(:,1) = round(c1_int, 2);
ExportData(:,2) = round(c2_int, 2);

function D_matrix = calculate_D(c1, c2)

    % Define the known constants.
    chi = 0;
    D1 = 1.33E-9; % [=] m2/s solution phase diffusion coefficient for Na or Ca
    D2 = 6.26E-10; % [=] m2/s solution phase diffusion coefficient for Ca or La
    D3 = 2.03E-9; % [=] m2/s solution phase diffusion coefficient for Cl
    z1 = abs(1); z2 = abs(3); z3 = abs(-1); % [=] unitless

    % Write out the equation for the diffusion coefficient denominator.
    % These equations are derived from solving the extended nerst-planck
    % flux equation assuming two cations and one common anion.
    denom = c1*(z1^2*D1 + z1*z3*D3) + c2*(z2^2*D2 + z2*z3*D3); % flux in solution

    % Compute Dij. These are the full diffusion coefficients derived from
    % the explanantion above. 
    D11 = (c1*(-z1^2*D1*D3 - z1*z3*D1*D3) + c2*(-z2^2*D1*D2 - z2*z3*D1*D3)) / denom; % [=] m2/s
    D12 = (c1*(z1*z2*D1*D2 - z1*z2*D1*D3)) / denom; % [=] m2/s

    D21 = (c2*(z1*z2*D1*D2 - z1*z2*D2*D3)) / denom; % [=] m2/s
    D22 = (c1*(-z1^2*D1*D2 - z1*z3*D2*D3) + c2*(-z2*z3*D2*D3 - z2^2*D2*D3)) / denom; % [=] m2/s

    % Put the diffusion coefficients into a 2x2 matrix.
    D_matrix = [D11, D12; D21, D22];
end

function MassTransferCoefficient = k(D_coeff) % units are m/s
        b = 0.0254; % Units [=] m
        KinematicViscosity = 1.003*10^-6; % [=] m2/s
        AngularVelocity = 2*pi*350/60; % rotations per minute (i.e., RPM)
        v0 = AngularVelocity*b/2;
        MassTransferCoefficient = 0.23*(v0^0.57)*(D_coeff^0.67)/(KinematicViscosity^0.24)/(b^0.43);
end
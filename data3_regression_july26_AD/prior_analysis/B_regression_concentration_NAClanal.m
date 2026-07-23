clear, clc, clf, close all

% Define the experiment you are analyzing.
FilePath = '/Users/wphillip/Documents/Students/Laurianne/multicomponent study/files from LLE/';
ExcelFile = 'Rejection_Analysis.xlsx';
% ExcelSheet = 'NF270_MC2 05.08.24_E1';
ExcelFile_Data = fullfile(FilePath, ExcelFile);

% Read in the data from the excel file.
NaCl_ss_c_int = readmatrix(ExcelFile_Data, 'Sheet', 'NF270_MC5 05.27.26_NaCl (2)', 'Range', 'B57:B748'); 
NaCl_ss_I_S = readmatrix(ExcelFile_Data, 'Sheet', 'NF270_MC5 05.27.26_NaCl (2)', 'Range', 'D57:D748');
NaCl_ss_c_p = readmatrix(ExcelFile_Data, 'Sheet', 'NF270_MC5 05.27.26_NaCl (2)', 'Range', 'H57:H748');
NaCl_ss_J_w = readmatrix(ExcelFile_Data, 'Sheet', 'NF270_MC5 05.27.26_NaCl (2)', 'Range', 'J57:J748');
NaCl_ss_B = readmatrix(ExcelFile_Data, 'Sheet', 'NF270_MC5 05.27.26_NaCl (2)', 'Range', 'K57:K748');

% Define known parameters and relationships.
DNa_m = 1.33E-12; % [=] m2/s
DCl_m = 2.03E-12; % [=] m2/s
l = 80E-9; % m
chi = abs(-44); % mM
D_NaCl_m = ((2*DNa_m*DCl_m)/(DNa_m+DCl_m)); % m2/s

perm_NaCl = D_NaCl_m/l*10^6; % units of um s-1 to compare with DATA2 result

% Work into form for regression
NaCl_ss_J_s = NaCl_ss_c_p.*NaCl_ss_J_w; %solute flux in units of mol m-2 s-1

deltaC_m_NaCl = (NaCl_ss_J_s*l)/(D_NaCl_m); %factor of 1000 in the denominator
plot(NaCl_ss_c_int, NaCl_ss_J_s)

figure(2)
plot(NaCl_ss_c_int, deltaC_m_NaCl)
hold on
plot(NaCl_ss_c_int, (NaCl_ss_c_int-NaCl_ss_c_p))

figure(3)
H_avg_NaCl = deltaC_m_NaCl./(NaCl_ss_c_int-NaCl_ss_c_p);
plot(NaCl_ss_c_int, H_avg_NaCl)

NaCl_ss_conc = [NaCl_ss_c_int, NaCl_ss_c_p];

% ((sqrt((beta_NaCl(1))^2+4*beta_NaCl(2).*NaCl_ss_conc(:,1).^2)-sqrt((beta_NaCl(1))^2+4*beta_NaCl(2).*NaCl_ss_conc(:,2).^2))/2);

rootfun = @(beta,c) arrayfun(@(ci) ...
    (sqrt((beta(1))^2+4.*beta(2)*ci.^2)), c);

fun_NaCl = @(beta_NaCl,NaCl_ss_conc) ...
    rootfun(beta_NaCl,NaCl_ss_conc(:,1)) - ...
    rootfun(beta_NaCl,NaCl_ss_conc(:,2));

fun_Hf_NaCl = @(beta_NaCl,NaCl_ss_conc) ...
    rootfun(beta_NaCl,NaCl_ss_conc(:,1))/NaCl_ss_conc(:,1);

 fun_Hp_NaCl = @(beta_NaCl,NaCl_ss_conc) ...
     rootfun(beta_NaCl,NaCl_ss_conc(:,1))/NaCl_ss_conc(:,1);

beta_guess_NaCl = [22.73, 0.09];

% 
% lb = [18.5 0.000000001];
% ub = [18.6 1];

[beta_estimate_NaCl, resnorm_NaCl, residual_NaCl, exitflag_NaCl, output_NaCl] = lsqcurvefit(fun_NaCl, beta_guess_NaCl, NaCl_ss_conc, deltaC_m_NaCl);

% % Plot the experimental B with the fitted B.
% label_NaCl_exp  = 'B_{NaCl} Experimental';
% label_NaCl_fit  = ['B_{NaCl} Regressed = ', num2str(beta_estimate_NaCl(1), '%5.3f'), '+', num2str(beta_estimate_NaCl(2), '%5.3f'), ' / c_{int}'];
% 
figure(4)
plot(NaCl_ss_c_int, deltaC_m_NaCl, 'r');
hold on
plot(NaCl_ss_c_int, fun_NaCl(beta_estimate_NaCl, NaCl_ss_conc), '--r');
% hold off
% xlabel('Interfacial Concentration (mM)');
% ylabel('B (\mum s^{-1})');
% legend('B_{NaCl} Experimental', 'B_{NaCl} Regressed', 'Location', 'best');
% 

figure(5)
hold on
plot(NaCl_ss_c_int, fun_Hf_NaCl(beta_estimate_NaCl, NaCl_ss_conc(:,1)), 'g');
plot(NaCl_ss_c_int, fun_Hp_NaCl(beta_estimate_NaCl, NaCl_ss_conc(:,2)), '--g');
plot(NaCl_ss_c_int, H_avg_NaCl)


chi_range = linspace(20,30,110);
del_range = logspace(-1.6,-1,110);

for i = 1:length(chi_range)
    for j = 1:length(del_range)
        calcdelC = fun_NaCl([chi_range(i), del_range(j)], NaCl_ss_conc);
        residual(i,j) = sum((calcdelC-deltaC_m_NaCl).^2);
%         figure(6)
%         plot(NaCl_ss_c_int, calcdelC, '--r');        
    end
end

figure(7)
contourf(chi_range, log10(del_range), log10(residual))
colorbar
% figure(2)
% plot(NaCl_ss_c_int, NaCl_ss_B, 'r');
% hold on
% plot(linspace(min(NaCl_ss_c_int), max(NaCl_ss_c_int), 200), fun_NaCl(beta_estimate_NaCl, linspace(min(NaCl_ss_c_int), max(NaCl_ss_c_int), 200)), '--r');
% hold off
% xlabel('Interfacial Concentration (mM)'); ylabel('B (\mum s^{-1})');
% legend(label_NaCl_exp, label_NaCl_fit, 'Location', 'best');
% 
% NaCl_ss_delC = NaCl_ss_c_int-NaCl_ss_c_p;
% NaCl_ss_Js = NaCl_ss_J_w.*NaCl_ss_c_p;
% NaCl_ss_Js_fit = fun_NaCl(beta_estimate_NaCl, NaCl_ss_c_int).*NaCl_ss_delC;
% 
% figure(3)
% plot(NaCl_ss_delC, NaCl_ss_Js, 'r');
% hold on
% plot(NaCl_ss_delC, NaCl_ss_Js_fit, '--r');
% hold off
% xlabel('c_{int} - c_{p} (mM)'); ylabel('J_{w} \cdot c_{p} (\mum^{3}/m^{2} \cdot s^{-1})');
% xlim([0 80]); ylim([-50 850]);
% legend(label_NaCl_exp, label_NaCl_fit, 'Location', 'best', 'FontSize', 8);
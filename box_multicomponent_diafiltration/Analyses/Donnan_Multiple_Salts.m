clear, clc, close all
%% User defined variables. 

% 1 and 2 are the valencies of the cations. 3 is the anion. 
% z corresponds to the valencies. 
z_1 = 1;
z_2 = 3;
z_3 = -1;
chi = -141; 

%% We will iterate through a range of concentrations. 

    % c corresponds to concentrations. The number is the ion, the subscript 
    % s is for the solution phase. 
index = 1;
for c_1_s = 10:5:200
    for c_2_s = 10:5:200    % Iterate through changing the bulk concentration of ion 2. 
        % The concentration of the anion is defined by the concentraion of the
        % cations. 
        c_3_s = (abs(z_1)*c_1_s+abs(z_2)*c_2_s)/abs(z_3);
        
        % Set up equations. 
        syms c_1_m c_2_m c_3_m
        % Simplify by taking absolute values (i.e., the equtions account for it)
        z1 = abs(z_1); z2 = abs(z_2); z3 = abs(z_3); chi = abs(chi);
        
        % Write out the right hand side and left hand side of the relevant equations. 
        c_1_left_HS = (c_1_s^(1/z1))*(((z1/z3)*c_1_s+(z2/z3)*(c_2_s))^(1/z3));
        c_1_right_HS = (c_1_m^(1/z1))*(((z1*c_1_m/z3)+(z2*c_2_s/z3)*((c_1_m/c_1_s)^(z2/z1))-(chi/z3))^(1/z3));
        c_2_left_HS = (c_2_s^(1/z2))*(((z1/z3)*c_1_s+(z2/z3)*(c_2_s))^(1/z3));
        c_2_right_HS = (c_2_m^(1/z2))*(((z2*c_2_m/z3)+(z1*c_1_s/z3)*((c_2_m/c_2_s)^(z1/z2))-(chi/z3))^(1/z3));
        c_3_left_HS = (c_3_s^(1/z3))*(((z3/z2)*c_3_s+(z1/z2)*(c_1_s))^(1/z2));
        c_3_right_HS = (c_3_m^(1/z3))*(((z3*c_3_m/z2)+(z1*c_1_s/z2)*((c_3_s/c_3_m)^(z1/z3))+(chi/z2))^(1/z2));

        % Use vpasolve to obtain the solution, turn them into numeric values using
        % 'double'. The only physically real values are the positive ones. 
        Conc_1_Mem = double(vpasolve(c_1_left_HS == c_1_right_HS,c_1_m));
        Conc_1_Mem = Conc_1_Mem(Conc_1_Mem>=0);
        Conc_2_Mem = double(vpasolve(c_2_left_HS == c_2_right_HS,c_2_m));
        Conc_2_Mem = Conc_2_Mem(Conc_2_Mem>=0);
        Conc_3_Mem = double(vpasolve(c_3_left_HS == c_3_right_HS,c_3_m));
        Conc_3_Mem = Conc_3_Mem(Conc_3_Mem>=0);
        
        % Obtain the partition coefficient of the ions relative to the bulk
        % solution concentration. 
        H_1 = Conc_1_Mem/c_1_s;
        H_2 = Conc_2_Mem/c_2_s;
        H_3 = Conc_3_Mem/c_3_s;

        Ionic_Strength = 0.5*(c_1_s*z1^2+c_2_s*z2^2+c_3_s*z3^2);
        
        InformationMatrix(index,:) = [c_1_s,c_2_s,c_3_s,Ionic_Strength,H_1,H_2,H_3];
        index = index+1;
    end
end

hold on, axis square, box on
xlabel('[+1^S]'),ylabel('Ionic Strength'),zlabel('Partition Coefficient')
plot3(InformationMatrix(:,1),InformationMatrix(:,4),InformationMatrix(:,5),'k.')
plot3(InformationMatrix(:,1),InformationMatrix(:,4),InformationMatrix(:,6),'b.')
plot3(InformationMatrix(:,1),InformationMatrix(:,4),InformationMatrix(:,7),'r.')

% hold on, axis square, box on
% xlabel('Concentration of 3^+ Cation'), ylabel('Partition Coefficient')
% plot(InformationMatrix(:,2),InformationMatrix(:,4),'b-')
% plot(InformationMatrix(:,2),InformationMatrix(:,5),'k-')
% legend('Cation +1', 'Cation +3')
% 




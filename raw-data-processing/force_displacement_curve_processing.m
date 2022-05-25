%Force displacement curve processing for tests with three circular plates,
%conducted 1/14/21 in the large sandbox 


% Reads the data for getting k-n power curves and plots all the trials

% folderpath = '../f-d-large-box/';
% folderpath = '';

%Grousered wheel
Aa = (60/1000/2)^2*3.14159;
Ab = (80/1000/2)^2*3.14159;
Ac = (100/1000/2)^2*3.14159;
[~, ~, Zs1, ~, Fzs1, ~] = read_force_displacement(folderpath, 'force-displacement_1', 11.03, 22, 1, 0, -9.62, 0.75); %60mm plate
[~, ~, Zs2, ~, Fzs2, ~] = read_force_displacement(folderpath, 'force-displacement_2', 13.73, 23.7, 1, 0, -8.55, 0); %60mm plate
[~, ~, Zs3, ~, Fzs3, ~] = read_force_displacement(folderpath, 'force-displacement_3', 10.97, 22.13, 1, 0, -9.58, .63); %60mm plate
[~, ~, Zs4, ~, Fzs4, ~] = read_force_displacement(folderpath, 'force-displacement_4', 11.67, 21.54, 1, 0, -9.0, .4); %100mm plate
[~, ~, Zs5, ~, Fzs5, ~] = read_force_displacement(folderpath, 'force-displacement_5', 9.39, 24.58, 1, 0, -9.68, 0); %100mm plate
[~, ~, Zs6, ~, Fzs6, ~] = read_force_displacement(folderpath, 'force-displacement_6', 9.6, 20, 1, 0, -9.42, 0); %100mm plate
[~, ~, Zs7, ~, Fzs7, ~] = read_force_displacement(folderpath, 'force-displacement_7', 12.65, 28, 1, 0, -8.72, 0); %80mm plate
[~, ~, Zs8, ~, Fzs8, ~] = read_force_displacement(folderpath, 'force-displacement_8', 9.9, 19.6, 1, 0, -9.33, .68); %80mm plate
[~, ~, Zs9, ~, Fzs9, ~] = read_force_displacement(folderpath, 'force-displacement_9', 12.08, 23.6, 1, 0, -9.45, .45); %80mm plate
Pzs1 = Fzs1/Aa;
Pzs2 = Fzs2/Aa; 
Pzs3 = Fzs3/Aa; 
Pzs4 = Fzs4/Ac; 
Pzs5 = Fzs5/Ac; 
Pzs6 = Fzs6/Ac; 
Pzs7 = Fzs7/Ab; 
Pzs8 = Fzs8/Ab; 
Pzs9 = Fzs9/Ab; 
expfit1 = fit(double(Zs1), double(Pzs1), 'power1');
expfit2 = fit(double(Zs2), double(Pzs2), 'power1');
expfit3 = fit(double(Zs3), double(Pzs3), 'power1');
expfit4 = fit(double(Zs4), double(Pzs4), 'power1');
expfit5 = fit(double(Zs5), double(Pzs5), 'power1');
expfit6 = fit(double(Zs6), double(Pzs6), 'power1');
expfit7 = fit(double(Zs7), double(Pzs7), 'power1');
expfit8 = fit(double(Zs8), double(Pzs8), 'power1');
expfit9 = fit(double(Zs9), double(Pzs9), 'power1');
coeffs1 = coeffvalues(expfit1);
k1 = coeffs1(1);
n1 = coeffs1(2);
coeffs2 = coeffvalues(expfit2);
k2 = coeffs2(1);
n2 = coeffs2(2);
coeffs3 = coeffvalues(expfit3);
k3 = coeffs3(1);
n3 = coeffs3(2);
coeffs4 = coeffvalues(expfit4);
k4 = coeffs4(1);
n4 = coeffs4(2);
coeffs5 = coeffvalues(expfit5);
k5 = coeffs5(1);
n5 = coeffs5(2);
coeffs6 = coeffvalues(expfit6);
k6 = coeffs6(1);
n6 = coeffs6(2);
coeffs7 = coeffvalues(expfit7);
k7 = coeffs7(1);
n7 = coeffs7(2);
coeffs8 = coeffvalues(expfit8);
k8 = coeffs8(1);
n8 = coeffs8(2);
coeffs9 = coeffvalues(expfit9);
k9 = coeffs9(1);
n9 = coeffs9(2);

n_p = median([n1 n2 n3 n4 n5 n6 n7 n8 n9]);
k_p = median([k1 k2 k3 k4 k5 k6 n7 k8 k9]);

figure()
title('Circular Plate (Large Sandbox)')
hold on
plot(expfit1, Zs1, Pzs1)
plot(expfit2, Zs2, Pzs2)
plot(expfit3, Zs3, Pzs3)
plot(expfit4, Zs4, Pzs4)
plot(expfit5, Zs5, Pzs5)
plot(expfit6, Zs6, Pzs6)
plot(expfit7, Zs7, Pzs7)
plot(expfit8, Zs8, Pzs8)
plot(expfit9, Zs9, Pzs9)
plot(Zs1, k_p*Zs1.^n_p, '--y')
xlabel('sinkage (m)');
ylabel('pressure (kPa)');
hold off

%% Force displacement curve plotting for tests with both wheels in large sandbox
% Reads the data for getting k-n power curves and plots all the trials
folderpath = 'C:\Users\cathe\Documents\MATLAB\Research_Repos\terramech-testbed\nptm_control\process_data\f-d-large-box\';

%Grousered wheel

[k1, n1, Zs1, Pzs1, Fzs1, expfit1] = read_force_displacement(folderpath, 'force-displacement_grouser_1', 12.78, 20.72, 0, 0, 10.17, 0);
[k2, n2, Zs2, Pzs2, Fzs2, expfit2] = read_force_displacement(folderpath, 'force-displacement_grouser_2', 6.59, 18.59, 0, 0, 9.15, 0);
[k3, n3, Zs3, Pzs3, Fzs3, expfit3] = read_force_displacement(folderpath, 'force-displacement_grouser_3', 10.52, 21.1, 0, 0, 8.25, -.76);

n_g = median([n1 n2 n3]);
k_g = median([k1 k2 k3]);

figure()
title('Grousered Wheel (Large Sandbox)')
hold on
plot(expfit1, Zs1, Pzs1)
plot(expfit2, Zs2, Pzs2)
plot(expfit3, Zs3, Pzs3)
% plot(Zs1, 745*Zs1.^.8, ':') %griffith 2011 values
plot(Zs1, k_g*Zs1.^n_g, '--y')
xlabel('sinkage (m)');
ylabel('pressure (kPa)');
hold off

figure()
hold on
plot(k1,n1, 'o')
plot(k2,n2, 'o')
plot(k3,n3, 'o')
plot(k_g, n_g, '*')
% plot(745, .8, '*') %griffith 2011 values
xlabel('k')
ylabel('n')
title('Grousered Wheel (Large Sandbox)')
hold off


% %Smooth wheel

[k6, n6, Zs6, Pzs6, Fzs6, expfit6] = read_force_displacement(folderpath, 'force-displacement_smooth_1', 10.2, 19.8, 0, 0, 9.98, -.71);
[k7, n7, Zs7, Pzs7, Fzs7, expfit7] = read_force_displacement(folderpath, 'force-displacement_smooth_2', 9.89, 17.7, 0, 0, 2.15, -1);
[k8, n8, Zs8, Pzs8, Fzs8, expfit8] = read_force_displacement(folderpath, 'force-displacement_smooth_3', 9.49, 17.18, 0, 0, 3.67, -.66);
n_s = median([n6 n7 n8]);
k_s = median([k6 k7 k8]);

figure()
hold on
title('Smooth Wheel (Large Sandbox)')
plot(expfit6, Zs6, Pzs6)
plot(expfit7, Zs7, Pzs7)
plot(expfit8, Zs8, Pzs8)
% %plot(Zs6, 745*Zs6.^.8, ':') %griffith 2011 values
plot(Zs6, k_s*Zs6.^n_s, '--y')
xlabel('sinkage (m)');
ylabel('pressure (kPa)');
hold off

figure()
hold on
plot(k6,n6, 'o')
plot(k7,n7, 'o')
plot(k8,n8, 'o')
plot(k_s, n_s, '*')
% plot(745, .8, '*') %griffith 2011 values
xlabel('k')
ylabel('n')
title('Smooth Wheel (Large Sandbox)')
hold off

%%
%Force displacement curve processing for tests with three circular plates,
%conducted 12/21/21 in the small sandbox 

% Reads the data for getting k-n power curves and plots all the trials
folderpath = 'C:\Users\cathe\Documents\MATLAB\Research_Repos\terramech-testbed\nptm_control\process_data\f-d-small-box\';

%Grousered wheel
Aa = (60/1000/2)^2*3.14159;
Ab = (80/1000/2)^2*3.14159;
Ac = (100/1000/2)^2*3.14159;
[~, ~, Zs11, ~, Fzs11, ~] = read_force_displacement(folderpath, 'force-displacement_1', 13.7, 20.5, 0, 0, -2.02, 0); %60mm plate
[~, ~, Zs12, ~, Fzs12, ~] = read_force_displacement(folderpath, 'force-displacement_2', 12.8, 18, 0, 0, -2.02, 0); %60mm plate
[~, ~, Zs13, ~, Fzs13, ~] = read_force_displacement(folderpath, 'force-displacement_3', 13.16, 17.5, 0, 0, -2.02, 0); %60mm plate
[~, ~, Zs14, ~, Fzs14, ~] = read_force_displacement(folderpath, 'force-displacement_4', 14.7, 22, 0, 0, -3.15, -.31); %100mm plate
[~, ~, Zs15, ~, Fzs15, ~] = read_force_displacement(folderpath, 'force-displacement_5', 10.88, 15.2, 0, 0, -3.65, 1.7); %100mm plate
[~, ~, Zs16, ~, Fzs16, ~] = read_force_displacement(folderpath, 'force-displacement_6', 10.15, 15.2, 0, 0, -3.15, 2.1); %100mm plate
[~, ~, Zs17, ~, Fzs17, ~] = read_force_displacement(folderpath, 'force-displacement_7', 13.32, 22, 0, 0, -4.18, -.6); %80mm plate
[~, ~, Zs18, ~, Fzs18, ~] = read_force_displacement(folderpath, 'force-displacement_8', 11.67, 22, 0, 0, -3.28, 0); %80mm plate
[~, ~, Zs19, ~, Fzs19, ~] = read_force_displacement(folderpath, 'force-displacement_9', 9.92, 15, 0, 0, -3.05, .39); %80mm plate
Pzs11 = Fzs11/Aa;
Pzs12 = Fzs12/Aa; 
Pzs13 = Fzs13/Aa; 
Pzs14 = Fzs14/Ac; 
Pzs15 = Fzs15/Ac; 
Pzs16 = Fzs16/Ac; 
Pzs17 = Fzs17/Ab; 
Pzs18 = Fzs18/Ab; 
Pzs19 = Fzs19/Ab; 
expfit11 = fit(double(Zs11), double(Pzs11), 'power1');
expfit12 = fit(double(Zs12), double(Pzs12), 'power1');
expfit13 = fit(double(Zs13), double(Pzs13), 'power1');
expfit14 = fit(double(Zs14), double(Pzs14), 'power1');
expfit15 = fit(double(Zs15), double(Pzs15), 'power1');
expfit16 = fit(double(Zs16), double(Pzs16), 'power1');
expfit17 = fit(double(Zs17), double(Pzs17), 'power1');
expfit18 = fit(double(Zs18), double(Pzs18), 'power1');
expfit19 = fit(double(Zs19), double(Pzs19), 'power1');
coeffs11 = coeffvalues(expfit11);
k11 = coeffs11(1);
n11 = coeffs11(2);
coeffs12 = coeffvalues(expfit12);
k12 = coeffs12(1);
n12 = coeffs12(2);
coeffs13 = coeffvalues(expfit13);
k13 = coeffs13(1);
n13 = coeffs13(2);
coeffs14 = coeffvalues(expfit14);
k14 = coeffs14(1);
n14 = coeffs14(2);
coeffs15 = coeffvalues(expfit15);
k15 = coeffs15(1);
n15 = coeffs15(2);
coeffs16 = coeffvalues(expfit16);
k16 = coeffs16(1);
n16 = coeffs16(2);
coeffs17 = coeffvalues(expfit17);
k17 = coeffs17(1);
n17 = coeffs17(2);
coeffs18 = coeffvalues(expfit18);
k18 = coeffs18(1);
n18 = coeffs18(2);
coeffs19 = coeffvalues(expfit19);
k19 = coeffs19(1);
n19 = coeffs19(2);

n_p_small = median([n11 n12 n13 n14 n15 n16 n17 n18 n19]);
k_p_small = median([k11 k12 k13 k14 k15 k16 k17 k18 k19]);

figure()
title('Circular Plate (Small Sandbox)')
hold on
plot(expfit11, Zs11, Pzs11)
plot(expfit12, Zs12, Pzs12)
plot(expfit13, Zs13, Pzs13)
plot(expfit14, Zs14, Pzs14)
plot(expfit15, Zs15, Pzs15)
plot(expfit16, Zs16, Pzs16)
plot(expfit17, Zs17, Pzs17)
plot(expfit18, Zs18, Pzs18)
plot(expfit19, Zs19, Pzs19)
plot(Zs11, k_p_small*Zs11.^n_p_small, '--y')
xlabel('sinkage (m)');
ylabel('pressure (kPa)');
hold off

%% Force displacement curve plotting for tests with both wheels in the small box
%Conducted 12/?/21 in the small sandbox
% Reads the data for getting k-n power curves and plots all the trials
folderpath = 'C:\Users\cathe\Documents\MATLAB\Research_Repos\terramech-testbed\nptm_control\process_data\f-d-small-box\';

%Grousered wheel

%[k4, n4, Zs4, Pzs4, Fzs4, expfit4] = read_force_displacement(folderpath,
%'force-displacement_grouser_4', 14.3, 23, 0, 0, -3.13, .3); %Displacement
%did not record
[k5, n5, Zs5, Pzs5, Fzs5, expfit5] = read_force_displacement(folderpath, 'force-displacement_grouser_5', 12.71, 28, 0, 0, -.18, -6.5);
%[k6, n6, Zs6, Pzs6, Fzs6, expfit6] = read_force_displacement(folderpath,
%'force-displacement_grouser_6', 10.64, 15.68, 0, 0, -3.64, 1);%
%Displacement did not record

n_g = median([n5]);
k_g = median([k5]);

figure()
title('Grousered Wheel (Small Sandbox)')
hold on
plot(expfit5, Zs5, Pzs5)
%plot(Zs5, 745*Zs5.^.8, ':') %griffith 2011 values
plot(Zs5, k_g*Zs5.^n_g, '--y')
xlabel('sinkage (m)');
ylabel('pressure (kPa)');
hold off

figure()
hold on
plot(k5,n5, 'o')
plot(k_g, n_g, '*')
% plot(745, .8, '*') %griffith 2011 values
xlabel('k')
ylabel('n')
title('Grousered Wheel (Small Sandbox)')
hold off


%Smooth wheel

%None of these trials recorded displacement data, so they can't be used.
%Have to compare to old trials below for wheels, it should be a valid
%comparison though as none of the testbed changes affected the wheel tests.
%[k1, n1, Zs1, Pzs1, Fzs1, expfit1] = read_force_displacement(folderpath, 'force-displacement_smooth_1', 6.42, 19.8, 0, 0, [], []);
% [k2, n2, Zs2, Pzs2, Fzs2, expfit2] = read_force_displacement(folderpath, 'force-displacement_smooth_2', 6.47, 17.7, 0, 0, [], []);
% [k3, n3, Zs3, Pzs3, Fzs3, expfit3] = read_force_displacement(folderpath, 'force-displacement_smooth_3', 7.92, 21.3, 0, 0, [], []);
% 
% n_s = median([n1 n2 n3]);
% k_s = median([k1 k2 k3]);

% figure()
% hold on
% title('Smooth Wheel')
% % plot(expfit1, Zs1, Pzs1)
% % plot(expfit2, Zs2, Pzs2)
% % plot(expfit3, Zs3, Pzs3)
% %plot(Zs6, 745*Zs6.^.8, ':') %griffith 2011 values
% plot(Zs6, k_s*Zs6.^n_s, '--y')
% xlabel('sinkage (m)');
% ylabel('pressure (kPa)');
% hold off
% 
% figure()
% hold on
% % plot(k1,n1, 'o')
% % plot(k2,n2, 'o')
% % plot(k3,n3, 'o')
% plot(k_s, n_s, '*')
% plot(745, .8, '*') %griffith 2011 values
% xlabel('k')
% ylabel('n')
% title('Smooth Wheel')
% hold off

%%
%Force displacement curve processing for old tests with a circular plate

% Reads the data for getting k-n power curves and plots all the trials
folderpath = 'C:\Users\cathe\Documents\MATLAB\Research_Repos\terramech-testbed\nptm_control\process_data\insertion_Bill\';

%Grousered wheel
Aa = (60/100/2)^2*3.14159;
[~, ~, Zs1, ~, Fzs1, ~] = read_force_displacement(folderpath, 'force-displacement_1', 24.4, 33.28, 0, 1, [], []);
[~, ~, Zs13, ~, Fzs13, ~] = read_force_displacement(folderpath, 'force-displacement_2', 15.38, 19.5, 0, 1, [], []);
Pzs1 = Fzs1/Aa;
Pzs13 = Fzs13/Aa;
expfit1 = fit(double(Zs1), double(Pzs1), 'power1');
expfit13 = fit(double(Zs13), double(Pzs13), 'power1');
coeffs1 = coeffvalues(expfit1);
k12 = coeffs1(1);
n12 = coeffs1(2);
coeffs13 = coeffvalues(expfit13);
k13 = coeffs13(1);
n13 = coeffs13(2);

n_p = median([n12 n13]);
k_p = median([k12 k13]);

figure()
title('Circular Plate (Old Tests)')
hold on
plot(expfit1, Zs1, Pzs1)
plot(expfit13, Zs13, Pzs13)
%plot(Zs1, 745*Zs1.^.8, ':') %griffith 2011 values
plot(Zs1, k_p*Zs1.^n_p, '--y')
xlabel('sinkage (m)');
ylabel('pressure (kPa)');
hold off

%% Force displacement curve plotting for old tests with the wheels
% Reads the data for getting k-n power curves and plots all the trials
folderpath = 'C:\Users\cathe\Documents\MATLAB\Research_Repos\terramech-testbed\nptm_control\process_data\force-displacement\';

%Grousered wheel

[k21, n21, Zs21, Pzs21, Fzs21, expfit21] = read_force_displacement(folderpath, 'force-displacement_1', 9.73, 27.29, 0, 1, [], []);
[k22, n22, Zs22, Pzs22, Fzs22, expfit22] = read_force_displacement(folderpath, 'force-displacement_2', 7.93, 19.5, 0, 1, [], []);
[k23, n23, Zs23, Pzs23, Fzs23, expfit23] = read_force_displacement(folderpath, 'force-displacement_3', 7.90, 21.1, 0, 1, [], []);
[k24, n24, Zs24, Pzs24, Fzs24, expfit24] = read_force_displacement(folderpath, 'force-displacement_4', 11.5, 24.7, 0, 1, [], []);
[k25, n25, Zs25, Pzs25, Fzs25, expfit25] = read_force_displacement(folderpath, 'force-displacement_5', 10.64, 21.9, 0, 1, [], []);

n_g_old = median([n21 n22 n23 n24 n25]);
k_g_old = median([k21 k22 k23 k24 k25]);

figure()
title('Grousered Wheel (Old Tests)')
hold on
plot(expfit21, Zs21, Pzs21)
plot(expfit22, Zs22, Pzs22)
plot(expfit23, Zs23, Pzs23)
plot(expfit24, Zs24, Pzs24)
plot(expfit25, Zs25, Pzs25)
%plot(Zs1, 745*Zs1.^.8, ':') %griffith 2011 values
plot(Zs21, k_g_old*Zs21.^n_g_old, '--y')
xlabel('sinkage (m)');
ylabel('pressure (kPa)');
hold off

figure()
hold on
plot(k21,n21, 'o')
plot(k22,n22, 'o')
plot(k23,n23, 'o')
plot(k24,n24, 'o')
plot(k25,n25, 'o')
plot(k_g_old, n_g_old, '*')
plot(745, .8, '*') %griffith 2011 values
xlabel('k')
ylabel('n')
title('Grousered Wheel (Old Tests)')
hold off


%Smooth wheel

[k26, n26, Zs26, Pzs26, Fzs26, expfit26] = read_force_displacement(folderpath, 'force-displacement_6', 6.42, 19.8, 0, 1, [], []);
[k27, n27, Zs27, Pzs27, Fzs27, expfit27] = read_force_displacement(folderpath, 'force-displacement_7', 6.47, 17.7, 0, 1, [], []);
[k28, n28, Zs28, Pzs28, Fzs28, expfit28] = read_force_displacement(folderpath, 'force-displacement_8', 7.92, 21.3, 0, 1, [], []);
[k29, n29, Zs29, Pzs29, Fzs29, expfit29] = read_force_displacement(folderpath, 'force-displacement_9', 13.16, 25.1, 0, 1, [], []);
[k30, n30, Zs30, Pzs30, Fzs30, expfit30] = read_force_displacement(folderpath, 'force-displacement_10', 7.75, 21.9, 0, 1, [], []);
[k31, n31, Zs31, Pzs31, Fzs31, expfit31] = read_force_displacement(folderpath, 'force-displacement_11', 4.5, 17.7, 0, 1, [], []);

n_s = median([n26 n27 n28 n29 n30 n31]);
k_s = median([k26 k27 k28 k29 k30 k31]);

figure()
hold on
title('Smooth Wheel (Old Tests)')
plot(expfit26, Zs26, Pzs26)
plot(expfit27, Zs27, Pzs27)
plot(expfit28, Zs28, Pzs28)
plot(expfit29, Zs29, Pzs29)
plot(expfit30, Zs30, Pzs30)
plot(expfit31, Zs31, Pzs31)
%plot(Zs6, 745*Zs6.^.8, ':') %griffith 2011 values
plot(Zs26, k_s*Zs26.^n_s, '--y')
xlabel('sinkage (m)');
ylabel('pressure (kPa)');
hold off

figure()
hold on
plot(k26,n26, 'o')
plot(k27,n27, 'o')
plot(k28,n28, 'o')
plot(k29,n29, 'o')
plot(k30,n30, 'o')
plot(k31,n31, 'o')
plot(k_s, n_s, '*')
plot(745, .8, '*') %griffith 2011 values
xlabel('k')
ylabel('n')
title('Smooth Wheel (Old Tests)')
hold off
%% Plot flat plate tests and wheel tests for both sandboxes
figure()
title('Circular Plates, both sand boxes')
hold on
plot(expfit11, Zs11, Pzs11, 'c')
plot(expfit1, Zs1, Pzs1, 'b')
plot(expfit2, Zs2, Pzs2, 'b')
plot(expfit3, Zs3, Pzs3, 'b')
plot(expfit4, Zs4, Pzs4, 'b')
plot(expfit5, Zs5, Pzs5, 'b')
plot(expfit6, Zs6, Pzs6, 'b')
plot(expfit7, Zs7, Pzs7, 'b')
plot(expfit8, Zs8, Pzs8, 'b')
plot(expfit9, Zs9, Pzs9, 'b')
plot(expfit12, Zs12, Pzs12, 'c')
plot(expfit13, Zs13, Pzs13, 'c')
plot(expfit14, Zs14, Pzs14, 'c')
plot(expfit15, Zs15, Pzs15, 'c')
plot(expfit16, Zs16, Pzs16, 'c')
plot(expfit17, Zs17, Pzs17, 'c')
plot(expfit18, Zs18, Pzs18, 'c')
plot(expfit19, Zs19, Pzs19, 'c')
plot(Zs11, k_p_small*Zs11.^n_p_small, '--y')
legend('Small Sandbox', 'Expfit', 'Large Sandbox')
xlabel('sinkage (m)');
ylabel('pressure (kPa)');
hold off


figure()
hold on
title('Smooth Wheel (Both sandboxes)')
plot(expfit26, Zs26, Pzs26, 'c')
plot(expfit6, Zs6, Pzs6, 'b')
plot(expfit7, Zs7, Pzs7, 'b')
plot(expfit8, Zs8, Pzs8, 'b')
plot(expfit27, Zs27, Pzs27, 'c')
plot(expfit28, Zs28, Pzs28, 'c')
plot(expfit29, Zs29, Pzs29, 'c')
plot(expfit30, Zs30, Pzs30, 'c')
plot(expfit31, Zs31, Pzs31, 'c')
%plot(Zs6, 745*Zs6.^.8, ':') %griffith 2011 values
plot(Zs26, k_s*Zs26.^n_s, '--y')
legend('Small Sandbox', 'Expfit', 'Large Sandbox')
xlabel('sinkage (m)');
ylabel('pressure (kPa)');
hold off

figure()
hold on
plot(expfit21, Zs21, Pzs21, 'c')
plot(expfit1, Zs1, Pzs1, 'b')
plot(expfit2, Zs2, Pzs2, 'b')
plot(expfit3, Zs3, Pzs3, 'b')
plot(expfit22, Zs22, Pzs22, 'c')
plot(expfit23, Zs23, Pzs23, 'c')
plot(expfit24, Zs24, Pzs24, 'c')
plot(expfit25, Zs25, Pzs25, 'c')
plot(745, .8, '*') %griffith 2011 values
legend('Small Sandbox', 'Expfit', 'Large Sandbox')
xlabel('k')
ylabel('n')
title('Grousered Wheel (Both sandboxes)')
hold off
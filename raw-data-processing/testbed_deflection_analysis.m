function linfit = testbed_deflection_analysis(disp_plot)
% Reads and plots the deflection of the testbed wrt force
folderpath = 'C:\Users\cathe\Documents\MATLAB\Research_Repos\terramech-testbed\nptm_control\process_data\testbed-deflection\';

[Zs1, Fzs1] = read_testbed_deflection(folderpath, 'force-displacement_1', 8, 19.5, 0, 1);
[Zs2, Fzs2] = read_testbed_deflection(folderpath, 'force-displacement_2', 13.9, 19.5, 0, 1);
[Zs3, Fzs3] = read_testbed_deflection(folderpath, 'force-displacement_3', 14.2, 23, 0, 1);
[Zs4, Fzs4] = read_testbed_deflection(folderpath, 'force-displacement_4', 11.6, 17, 0, 1);
[Zs5, Fzs5] = read_testbed_deflection(folderpath, 'force-displacement_5', 9.45, 13.85, 0, 1);
[Zs6, Fzs6] = read_testbed_deflection(folderpath, 'force-displacement_6', 5.9, 10.7, 0, 1);

allZs = [Zs1; Zs2; Zs3; Zs4; Zs5; Zs6];
allFzs = [Fzs1; Fzs2; Fzs3; Fzs4; Fzs5; Fzs6];

ft = fittype({'x'}); %Force linear fit to have a y-intercept of 0
linfit = fit(double(allFzs), double(allZs), ft);

if disp_plot
    figure()
    hold on
    plot(Fzs1, Zs1)
    plot(Fzs2, Zs2)
    plot(Fzs3, Zs3)
    plot(Fzs4, Zs4)
    plot(Fzs5, Zs5)
    plot(Fzs6, Zs6)
    plot(linfit, allFzs, allZs)
    ylabel('Sinkage (m)');
    xlabel('Force (kN)');
end
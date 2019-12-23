% This script plots the p-value matrices obtained from the R-script
% MultivariateAOV. Richard Ott, 2019
clc
clear
close all
printfig = 1;

% load data
sigL = [0.1,0.05];                % significance levels
[~,sheets] = xlsfinfo('tukey.xlsx');
data = cell(length(sheets),1);
raw = cell(length(sheets),1);
for i = 1:length(sheets)
    [data{i},~,raw{i}] = xlsread('tukey.xlsx',i);
end
labs = {'su','vc','pl','ss','sm','sc','mt'};


% PLOT P-VALUE MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
for i = 1:length(sheets)
    % set up matrix to plot
    pplot = nan(7,7);
    pplot(2:7,1) = data{i}(1:6,5);
    pplot(3:7,2) = data{i}(7:11,5);
    pplot(4:7,3) = data{i}(12:15,5);
    pplot(5:7,4) = data{i}(16:18,5);
    pplot(6:7,5) = data{i}(19:20,5);
    pplot(7,6) = data{i}(21,5);
    % color pplot by significance
    pplot(pplot < 1 & pplot > sigL(1)) = 1;
    pplot(pplot < sigL(1) & pplot > sigL(2)) = 2;
    pplot(pplot < sigL(2)) = 3;
    pplot(triu(pplot) ~= 0) = 2;
    
    subplot(3,3,i)
    imagesc(pplot)
    % imagesc(plotCorr)
    colormap redblue
    hold on;
    for j = 1:length(pplot)+1
       plot([.5,length(pplot)+.5],[j-.5,j-.5],'k-');
       plot([j-.5,j-.5],[.5,length(pplot)+.5],'k-');
    end
    caxis([1,3])
    title(sheets{i});
    axis tight
    textStrings = num2str(pplot(tril(pplot)~=0), '%0.2f');       % Create strings from the matrix values
    textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
    [x, y] = meshgrid(1:length(pplot));  % Create x and y coordinates for the strings
    set(gca, 'XTick', 1:length(pplot), ...                             % Change the axes tick marks
             'YTick', 1:length(pplot), ...
             'XTickLabel',labs,...
             'YTickLabel',labs,...
             'TickLength', [0 0],...
             'FontSize',8);
end

if printfig
    print('Tukey_pvals','-depsc','-painters')
end

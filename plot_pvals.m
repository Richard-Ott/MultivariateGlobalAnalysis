% This script plots the p-value matrices obtained from the R-script
% MultivariateAOV. Richard Ott, 2019
clc
clear
close all

% load data
sigL = [0.1,0.05];                % significance levels
[raw,~,~] = xlsread('pvalues2.xlsx');
n = length(find(isnan(raw(:,1))));
vars = {"Elev","Locrel","slope","P","NDVI","ksn","tet_corr","Amph_corr","Mam_corr"};
for j = 1:n
    pvals{j} = raw(2+(j-1)*8:2+(j-1)*8+6,2:end);
end

% PLOT P-VALUE MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
for i = 1:n
    subplot(3,3,i)
    
    p = pvals{i};
    pplot = nan(size(p));
    pplot(p < 1 & p > sigL(1)) = 2;
    pplot(p < sigL(1) & p > sigL(2)) = 1;
    pplot(p < sigL(2)) = 3;
    pplot(triu(pplot) ~= 0) = 2;

    imagesc(pplot)
    % imagesc(plotCorr)
    colormap redblue
    hold on;
    for j = 1:length(p)+1
       plot([.5,length(p)+.5],[j-.5,j-.5],'k-');
       plot([j-.5,j-.5],[.5,length(p)+.5],'k-');
    end
    caxis([1,3])
    title(vars{i});
    axis tight
    textStrings = num2str(p(tril(p)~=0), '%0.2f');       % Create strings from the matrix values
    textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
    [x, y] = meshgrid(1:length(p));  % Create x and y coordinates for the strings
    set(gca, 'XTick', 1:length(p), ...                             % Change the axes tick marks
             'YTick', 1:length(p), ...
             'TickLength', [0 0],...
             'FontSize',8);
end

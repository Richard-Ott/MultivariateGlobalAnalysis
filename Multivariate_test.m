clear
clc

% Data Input
printfig = 0;

load testData.mat
data = testData;

depVar = 7;             % dependent variable
inds = logical([1 1 1 0 1 1 0 1 1]); % which parameters to start with
vars = {'elev','locrel','slope','P','NDVI','Lat','PSI'};

% lithologies have to be grouped by major rocktypes -----------------------
Geo = nan(length(data),1);
Geo(data(:,depVar) == 1) = 1;                   % su
Geo(ismember(data(:,depVar),[3,7,8,9])) = 2;    % vc
Geo(ismember(data(:,depVar),[10,11,12])) = 3;   % pl
Geo(data(:,depVar) == 2) = 7;                   % ss
Geo(data(:,depVar) == 4) = 5;                   % sm 
Geo(data(:,depVar) == 5) = 6;                   % sc
Geo(data(:,depVar) == 13) = 4;                  % mt
ids = {'su','vc','pl','mt','sm','sc','ss'};
% -------------------------------------------------------------------------

% Do the regression
X = data(:,inds);
Y = categorical(Geo,1:7,ids);
[B,dev,stats] = mnrfit(X,Y);

% find likelihood values


%% Plot results --------------------------------------------------------- %
% colorcode t-values by red, yellow green for significance
cind = nan(size(stats.t));
cind(stats.t < 1) = 1;
cind(stats.t > 1 & stats.t < 5) = 2;
cind(stats.t > 5) = 3;

[~,n] = size(X);
cc = [1,0,0;1,1,0;0,1,0];   % red,yellow,green   
% plot coefficients
for i = 1:n                 % loop through predictor variables    
    subplot(3,3,i)
    scatter(B(i+1,:),1:length(ids)-1,20,cc(cind(i+1,:),:),'filled','MarkerEdgeColor','k');  
    set(gca,'Ytick',1:length(ids)-1,'YTickLabels',ids)
    hold on
    vline(0,'k-')
    title(vars(i))
end


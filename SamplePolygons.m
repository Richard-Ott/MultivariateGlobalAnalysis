% THis script samples from raster with different polygon ID's. For every polygon the mean
% parameters from different grids will de determined.
% Richard Ott 2019

clc
clear 
close all

addpath('E:\Richard\Global_Hypsometry\matlab\Seans_functions')
addpath('E:\Richard\Global_Hypsometry\matlab\river_chi_profiles_with_geology')
addpath('E:\Richard\Global_Hypsometry\GIS\sinusoidal') 

% Input data
DEM = GRIDobj('DEM_sinus1km_50.tif');
DEM.Z(DEM.Z > 10000) = nan;
DEM.Z(DEM.Z == 0)= nan; 
[xMinP,xMaxP,yMinP,yMaxP] = findCorners(DEM);

% load Polygon raster
Polygons = GRIDobj('Geo_sin.tif');
Polygons = gridReSize(Polygons,xMinP,xMaxP,yMinP,yMaxP);
Polygons = resample(Polygons,DEM,'nearest');

% keep only the Polygons areas that are covered by DEM
ind = isnan(DEM.Z);
Polygons.Z(ind) = nan;
n = unique(Polygons.Z(:),'sorted');
n(isnan(n)) = [];

% sample FIDs in shapefile
data = nan(length(n),15);                     % set up empty data matrix
pinds = cell(length(n),1);
for i = 1:length(n)                           % loop through polygons
    in = find(Polygons.Z == n(i));       % get indices of this polygon
    pinds{i} = in;
    data(i,end) = length(in);         % weighting factor for polygon by pixel count
end
clear Polygons

meanVals = AreaAverage(DEM,pinds,'numerical');                  % Average elevation
meanVals(isnan(meanVals)) = -9999;
data(:,1) = int16(meanVals); 

Locrel = GRIDobj('Locrel_sin.tif');                             % Local relief
Locrel = gridReSize(Locrel,xMinP,xMaxP,yMinP,yMaxP);
Locrel = resample(Locrel,DEM,'bilinear');
meanVals = AreaAverage(Locrel,pinds,'numerical');
meanVals(isnan(meanVals)) = -9999;
data(:,2) = int16(meanVals);       clear Locrel

slope = GRIDobj('slope_sin.tif');                               % slope value times 10 to make integer but save one decimal
slope = gridReSize(slope,xMinP,xMaxP,yMinP,yMaxP);
slope.Z(slope.Z == -32768) = nan;
slope = resample(slope,DEM,'bilinear');
meanVals = AreaAverage(slope,pinds,'numerical');
meanVals(isnan(meanVals)) = -9999;  
meanVals(meanVals ~= -9999) = meanVals(meanVals ~= -9999).* 10;
data(:,3) = int16(meanVals);       clear slope

P = GRIDobj('CHELSA_P_sin.tif');                               % Precipitaiton
P = gridReSize(P,xMinP,xMaxP,yMinP,yMaxP);
P = resample(P,DEM,'bilinear');
meanVals = AreaAverage(P,pinds,'numerical');
meanVals(isnan(meanVals)) = -9999;  
data(:,4) = int16(meanVals);       clear P

T = GRIDobj('CHELSA_T_sin.tif');                               % Temperature
T = gridReSize(T,xMinP,xMaxP,yMinP,yMaxP);
T = resample(T,DEM,'bilinear');
meanVals = AreaAverage(T,pinds,'numerical');
meanVals(isnan(meanVals)) = -9999;  
data(:,5) = int16(meanVals);       clear T

NDVI = GRIDobj('NDVI_avg_sin.tif');                       % NDVI avg, leave values as 0 to 255 so that I can store integer 
% NDVI.Z = NDVI.Z./255;                                   % original values are from zero to 255, normalize again so that values are from 0 to 1
NDVI = gridReSize(NDVI,xMinP,xMaxP,yMinP,yMaxP);
NDVI = resample(NDVI,DEM,'bilinear');
meanVals = AreaAverage(NDVI,pinds,'numerical');
meanVals(isnan(meanVals)) = -9999;  
data(:,6) = int16(meanVals);    clear NDVI

Ksn = GRIDobj('Ksn_sin.tif');                               % Ksn
Ksn = gridReSize(Ksn,xMinP,xMaxP,yMinP,yMaxP);
Ksn = resample(Ksn,DEM,'bilinear');
meanVals = AreaAverage(Ksn,pinds,'numerical');
meanVals(isnan(meanVals)) = -9999;  
data(:,7) = int16(meanVals);       clear Ksn

Tetrapod = GRIDobj('tetrapod_sin.tif');                               % Tetrapod
Tetrapod = gridReSize(Tetrapod,xMinP,xMaxP,yMinP,yMaxP);
Tetrapod = resample(Tetrapod,DEM,'bilinear');
meanVals = AreaAverage(Tetrapod,pinds,'numerical');
meanVals(isnan(meanVals)) = -9999;  
data(:,8) = int16(meanVals);       clear Tetrapod

Amph = GRIDobj('amph_sin.tif');                               % Amphibians
Amph = gridReSize(Amph,xMinP,xMaxP,yMinP,yMaxP);
Amph = resample(Amph,DEM,'bilinear');
meanVals = AreaAverage(Amph,pinds,'numerical');
meanVals(isnan(meanVals)) = -9999;  
data(:,9) = int16(meanVals);       clear Amph

Mammal = GRIDobj('mammal_sin.tif');                               % Mammals
Mammal = gridReSize(Mammal,xMinP,xMaxP,yMinP,yMaxP);
Mammal = resample(Mammal,DEM,'bilinear');
meanVals = AreaAverage(Mammal,pinds,'numerical');
meanVals(isnan(meanVals)) = -9999;  
data(:,10) = int16(meanVals);       clear Mammal

% The next 2 should be resampled to the same resolution beforehand in
% ArcGIS directly from their polygons...
% decide which Koeppen Geiger classes to group 
idkg = nan(5,28-16);                              % Group KG classes
idkg(1,1:3) = [1,2,3];                            % Tropical
idkg(2,1:4) = [4,5,6,7];                          % Dry
idkg(3,1:length(8:1:16)) = 8:1:16;                % Temperate
idkg(4,:) = 17:1:28;                              % Continental
idkg(5,1:4) = [29,30,31,32];                      % Polar
idkg_text = {'tropical','dry','temperate','continental','polar'};
KG = GRIDobj('koeppengeiger_sin.tif');                      % Köppen Geiger classification
% group values from KG
[row,~] = size(idkg);
for i = 1:row
    KG.Z(ismember(KG.Z,idkg(i,:))) = i;
end
KG = gridReSize(KG,xMinP,xMaxP,yMinP,yMaxP);
KG = resample(KG,DEM,'nearest');
meanVals = AreaAverage(KG,pinds,'categorical');
meanVals(isnan(meanVals)) = -9999;  
data(:,11) = int16(meanVals);       clear KG

GEO = GRIDobj('globalgeo1km_sin.tif');                     % main lithology from Hartmann 2012
GEO = gridReSize(GEO,xMinP,xMaxP,yMinP,yMaxP);
GEO = resample(GEO,DEM,'nearest');
meanVals = AreaAverage(GEO,pinds,'categorical');
meanVals(isnan(meanVals)) = -9999;  
data(:,12) = int16(meanVals);       clear GEO

LC = GRIDobj('landcover_sin.tif');                     % Landcover data from GLCNMO
LC = gridReSize(LC,xMinP,xMaxP,yMinP,yMaxP);
LC = resample(LC,DEM,'nearest');
meanVals = AreaAverage(LC,pinds,'categorical');
meanVals(isnan(meanVals)) = -9999;  
data(:,13) = int16(meanVals);       clear LC

SR = GRIDobj('GSRM_sin.tif');                     % Global strain rate data 
SR = gridReSize(SR,xMinP,xMaxP,yMinP,yMaxP);
SR = resample(SR,DEM,'bilinear');
meanVals = AreaAverage(SR,pinds,'numerical');
meanVals(isnan(meanVals)) = -9999;  
data(:,14) = int16(meanVals);       clear SR

% only keep rows where all variables are defined
% HERE EVERYTHING > 60 DEGREE GETS REMOVED BECAUSE OF MISSING KSN!!!!!
data(any([data(:,1:8),data(:,10:end)] == -9999, 2), :) = []; % leave out amphibians for this, because they have large undefinded areas
save('global_poly_50mLC.mat','data','-v7.3')
dlmwrite('global_poly_50mLC.txt',data)
load handel
sound(y,Fs)
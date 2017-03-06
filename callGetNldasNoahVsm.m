% This script will call the getNldasVsm function
% P. Shellito
% 7/29/16

clear all
close all

% -------------------------------------------------------------------------
% In this example, there are sites with latitude and longitude in the
% following input file:
inFile = './inFile_test.txt';
% Date range requested
qStart = [2006,1,1];
% qStart = 'apnd';
qEnd = [2017,1,27];

% -------------------------------------------------------------------------
% What variables to record

% Rain is the liquid precipitation, accumulated over the previous hour (kg/m2)
% Snow is the frozen precipitation, accumulated over the previous hour (kg/m2)
% LSM_0_10 is the liquid soil moisture content between 0 and 10 cm, instantaneous (kg/m2)
% LSM_10_40 is the liquid soil moisture content between 10 and 40 cm, instantaneous (kg/m2)
% LSM_40_100 is the liquid soil moisture content between 40 and 100 cm, instantaneous (kg/m2)
% LSM_100_200 is the liquid soil moisture content between 100 and 200 cm, instantaneous (kg/m2)
% TSM_0_10 is the total soil moisture content between 0 and 10 cm, instantaneous (kg/m2)
% TSM_10_40 is the total soil moisture content between 10 and 40 cm, instantaneous (kg/m2)
% TSM_40_100 is the total soil moisture content between 40 and 100 cm, instantaneous (kg/m2)
% TSM_100_200 is the total soil moisture content between 100 and 200 cm, instantaneous (kg/m2)
% LHF is the latent heat flux, averaged over the previous hour (W/m2)
% PotLHF is the potential latent heat flux, averaged over the previous hour (W/m2)
% SHF is the sensible heat flux, averaged over the previous hour
% G is the ground heat flux, averaged over the previous hour (W/m2)
% ET is the total evapotranspiration, accumulated over the previous hour (same as LHF but different units) (kg/m2)
% Transp is the transpiration, averaged over the previous hour (W/m2)
% EDir is the direct evaporation from bare soil, averaged over the previous hour (W/m2)
% ECanopy is the canopy water evaporation, averaged over the previous hour (W/m2)
% Sublim is the sublimation (evaporation from snow), averaged over the previous hour (W/m2)
% CanopyH2O is the plant canopy surface water, instantaneous (kg/m2)
% SnoDepth is the snow depth, measured instantaneously (m)
% SWE is the water equivalent of accumulated snow depth, instantaneous (kg/m2)
% SnoFrac is the snow cover, instantaneous (fraction)
% SnoMelt is the snow melt, accumulated over the previous hour (kg/m2)
% QSub is the subsurface runoff, accumulated over the previous hour (kg/m2)
% QSurf is the surface runoff (non-infiltrating), accumulated over the previous hour (kg/m2)
% SolDn is the shortwave radiation down, averaged over the previous hour (W/m2)
% LWDn is the long wave radiation down, averaged over the previous hour (W/m2)
% SolNet is the net shortwave radiation, averaged over the previous hour (W/m2)
% LWNet is the net long wave radiation, averaged over the previous hour (W/m2)
% Albedo is the albedo, instantaneous (percent)
% LAI is the leaf area index, instantaneous (unitless)
% Veg is the vegetation cover, instantaneous (fraction)

qVars = {'Rain', 'Snow', ...
    'LSM_0_10', 'LSM_10_40', 'LSM_40_100', 'LSM_100_200', ...
    'TSM_0_10', 'TSM_10_40', 'TSM_40_100', 'TSM_100_200', ...
    'PotLHF', 'LHF', 'ET', 'EDir', 'ECanopy', 'Transp', 'Sublim', ...
    'SolDn', 'LWDn', 'SolNet', 'LWNet', 'Albedo'};
qVars = {'Rain', 'Snow', ...
    ... 'LSM_0_10', 'LSM_10_40', 'LSM_40_100', 'LSM_100_200', ...
    'TSM_0_10', 'TSM_10_40', 'TSM_40_100', 'TSM_100_200'};

% -------------------------------------------------------------------------
% Read the input data from the text file
% Open the input file
fid = fopen(inFile);
% Read the data in the input file
data = textscan(fid,'%s\t%f\t%f', 'headerlines', 1);
% Close the input file
fclose(fid);

% -------------------------------------------------------------------------
% Organize input data and create the needed vectors to pass into the function
% A cell array of strings
qNames = data{1,1};
% Latitude of the sites in qNames
qLat = data{1,2};
% Longitude of the sites in qNames
qLon = data{1,3};
% Directory to hold the output text files
outDir = './nldasVsm';

% -------------------------------------------------------------------------
% Record what time is is before the function is called
disp('Starting the script at')
disp(datetime)
startTime = datetime;

fid = fopen('startStopTime.txt','w');
fprintf(fid, ['script started at ' datestr(startTime) '\n']);
fclose(fid);

% -------------------------------------------------------------------------
% Call the function
outDirectory = getNldasNoahVsm(qNames, qLat, qLon, qStart, qEnd, qVars, outDir);

% -------------------------------------------------------------------------
% Report where the data are held and how long the script took to run

fid = fopen('startStopTime.txt', 'a');
fprintf(fid, ['script ended at   ' datestr(datetime)]);
fclose(fid);

disp('Finished! Site data can be found here:')
disp(outDirectory)
disp('Start and finish times were:')
disp(startTime)
disp(datetime)

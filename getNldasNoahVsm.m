function [ outDir ] = getNldasNoahVsm(qNames, qLat, qLon, qStart, qEnd, qVars, outDir)

% GETNLDASFORCING This script will download NLDAS primary forcing data from
%       Nasa's servers.
% Created by Peter J. Shellito 2/16/16
% 
% qNames: cell array of site names
% qLat: latitude corresponding to each site
% qLon: longitude corresponding to each site
% qStart: a vector [yyyy, mm, dd] specifying the day to start downloading.
%       Start date must be [1979, 1, 2] or later. If set to 'apnd,'
%       script will look for an existing file with the
%       same site name and start at the last date in that file and
%       append new data to it.
% qEnd: a vector [yyyy, mm, dd] specifying the day to stop downloading. End
%       date must be 4 days before today or earlier. Neither qStart nor 
%       qEnd support a starting hour and minute.
% qVars: a cell array specifying which variables to record. Choose between:
%       
% outDir: Directory to place the output files. If nothing is provided,
%       default is to create a directory in the present directory titled
%       './outFiles/'
% 
% For your own reference, or if something goes wrong, the directory where
% you can find all the NLDAS data:
% http://disc.sci.gsfc.nasa.gov/hydrology/data-holdings
% And the location of a sample file:
% fileDir = 'ftp://hydro1.sci.gsfc.nasa.gov/data/s4pa/NLDAS/NLDAS_NOAH0125_H.002/2016/001/NLDAS_NOAH0125_H.A20160101.0000.002.grb';
% 
% Initial testing revealed that processing one month took 20 minutes using
% a CU internet connection.
% 
% ========================================================================
% NB: This script requires the user to have downloaded and installed the
% nctoolbox: http://nctoolbox.github.io/nctoolbox/
% Download the stable release zip file nctoolbox-20130305 (NOT THE MOST 
% RECENT ONE, which has a bug. see 
% https://github.com/nctoolbox/nctoolbox/issues/61). 
% Extract nctoolbox-20130305, and then in the matlab command line, navigate 
% to the nctoolbox: cd('Path/to/toolbox/nctoolbox-1.1.0')
% Then type "setup_nctoolbox" and you should be able to use the toolbox.
% It may be necessary to run setup_nctoolbox every time you start a new 
% matlab session. To automate this, add the '/Path/to/toolbox/ to your 
% matlab path. Then edit (or create) a startup.m file that is also in your
% matlab path. Add the following to startup.m: setup_nctoolbox;
% ========================================================================
% 
% Copyright (c) 2016, Peter J. Shellito, University of Colorado Boulder
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
% 
% 2. Redistributions in binary form must reproduce the above copyright 
% notice, this list of conditions and the following disclaimer in the 
% documentation and/or other materials provided with the distribution.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS 
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR 
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% -------------------------------------------------------------------------
% If no output directory was provided
if nargin<7
    outDir = './outFiles';
end

% -------------------------------------------------------------------------
% Variables available to request
varsAvail = {'Rain', 'Snow',...
    'LSM_0_10', 'LSM_10_40', 'LSM_40_100', 'LSM_100_200', ...
    'TSM_0_10', 'TSM_10_40', 'TSM_40_100', 'TSM_100_200', ...
    'LHF', ...
    'PotLHF', 'SHF', 'G', 'Transp', 'EDir', 'ET', 'ECanopy', ...
    'Sublim', 'CanopyH2O', 'SnoDepth', 'SWE', 'SnoFrac', 'SnoMelt', ...
    'QSub', 'QSurf', 'SolDn', 'LWDn', 'SolNet', 'LWNet', 'Albedo', ...
    'LAI', 'Veg'};
% The strings found in the grib file that correspond to these variables.
% All variable descriptions are taken from the NDLAS2 readme file included
%   in this directory.
% Note that LSM and TSM are repeated because all levels of those variables
% come from the same matrix in the grib file.
% "Instantaneous" means the data values are "at exactly 00 minute of every 
%   hour"
% "Averaged" means the data values are the average over the previous hour 
%   of the time listed in the file. For example, for the 03Z files, the data 
%   values are the average over the time from 02Z to 03Z.
% "Accumulated" means the data values are the accumulation over the previous 
%   hour of the time listed in the file. For example, for the 03Z files, the 
%   data values are the accumulation over the time from 02Z to 03Z.
varsStrings = {'Liquid_precipitation_rainfall_surface_1_Hour_Average', ... % 1x224x464. Rainfall. kg/m2 accumulated.
    'Frozen_precipitation_eg_snowfall_surface_1_Hour_Average', ... % 1x224x464. Snowfall. kg/m2 accumulated.
    'Liquid_soil_moisture_content_non-frozen_layer_between_two_depths_below_surface_layer', ... % a 1 by 4 by 224 by 464 matrix. Nonfrozen moisture. Depths are: 0 to 10, 10 to 40, 40 to 100, 100 to 200 cm layer depths. (kg/m2) Instantaneous.
    'Liquid_soil_moisture_content_non-frozen_layer_between_two_depths_below_surface_layer', ... % a 1 by 4 by 224 by 464 matrix. Nonfrozen moisture. Depths are: 0 to 10, 10 to 40, 40 to 100, 100 to 200 cm layer depths. (kg/m2) Instantaneous.
    'Liquid_soil_moisture_content_non-frozen_layer_between_two_depths_below_surface_layer', ... % a 1 by 4 by 224 by 464 matrix. Nonfrozen moisture. Depths are: 0 to 10, 10 to 40, 40 to 100, 100 to 200 cm layer depths. (kg/m2) Instantaneous.
    'Liquid_soil_moisture_content_non-frozen_layer_between_two_depths_below_surface_layer', ... % a 1 by 4 by 224 by 464 matrix. Nonfrozen moisture. Depths are: 0 to 10, 10 to 40, 40 to 100, 100 to 200 cm layer depths. (kg/m2) Instantaneous.
    'Soil_moisture_content_layer_between_two_depths_below_surface_layer', ... % a 1 by 6 by 224 by 464 matrix. Depths are 0 to 10, 10 to 40, 0 to 100, 40 to 100, 0 to 200, 100 to 200. (kg/m2) Instantaneous.
    'Soil_moisture_content_layer_between_two_depths_below_surface_layer', ... % a 1 by 6 by 224 by 464 matrix. Depths are 0 to 10, 10 to 40, 0 to 100, 40 to 100, 0 to 200, 100 to 200. (kg/m2) Instantaneous.
    'Soil_moisture_content_layer_between_two_depths_below_surface_layer', ... % a 1 by 6 by 224 by 464 matrix. Depths are 0 to 10, 10 to 40, 0 to 100, 40 to 100, 0 to 200, 100 to 200. (kg/m2) Instantaneous.
    'Soil_moisture_content_layer_between_two_depths_below_surface_layer', ... % a 1 by 6 by 224 by 464 matrix. Depths are 0 to 10, 10 to 40, 0 to 100, 40 to 100, 0 to 200, 100 to 200. (kg/m2) Instantaneous.
    'Latent_heat_flux_surface_1_Hour_Average', ... % a 1 by 224 by 264. Latent heat flux. W/m2 averaged.
    ... % new ones below
    'Potential_latent_heat_flux_potential_evaporation_surface_1_Hour_Average', ... % Potential latent heat flux (potential evaporation). W/m2 averaged.
    'Sensible_heat_flux_surface_1_Hour_Average', ... % Sensible heat flux. W/m2 averaged.
    'Ground_Heat_Flux_surface_1_Hour_Average', ... % Ground heat flux. W/m2 averaged.
    'Transpiration_surface_1_Hour_Average', ... % Transpiration. W/m2 averaged.
    'Direct_evaporation_from_bare_soil_surface_1_Hour_Average', ... % Direct evaporation from bare soil. W/m2 averaged.
    'Evaporation_surface_1_Hour_Average', ... % Evaporation. kg/m2 accumulated. This is actually evapotranspiration.
    'Canopy_water_evaporation_surface_1_Hour_Average', ... % Canopy water evaporation. W/m2 averaged.
    'Sublimation_evaporation_from_snow_surface_1_Hour_Average', ... % Sublimation (evaporation from snow). W/m2 averaged.
    'Plant_canopy_surface_water_surface', ... % Plant canopy surface water. kg/m2 instantaneous.
    'Snow_depth_surface', ... % Snow depth. m instantaneous.
    'Water_equivalent_of_accumulated_snow_depth_surface', ... % Water equivalent of accumulated snow depth. kg/m2 instantaneous.
    'Snow_cover_surface', ... % Snow cover. fraction instantaneous.
    'Snow_melt_surface_1_Hour_Average', ... % Snow melt. kg/m2 accumulated.
    'Subsurface_runoff_baseflow_surface_1_Hour_Average', ... % Subsurface runoff. kg/m2 accumulated.
    'Surface_runoff_non-infiltrating_surface_1_Hour_Average', ... % Surface runoff (non-infiltrating). kg/m2 accumulated.
    'Downward_shortwave_radiation_flux_surface_1_Hour_Average', ... % Shortwave radiation down. W/m2 averaged.
    'Downward_longwave_radiation_flux_surface_1_Hour_Average', ... % Long wave radiation down. W/m2 averaged.
    'Net_short-wave_radiation_flux_surface_surface_1_Hour_Average', ... % Net shortwave radiation. W/m2 averaged.
    'Net_longwave_radiation_flux_surface_surface_1_Hour_Average', ... % Net long wave radiation. W/m2 averaged.
    'Albedo_surface', ... % Albedo. Percent instantaneous.
    'Leaf_area_index_0-9_surface', ... % Leaf area index (0-9). Unitless instantaneous.
    'Vegetation_surface' ... % Vegetation. Fraction instantaneous.
    }; 
% The units of the above variables (must be in the same order) as varsAvail
varUnits = {'[kg/m2]', '[kg/m2]', ...
    '[kg/m2]', '[kg/m2]', '[kg/m2]', '[kg/m2]',...
    '[kg/m2]', '[kg/m2]', '[kg/m2]', '[kg/m2]',...
    '[W/m2]', ...
    '[W/m2]', '[W/m2]', '[W/m2]', '[W/m2]', '[W/m2]', ...
    '[kg/m2]', '[W/m2]', '[W/m2]', '[kg/m2]', '[m]', ...
    '[kg/m2]', '[frac]', '[kg/m2]', '[kg/m2]', '[kg/m2]', ...
    '[W/m2]', '[W/m2]', '[W/m2]', '[W/m2]', '[percent]', '[-]', '[frac]'};
% The length each variable name will take up in the output file
varLength = [10 10 ...
    9 10 11 12 ...
    9 10 11 12 ...
    9 9 9 9 9 9 ...
    11 8 8 10 9 9 ...
    9 9 9 9 9 9 9 10 10 5 7];
% The format to print each variable
varFmt = {'%5.3e', '%5.3e', ...
    '%7.4f ', '%8.4f ', '%8.4f  ', '%8.4f   ', ... % Liquid soil moisture
    '%7.4f ', '%8.4f ', '%8.4f  ', '%8.4f   ', ... % Total soil moisture
    '%7.2f ', ... % LHF
    '%7.2f ', '%7.2f ', '%+7.2f ', '%7.2f ', '%7.2f ', ... % PotLHF, SHF, G, Transp, EDir
    '%9.2e ', '%6.2f ', '%6.2f ', '%8.4f ', '%7.4f ', ... % ET, ECanopy, Sublim, CanopyH2O, SnowD
    '%7.2f ', '%5.4f  ', '%7.2f ', '%7.2f ', '%7.2f ', ... % SWE, SnowFrac, SnowMelt, Qsub, Qsurf
    '%7.2f ', '%7.2f ', '%7.2f ', '%+8.2f ', ... % SolDn, LWDn, SolNet, LWNet
    '%4.1f     ', '%3.1f ', '%4.3f ', ... % Albedo, LAI, Veg
    };
% -------------------------------------------------------------------------
% Some initial checks

% Make sure data are of the proper type and dimensions
if ~iscellstr(qNames)
    error('Input names must be a cell array of strings.')
end

% If a start date is provided
if isfloat(qStart)
    % Do not append files
    apnd = false;
    % The start datenum is provided in qStart
    dnStart = floor(datenum(qStart));
    if dnStart < datenum(1979,1,3)
        dnStart = datenum(1979,1,3);
    end
% If the string 'apnd' is provided, use the date of the last record
elseif strcmp(qStart,'apnd')
    % Append these files, starting with the datenum after the last line,
    % which should be at hour 2300.
    apnd = true;
    % Open the last file
    fid = fopen([outDir '/' qNames{end} '.txt'], 'r');
    if fid == -1
        error(['There is no file by name ' outDir '/' qNames{end} '.txt to which to append data'])
    end
    % Count the lines in that file
    nLines = 0;
    tline = fgetl(fid);
    while ischar(tline)
        tline = fgetl(fid);
        nLines = nLines+1;
    end
    % Read the last line to get the last date
    frewind(fid)
    lastLine = textscan(fid, '%f', 'headerlines', nLines-1); 
    % Assign the last date in the file as the date to start with
    qStart = [lastLine{1}(1:5)' 0];
    dnStart = ceil(datenum(qStart));
else
    error('qStart must be either a start date or the string ''apnd'' if you wish to append to an existing record')
end
% Make sure start date is not after end date
dnEnd = floor(datenum(qEnd));
if dnEnd > floor(datenum(now))-4
    dnEnd = floor(datenum(now))-4;
end
if dnEnd < dnStart
    error('Start date cannot be after end date')
end
% If longitude is over 180, express as a negative
for qq = 1:length(qLon)
    if qLon(qq)>180
        qLon(qq) = qLon(qq)-360;
    end
end
% Check to be sure each variable requested is available
for vv = 1:length(qVars)
    isPresent = any(strcmp(qVars{vv}, varsAvail));
    if ~isPresent
        error(['You have reqeseted ' qVars{vv} ', a variable that is either misspelled or not yet available via this code.'])
    end
end
% -------------------------------------------------------------------------
% Set up some variables
% Number of sites requested
nSites = length(qNames);
% Number of variables requested
nVars = length(qVars);
% Initialize vectors for lat/lon idcs
latIdcs = nan(nSites,1);
lonIdcs = nan(nSites,1);
% The years and days of years of the query
qDatenums = dnStart:dnEnd;
% Initialize a day of year vector
qDoy = nan(1,length(qDatenums));
% The years, months, and days of the query
[qYears, qMonths, qDays] = datevec(qDatenums);
% The unique years in the query
qUniqueYears = unique(qYears);
% For each unique year
for yy = 1:length(qUniqueYears)
    % Find all the indices where qYears matches this unique year
    sameYrIdcs = find(qYears == qUniqueYears(yy));
    % Subtract the datenum of day zero of the year from the datenums in the
    % year to get day of year
    qDoy(sameYrIdcs) = datenum(qDatenums(sameYrIdcs)) - datenum(qUniqueYears(yy),0,0);
end % Loop through each unique year
% Convert qYears to string
qYearStr = num2str(qYears');
% Convert qMonths to string
qMonthStr = num2str(qMonths', '%02d');
% Convert qDays to string
qDayStr = num2str(qDays', '%02d');
% Convert qDoy to string
qDoyStr = num2str(qDoy', '%03d');
% Create hour strings
qHourStr = num2str((0:100:2300)','%04d');

% Date formats (for year month day hour minute)
dateFmt = '%04d   %02d    %02d  %02d   %02d     ';

% The directory where nldas forcings are held
ftpBaseDir = '/data/s4pa/NLDAS/NLDAS_NOAH0125_H.002/';
% The local directory where nldas forcings will be placed
localBaseDir = [pwd '/data'];
% If there is already a directory named data in the working directory, do not continue because
% it will be deleted at the end of this script and I don't want to delete
% anything that this script itself did not create.
if exist(localBaseDir, 'dir') == 7
    error('%s\n%s\n%s','This function requires use of the local directory named ''./data.'',', ...
        'and you already have a directory with that name. We are stopping here',...
        'because the end of this function will DELETE ./data and all files within it.')
end
% The beginning of the file name of forcings
ftpBaseFn = 'NLDAS_NOAH0125_H.A';
% The end of the file name of forcings
ftpEndFn = '.002.grb';

% Strings for lat and lon
% These are 1-d matrices (1 by 224 or 464).
latStr = 'lat'; % 224x1. Center of 1/8 degree pixel
lonStr = 'lon'; % 464x1. Center of 1/8 degree pixel

% -------------------------------------------------------------------------
% Create a directory to hold output data
if exist(outDir, 'dir') ~= 7
    disp(['Making an output directory here: ' outDir])
    mkdir(outDir);
end

% -------------------------------------------------------------------------
% Set up ftp connection and download the first file to get its lat/lon data
% The Nasa host that holds nldas forcings
nasaHost = 'hydro1.sci.gsfc.nasa.gov';
% Create an ftp object (open the connection)
ftpObj = ftp(nasaHost);
% The file name of the first date requested
qFileName = [ftpBaseDir qYearStr(1,:) '/' qDoyStr(1,:) '/' ftpBaseFn qYearStr(1,:) qMonthStr(1,:) qDayStr(1,:) '.' qHourStr(1,:) ftpEndFn];
% Get that first file
disp(['Getting ' qFileName '...'])
localFileName = mget(ftpObj,qFileName);
% Create ncgeodataset object
geo = ncdataset(localFileName{1});
% To list the variables available: geo.variables
% Extract lat and lon from the grib file
lat = geo.data(latStr);
lon = geo.data(lonStr);
% Size of these vectors
nLat = length(lat);
nLon = length(lon);

% -------------------------------------------------------------------------
% Set up output files. Write headers with lat/lon data in them. Get lat/lon
% idcs.

% Initialize a vector of file IDs
fid = nan(1,nSites);

% Set up header strings
% Initialize a string to hold variable names
namesStrAll = [];
% Initialize a string to hold the units
unitsStrAll = [];
% Initialize a string to hold the variable's format
fmtStrAll = [];
% Initialize a cell array to hold the varaibles' names
varStrAll = cell(1,nVars);
% Initialize a vector to hold the indices of each variable
idcsAll = nan(1,nVars);

% Loop through the variables requested
for vv = 1:nVars
    % This variable is
    thisVar = qVars{vv};
    % The index of this requested variable (refers to all that are
    % available)
    idcsAll(vv) = find(strcmp(varsAvail, thisVar));
    % This variable is referred to in the grib file as
    varStrAll{vv} = varsStrings{idcsAll(vv)};
    % The units of this variable
    thisUnits = varUnits{idcsAll(vv)};
    % The format for this varaible's data
    fmtStrAll = [fmtStrAll varFmt{idcsAll(vv)} ' '];
    % The number of spaces to add after the var name
    nSpaces(1) = varLength(idcsAll(vv)) - length(thisVar);
    % The number of spaces to add after the units
    nSpaces(2) = varLength(idcsAll(vv)) - length(thisUnits);
    % Fill in the names string
    namesStrAll = [namesStrAll thisVar blanks(nSpaces(1))];
    % Fill in the units string
    unitsStrAll = [unitsStrAll thisUnits blanks(nSpaces(2))];
end
% Add a line return to the end of the format string
fmtStrAll = [fmtStrAll '\n'];

% Loop through each site
for ss = 1:nSites
    % Get lat/lon idcs
    [latDiff latIdcs(ss)] = min(abs(lat-qLat(ss)));
    [lonDiff lonIdcs(ss)] = min(abs(lon-qLon(ss)));
    % Get the rounded lat and lon. This is the center point of the NLDAS
    % cell.
    nearestLat = lat(latIdcs(ss));
    nearestLon = lon(lonIdcs(ss));
    % If the difference between the queried and actual lat or lon is too
    % big, display a warning
    if latDiff>0.125 || lonDiff>0.125
        warning('%s site does not have a nearby nldas cell. \nQueried lat/lon: %f %f\nClosest NLDAS lat/lon: %f %f',...
            qNames{ss}, qLat(ss), qLon(ss), nearestLat, nearestLon)
    end
    % The name for this output file
    outFile = [outDir '/' qNames{ss} '.txt'];
    % If we are not appending data, write a header
    if ~apnd
        % Open the output file for the first time
        fid(ss) = fopen(outFile,'w');
        % Print header lines to the file
        fprintf(fid(ss),['%% Site: ' qNames{ss} '\n']);
        fprintf(fid(ss),['%% Site lat/lon: ' num2str([qLat(ss) qLon(ss)]) '\n']);
        fprintf(fid(ss),['%% Closest NLDAS pixel (1/8 degree) center: ' num2str([nearestLat nearestLon]) '\n']);
        fprintf(fid(ss),['%% File created on ' date '.\n']);
        fprintf(fid(ss),['%% Date (UTC)                 ' namesStrAll '\n']);
        fprintf(fid(ss),['%% year month day hour minute ' unitsStrAll '\n']);
    else 
        % If we are appending data, open the output file with append
        % permission only
        fid(ss) = fopen(outFile, 'a');
    end % If we are appending data
end % Loop through each site

% -------------------------------------------------------------------------
% Loop through each day in the record
for dd = 1:length(qDatenums)
    % Location of this day's data on the local machine
    localDir = [pwd ftpBaseDir qYearStr(dd,:) '/' qDoyStr(dd,:)];
    % Loop through each hour of the day
    for hh = 1:24
        % Create strings and such to define where the files are located on the host
        qFileName = [ftpBaseDir qYearStr(dd,:) '/' qDoyStr(dd,:) '/' ftpBaseFn qYearStr(dd,:) qMonthStr(dd,:) qDayStr(dd,:) '.' qHourStr(hh,:) ftpEndFn];
        % Get the file from Nasa's server
        disp(['Getting ' qFileName '...'])
        localFileName = mget(ftpObj,qFileName);
        % Create ncgeodataset object
        geo = ncdataset(localFileName{1});
        % Initialize a vector to hold the timestep's data (all variables) for entire domain
        domainData = nan(nLat, nLon, nVars);
        % Get the data for each variable
        for vv = 1:nVars
            varData = squeeze(geo.data(varStrAll{(vv)}));
            % If these data are liquid soil moisture content or total soil
            % moisture content, the varaible will be 3-dimentional
            if ndims(varData) == 3
                % Specify which index corresponds to the requested layer
                switch qVars{vv}
                    case 'LSM_0_10'
                        dimIdx = 1;
                    case 'LSM_10_40'
                        dimIdx = 2;
                    case 'LSM_40_100'
                        dimIdx = 3;
                    case 'LSM_100_200'
                        dimIdx = 4;
                    case 'TSM_0_10'
                        dimIdx = 1;
                    case 'TSM_10_40'
                        dimIdx = 2;
                    case 'TSM_40_100'
                        dimIdx = 4;
                    case 'TSM_100_200'
                        dimIdx = 6;
                end
%                 qVars{vv}
%                 dimIdx
                % Save only the requested layer
                domainData(:,:,(vv)) = varData(dimIdx,:,:);
            else % This is just a 2-d matrix
                domainData(:,:,(vv)) = varData;
            end % If this is a 3-d matrix
        end % Loop through each variable
        % Loop through each site 
        for ss = 1:nSites
            % Extract the point data at each site from the domain data
            siteData = squeeze(domainData(latIdcs(ss), lonIdcs(ss), :));
            % Print the data. Hour is hh-1 because hours are listed from
            % 00:00 to 23:00.
            fprintf(fid(ss), [dateFmt fmtStrAll], [qYears(dd) qMonths(dd) qDays(dd) (hh-1) 0 siteData']);
        end % loop through each site
    end % Loop through each hour of the day
    % Delete this day's directory and all data within
    rmdir(localDir, 's')
end % Loop through each day in the record
% -------------------------------------------------------------------------
% Clean up and close files
% Delete all data stored in this session
disp(['Cleaning up...'])
rmdir(localBaseDir, 's')
% Loop through each site and close the files
for ss = 1:nSites
    outFile = qNames{ss};
    fclose(fid(ss));
end
% Close the ftp connection
close(ftpObj)

end
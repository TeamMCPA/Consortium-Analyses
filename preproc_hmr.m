function new_filename = preproc_hmr(filename,params_file)
%PREPROC_HMR  Use Homer2 functions to preprocess .nirs data file according
%to standard Princeton Baby Lab pipeline
%
%   NEW_FILENAME = PREPROC_HMR( FILENAME, PARAMS_FILE ) applys the 
%   parameters defined in PARAMS_FILE to the .nirs file specified in 
%   FILENAME. Name of the new file is returned as a string.
%
%   The following preprocessing procedures are performed (in order):
%   1. Convert intensity (d) to optical density (dod)
%   2. Detect motion artifacts in dod and apply spline correction
%   3. Apply wavelet correction to dod for motion artifacts
%   4. Bandpass filter dod data
%   5. Convert optical density (dod) to concentrations (dc)
%   6. Remove bad channels based on max/min values for dc
%   7. Perform block averaging on all channels (except bad channels)
%
%   bdz 05 oct 2018

%% TODO:
% - Need to revise to pull command sequence and parameters from .cfg file
% - Gracefully handle exceptions for missing files, missing parameters
% - Test batching

%% 
% these were originally hardcoded for adult Shimadzu
% nChannels = 139;
% fs = 8.1301;

%% Handle batching of filename with cell array
if iscell(filename)
    new_filename = cell(length(filename),1);
    for filenum = 1:length(filename)
        new_filename(filenum) = preproc_hmr(filename{filenum},params_file);
    end
end

%% Import nirs file and preprocessing parameters file
nirs_dat = load(filename,'-mat');
params = load(params_file,'-mat');

%% Preapre some additional needed parameters based on provided data
params.nWav = length(params.ppm);
params.nChannels = size(nirs_dat.d,2) / nWav;
params.fs = mean(diff(nirs_dat.t));
if ~exist('nirs_dat.tIncMan','var')
    nirs_dat.tIncMan = ones(size(t));
end

%% Perform preprocessing

% Convert intensity data to optical density
procResult.dod = hmrIntensity2OD(nirs_data.d);

% Perform motion artifact detection
[tInc,tIncCh] = hmrMotionArtifactByChannel(...
    nirs_dat.d, params.fs, nirs_dat.SD, nirs_dat.tIncMan, params.tMotion, params.tMask, params.STDEVthresh, params.AMPthresh);

% Apply spline correction to motion artifacts
nirs_dat.procResult.dodSpline = hmrMotionCorrectSpline(...
    nirs_dat.procResult.dod,nirs_dat.t,nirs_dat.SD,tIncCh,params.p);

% Apply wavelet correction for motion artifacts
nirs_dat.procResult.dodWavelet = hmrMotionCorrectWavelet(...
    nirs_dat.procResult.dodSpline,nirs_dat.SD,params.IQR);

% Bandpass filter the optical density data
nirs_dat.procResult.dodBP = hmrBandpassFilt(...
    nirs_dat.procResult.dodWavelet,params.fs,params.hpf,params.lpf);

% Convert corrected/filtered optical density data to concentration data
nirs_dat.procResult.dc = hmrOD2Conc(...
    nirs_dat.procResult.dodBP,nirs_dat.SD,params.ppf);

% Generate a list of bad channels and NaN-out their data
nirs_dat.badCh = any(squeeze(any(nirs_dat.procResult.dc < params.dRange(1) | nirs_dat.procResult.dc > params.dRange(2))));
nirs_dat.procResult.dcFix(:,:,badCh) = nan;

% Perform block-averaging on the dcFix data (motion-corrected, filtered,
% and bad channels removed)
[...
    nirs_dat.procResult.dcAvg,...
    nirs_dat.procResult.dcStd,...
    nirs_dat.procResult.tHRF,...
    nirs_dat.procResult.nTrials,...
    nirs_dat.procResult.dcSum2,...
    nirs_dat.procResult.dcTrials...
    ]...
    = hmrBlockAvg( procResult.dcFix, nirs_dat.sFix, nirs_dat.tFix, nirs_dat.tRange );

%% Save the preprocessed data out to a new file and return the filename
[filePath,fileName] = fileparts(filename);
new_filename = fullfile(filePath,strcat(fileName,'_proc.nirs'));
save(new_filename,'-struct','nirs_dat');
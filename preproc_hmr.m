function new_filename = preproc_hmr(filename,params_file,overwrite)
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

if ~exist('overwrite','var')
    overwrite = false;
    elseif ~islogical(overwrite)
        if isnumeric(overwrite), overwrite = logical(overwrite);
        else overwrite = false;
        end
end

%% Handle batching of filename with cell array
if iscell(filename)
    new_filename = cell(length(filename),1);
    for filenum = 1:length(filename)
        new_filename{filenum} = preproc_hmr(filename{filenum},params_file,overwrite);
    end
    return
end

%% Import nirs file and preprocessing parameters file
nirs_dat = load(filename,'-mat');
fprintf('Loaded file: %s... ',filename);
params = load(params_file,'-mat');

%% Preapre some additional needed parameters based on provided data
params.nWav = length(params.ppf);
params.nChannels = size(nirs_dat.d,2) / params.nWav;
params.fs = 1/mean(diff(nirs_dat.t));
if ~exist('nirs_dat.tIncMan','var')
    nirs_dat.tIncMan = ones(size(nirs_dat.t));
end

%% Perform preprocessing

if ~exist('nirs_dat.procResult','var')
    nirs_dat.procResult = struct();
end

% 1. Convert intensity data to optical density
nirs_dat.procResult.dod = hmrIntensity2OD(nirs_dat.d);

% 2. Perform motion artifact detection
[tInc,tIncCh] = hmrMotionArtifactByChannel(...
    nirs_dat.d, params.fs, nirs_dat.SD, nirs_dat.tIncMan, params.tMotion, params.tMask, params.STDEVthresh, params.AMPthresh);

% Apply spline correction to motion artifacts
if isfield(params,'p')
nirs_dat.procResult.dodSpline = hmrMotionCorrectSpline(...
    nirs_dat.procResult.dod,nirs_dat.t,nirs_dat.SD,tIncCh,params.p);
else
    nirs_dat.procResult.dodSpline = nirs_dat.procResult.dod;
end
% 3. Apply wavelet correction for motion artifacts
if isfield(params,'IQR')
    try
        nirs_dat.procResult.dodWavelet = hmrMotionCorrectWavelet(...
            nirs_dat.procResult.dodSpline,nirs_dat.SD,params.IQR);
    catch
        Just move the Spline output to Wavelet
        disp('Warning Wavelet Motion Correction failed.')
        nirs_dat.procResult.dodWavelet = nirs_dat.procResult.dodSpline;
    end
else
    nirs_dat.procResult.dodWavelet = nirs_dat.procResult.dodSpline;
end

% 4. Bandpass filter the optical density data
nirs_dat.procResult.dodBP = hmrBandpassFilt(...
    nirs_dat.procResult.dodWavelet,params.fs,params.hpf,params.lpf);

% 5. Convert corrected/filtered optical density data to concentration data
nirs_dat.procResult.dc = hmrOD2Conc(...
    nirs_dat.procResult.dodBP,nirs_dat.SD,params.ppf);

% Generate a list of bad channels and NaN-out their data
nirs_dat.badCh = any(squeeze(any(nirs_dat.procResult.dc < params.dRange(1) | nirs_dat.procResult.dc > params.dRange(2))));
nirs_dat.procResult.dcFix = nirs_dat.procResult.dc;
nirs_dat.procResult.dcFix(:,:,nirs_dat.badCh) = nan;
if isfield(params,'dropChans')
    if params.dropChans
        nirs_dat.procResult.dc = nirs_dat.procResult.dcFix;
    end
end

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
    = hmrBlockAvg( nirs_dat.procResult.dcFix, nirs_dat.s, nirs_dat.t, params.tRange );

fprintf('Finished preproc... ')

%% Clean up the struct before writing it out.

% aux field is obligatory in Homer2, so create it if needed
if ~isfield(nirs_dat,'aux')
    % This will take the binary s matrix and make it a single vector of
    % length t (time) with integer values for the triggers. If two triggers
    % co-occur, the aux value will be their sum, so it's still better to
    % use the s matrix when possible.
    nirs_dat.aux = nirs_dat.s * [1:size(nirs_dat.s,2)]';
end

% remove extraneous processing steps
nirs_dat.procResult= rmfield(nirs_dat.procResult,'dodSpline');
nirs_dat.procResult= rmfield(nirs_dat.procResult,'dodWavelet');
nirs_dat.procResult= rmfield(nirs_dat.procResult,'dodBP');
nirs_dat.procResult= rmfield(nirs_dat.procResult,'dcFix');

%% Save the preprocessed data out to a new file and return the filename
[filePath,fileName] = fileparts(filename);
if overwrite, new_filename = filename;
else, new_filename = fullfile(filePath,strcat(fileName,'_proc.nirs'));
end
fprintf('Saved.\n')
save(new_filename,'-struct','nirs_dat');
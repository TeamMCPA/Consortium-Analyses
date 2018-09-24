function decodingPreproc(fullFileName)
nChannels = 139;
fs = 8.1301;

load(fullFileName,'-mat')

procResult.dod = hmrIntensity2OD(d);

tIncMan = ones(size(t));
tMotion = 1;
tMask = 1;
STDEVthresh = 50;
AMPthresh = .5;
[tInc,tIncCh] = hmrMotionArtifactByChannel(d, fs, SD, tIncMan, tMotion, tMask, STDEVthresh, AMPthresh);
p = .99;
procResult.dodSpline = hmrMotionCorrectSpline(procResult.dod,t,SD,tIncCh,p);
IQR = 1;
procResult.dodWavelet = hmrMotionCorrectWavelet(procResult.dodSpline,SD,IQR);
hpf = 0.01;
lpf = 1;
procResult.dodBP = hmrBandpassFilt(procResult.dodWavelet,fs,hpf,lpf);
ppf = [6,6,6];
procResult.dc = hmrOD2Conc(procResult.dodBP,SD,ppf);

dRange = [-5e-5,5e-5];
badCh = any(squeeze(any(procResult.dc < dRange(1) | procResult.dc > dRange(2))));
procResult.dcFix(:,:,badCh) = nan;

[procResult.dcAvg, procResult.dcStd, procResult.tHRF, procResult.nTrials, procResult.dcSum2, procResult.dcTrials] = hmrBlockAvg( procResult.dcFix, sFix, tFix, trange );

[filePath,fileName] = fileparts(fullFileName);
save(fullfile(filePath,strcat(fileName,'_proc.nirs')),'d','fs','s','SD','t','procResult','badCh','-mat')




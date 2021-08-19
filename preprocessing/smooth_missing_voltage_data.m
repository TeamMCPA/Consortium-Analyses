
function mask = smooth_missing_voltage_data(tMask, d, fs)
%% Finds times where probes were not in contact with head and sets the 
% half second before and after to 0

buffer = round(fs*tMask);
halfBuffer = round(buffer/2);
mask = [];

for chan = 1:size(d,2)
    
    bad_idx = find(d(:,chan) == 0);
    correctDat = ones(length(d),1);
    
    if ~isempty(bad_idx)
        for sample = 1:length(bad_idx)
            if bad_idx(sample)-halfBuffer <= 0
                start_correction = 1;
            else
                start_correction = bad_idx(sample)-halfBuffer;
            end

            if bad_idx(sample)+halfBuffer > length(d)
                end_correction = correctDat(end);
            else
                end_correction = bad_idx(sample)+halfBuffer;
            end
            correctDat(start_correction:end_correction) = 0;
        end
    end
    mask = [mask correctDat];
end
end


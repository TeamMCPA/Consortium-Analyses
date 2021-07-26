function quick_plot_chans(MCP_data,show_grid)
%% simple plotter for probe layouts in MCP data
% Input one (singleton) MCP struct (e.g., one subject), and click on any
% detector in the figure to view the Oxy fNIRS data for the channels
% associated with that detector. Set show_grid to 0 to prevent channel
% lines from being drawn. Lines will appear for each detector you click
% instead.

if ~exist('show_grid','var'), show_grid=1; end

%% Set up the SD struct that you're going to work with
if length(MCP_data)>1
    warning('%g elements found in MCP struct. Visualizing only the first.',length(MCP_data));
    MCP_data = MCP_data(1);
end
% For MCP
SD = MCP_data.Experiment.Probe_arrays.Geometry;
% % For *.nirs
% SD = SD;

%% DELETE THIS CODE BLOCK AFTER REPAIRING THE MCP DATA
% It manually switches the sources and detectors, whose positions were
% incorrectly coded in the SD struct originally used to import the data.
warning('HARD-CODED SWAP OF SOURCE AND DETECTOR POSITIONS');
warning('This swap is a correction for incorrect SD data used to build the Consortium data. Remove for new data.');
SD2 = SD;
SD2.SrcPos = SD.DetPos;
SD2.DetPos = SD.SrcPos;
SD = SD2;

%% Initiate the figure
chanfig = figure;
chanplot = axes;
hold(chanplot,'on');
xlabel('X axis');
ylabel('Y axis');
zlabel('Z axis');

%% Plot channel positions
% MeasList contains:
% [ source_number, detector_number, placeholder, wavelength]
chan_locs = SD.MeasList;
% Since wavelengths repeat the source-detector pairings, just look at the
% pairings for the first wavelength
chan_locs = chan_locs(chan_locs(:,4)==1,:);

% Cycle through the list of channels and plot
for chan = 1:size(chan_locs,1)
    
    % First column is source number
    src=chan_locs(chan,1);
    % Second column is detector number
    det=chan_locs(chan,2);
    
    % If the show_grid feature is active, preprint the channel lines
    if show_grid
        % Each line is drawn from the source position to the detector position.
        % Recall that plot3 requires three separate inputs for the three
        % dimensions, so src-det pair are listed in x, then pair listed in y,
        % and then pair listed in z. Draw black line.
        p(1) = plot3(chanplot,...
            [SD.SrcPos(src,1);SD.DetPos(det,1)],...
            [SD.SrcPos(src,2);SD.DetPos(det,2)],...
            [SD.SrcPos(src,3);SD.DetPos(det,3)],...
            'k--','LineWidth',1);
    end
    
    % Add channel number labels. These are the simply the index of the
    % channel's position in the MeasList, which should also correspond to
    % the channel number in the Transformation Matrix. The text is printed
    % at the mean point between the source & detector.
    p(2) = text(chanplot,...
        mean([SD.SrcPos(src,1);SD.DetPos(det,1)]),...
        mean([SD.SrcPos(src,2);SD.DetPos(det,2)]),...
        mean([SD.SrcPos(src,3);SD.DetPos(det,3)]),...
        num2str(chan));
end

%% Plot probe positions
% Plot the positions of the sources in 3d space
p(3) = plot3(chanplot,...
    SD.SrcPos(:,1),...
    SD.SrcPos(:,2),...
    SD.SrcPos(:,3),...
    'r*','LineWidth',2,'MarkerSize',8,'MarkerFaceColor','r');
% Plot the positions of the detectors in 3d space
p(4) = plot3(chanplot,...
    SD.DetPos(:,1),...
    SD.DetPos(:,2),...
    SD.DetPos(:,3),...
    'bo','LineWidth',2,'MarkerSize',8,'MarkerFaceColor','b');
grid on;
title('Click on Detectors to view time series data')


% Hold aspect ratio
daspect([1 1 1])

% Initialize a global variable to store the detector number
global clickDetectorID
clickDetectorID = NaN;

% Create user interface for detectors
p(4).PickableParts = 'all';
p(4).ButtonDownFcn = @showChanData;

% Monitor for feedback as long as the figure is open
while ishghandle(chanfig)
    pause(.01)
    if ~isnan(clickDetectorID)
        chansToPlot = find(chan_locs(:,2)==clickDetectorID);
        
        % If the lines aren't already printed, show them now
        if ~show_grid
            for c = 1:size(chansToPlot,1)
                p(1) = plot3(chanplot,...
                    [SD.SrcPos(chan_locs(chansToPlot(c),1),1);SD.DetPos(clickDetectorID,1)],...
                    [SD.SrcPos(chan_locs(chansToPlot(c),1),2);SD.DetPos(clickDetectorID,2)],...
                    [SD.SrcPos(chan_locs(chansToPlot(c),1),3);SD.DetPos(clickDetectorID,3)],...
                    'k-','LineWidth',1);
            end
        end
        
        tsfig = figure;
        hold on;
        time_vec = cat(1,MCP_data.Experiment.Runs(:).Time);
        plot(MCP_data.fNIRS_Data.Hb_data.Oxy(:,chansToPlot));
        %Add lines for session cut-offs. Doesn't work yet.
        %plot(find(diff(time_vec)<0),ones(size(find(diff(time_vec)<0))));
        set(gca,'XTick',1:1000:length(time_vec));
        set(gca,'XTickLabel',round(time_vec(1:1000:end)));
        title(['Detector ' num2str(clickDetectorID)]);
        chanleg = legend(num2str(chansToPlot));
        chanleg.Title.String = 'Channels';
        drawnow update;
        hold off;
        clickDetectorID=NaN;
    end
end

% Clear the global variable before leaving
clear clickDetectorID;




    function showChanData(hObj, event)
        
        x = hObj.XData;
        y = hObj.YData;
        z = hObj.ZData;
        
        pt = event.IntersectionPoint;
        coords = [x(:),y(:),z(:)];
        dist = pdist2(pt,coords);
        [~, minIdx] = min(dist);
        coordsSelect = coords(minIdx,:);
        
        fprintf('Det %g: X: %0.1f, Y: %0.1f, Z: %0.1f\n', minIdx, coordsSelect);
        
        %global clickDetectorID
        clickDetectorID = minIdx;
        
    end

end
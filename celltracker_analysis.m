clear all
close all
frames=60;
dt=2*15*60; %time interval in seconds
pixelsize=0.55; %pix/micron
path="C:\Users\G-mo10\Desktop\Test images nuclei Agata\CellTracker results\";
filename = "TracksCoordinates_C0506A1R3_007.mat";
filename=path+filename;
load(filename);
CELLS=pos(end:end);
matrix_x=nan(frames,CELLS); matrix_y=matrix_x;
%%%%%%%Fills matrices of positions: cell 0 | cell 1 | cell 2 | etc.:%%%%%%
for jj=1:CELLS
    kk=1;
for ii=1:2:length(pos)
    if pos(ii,4)==jj
        matrix_x(kk,jj)=pos(ii,1)/pixelsize; %in microns
        matrix_y(kk,jj)=pos(ii,2)/pixelsize; %in microns
        kk=kk+1;
    end
end
end
clear ii jj kk
%%%%%%%Plots positions:%%%%%%%%%
% sz=100; a = (1:frames)'; b = num2str(flipud(a)); c = cellstr(b);
% for ll=1:CELLS
% scatter(matrix_x(:,ll),matrix_y(:,ll),sz,'filled')
% text(matrix_x(:,ll)-3.5, matrix_y(:,ll)+.5, c,'Color','white','FontWeight','bold')
% hold on
% plot(matrix_x(:,ll),matrix_y(:,ll))
% hold on
% %axis([0 200 0 120])
% end
% title('Cell trajectories')
% xlabel('Microns'); ylabel('Microns')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%Separate tracks between short and long:
%Get the length of each track:
track_lengths=[];
for ll=1:CELLS
single_track=length(matrix_x(:,ll))-sum(isnan(matrix_x(:,ll)));
track_lengths=[track_lengths single_track];
end
clear sincle_track ll
mean_track=mean(track_lengths);
% histogram(track_lengths)
% title('Length in frames of tracks')
% xlabel('Length (# of frames)'); ylabel('Count')

%%%%%%%%%%%%%%%%%%%% Velocities: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vel_x=diff(matrix_x)/dt;vel_y=diff(matrix_y)/dt; %microns/s
%vel_x=matrix_x*3600;vel_y=matrix_y*3600;        %microns/h
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%% MSD: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SPACE_UNITS = 'µm';
TIME_UNITS = 's';
tracks = cell(CELLS, 1);
for mm = 1 : CELLS
    % Time
    time = (0 : frames-1)' * dt;
    % Positions
    X=[matrix_x(:,mm) matrix_y(:,mm)];
    % Store
    tracks{mm} = [time X];
end
clear mm X time
%Initiate msd analyzer (2 means 2D):
ma = msdanalyzer(2, SPACE_UNITS, TIME_UNITS);
ma = ma.addAll(tracks);
%Calculates msd:
ma = ma.computeMSD;
meansqdis=ma.msd;
%Plot weighted average over all MSD curves with errorbar:
figure
ma.plotMeanMSD(gca, true)
mmsd = ma.getMeanMSD;
t = mmsd(:,1);
x = mmsd(:,2);
dx = mmsd(:,3) ./ sqrt(mmsd(:,4));
errorbar(t, x, dx, 'k')
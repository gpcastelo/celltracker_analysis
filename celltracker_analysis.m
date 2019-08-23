clear all
close all
frames=60;
dt=2*15*60; %time interval in seconds
pixelsize=0.55; %pix/micron
path="C:\Users\G-mo10\Desktop\Test images nuclei Agata\CellTracker results\";
filename = "TracksCoordinates_C0506A1R3_007.mat";
filename=path+filename;
load(filename); save_pos=pos;
CELLS=pos(end:end);
matrix_x=nan(frames,CELLS); matrix_y=matrix_x; dummy_frames=nan(frames,CELLS);
%%%%%%%Fills matrices of positions:         cell 1 | cell 2 | cell 3 | etc.:%%%%%%
%%%%%%% and frames                  frame 1
%                                   frame 2
%                                   frame 3 etc.
for jj=1:CELLS
    kk=1;
    for ii=1:length(pos)
        if pos(ii,4)==jj && mod(pos(ii,3),2)==1
            matrix_x(ceil(pos(ii,3)/2),jj)=pos(ii,1)/pixelsize; %in microns
            matrix_y(ceil(pos(ii,3)/2),jj)=pos(ii,2)/pixelsize; %in microns
            dummy_frames(kk,jj)=pos(ii,3);
            kk=kk+1;
        elseif pos(ii,4)==jj+1
            break
        end
    end
end
clear ii jj kk
%%%%%%%%%%Now for the primary and secondary cells matrix:
matrix_x_primary=[]; matrix_y_primary=[];
matrix_x_daughters=[]; matrix_y_daughters=[];

for jj=1:CELLS
    if dummy_frames(1,jj)<=16                               %Primary cells 
                                                            %(first frame 
                                                            %<=8*2)
        matrix_x_primary=[matrix_x_primary matrix_x(:,jj)];
        matrix_y_primary=[matrix_y_primary matrix_y(:,jj)];
    else                                                    %Daughter cells
        matrix_x_daughters=[matrix_x_daughters matrix_x(:,jj)];
        matrix_y_daughters=[matrix_y_daughters matrix_y(:,jj)];
     end
end
clear jj
%%%%%%%Polar coordinates: rho=sqrt(x^2+y^2); phi=atan2(y,x)
matrix_rho=sqrt(matrix_x.^2+matrix_y.^2); %microns
[len_x,wid_x]=size(matrix_rho); 
matrix_phi=nan(frames,CELLS);
for ii=1:len_x
    for jj=1:wid_x
        matrix_phi(ii,jj)=atan2(matrix_y(ii,jj),matrix_x(ii,jj)); %rad
    end
end
clear ii jj len_x wid_x
%%%%%%%Fills matrices of positions: frame 1 | frame 2 | frame 3 | etc.:%%%%%%
% matrix_x_for_frameav=nan(CELLS,2*frames); matrix_y_for_frameav=matrix_x_for_frameav;
% for jj=1:2:frames
%     kk=1;
% for ii=1:length(pos)
%     if pos(ii,3)==jj && mod(pos(ii,3),2)==1
%         matrix_x_for_frameav(kk,jj)=pos(ii,1)/pixelsize; %in microns
%         matrix_y_for_frameav(kk,jj)=pos(ii,2)/pixelsize; %in microns
%         kk=kk+1;
%     end
% end
% end
% clear ii jj kk

%%%%%%%Separate tracks between short and long:
%Get the length of each track:
track_lengths=[];
for ll=1:CELLS
single_track=length(matrix_x(:,ll))-sum(isnan(matrix_x(:,ll)));
track_lengths=[track_lengths single_track];
end
clear single_track ll
disp('Average length of track in number of frames:')
mean_track=mean(track_lengths)
disp('Average length of track in seconds:')
mean_track_seconds=mean_track*dt
figure
histogram(track_lengths)
title('Length in frames of tracks')
xlabel('Length (# of frames)'); ylabel('Count')

%%%%%%%%%%%%%%%%%%%% Velocities: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%disp('Velocities in microns/s.')
vel_x=diff(matrix_x)/dt;vel_y=diff(matrix_y)/dt; %microns/s
%vel_x=matrix_x*3600;vel_y=matrix_y*3600;        %microns/h
% % % % % % %Magnitude:
vel_mag=sqrt(vel_x.^2+vel_y.^2);
% % % % % % %Average:
%disp('Average velocity for each cell:')
vel_cell_avg=nanmean(vel_mag); %average velocity for each cell
% figure
% plot(vel_cell_avg,'b')
% title('Average velocity for each cell')
% xlabel('Cell'); ylabel('Velocity (\mu m/s)')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % %For each frame of all cells:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% vel_x_frame=diff(matrix_x_for_frameav)/dt; vel_y_frame=diff(matrix_y_for_frameav)/dt;
% vel_frame=sqrt(vel_x_frame.^2+vel_y_frame.^2);
% [~,widmat]=size(vel_frame);
% vel_frame_avg=[];
% for aa=1:2:widmat
%     vel_avg_single_frame=nanmean(vel_frame(:,aa));
%     vel_frame_avg=[vel_frame_avg vel_avg_single_frame];
% end
% clear aa widmat vel_avg_single_frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Create matrix of instant velocities:
vel_x_frame=[]; vel_y_frame=[];
[len_x,~]=size(matrix_x);
for aa=2:len_x
    instant_vel_x=(matrix_x(aa,:)-matrix_x(aa-1,:))/dt;
    vel_x_frame=[vel_x_frame; instant_vel_x];
    instant_vel_y=(matrix_y(aa,:)-matrix_y(aa-1,:))/dt;
    vel_y_frame=[vel_y_frame; instant_vel_y];
end
clear aa len_x instant_vel_x instant_vel_y
vel_all=sqrt(vel_x_frame.^2+vel_y_frame.^2);
%%From the matrix of instant velocities, create vector of avg velocity:
avg_vel_all=[];
[len_x,~]=size(vel_all);
for aa=1:len_x
    dummy_vel=nanmean(vel_all(aa,:));
    avg_vel_all=[avg_vel_all dummy_vel];
end
clear aa len_x dummy_vel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot:
err=std(avg_vel_all)*ones(size(avg_vel_all));
figure
ebar=errorbar(1:59, avg_vel_all, err, 'LineStyle','none');
ebar.Color = 'k';
hold on
scatter(1:59,avg_vel_all,'r','filled')
title('Average velocity for each frame: All cells')
xlabel('Frame'); ylabel('Velocity (\mu m/s)')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% For each frame of primary cells:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% vel_x_primary=diff(matrix_x_primary')/dt; vel_y_primary=diff(matrix_y_primary')/dt;
% vel_x_primary=vel_x_primary'; vel_y_primary=vel_y_primary';
% vel_frame_primary=sqrt(vel_x_frame.^2+vel_y_frame.^2);
% [~,widmat]=size(vel_frame_primary); 
% vel_primary_avg=[];
% for aa=1:2:widmat
%     vel_avg_single_frame=nanmean(vel_frame_primary(:,aa));
%     vel_primary_avg=[vel_primary_avg vel_avg_single_frame];
% end
% clear aa widmat vel_avg_single_frame
% figure
% scatter(1:60,vel_primary_avg,'b')
% title('Average velocity for each frame: Primary cells')
% xlabel('Frame'); ylabel('Velocity (\mu m/s)')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dummy_vector_for_avg_vel=[];
% for ii=1:frames-1
%     for jj=1:CELLS
%         dummy_vector_for_avg_vel=[dummy_vector_for_avg_vel vel_mag(ii,jj)];
%     end
% end
% disp('Average velocity of all cells:')
% vel_total_avg=nanmean(dummy_vector_for_avg_vel) %average velocity of all cells
clear ii jj dummy_vector_for_avg_vel
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
%Plot tracks:
figure
ma.plotTracks;
ma.labelPlotTracks;
title('Cell trajectories')
xlabel('X (\mu m)'); ylabel('Y (\mu m)')
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

%%%%%Velocities autocorrelation:
ma = ma.computeVCorr;
ma.vcorr
figure
ma.plotMeanVCorr
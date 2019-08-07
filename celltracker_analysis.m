clear all
close all
frames=120;
dt=2*15*60; %time interval in seconds
pixelsize=0.55; %pix/micron
path="C:\Users\G-mo10\Desktop\Test images nuclei Agata\CellTracker results\";
filename = "TracksCoordinates_C0506A1R3_007.mat";
filename=path+filename;
load(filename);
cells=pos(end:end);
matrix_x=nan(frames,cells); matrix_y=matrix_x;
%%%%%%%Fills matrices of positions: cell 0 | cell 1 | cell 2 | etc.:%%%%%%
for jj=1:cells
    kk=1;
for ii=1:length(pos)
    if pos(ii,4)==jj
        matrix_x(kk,jj)=pos(ii,1)/pixelsize; %in microns
        matrix_y(kk,jj)=pos(ii,2)/pixelsize; %in microns
        kk=kk+1;
    end
end
end
%%%%%%%Plots positions:%%%%%%%%%
% sz=100; a = (1:frames)'; b = num2str(flipud(a)); c = cellstr(b);
% for ll=1:cells
% scatter(matrix_x(:,ll),matrix_y(:,ll),sz,'filled')
% text(matrix_x(:,ll)-3.5, matrix_y(:,ll)+.5, c,'Color','white','FontWeight','bold')
% hold on
% plot(matrix_x(:,ll),matrix_y(:,ll))
% hold on
% %axis([0 200 0 120])
% end
% title('Cell trajectories')
% xlabel('Microns'); ylabel('Microns')

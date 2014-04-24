%% load data

clear all
close all
clc

load data/data0305_7_19.mat
%% visualize measurements 

measurements = route.densityMeasured(route.sensorCellMap,:);
surf(measurements,'Linestyle','None');
%surf(route.densityMeasured,'Linestyle','None');
view(2)
colorbar
%set(gca,'fontsize',14)
xlabel('time')
ylabel('position')

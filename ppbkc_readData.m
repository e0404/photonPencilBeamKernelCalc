clc
close all
clear

%% tpr

% read data from file
tmpTable = readtable('tpr.dat');
tprTmp = table2array(tmpTable);

tprFieldSizes = tprTmp(1,2:end);
tprDepths     = tprTmp(2:end,1);

tpr           = tprTmp(2:end,2:end);

% plot
figure
hold on
plot(tprDepths,tpr);
title('tpr')
xlabel('depth [mm]')
ylabel('a.u.')
grid minor
box on
legendStringsFieldSizes = mat2cell(tprFieldSizes',ones(length(tprFieldSizes),1));
legendStringsFieldSizes = cellfun(@num2str,legendStringsFieldSizes,'UniformOutput',false);
legend(legendStringsFieldSizes);

%% output factor

% read data from file
fileHandle = fopen('of.dat');
outputFactor = cell2mat(textscan(fileHandle,'%f %f','CommentStyle',{'#'}));
fclose(fileHandle);

% plot
figure
plot(outputFactor(:,1),outputFactor(:,2))
title('output factor')
xlabel('field size [mm]')
ylabel('a.u.')
grid minor
box on

%% primary fluence

% read data from file
fileHandle = fopen('primflu.dat');
primaryFluence = cell2mat(textscan(fileHandle,'%f %f','CommentStyle',{'#'}));
fclose(fileHandle);

% plot
figure
plot(primaryFluence(:,1),primaryFluence(:,2))
title('primary fluence')
xlabel('r [mm]')
ylabel('a.u.')
grid minor
box on

%% read parameters
fileHandle = fopen('params.dat');
tmp = textscan(fileHandle,'%s %f','CommentStyle',{'#'});
fclose(fileHandle);

params = containers.Map(tmp{1},tmp{2});

%% clear up
clear tprTmp tmp ans fileHandle;
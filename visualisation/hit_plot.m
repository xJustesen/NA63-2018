clear all; clc; close all
%% load data

filename = 'hits_coord_data_parallel_test';
filename1 = 'hits_coord_data_4'; % fake hits

hitdat   = loaddat(filename);
fakehitdat  = loaddat(filename1);
fakeevents = 130042
events = 201678


%% load data
input_ang   = load('/home/christian/Dropbox/Cern2018Experiment/kode/beam_direction_4mm.txt'); 

x = load('../beamParameters/xdat_80GeV_beam_params_1.5mm.txt');
y = load('../beamParameters/ydat_80GeV_beam_params_1.5mm.txt');
xw = load('../beamParameters/xweight_80GeV_beam_params_1.5mm.txt');
yw = load('../beamParameters/yweight_80GeV_beam_params_1.5mm.txt');

%% plot data
for i = 1:numel(fieldnames(hitdat))
    field = strcat('plane',num2str(i-1));
    fprintf(strcat(field,'\n'));

    titlestr = strcat('Plane ',num2str(i-1),' hits');
       
    [xcounts, xbins] = hist(hitdat.(field)(:,1),100);
    [ycounts, ybins] = hist(hitdat.(field)(:,2),100);
    
    f = figure;
    title(titlestr)
    plot(hitdat.(field)(:,1), hitdat.(field)(:,2),'.');
    axis equal
    box on
    
    if i == 1             
        f = figure;
        title('x')
        hold on
        plot(xbins, xcounts/max(xcounts))
        plot(x, xw/max(xw))
        
        figure
        title('y')
        hold on
        plot(ybins, ycounts/max(ycounts))
        plot(y, yw/max(yw))
    end
end

%% functions
function hitdat = loaddat(filename)
    datpath = '/home/christian/Documents/cern2018/simdata/';
    filepath = strcat(datpath,filename,'.txt');
    formatSpec = '%f %f';
    fileID = fopen(filepath,'r');
    datblocks = 6;
    for i = 1:datblocks
        field = strcat('plane',num2str(i-1));
        data = textscan(fileID,formatSpec,'HeaderLines',1);
        hitdat.(field) = [data{1},data{2}];
    end
end

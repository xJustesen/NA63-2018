clc;clear all;

% Combine data form multiple files in single file
number_of_files = 20;
parrent_file_name_energy = '/home/christian/Dropbox/Cern2018Experiment/grendel/spectre/cuts/energy_sim_background_20GeV';
parrent_file_name_events = '/home/christian/Dropbox/Cern2018Experiment/grendel/spectre/cuts/events_run_sim_background_20GeV';
ConcatenateFiles(parrent_file_name_energy, number_of_files);
SumFiles(parrent_file_name_events, number_of_files);

fprintf("Success!\n");

function ConcatenateFiles(parrent_file_name, number_of_files)
    complete_data_array = [];
    for i = 1:number_of_files
        partial_data_array = load(strcat(parrent_file_name,'_',num2str(i),'_finite_size.txt'));
        complete_data_array = [complete_data_array; partial_data_array];
    end
    dlmwrite(strcat(parrent_file_name,'.txt'),complete_data_array);
end

function SumFiles(parrent_file_name, number_of_files)
    sum = 0;
    for i = 1:number_of_files
        data = load(strcat(parrent_file_name,'_',num2str(i),'_finite_size.txt'));
        sum = sum + data;
    end
    dlmwrite(strcat(parrent_file_name,'.txt'),sum,'precision','%d');
end
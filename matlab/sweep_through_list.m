clear all; clc; close all

EMAIL = false;
%%% Send Email Setup
%%% This will email you if there is an error as well as the files when they
%%% are completed.
%%%
%%% 1. Add email + password to start email file
%%% 2. Uncomment the below lines
% start_email();
% EMAIL = true;
% EMAIL_ADDRESS = '...@....com';

% Setup
import com.comsol.model.*
import com.comsol.model.util.*
ModelUtil.showProgress(true);
COMSOL_MODEL_FILE = '../comsol/P2D-NLEIS.mph';

% Any .csv file in the SWEEP_DIR directory will be read in and run
SWEEP_DIR = '../supplementary-files/model-inputs/';
files = dir([SWEEP_DIR '*.csv']);

% Set the frequencies for COMSOL to simulate
FREQUENCIES = '10^range(5,-(1/3),-3)';

names_for_sending = {};
for file = files'
    tic
    RUN_LIST = readtable([SWEEP_DIR file.name], 'ReadVariableNames', false, 'HeaderLines', 1);

    file_name = strsplit(file.name,'.csv');
    output_folder = ['../data/simulations/' file_name{1} '/'];
    mkdir(output_folder);
    
    n_parameters = size(RUN_LIST,2);
    
    fid = fopen([SWEEP_DIR file.name], 'r');
    header = textscan(fid, repmat('%[^,],', [1,n_parameters]), 1);
    fclose(fid);
    
    names = cell(n_parameters,1);
    units = cell(size(names));

    for i=2:n_parameters
        string = cell2mat(header{i});
        split = strsplit(string,'[');
        names{i} = cell2mat(split(1));
        string2 = cell2mat(split(2));
        split2 = strsplit(string2,']');
        units{i} = split2(1);
    end
    
    for j = 1:height(RUN_LIST)
        tic
        run = table2array(RUN_LIST(j,1));
        disp(['Run:' cell2mat(run)]);
    
        model = mphload(COMSOL_MODEL_FILE);
        model.hist.disable;
    
        % Load parameter values for run
        parameter_values = table2array(RUN_LIST(j,2:end));
        for k = 1:size(parameter_values,2)
            variable = names(k+1);
            value = [num2str(parameter_values(k)) '[' cell2mat(units{k+1}) ']'];
            model.param.set(variable, value);
        end
    
        % Update location of probe
        L = parameter_values(1)+ parameter_values(2) ...
                        + parameter_values(3);
        disp(['Length = ' num2str(L)])
        model.probe('pdom1').setIndex('coords1', num2str(L),0,0);
    
        % Set frequency
        model.study('std1').feature('param').set('plistarr', FREQUENCIES);
        model.study('std1').feature('param').set('pname', 'f');
    
        model.study('std1').run;
            
        % Extract the harmonics
        try
            str = mphtable(model,'tbl5');
            tbl_data = str.data;
            f = tbl_data(:,1);
        catch
            warning('tbl_data does not exist for these parameters');
            tbl_data = NaN(1,7);
        end
        harmonics = tbl_data(:,:)';

        % Save the harmonics to file
        OUTPUT_FILE_NAME = [output_folder 'spectra-' cell2mat(run) '.txt'];
        names_for_sending = [names_for_sending OUTPUT_FILE_NAME];

        [fileID2, msg] = fopen(OUTPUT_FILE_NAME, 'a');
        if fileID2 == -1
            if EMAIL
                sendmail(EMAIL_ADDRESS, ['COMSOL FAILURE - Parameter Sweep (COMP ' COMPUTER ')'],...
                    ['fopen failed. harmonics data from run ' cell2mat(run) ' failed to save. \n Failure Message: ' msg]);
            end
            fclose('all');     
        else
            for line = 1:size(harmonics,1)
                formatSpec = ['%e' repmat(' %e', 1, size(harmonics,2)-1) '\n'];
                fprintf(fileID2, formatSpec, harmonics(line, :));
            end
            fclose('all'); 
        end 
        toc
    end
    if EMAIL
        sendmail(EMAIL_ADDRESS, ['Finished Sweep Step - ' file_name{1}],...
                ['Finished sweep section in ' num2str(toc) ' seconds'],...
                names_for_sending);
    end
end
% This test code compares bounded accumulation between
% Monte-Carlo simulation and Fokker-Plank propagation.
% There is a good agreement between them.
clear
close all

id = 1;

root_dir = '~/Documents/MATLAB/';
data_dir = [root_dir,'MonkeyPhys/790_sk/data/'];
% load([data_dir,'monkInfo_N-1_only']);
% load([data_dir,'monkInfo_N-1']);
file_name = 'monkInfo_N-1_N-2_N-3';
load([data_dir,file_name]);

num_param = [2,4,2,4];

if 0
    for i = 1:2
        switch i
            case 1
                load('popInfoEliSelective_200_131209');
            case 2
                load('popInfoJoeyAll_130_131209');
        end
        info = rmfield(info,'elMtrx');
        if isfield(info,'adMtrx')
            info = rmfield(info,'adMtrx');
        end
        info = rmfield(info,'SM');
        info = rmfield(info,'RM');
        info = rmfield(info,'FM');
        info = rmfield(info,'KM');
        info = rmfield(info,'LM');
        info = rmfield(info,'SpM');
        RT = info.TM(:,3);
        monkInfo{i}.RT = RT;
        save([data_dir,file_name],'monkInfo','-append')
    end
end

model_list = [{'simple'},{'urgency'},{'gain'},{'xjw'}];

for mi = 1
    fit_info = monkInfo{id}; % Eli

    fit_param = [];
    USfunc = 'Linear';

    % 'shared' or 'observed'
    fit_info.noise_property = 'shared';
    % 'ind' (independent) or 'anti' (anti-correlated)
    fit_info.noise_correlation_type = 'ind';
    fit_info.model = model_list{mi};
    fit_info.id = id;

    file_path = '~/Documents/MATLAB/MonkeyPhys/790_sk/model/frWoe2dSimFitDir22/';
    temp = fit_info;

    for rep = 1:30;

        % Deal with the file name
        noise_name = ['_',fit_info.noise_property,'_',fit_info.noise_correlation_type,'_noise'];
        id_name = sprintf('_%d',id);
        model_name = ['_',model_list{mi}];
        fit_num_name = sprintf('_%03d',rep);

        D = dir([file_path,'frWoe2dFit',noise_name,id_name,model_name,fit_num_name,'*']);

        % Check whether the file already exists
        if isempty(D)
            continue
        end

        load([file_path,D(1).name]);
        if Err(end)==0
            continue
        end
        fit_param = [fit_param;Prmts(end,:),Err(end)];

    end
    [fit_param,order] = sortrows(fit_param,size(fit_param,2));
    model_list{mi}
    best_param = fit_param(1,1:end-1)
    best_llk = fit_param(1,end)
    BIC = 2*best_llk + num_param(mi)*log(sum(fit_info.n_accum))

    if 0
        % converting from 4 parameter model to 3 parameter model
        best_param(2) = best_param(2)*exp(best_param(3)*best_param(4));
        best_param(1) = best_param(1)*26/18;
        best_param = best_param(1:3);
        fit_info.num_param = 3;
    end
    % best_param(2) = 60;
    % out = frWoe2dSimCalc(best_param,temp,'fig_switch',1);
    out = frWoe2dDynSimCalc(best_param,temp,'fig_switch',1);
end

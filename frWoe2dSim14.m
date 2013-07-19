% This test code compares bounded accumulation between
% Monte-Carlo simulation and Fokker-Plank propagation.
% There is a good agreement between them.

clear
close all
id = 1;
num_step = 10;
% Best fitted parameters for singly stochastic model.
switch id
    case 1
        sqrt_r = 1;% 0.07;
        fr_offset = 20;
        B = 55; % threshold
        min_fr = 0; % reflecting lower bound
        
        % all_epoch_gain = 8.1;
        % all_epoch_pos_gain = 13.7;
        % all_epoch_neg_gain = 2.3;
        
        % all_gain = 3.1 + 1.4*(1:num_step);
        % relative_gain = all_gain/all_epoch_gain;
        
        % pos_dfr_gain = all_epoch_pos_gain * relative_gain;
        % neg_dfr_gain = all_epoch_neg_gain * relative_gain;
        
        % pos_dfr_offset = 0*ones(1,num_step); % shift in the baseline FR in each epoch (i.e. urgency)
        % neg_dfr_offset = 0*ones(1,num_step);
        
        if 0
            % without offset
            pos_dfr_gain = 6.1 + 2.1*(1:num_step);
            neg_dfr_gain = 0.2 + 0.6*(1:num_step);

            pos_dfr_offset = 0*ones(1,num_step); % shift in the baseline FR in each epoch (i.e. urgency)
            neg_dfr_offset = 0*ones(1,num_step);
        else
            % with offset
            pos_dfr_gain = 9.3 + 1.6*(1:num_step);
            neg_dfr_gain = 0.2 + 0.6*(1:num_step);

            pos_dfr_offset = -1.7 + 0.3*(1:num_step); % shift in the baseline FR in each epoch (i.e. urgency)
            neg_dfr_offset = 1.0 - 0.1*(1:num_step);
        end
        relative_gain = ones(1,10);
        dfr_sd = 4;
        init_dfr_sd = 10;
        
    case 2
        sqrt_r = 0.13;
        fr_offset = 30;
        B = 45; % threshold
        min_fr = 10; % reflecting lower bound
        
        all_epoch_gain = 18.0;
        all_epoch_pos_gain = 24.2;% 30.6;%
        all_epoch_neg_gain = 12.8;
        
        all_gain = 21.3 - 1.4*(1:num_step);
        % all_gain = [21.3346   15.2767   18.3673   15.7779   17.7635     8.7154   13.6475   16.9276    9.4715    12.7828];
        % relative_gain = all_gain/all_epoch_gain;
        relative_gain = ones(1,10);
        
        pos_dfr_gain = all_epoch_pos_gain * relative_gain;
        neg_dfr_gain = all_epoch_neg_gain * relative_gain;

        pos_dfr_offset = 0*ones(1,num_step); % shift in the baseline FR in each epoch (i.e. urgency)
        neg_dfr_offset = 0*ones(1,num_step);
        % relative_gain = (pos_dfr_gain + pos_dfr_gain)/(pos_dfr_gain(1) + pos_dfr_gain(1));
        
        dfr_sd = 5;
        init_dfr_sd = 5;
end

if ismac
    root_dir = '~/Documents/MATLAB/';
else
    root_dir = '~/MATLAB/';
end
data_dir = [root_dir,'MonkeyPhys/790_sk/data/'];
sim_dir = [root_dir,'MonkeyPhys/790_sk/model/sim/'];
    % If the directory does not exist, create one.
if ~isdir(sim_dir)
    mkdir(sim_dir)
end

load([data_dir,'monkInfo']);
info = monkInfo{id}; % Eli
% Quick calculation of s.d. in baseline
% base_sd = nanstd(info.TM(:,51));

sim_switch = [1 0];
fig_switch = 1;
compute_flag = 1;
fig = 1;

    % Choose the noise type:
    % either 'anti' (anti-correlated) or 'ind' (independent)
noise_type = 'ind';
switch noise_type
    case 'anti'
        noise_name = 'Anti';
    case 'ind'
        noise_name = 'Ind';
end

shape_prob = make790StimProb(1);
cum_shape_prob = cumsum(shape_prob);

WOE = -0.9:0.2:0.9;
deci_WOE = -9:2:9;
num_sim = 1;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Numerical method (Monte-Carlo simulation)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sim_switch(1)
    for rep = 1:num_sim
        save_name = ['numericalSim',noise_name,sprintf('%03d',rep)];
        date_str = datestr(now,31);
        if ~exist([sim_dir,save_name,'.mat'],'file')
            save([sim_dir,save_name],'save_name','date_str');
        else
            % continue
        end
        
        % preallocation
        num_trial = 1e4;
        x1_hfr = nan(num_trial,num_step+1);
        x2_hfr = nan(num_trial,num_step+1);
        x1_fr = nan(num_trial,num_step+1);
        x2_fr = nan(num_trial,num_step+1);
        
%         x1_hfr(:,1) = 0;
%         x2_hfr(:,1) = 0;
%         x1_fr(:,1) = 0;
%         x2_fr(:,1) = 0;
        
        cum_woe = nan(num_trial,num_step+1);
        x1_correct_fr = nan(num_trial,1);
        x2_correct_fr = nan(num_trial,1);
        x1_error_fr = nan(num_trial,1);
        x2_error_fr = nan(num_trial,1);
        result = nan(num_trial,1);
        num_accum = nan(num_trial,1);
        cum_woe_correct = nan(num_trial,1);
        cum_woe_error = nan(num_trial,1);
        
    %%
        for k = 1:num_trial
            
            if mod(k,1e3)==0
                k
            end
            
            x1_hfr(k,1) = fr_offset + init_dfr_sd*randn(1);
            x2_hfr(k,1) = fr_offset + init_dfr_sd*randn(1);
            x1_fr(k,1) = x1_hfr(k,1);
            x1_fr(k,1) = x1_hfr(k,1);
            
            cum_woe(k,1) = 0;
            for ei = 1:num_step
                    % WOE
                shape_ind = find(rand(1)<=cum_shape_prob,1,'first');
                shape_woe = deci_WOE(shape_ind);
                cum_woe(k,ei+1) = cum_woe(k,ei) + shape_woe;
                    % FR
                noise1 = dfr_sd * randn(1);
                switch noise_type
                    case 'anti'
                        noise2 = -noise1;
                    case 'ind'
                        noise2 = dfr_sd * randn(1);
                end
                pos_ind = heaviside(shape_woe);
                if pos_ind
                    sig1 = shape_woe * 0.1 * pos_dfr_gain(ei) + noise1;
                    sig2 = -shape_woe * 0.1 * neg_dfr_gain(ei) + noise2;
                else
                    sig1 = shape_woe * 0.1 * neg_dfr_gain(ei) + noise1;
                    sig2 = -shape_woe * 0.1 * pos_dfr_gain(ei) + noise2;
                end
                
                x1_hfr(k,ei+1) = x1_hfr(k,ei) + sig1;
                x2_hfr(k,ei+1) = x2_hfr(k,ei) + sig2;
                
                x1_fr(k,ei+1) = relative_gain(ei)*x1_hfr(k,ei+1);
                x2_fr(k,ei+1) = relative_gain(ei)*x2_hfr(k,ei+1);

                if x1_fr(k,ei+1)<min_fr
                    x1_fr(k,ei+1)=min_fr;
                end
                if x2_fr(k,ei+1)<min_fr
                    x2_fr(k,ei+1)=min_fr;
                end

                % if FR of either neuron reaches the bound...
                if (x1_fr(k,ei+1)>B) || (x2_fr(k,ei+1)>B)
                    num_accum(k) = ei;
                    
                    if x1_fr(k,ei+1)>x2_fr(k,ei+1)
                        x1_correct_fr(k) = x1_fr(k,ei+1);
                        x2_correct_fr(k) = x2_fr(k,ei+1);
                        cum_woe_correct(k) = cum_woe(k,ei+1);
                        % x1_fr(k,ei+1) = nan;
                        % x2_fr(k,ei+1) = nan;
                        result(k) = 1;
                    elseif x1_fr(k,ei+1)<x2_fr(k,ei+1)
                        x1_error_fr(k) = x1_fr(k,ei+1);
                        x2_error_fr(k) = x2_fr(k,ei+1);
                        cum_woe_error(k) = cum_woe(k,ei+1);
                        % x1_fr(k,ei+1) = nan;
                        % x2_fr(k,ei+1) = nan;
                        result(k) = 0;
                    else % Tied case: this heppens very rarely, if any.
                        if rand(1)>0.5
                            x1_correct_fr(k) = x1_fr(k,ei+1);
                            x2_correct_fr(k) = x2_fr(k,ei+1);
                            cum_woe_correct(k) = cum_woe(k,ei+1);
                            % x1_fr(k,ei+1) = nan;
                            % x2_fr(k,ei+1) = nan;
                            result(k) = 1;
                        else
                            x1_correct_fr(k) = x1_fr(k,ei+1);
                            x2_correct_fr(k) = x2_fr(k,ei+1);
                            cum_woe_error(k) = cum_woe(k,ei+1);
                            % x1_fr(k,ei+1) = nan;
                            % x2_fr(k,ei+1) = nan;
                            result(k) = 0;
                        end
                    end
                    if x1_fr(k,ei+1)>B
                        x1_fr(k,ei+1) = B;
                    end
                    if x2_fr(k,ei+1)>B
                        x2_fr(k,ei+1) = B;
                    end
                    break
                end
            end
            
        end
        
        sim_struct = struct('B',B,...
                            'min_fr',min_fr,...
                            'fr_offset',fr_offset,...
                            'pos_dfr_gain',pos_dfr_gain,...
                            'neg_dfr_gain',neg_dfr_gain,...
                            'pos_dfr_offset',pos_dfr_offset,...
                            'neg_dfr_offset',neg_dfr_offset,...
                            'dfr_sd',dfr_sd,...
                            'WOE',WOE,...
                            'shape_prob',shape_prob,...
                            'num_trial',num_trial,...
                            'num_step',num_step,...
                            'noise_type',noise_type,...
                            'x1_fr',x1_fr,...
                            'x2_fr',x2_fr,...
                            'x1_hfr',x1_hfr,...
                            'x2_hfr',x2_hfr,...
                            'result',result,...
                            'num_accum',num_accum,...
                            'x1_correct_fr',x1_correct_fr,...
                            'x2_correct_fr',x2_correct_fr,...
                            'x1_error_fr',x1_error_fr,...
                            'x2_error_fr',x2_error_fr,...
                            'cum_woe',cum_woe,...
                            'cum_woe_correct',cum_woe_correct,...
                            'cum_woe_error',cum_woe_error);
                            
        save([sim_dir,save_name],'sim_struct','-append');
    end
    
    %% Format figures
    if fig_switch
        FH1 = figure(fig);hold on
        set(FH1,'position',[100 100 1200 400],'color','w');
        FH2 = figure(fig+1);hold on
        set(FH2,'position',[300 50 1200 400],'color','w');
        h_margin = 0.1;
        v_margin = 0.1;
        h = (1-h_margin)/5;
        v = (1-v_margin)/2;

        for fig_ind = 1:2
            switch fig_ind
                case 1
                    figure(fig);
                case 2
                    figure(fig+1);
            end
            for i = 1:3
                switch i
                    case 1
                        for j = 1:10
                            handles{fig_ind}(i,j) = subplot('position',[h_margin+h*mod(j-1,5)+h*0.1,1-v*(ceil(j/5))+v*0.1,h*0.6,v*0.6]);
                        end
                    case 2
                        for j = 1:10
                            handles{fig_ind}(i,j) = subplot('position',[h_margin+h*mod(j-1,5)+h*0.75,1-v*(ceil(j/5))+v*0.1,h*0.15,v*0.6]);
                        end
                    case 3
                        for j = 1:10
                            handles{fig_ind}(i,j) = subplot('position',[h_margin+h*mod(j-1,5)+h*0.1,1-v*(ceil(j/5))+v*0.75,h*0.6,v*0.15]);
                        end
                end
            end
        end
        
        %% 2D plot
        
            % Load the data first.
        if 1
            % Assembling simulation files
            x1_fr = [];
            x2_fr = [];
            cum_woe = [];
            result = [];
            num_accum = [];
            x1_correct_fr = [];
            x2_correct_fr = [];
            x1_error_fr = [];
            x2_error_fr = [];
            cum_woe_correct = [];
            cum_woe_error = [];
            num_trial = 0;
            for rep = 1:num_sim
                save_name = ['numericalSim',noise_name,sprintf('%03d',rep)];
                load([sim_dir,save_name]);
                x1_fr = [x1_fr; sim_struct.x1_fr];
                x2_fr = [x2_fr; sim_struct.x2_fr];
                x1_hfr = [x1_fr; sim_struct.x1_hfr];
                x2_hfr = [x2_fr; sim_struct.x2_hfr];
                cum_woe = [cum_woe; sim_struct.cum_woe];
                result = [result; sim_struct.result];
                num_accum = [num_accum; sim_struct.num_accum];
                x1_correct_fr = [x1_correct_fr; sim_struct.x1_correct_fr];
                x2_correct_fr = [x2_correct_fr; sim_struct.x2_correct_fr];
                x1_error_fr = [x1_error_fr; sim_struct.x1_error_fr];
                x2_error_fr = [x2_error_fr; sim_struct.x2_error_fr];
                cum_woe_correct = [cum_woe_correct; sim_struct.cum_woe_correct];
                cum_woe_error = [cum_woe_error; sim_struct.cum_woe_error];
                num_trial = num_trial + sim_struct.num_trial;
            end
            save([sim_dir,'numericalSim',noise_name,'All'],'sim_struct','x1_fr','x2_fr',...
                'cum_woe','result','num_accum','x1_correct_fr','x2_correct_fr',...
                'x1_error_fr','x2_error_fr','cum_woe_correct','cum_woe_error','num_trial');
        else
            load([sim_dir,'numericalSim',noise_type,'All'])
        end
        
        for fig_ind = 1:2
            for ei = 1:10
                woe_conv = -9*ei:2:9*ei;
                fr_conv = (-100*ei):(100*ei);

                figure(fig+fig_ind-1)
                subplot(handles{fig_ind}(1,ei))
                switch fig_ind
                    case 1
                        histmat_x1 = hist2(cum_woe(:,ei+1),x1_fr(:,ei+1),woe_conv,fr_conv);
                        imagesc(woe_conv/10,fr_conv,histmat_x1)
                    case 2
                        histmat_x2 = hist2(cum_woe(:,ei+1),x2_fr(:,ei+1),woe_conv,fr_conv);
                        imagesc(woe_conv/10,fr_conv,histmat_x2)
                end
                axis xy
                xlabel('WOE')
                ylabel('FR (subjective WOE)')
                xlim([-5 5])
                ylim([-10 90])
                hold off
                
                subplot(handles{fig_ind}(2,ei))
                switch fig_ind
                    case 1
                        x1_correct_fr_hist = histc(x1_correct_fr(num_accum==ei),fr_conv);
                        plot(fr_conv,x1_correct_fr_hist)
                    case 2
                        x2_error_fr_hist = histc(x2_error_fr(num_accum==ei),fr_conv);
                        plot(fr_conv,x2_error_fr_hist)
                end
                view(90,-90);
                xlim([-10 90])
                set(handles{fig_ind}(2,ei),'xTickLabel',[]);

                subplot(handles{fig_ind}(3,ei))
                switch fig_ind
                    case 1
                        cum_woe_correct_hist = histc(cum_woe_correct(num_accum==ei),woe_conv);
                        bar(woe_conv/10,cum_woe_correct_hist);
                    case 2
                        cum_woe_error_hist = histc(cum_woe_error(num_accum==ei),woe_conv);
                        bar(woe_conv/10,cum_woe_error_hist)
                end
                xlim([-5 5])
                set(handles{fig_ind}(3,ei),'xTickLabel',[]);
            end
        end
        
        %% Compute descriptive statistics of numerical solution

            % RT histogram
        figure(fig+2);clf
        set(gcf,'position',[0 550,500,700],'color','w')
        
        bins = 1:10;
        num_accum_hist = histc(num_accum,bins);
        num_accum_hist_correct = histc(num_accum(result==1),bins);
        num_accum_hist_error = histc(num_accum(result==0),bins);
        
        n_accum_pdf = info.n_accum./sum(info.n_accum);
        n_accum_correct_pdf = info.n_accum_correct./sum(info.n_accum);
        n_accum_error_pdf = info.n_accum_wrong./sum(info.n_accum);
        
        shift = 0.2;
        bar_width = 0.4;
        
        subplot(3,1,1);hold on;
        bar((1:10)-shift,n_accum_pdf(1:10),bar_width,'k')
        bar((1:10)+shift,num_accum_hist/sum(num_accum_hist),bar_width,'r')
        ylabel('Prob','FontSize',24,'FontWeight','bold')
        set(gca,'XTickMode','manual','XTick',1:10);
        set(gca,'FontSize',18,'FontWeight','bold','Box','OFF','TickDir','out');
        y = get(gca,'Ylim');
        xlim([0 11])
        
        subplot(3,1,2);hold on;
        bar((1:10)-shift,n_accum_correct_pdf(1:10),bar_width,'k')
        bar((1:10)+shift,num_accum_hist_correct/sum(num_accum_hist),bar_width,'r')
        ylabel('Prob','FontSize',24,'FontWeight','bold')
        set(gca,'XTickMode','manual','XTick',1:10);
        set(gca,'FontSize',18,'FontWeight','bold','Box','OFF','TickDir','out');
        xlim([0 11])
        ylim(y);
        
        subplot(3,1,3);hold on;
        bar((1:10)-shift,n_accum_error_pdf(1:10),bar_width,'k')
        bar((1:10)+shift,num_accum_hist_error/sum(num_accum_hist),bar_width,'r')
        xlabel('Number of shapes used for decision','FontSize',24,'FontWeight','bold')
        ylabel('Prob','FontSize',24,'FontWeight','bold')
        set(gca,'XTickMode','manual','XTick',1:10);
        set(gca,'FontSize',18,'FontWeight','bold','Box','OFF','TickDir','out');
        xlim([0 11])
        ylim(y);
        
%%
            % Survivor function
        figure(fig+3)
        set(gcf,'position',[500 550,500,400],'color','w')
        res_trial = num_trial-cumsum(num_accum_hist);
        bar(1:10,res_trial)
%%
            % mean and s.d. of bound reaching WOE
        figure(fig+4);hold on
        set(gcf,'position',[1000,550,800,500],'color','w')
        for ei = 1:10
            mean_cum_woe_correct(ei) = nanmean(cum_woe_correct(num_accum==ei));
            sd_cum_woe_correct(ei) = nanstd(cum_woe_correct(num_accum==ei));
            mean_cum_woe_error(ei) = nanmean(-cum_woe_error(num_accum==ei));
            sd_cum_woe_error(ei) = nanstd(-cum_woe_error(num_accum==ei));
            mean_cum_woe_all(ei) = nanmean([cum_woe_correct(num_accum==ei);-cum_woe_error(num_accum==ei)]);
            sd_cum_woe_all(ei) = nanstd([cum_woe_correct(num_accum==ei);-cum_woe_error(num_accum==ei)]);
        end
        
        % ploterr([1:10]-0.2,mean_cum_woe_correct/10,[],sd_cum_woe_correct/10,1,'m-.','abshhy',0.1);
        % ploterr([1:10]-0.1,mean_cum_woe_error/10,[],sd_cum_woe_error/10,1,'c-.','abshhy',0.1);
        errorbar([1:10]-0.2,mean_cum_woe_correct/10,sd_cum_woe_correct/10,'m-','LineWidth',2);
        errorbar([1:10]-0.1,mean_cum_woe_error/10,sd_cum_woe_error/10,'c-','LineWidth',2);
        errorbar([1:10]-0.1,mean_cum_woe_all/10,sd_cum_woe_all/10,'k-','LineWidth',2);
        
        for ei = 1:10;
            if size(info.cumLLR_all,1)>1
                meanCumLLR(ei) = mean([-info.cumLLR_all{1,1,ei};info.cumLLR_all{1,2,ei}]);
                stdCumLLR(ei) = std([-info.cumLLR_all{1,1,ei};info.cumLLR_all{1,2,ei}]);
            else
                meanCumLLR_all(ei) = mean(info.cumLLR_all{1,ei});
                stdCumLLR_all(ei) = std(info.cumLLR_all{1,ei});
                meanCumLLR_correct(ei) = mean(info.cumLLR_correct{1,ei});
                stdCumLLR_correct(ei) = std(info.cumLLR_correct{1,ei});
                meanCumLLR_error(ei) = mean(info.cumLLR_error{1,ei});
                stdCumLLR_error(ei) = std(info.cumLLR_error{1,ei});
            end
        end
    
        errorbar([1:10],meanCumLLR_all,stdCumLLR_all,'--k')
        errorbar([1:10],meanCumLLR_correct,stdCumLLR_correct,'--r')
        errorbar([1:10],meanCumLLR_error,stdCumLLR_error,'--b')
        
        % errorbar([1:10]-0.1,mean_bound_woe{3},sd_bound_woe{3},'-r')

        % legend('data','model','location','NorthWest')
        % legend('correct','error','location','NorthWest')
        ylim([-3 3]);
        xlabel('Shape epoch','FontSize',24,'FontWeight','bold')
        ylabel('Cumulative logLR','FontSize',24,'FontWeight','bold')
        set(gca,'FontSize',24,'FontWeight','bold','Box','OFF','TickDir','out');
%%        
            % Time-dependent accuracy (TDA)
        figure(fig+5); hold on
        set(gcf,'position',[0 50,500,400],'color','w')
        for ei = 1:10
            tda(ei) = sum(result(num_accum==ei))./sum(num_accum==ei);
        end
        plot(1:10,info.LR_p_correct(:,1),'-b');
        plot(1:10,info.LR_p_correct(:,2),'-b');
        plot(1:10,tda,'-r')
        set(gca,'XTickMode','manual','XTick',1:10);
        xlabel('Number of shapes used for decision','FontSize',24,'FontWeight','bold')
        ylabel('Proportion correct','FontSize',24,'FontWeight','bold')
        title('Time Dependent Accuracy','FontSize',24,'FontWeight','bold')
        xlim([0 11])
        ylim([0.5 1])
        set(gca,'FontSize',24,'FontWeight','bold','Box','OFF','TickDir','out');
%%        
            % fr vs total logLR
        figure(fig+6);clf;hold on;
        set(gcf,'position',[500 50,500,400],'color','w')
        for ei = 2:11
            unique_cum_woe = unique([cum_woe(:,ei);-cum_woe(:,ei)]);
            unique_cum_woe = unique_cum_woe(isfinite(unique_cum_woe));
            for wi = 1:length(unique_cum_woe)
                pick1 = logical(cum_woe(:,ei)==unique_cum_woe(wi));
                pick2 = logical(-cum_woe(:,ei)==unique_cum_woe(wi));
                mean_fr_given_cum_woe{ei}(wi) = nanmean([x1_fr(pick1,ei);x2_fr(pick2,ei)]);
                sd_fr_given_cum_woe{ei}(wi) = nanstd([x1_fr(pick1,ei);x2_fr(pick2,ei)]);
                se_fr_given_cum_woe{ei}(wi) = sd_fr_given_cum_woe{ei}(wi)/sqrt(sum(pick1)+sum(pick2));
            end
            if all(isnan([x1_fr(:,ei);x2_fr(:,ei)]))
                continue
            end
            [beta bint r rint stats] = regress([x1_fr(:,ei);x2_fr(:,ei)],[ones(size(cum_woe,1)*2,1),[cum_woe(:,ei);-cum_woe(:,ei)]],0.05);
            offset(ei) = beta(1);
            offset_ci(ei) = (bint(1,2)-bint(1,1))/2;
            slope(ei) = beta(2)*10;
            slope_ci(ei) = (bint(2,2)-bint(2,1))/2*10;
            p_values(ei) = stats(3);

            subplot(2,5,ei-1);hold on;
            ploterr(unique_cum_woe/10,mean_fr_given_cum_woe{ei},[],se_fr_given_cum_woe{ei},1,'bo-');
            plot(unique_cum_woe/10,slope(ei).*unique_cum_woe./10+offset(ei),'k--');
            axis([-2.5 2.5 0 60]);
            maxY = 60;
            text(-2,maxY*0.9,sprintf('Slope'),'FontSize',16,'FontWeight','bold');
            text(-2,maxY*0.8,sprintf('%.1f\\pm%.1f',slope(ei),slope_ci(ei)),'FontSize',16,'FontWeight','bold');
            text(-2,maxY*0.6,sprintf('Offset'),'FontSize',16,'FontWeight','bold');
            text(-2,maxY*0.5,sprintf('%.1f\\pm%.1f',offset(ei),offset_ci(ei)),'FontSize',16,'FontWeight','bold');
        end
  %%      
        % delta_fr vs delta_logLR
        
        figure(fig+7)
        set(gcf,'position',[1000 50,500,400],'color','w')
        for ei = 2:11
            all_delta_fr = [];
            all_delta_woe = [];
            for wi = 1:length(deci_WOE)
                delta_woe = cum_woe(:,ei) - cum_woe(:,ei-1);
                pick1 = logical(delta_woe==deci_WOE(wi));
                delta_x1_fr = x1_fr(pick1,ei)-x1_fr(pick1,ei-1);
                pick2 = logical(delta_woe==-deci_WOE(wi));
                delta_x2_fr = x2_fr(pick2,ei)-x2_fr(pick2,ei-1);
                mean_delta_fr(ei-1,wi) = nanmean([delta_x1_fr;delta_x2_fr]);
                sd_delta_fr(ei-1,wi) = nanstd([delta_x1_fr;delta_x2_fr]);
                se_delta_fr(ei-1,wi) = sd_delta_fr(ei-1,wi)/sqrt(sum(pick1)+sum(pick2));
                all_delta_fr = [all_delta_fr;delta_x1_fr;delta_x2_fr];
                all_delta_woe = [all_delta_woe;deci_WOE(wi)*ones(sum(pick1)+sum(pick2),1)];
            end
            
            [beta bint r rint stats] = regress(all_delta_fr,[ones(size(all_delta_fr)),all_delta_woe],0.05);
            offset(ei) = beta(1);
            offset_ci(ei) = (bint(1,2)-bint(1,1))/2;
            slope(ei) = beta(2)*10;
            slope_ci(ei) = (bint(2,2)-bint(2,1))/2*10;
            p_values(ei) = stats(3);
            
            subplot(2,5,ei-1);hold on;
            ploterr(WOE,mean_delta_fr(ei-1,:),[],se_delta_fr(ei-1,:),1,'ok','abshhy',0);
            plot(deci_WOE/10,slope(ei).*deci_WOE./10+offset(ei),'k--');
            axis([-1 1 -20 20])
            y = get(gca,'Ylim');
            text(-0.9,y(1)+diff(y)*0.9,sprintf('Slope'),'FontSize',16,'FontWeight','bold');
            text(-0.9,y(1)+diff(y)*0.8,sprintf('%.1f\\pm%.1f',slope(ei),slope_ci(ei)),'FontSize',16,'FontWeight','bold');
            text(-0.9,y(1)+diff(y)*0.6,sprintf('Offset'),'FontSize',16,'FontWeight','bold');
            text(-0.9,y(1)+diff(y)*0.5,sprintf('%.1f\\pm%.1f',offset(ei),offset_ci(ei)),'FontSize',16,'FontWeight','bold');
        end
%%        
        % Psychometric curve
        figure(fig+8);hold on;
        set(gcf,'position',[1500 50,500,400],'color','w')
        unique_cum_woe = unique([cum_woe_correct;cum_woe_error]);
        unique_cum_woe = unique_cum_woe(isfinite(unique_cum_woe));
        for wi = 1:length(unique_cum_woe)
            R_num1 = sum(cum_woe_correct==unique_cum_woe(wi));
            L_num1 = sum(cum_woe_error==unique_cum_woe(wi));
            
            R_num2 = sum(cum_woe_error==-unique_cum_woe(wi));
            L_num2 = sum(cum_woe_correct==-unique_cum_woe(wi));
            
            R_num = R_num1 + R_num2;
            L_num = L_num1 + L_num2;
            PR(wi) = R_num/(R_num + L_num);
        end
        plot(info.uniqueLLR/10,info.PR,'ko','MarkerFaceColor','k');
        plot(unique_cum_woe/10,PR,'ro','MarkerFaceColor','r');
        xlabel('total logLR','FontSize',24,'FontWeight','bold')
        ylabel('Proportion of R choice','FontSize',24,'FontWeight','bold')
        set(gca,'FontSize',18,'FontWeight','bold','Box','OFF','TickDir','out');
        axis([-3 3 0 1]);
        
    end
end
num_fig = 7;
fig = fig+num_fig;
return
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Analytical Method
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sim_switch(2)
    save_name = ['analytical',noise_name];
    if compute_flag % set to 0 to just load analytical solutions
        fr = -100:100;
        delta_fr_mean_temp = [epoch_delta_fr_mean{1}(1:4), 0, 0, epoch_delta_fr_mean{1}(5:8)];
        delta_fr_sd_temp = [epoch_delta_fr_sd{1}(1:4), 1, 1, epoch_delta_fr_sd{1}(5:8)];
        
        for i = 1:length(WOE)
            px{1,1}(:,i) = normpdf(fr,delta_fr_mean_temp(i)+fr_offset,delta_fr_sd_temp(i)+base_sd).*shape_prob(i);
        end
        for i = 1:length(WOE)
            px{2,1}(:,i) = normpdf(fr,delta_fr_mean_temp(11-i)+fr_offset,delta_fr_sd_temp(11-i)+base_sd).*shape_prob(i);
        end

        %% Format figures
        if fig_switch
            FH1 = figure(fig);hold on
            set(FH1,'position',[100 100 1200 400],'color','w');
            FH2 = figure(fig+1);hold on
            set(FH2,'position',[300 50 1200 400],'color','w');
            h_margin = 0.1;
            v_margin = 0.1;
            h = (1-h_margin)/5;
            v = (1-v_margin)/2;

            for fig_ind = 1:2
                switch fig_ind
                    case 1
                        figure(fig);
                    case 2
                        figure(fig+1);
                end
                for i = 1:3
                    switch i
                        case 1
                            for j = 1:10
                                handles{fig_ind}(i,j) = subplot('position',[h_margin+h*mod(j-1,5)+h*0.1,1-v*(ceil(j/5))+v*0.1,h*0.6,v*0.6]);
                            end
                        case 2
                            for j = 1:10
                                handles{fig_ind}(i,j) = subplot('position',[h_margin+h*mod(j-1,5)+h*0.75,1-v*(ceil(j/5))+v*0.1,h*0.15,v*0.6]);
                            end
                        case 3
                            for j = 1:10
                                handles{fig_ind}(i,j) = subplot('position',[h_margin+h*mod(j-1,5)+h*0.1,1-v*(ceil(j/5))+v*0.75,h*0.6,v*0.15]);
                            end
                    end
                end
            end
        end
    %%
        tic
        fr_conv = -100:100;
        woe_full = (-0.9*num_step):0.1:(0.9*num_step);
        bound_fr_woe_sum{1} = zeros(length(fr_conv),length(woe_full));
        bound_fr_woe_sum{2} = zeros(length(fr_conv),length(woe_full));
        for ei = 1:num_step
            woe_conv = (-0.9*ei):0.2:(0.9*ei);
            % fr_conv = (-100*ei):(100*ei);
            %% Anti-correlated noise          
            if strcmp(noise_type,'anti')
                
                %%%%%%%%%%%%%%%% deal with a positive bound %%%%%%%%%%%%%%%%%
                bound_fr_woe{1} = zeros(length(fr_conv),length(woe_conv));
                bound_fr_woe{1}(fr_conv>B,:) = px{1,ei}(fr_conv>B,:);
                bound_reach_prob{1}(ei) = sum(sum(bound_fr_woe{1}));

                % WOE
                bound_woe_prob{1} = sum(bound_fr_woe{1},1)./bound_reach_prob{1}(ei);
                mean_bound_woe{1}(ei) = sum(bound_woe_prob{1}.*woe_conv);
                sd_bound_woe{1}(ei) = sqrt(sum(bound_woe_prob{1}.*(woe_conv.^2)) - mean_bound_woe{1}(ei).^2);

                % FR
                bound_fr_prob{1} = sum(bound_fr_woe{1},2)'./bound_reach_prob{1}(ei);
                mean_bound_fr{1}(ei) = sum(bound_fr_prob{1}.*fr_conv);
                sd_bound_fr{1}(ei) = sqrt(sum(bound_fr_prob{1}.*(fr_conv.^2)) - mean_bound_fr{1}(ei).^2);
                px{1,ei}(fr_conv>B,:)=0;

                %%%%%%%%%%%%%%%% deal with a negative bound %%%%%%%%%%%%%%%%%%
                bound_fr_woe{2} = zeros(length(fr_conv),length(woe_conv));
                bound_fr_woe{2}(fr_conv>B,:) = px{2,ei}(fr_conv>B,:);
                bound_reach_prob{2}(ei) = sum(sum(bound_fr_woe{2}));

                % WOE
                bound_woe_prob{2} = sum(bound_fr_woe{2},1)./bound_reach_prob{2}(ei);
                mean_bound_woe{2}(ei) = sum(bound_woe_prob{2}.*woe_conv);
                sd_bound_woe{2}(ei) = sqrt(sum(bound_woe_prob{2}.*(woe_conv.^2)) - mean_bound_woe{2}(ei).^2);

                % FR
                bound_fr_prob{2} = sum(bound_fr_woe{2},2)'./bound_reach_prob{2}(ei);
                mean_bound_fr{2}(ei) = sum(bound_fr_prob{2}.*fr_conv);
                sd_bound_fr{2}(ei) = sqrt(sum(bound_fr_prob{2}.*(fr_conv.^2)) - mean_bound_fr{2}(ei).^2);
                px{2,ei}(fr_conv>B,:)=0;
                %% Eliminating the bound-reaching probability from the competing neuron
            
                B_ind = find(fr_conv>=B,1,'first');
                min_fr_ind = find(fr_conv>=min_fr,1,'first');

                %%%%%%%%%%%%%%%% Anti-correlated noise %%%%%%%%%%%%%%%%%%
                % Probability of reaching a Tout bound is subtracted from probability
                % of lowest FR in positive neuron.
                res_prob = cumsum(px{1,ei},1) - repmat(sum(bound_fr_woe{2},1),length(fr_conv),1);
                unbound_fr_woe{1} = zeros(length(fr_conv),length(woe_conv));
                for i = 1:length(woe_conv)
                    ind = find(res_prob(:,i)>0,1,'first');
                    unbound_fr_woe{1}(1:ind-1,i) = px{1,ei}(1:ind-1,i);
                    unbound_fr_woe{1}(ind,i) = px{1,ei}(ind,i)-res_prob(ind,i);
                    px{1,ei}(1:ind-1,i) = 0;
                    px{1,ei}(ind,i) = res_prob(ind,i);
                    % FR cannot go below min_fr
                    px{1,ei}(min_fr_ind,i) = px{1,ei}(min_fr_ind,i) + sum(px{1,ei}(1:min_fr_ind-1,i),1);
                    px{1,ei}(1:min_fr_ind-1,i) = 0;
                    unbound_fr_woe{1}(min_fr_ind,i) = unbound_fr_woe{1}(min_fr_ind,i) + sum(unbound_fr_woe{1}(1:min_fr_ind-1,i),1);
                    unbound_fr_woe{1}(1:min_fr_ind-1,i) = 0;
                end

                %%%%%%%%%%%%%%%% Anti-correlated noise %%%%%%%%%%%%%%%%%%
                % Probability of reaching a Tin bound is subtracted from probability
                % of lowest FR in negative neuron.
                res_prob = cumsum(px{2,ei},1) - repmat(sum(bound_fr_woe{1},1),length(fr_conv),1);
                unbound_fr_woe{2} = zeros(length(fr_conv),length(woe_conv));
                for i = 1:length(woe_conv)
                    ind = find(res_prob(:,i)>0,1,'first');
                    unbound_fr_woe{2}(1:ind-1,i) = px{2,ei}(1:ind-1,i);
                    unbound_fr_woe{2}(ind,i) = px{2,ei}(ind,i)-res_prob(ind,i);
                    px{2,ei}(1:ind-1,i) = 0;
                    px{2,ei}(ind,i) = res_prob(ind,i);
                    % FR cannot go below min_fr
                    px{2,ei}(min_fr_ind,i) = px{2,ei}(min_fr_ind,i) + sum(px{2,ei}(1:min_fr_ind-1,i),1);
                    px{2,ei}(1:min_fr_ind-1,i) = 0;
                    unbound_fr_woe{2}(min_fr_ind,i) = unbound_fr_woe{2}(min_fr_ind,i) + sum(unbound_fr_woe{2}(1:min_fr_ind-1,i),1);
                    unbound_fr_woe{2}(1:min_fr_ind-1,i) = 0;
                end
            end
            
            %% independent noise
            if strcmp(noise_type,'ind')
                % Preallocation
                bound_fr_woe{1} = zeros(length(fr_conv),length(woe_conv));
                bound_fr_woe{2} = zeros(length(fr_conv),length(woe_conv));
                unbound_fr_woe{1} = zeros(length(fr_conv),length(woe_conv));
                unbound_fr_woe{2} = zeros(length(fr_conv),length(woe_conv));
                
                % Computing only for [-100 100] Hz, given that the probability
                % mass exists only within that range.
                
                
                t3(ei)=toc;
                min_fr_conv = -20;
                max_fr_conv = 80;
                pick = find(fr_conv>=min_fr_conv & fr_conv<=max_fr_conv);
                B_ind = find(fr_conv(pick)==B,1,'first');
                for i = 1:length(woe_conv)
                    % pick = 1:length(fr_conv);
                    
                    
                    % [FR1_choice1,FR2_choice1,FR1_choice2,FR2_choice2,FR1_res,FR2_res]=indNoiseRace(px{1,ei}(pick,i),px{2,ei}(pick,i),B,fr_conv(pick));
                    % [fr1_choice1,fr2_choice1,fr1_choice2,fr2_choice2,fr1_res,fr2_res]=indNoiseRace(px{1,ei}(pick,i),px{2,ei}(pick,i),B,fr_conv(pick));
                    % t1(i)=toc;
                    % [A,FR2_choice2,FR1_choice2,FR1_choice1,FR2_choice1,FR2_res,FR1_res]=mexIndRace(px{1,ei}(pick,i),px{2,ei}(pick,i),B_ind);
                    [A,FR1_choice1,FR2_choice1,FR1_choice2,FR2_choice2,FR1_res,FR2_res]=mexIndRace(px{1,ei}(pick,i),px{2,ei}(pick,i),B_ind);
                    % t2(i)=toc;
                    % size(fr1_choice1)
                    % size(FR1_choice1')
                    % sum((fr1_choice1 - FR1_choice1').^2)
                    bound_fr_woe{1}(pick,i) = FR1_choice1'; % correct trial
                    bound_fr_woe{2}(pick,i) = FR2_choice2'; % error trial
                    unbound_fr_woe{1}(pick,i) = FR1_choice2'; % error trial
                    unbound_fr_woe{2}(pick,i) = FR2_choice1'; % correct trial
                    px{1,ei}(pick,i) = FR1_res';
                    px{2,ei}(pick,i) = FR2_res';
                end
                for i = 1:length(woe_conv)
                    pick = find(abs(woe_full-woe_conv(i))<1e-3);
                    bound_fr_woe_sum{1}(:,pick) = bound_fr_woe_sum{1}(:,pick) + bound_fr_woe{1}(:,i);
                    bound_fr_woe_sum{2}(:,pick) = bound_fr_woe_sum{2}(:,pick) + bound_fr_woe{2}(:,i);
                end
                t4(ei)=toc;
                % Positive-preferring neuron
                
                bound_reach_prob{1}(ei) = sum(sum(bound_fr_woe{1}));
                
                % WOE
                bound_woe_prob{1,ei} = sum(bound_fr_woe{1},1)./bound_reach_prob{1}(ei);
                mean_bound_woe{1}(ei) = sum(bound_woe_prob{1,ei}.*woe_conv);
                sd_bound_woe{1}(ei) = sqrt(sum(bound_woe_prob{1,ei}.*(woe_conv.^2)) - mean_bound_woe{1}(ei).^2);

                % FR
                bound_fr_prob{1,ei} = sum(bound_fr_woe{1},2)'./bound_reach_prob{1}(ei);
                mean_bound_fr{1}(ei) = sum(bound_fr_prob{1}.*fr_conv);
                sd_bound_fr{1}(ei) = sqrt(sum(bound_fr_prob{1}.*(fr_conv.^2)) - mean_bound_fr{1}(ei).^2);
                
                % Negative-preferring neuron
                
                bound_reach_prob{2}(ei) = sum(sum(bound_fr_woe{2}));
                
                % WOE
                bound_woe_prob{2,ei} = sum(bound_fr_woe{2},1)./bound_reach_prob{2}(ei);
                mean_bound_woe{2}(ei) = sum(bound_woe_prob{2,ei}.*woe_conv);
                sd_bound_woe{2}(ei) = sqrt(sum(bound_woe_prob{2,ei}.*(woe_conv.^2)) - mean_bound_woe{2}(ei).^2);

                % FR
                bound_fr_prob{2,ei} = sum(bound_fr_woe{2},2)'./bound_reach_prob{2}(ei);
                mean_bound_fr{2}(ei) = sum(bound_fr_prob{2}.*fr_conv);
                sd_bound_fr{2}(ei) = sqrt(sum(bound_fr_prob{2}.*(fr_conv.^2)) - mean_bound_fr{2}(ei).^2);
                
                % correct and error trials combined together
                
                % WOE
                bound_woe_prob{3,ei} = sum(bound_fr_woe{1} + fliplr(bound_fr_woe{2}),1)./(bound_reach_prob{1}(ei)+bound_reach_prob{2}(ei));
                mean_bound_woe{3}(ei) = sum(bound_woe_prob{3,ei}.*woe_conv);
                sd_bound_woe{3}(ei) = sqrt(sum(bound_woe_prob{3,ei}.*(woe_conv.^2)) - mean_bound_woe{3}(ei).^2);
                
                bound_fr_prob{3} = sum(bound_fr_woe{1} + bound_fr_woe{2},2)'./(bound_reach_prob{1}(ei)+bound_reach_prob{2}(ei));
                mean_bound_fr{3}(ei) = sum(bound_fr_prob{3}.*fr_conv);
                sd_bound_fr{3}(ei) = sqrt(sum(bound_fr_prob{3}.*(fr_conv.^2)) - mean_bound_fr{3}(ei).^2);
                
                %%%%%%% Reflecting lower bound %%%%%%%
                min_fr_ind = find(fr_conv>=min_fr,1,'first');
                
                px{1,ei}(min_fr_ind,:) = px{1,ei}(min_fr_ind,:) + sum(px{1,ei}(1:min_fr_ind-1,:),1);
                px{1,ei}(1:min_fr_ind-1,:) = 0;
                
                px{2,ei}(min_fr_ind,:) = px{2,ei}(min_fr_ind,:) + sum(px{2,ei}(1:min_fr_ind-1,:),1);
                px{2,ei}(1:min_fr_ind-1,:) = 0;

            end
            
        %% Computing the mean and standard deviation of FR across WOE
            for i = 1:2 % neuron index

                px_norm{i} = px{i,ei}./repmat(sum(px{i,ei},1),size(px{i,ei},1),1);

                mean_px{i,ei} = fr_conv*px_norm{i};
                sd_px{i,ei} = sqrt((fr_conv.^2)*px_norm{i}-(fr_conv*px_norm{i}).^2);

                bound_fr_woe_norm{i} = ...
                    bound_fr_woe{i}./repmat(sum(bound_fr_woe{i},1),size(bound_fr_woe{i},1),1);
                mean_bound_fr_woe{i,ei} = fr_conv*bound_fr_woe_norm{i};
                sd_bound_fr_woe{i,ei} = sqrt((fr_conv.^2)*bound_fr_woe_norm{i}-(fr_conv*bound_fr_woe_norm{i}).^2);

                unbound_fr_woe_norm{i} = ...
                    unbound_fr_woe{i}./repmat(sum(unbound_fr_woe{i},1),size(unbound_fr_woe{i},1),1);
                mean_unbound_fr_woe{i,ei} = fr_conv*unbound_fr_woe_norm{i};
                sd_unbound_fr_woe{i,ei} = sqrt((fr_conv.^2)*unbound_fr_woe_norm{i}-(fr_conv*unbound_fr_woe_norm{i}).^2);
            end

        %%  
            if fig_switch && ei<=10   
                for fig_ind = 1:2
                    figure(fig+fig_ind-1)
                    subplot(handles{fig_ind}(1,ei));hold on

                    % zero padding for visualization
                    if 0.9*ei<5
                        if mod(ei,2)
                            WOE_temp = (-0.9*7):0.2:(0.9*7);
                            px_temp{fig_ind} = zeros(size(px{fig_ind,ei},1),length(WOE_temp));
                            px_temp{fig_ind}(:,length(WOE_temp)/2-length(woe_conv)/2+1:length(WOE_temp)/2+length(woe_conv)/2)=px{fig_ind,ei};
                        else
                            WOE_temp = (-0.9*6):0.2:(0.9*6);
                            px_temp{fig_ind} = zeros(size(px{fig_ind,ei},1),length(WOE_temp));
                            px_temp{fig_ind}(:,length(WOE_temp)/2-length(woe_conv)/2+1:length(WOE_temp)/2+length(woe_conv)/2)=px{fig_ind,ei};
                        end
                        imagesc(WOE_temp,fr_conv,px_temp{fig_ind})
                    else
                        imagesc(woe_conv,fr_conv,px{fig_ind,ei})
                    end

                    % imagesc(woe_conv,fr_conv,px{1,ei},[0 max(max(maxP*0.5))])
                    axis xy
                    xlabel('WOE')
                    ylabel('FR (sbujective WOE)')
                    xlim([-5 5])
                    ylim([-10 90])
                    hold off

                    subplot(handles{fig_ind}(2,ei))
                    switch fig_ind
                        case 1
                            plot(fr_conv,bound_fr_prob{1})
                            % zero_padded_bound_fr_prob{1} = [zeros(1,length(fr_conv)-length(bound_fr_prob{1})),sum(bound_fr_woe{1},2)'];
                            % plot(fr_conv,zero_padded_bound_fr_prob{1})
                        case 2
                            plot(fr_conv,bound_fr_prob{2})
                            % zero_padded_bound_fr_prob{2} = [zeros(1,length(fr_conv)-length(bound_fr_prob{1})),sum(bound_fr_woe{2},2)'];
                            % plot(fr_conv,zero_padded_bound_fr_prob{2})
                    end
                    view(90,-90);
                    xlim([-10 90])
                    set(handles{fig_ind}(2,ei),'xTickLabel',[]);

                    subplot(handles{fig_ind}(3,ei))
                    switch fig_ind
                        case 1
                            bar(woe_conv,sum(bound_fr_woe{1},1))
                        case 2
                            bar(woe_conv,sum(bound_fr_woe{2},1))
                    end
                    xlim([-5 5])
                    set(handles{fig_ind}(3,ei),'xTickLabel',[]);
                    % imagesc(woe_conv,fr_conv,px{1,ei},[0 max(max(maxP*0.5))])
                end
            end

            if ei==num_step
                break
            end
            
            delta_fr_mean_temp = [epoch_delta_fr_mean{ei+1}(1:4), 0, 0, epoch_delta_fr_mean{ei+1}(5:8)];
            delta_fr_sd_temp = [epoch_delta_fr_sd{ei+1}(1:4), 1, 1, epoch_delta_fr_sd{ei+1}(5:8)];
            for i = 1:length(WOE)
                T1(:,i) = normpdf(fr,delta_fr_mean_temp(i),delta_fr_sd_temp(i)).*shape_prob(i);
            end
            for i = 1:length(WOE)
                T2(:,i) = normpdf(fr,delta_fr_mean_temp(11-i),delta_fr_sd_temp(11-i)).*shape_prob(i);
            end
            toc;
            px1_temp = conv2(px{1,ei},T1);
            px2_temp = conv2(px{2,ei},T2);
            
            px{1,ei+1} = px1_temp(101:301,:);
            px{2,ei+1} = px2_temp(101:301,:);
            
            % px{1,ei+1} = px1_temp;
            % px{2,ei+1} = px2_temp;
            
            % px{1,ei+1} = convolve2(px{1,ei},T1);
            % px{2,ei+1} = convolve2(px{2,ei},T2);
            toc;
        end
        % cd(sim_dir);
        % save([sim_dir,save_name]);
        toc
    else
        load([sim_dir,save_name]);
    end
    
        % Compute descriptive statistics of analytical solution
    
        % RT histogram
        
    n_accum_pdf = info.n_accum./sum(info.n_accum);
    n_accum_correct_pdf = info.n_accum_correct./sum(info.n_accum);
    n_accum_wrong_pdf = info.n_accum_wrong./sum(info.n_accum);
    
    figure(fig+2)
    subplot(3,1,1);
    bar(1:10,[n_accum_correct_pdf(1:10),bound_reach_prob{1}'])
    subplot(3,1,2);
    bar(1:10,[n_accum_wrong_pdf(1:10),bound_reach_prob{2}'])
    subplot(3,1,3);
    bar(1:10,[n_accum_pdf(1:10),bound_reach_prob{1}'+bound_reach_prob{2}'])
    
        % RT cdf
    figure(fig+4); hold on
    n_accum_cdf = cumsum(n_accum_pdf);
    bound_reach_cdf = cumsum(bound_reach_prob{1}+bound_reach_prob{2});

    plot(1:10,n_accum_cdf(1:10),'-b')
    plot(1:10,bound_reach_cdf,'-r')
    
    
        % Time-dependent accuracy (TDA)
    figure(fig+3); hold on
    tda = bound_reach_prob{1}./(bound_reach_prob{1}+bound_reach_prob{2});
    plot(1:10,info.LR_p_correct(:,1),'-b');
    plot(1:10,info.LR_p_correct(:,2),'-b');
    plot(1:10,tda,'-r')
    ylim([0.5 1])
    
        
    
        % mean ± s.d. of bound reaching WOE
    if all(sim_switch)
        figure(fig-1);hold on
    else
        figure(fig+5);hold on
    end
    % ploterr(1:10,mean_bound_woe{1},[],sd_bound_woe{1},1,'r-.','abshhy',0.1);
    % ploterr([1:10]+0.1,mean_bound_woe{2},[],sd_bound_woe{2},1,'b-.','abshhy',0.1);
    
    % errorbar(1:10,mean_bound_woe{1},sd_bound_woe{1},'-r')
    % errorbar([1:10]+0.1,-mean_bound_woe{2},sd_bound_woe{2},'-b')
    
    for ei = 1:10;
        if size(info.cumLLR_all,2)>1
            meanCumLLR(ei) = mean([-info.cumLLR_all{1,1,ei};info.cumLLR_all{1,2,ei}]);
            stdCumLLR(ei) = std([-info.cumLLR_all{1,1,ei};info.cumLLR_all{1,2,ei}]);
        else
            meanCumLLR(ei) = mean(info.cumLLR_all{1,1,ei});
            stdCumLLR(ei) = std(info.cumLLR_all{1,1,ei});
        end
    end
    
    errorbar([1:10],meanCumLLR,stdCumLLR,'-b')
    errorbar([1:10]-0.1,mean_bound_woe{3},sd_bound_woe{3},'-r')
    
    legend('data','model','location','NorthWest')
    xlabel('Shape epoch')
    ylabel('Cumulative logLR')
    
        % mean ± s.d of bound reaching FR as function of WOE in each epoch  
    figure(fig+6)
    ymin = 0;
    ymax = 70;
    for ei = 1:10
        woe_conv = (-0.9*ei):0.2:(0.9*ei);
        subplot(4,10,ei)
        pick = logical(bound_woe_prob{1,ei}>1e-3); % Do not trust if the prob is too small.
        errorbar(woe_conv(pick),mean_bound_fr_woe{1,ei}(pick),sd_bound_fr_woe{1,ei}(pick),'-r')
        xlim([-4 4])
        ylim([ymin ymax])
        subplot(4,10,10+ei)
        pick = logical(bound_woe_prob{2,ei}>1e-3); % Do not trust if the prob is too small.
        errorbar(woe_conv(pick),mean_bound_fr_woe{2,ei}(pick),sd_bound_fr_woe{2,ei}(pick),'-b')
        xlim([-4 4])
        ylim([ymin ymax])
        subplot(4,10,20+ei)
        pick = logical(bound_woe_prob{1,ei}>1e-3); % Do not trust if the prob is too small.
        errorbar(woe_conv(pick),mean_unbound_fr_woe{1,ei}(pick),sd_unbound_fr_woe{1,ei}(pick),'-r')
        xlim([-4 4])
        ylim([ymin ymax])
        subplot(4,10,30+ei)
        pick = logical(bound_woe_prob{2,ei}>1e-3); % Do not trust if the prob is too small.
        errorbar(woe_conv(pick),mean_unbound_fr_woe{2,ei}(pick),sd_unbound_fr_woe{2,ei}(pick),'-b')
        xlim([-4 4])
        ylim([ymin ymax])
    end
%%
    figure(fig+7);clf;hold on;
    % Obtaining the slope of a linear regression for fr vs total logLR by
    % Monte-Carlo method
    n_row = size(bound_fr_woe_sum{1},1);
    n_col = size(bound_fr_woe_sum{1},2);
    % Creating a giant cdf for two matrices
    % bound_fr_woe_cdf = cumsum([reshape(bound_fr_woe_sum{1},[],1);reshape(bound_fr_woe_sum{2},[],1)]);
    bound_fr_woe_combined = bound_fr_woe_sum{1}+fliplr(bound_fr_woe_sum{2});
    bound_fr_woe_cdf = cumsum(reshape(bound_fr_woe_combined,[],1));
    bound_fr_woe_cdf = bound_fr_woe_cdf./bound_fr_woe_cdf(end);
    num_pts = 1e4;
    choice_picked = nan(num_pts,1);
    woe_picked = nan(num_pts,1);
    fr_picked = nan(num_pts,1);
    for i = 1:num_pts
        pick = find(bound_fr_woe_cdf>rand(1),1,'first');
        choice_picked(i) = ceil(pick/(n_row*n_col));
        col_picked = ceil(pick/n_row);
        woe_picked(i) = woe_full(mod(col_picked-1,n_col)+1);
        fr_picked(i) = fr_conv(mod(pick-1,n_row)+1);
    end
    for i = 1:1
        pick = logical(choice_picked==i);
        [beta bint r rint stats]=regress(fr_picked(pick), [ones(sum(pick),1),woe_picked(pick)], 0.3173);
        slopes(i) = beta(2);
        slopes_ci(i) = (bint(2,2)-bint(2,1))/2;
        p_values(i) = stats(3);

        offsets(i) = beta(1);
        offsets_ci(i) = (bint(1,2)-bint(1,1))/2;
    end
    
    % Dividing trials into 5 quantiles
    M = sortrows([woe_picked,fr_picked],1);
    num_group = 5;
    num_per_group = floor(num_pts/num_group);
    for i = 1:num_group
        pick = (num_per_group*(i-1)+1):(num_per_group*i);
        mean_woe(i) = mean(M(pick,1));
        mean_fr(i) = mean(M(pick,2));
    end
            
    pick = logical(sum(bound_fr_woe_combined,1)>0.001);
    imagesc(woe_full,fr_conv,bound_fr_woe_sum{1}+fliplr(bound_fr_woe_sum{2}));
    plot(woe_full(pick),woe_full(pick)*slopes+offsets,'-r','LineWidth',5);
    plot(mean_woe,mean_fr,'ko','MarkerFaceColor','k');
    xlim([-2 4])
    ylim([30 60])
    axis xy
%%    
    figure(fig+8);clf;hold on
    for i = 1:length(woe_full)
        p_r(i) = sum(bound_fr_woe_combined(:,i),1)./...
            sum((bound_fr_woe_combined(:,i)+bound_fr_woe_combined(:,length(woe_full)-(i-1))),1);
    end
    pick = logical(sum((bound_fr_woe_combined+fliplr(bound_fr_woe_combined))./2,1)>0.001);
    plot(woe_full(pick),p_r(pick),'ro','MarkerSize',4,'MarkerFaceColor','r')
    xlabel('P(R)')
    ylabel('total logLR')
    
    
    p_right = info.LLR_num_right./info.LLR_num_trial;
    plot(info.uniqueLLR/10,p_right,'ko','MarkerSize',4,'MarkerFaceColor','k');
    
    xlim([-3 3])
    
end

return
%% Plotting the error between numerical and analytical solutions
if all(sim_switch)
    figure(fig)
    subplot(2,2,1);hold on
    plot(1:10,mean_x1(1:10),'-b')
    plot(1:10,mean_px1,'-r')
    xlabel('time step')
    ylabel('neuron 1 mean')
    subplot(2,2,2);hold on
    plot(1:10,var_x1(1:10),'-b')
    plot(1:10,var_px1,'-r')
    xlabel('time step')
    ylabel('neuron 1 var')
    subplot(2,2,3);hold on
    plot(1:10,mean_x2(1:10),'-b')
    plot(1:10,mean_px2,'-r')
    xlabel('time step')
    ylabel('neuron 2 mean')
    subplot(2,2,4);hold on
    plot(1:10,var_x2(1:10),'-b')
    plot(1:10,var_px2,'-r')
    xlabel('time step')
    ylabel('neuron 2 var')
    
    figure(fig+1)
    subplot(2,2,1)
    plot(mean_x1(1:10)-mean_px1)
    subplot(2,2,2)
    plot(var_x1(1:10)-var_px1)
    subplot(2,2,3)
    plot(mean_x2(1:10)-mean_px2)
    subplot(2,2,4)
    plot(var_x2(1:10)-var_px2)
    supertitle('Error between numerical and analytical solutions')
end
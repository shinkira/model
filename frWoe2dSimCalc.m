function out = frWoe2dSimCalc(theta,info,varargin)

% This test code compares bounded accumulation between
% Monte-Carlo simulation and Fokker-Plank propagation.
% There is a good agreement between them.

set(0,'defaultaxesfontsize',18);
set(0,'defaulttextfontsize',18);
set(0,'defaultaxesfontweight','bold');
set(0,'defaulttextfontweight','bold');
if 0
    set(0,'defaultaxesfontsize',12);
    set(0,'defaulttextfontsize',12);
end
set(0,'defaultaxestickdir','out');
set(0,'defaultaxesbox','off');
set(0,'defaultFigureColor','w');

% Reset the random number generator
rng(0)

fig_switch = 0;
varargin2V(varargin);
symmetric_flag = 0; % if 1, x1_dfr(ei) = x2_dfr(-ei); if 0, if 1, x1_dfr(ei) = -x2_dfr(ei)

num_step = 10;
% Best fitted parameters for singly stochastic model.
% num_param: 
% 1: dfr_sd (SD of noise in delta FR)
% 2: B_init (initial bound height)
% 3: a (rate of bound/gain change) 
% 4: b (rate of bound/gain change) 
% 5: d (time delay)
% 6: init_dfr_sd (SD of noise in initial state)
% 7: min_fr (reflection bound)
% 8: alpha (self-excitation)
% 9: beta (mutual-inhibition)
dfr_sd = theta(1);
B = theta(2); % threshold
a = theta(3);
b = theta(4);
d = theta(5);
init_dfr_sd = 0; % theta(6);
min_fr = theta(7); % reflecting lower bound
alpha = theta(8); % self excitation
beta = theta(9); % mutual inhibition

burst = 0; % saccadic burst

WOE = -0.9:0.2:0.9;
deci_WOE = -9:2:9;

shape_prob = make790StimProb(1);
cum_shape_prob = cumsum(shape_prob);
cum_shape_prob_flip = cumsum(fliplr(shape_prob));

% As in Fig. 5b, gain slopes were computed piecewise linearly
% for positive and negative logLR

% including N*th shape
% all_epoch_gain = 8.1;
% all_epoch_pos_gain = 13.7;
% all_epoch_neg_gain = 2.3;


% excluding N*th shape
% [derived from normalized firing rate]
% all_epoch_pos_gain = 12.7;
% all_epoch_neg_gain = 1.7;
% [derived from raw firng rate]

if 1
    switch info.id
        case 1
            fr_offset = 20;
            t_nd_sd = 60;
            t_nd_min = 270;
            all_epoch_pos_gain = 16.3;
            all_epoch_neg_gain = 2.3;
            % 1st to N* 
            % dfr_v = [-2.0724   -2.3793   -0.9555   -0.3809  3.7700    6.8718   10.5923   11.6678];
                % 1st to N*-1
            % dfr_v = [-1.7144   -1.7373   -0.5164   -0.1545  3.1733    6.2534    9.5183   10.4254];
            % dfr_v = [-2.7433   -3.3367   -1.6834   -1.2381    4.4063    9.7445   13.4855   13.8751];
                %% N*-1 only (200-250)
            % dfr_v = [-2.4993 -2.4030 -0.2303 0.1015 5.7875 10.3396 15.1972 15.9046];
            if 0
                dfr_v = [-2.8782   -3.4534   -0.9936   -0.2898    6.6742   13.3659   18.3843   19.3156];
                mean_epoch_dfr_v = [3.9, 6.3, 3.3, 4.4, 6.4, 5.1, 6.4, 6.9, -0.3, 5.7];
                % mean_epoch_dfr_v = [3.9, 3, 3.3, 4.4, 6.4, 5.1, 6.4, 6.9, -0.3, 5.7];
            end
                %% N*-1, N*-2
            % dfr_v = [-2.1876   -2.5669   -0.7550   -0.3116    4.3559    8.9964   12.3569   12.7768];
            if 1
                dfr_v = [-2.9878   -3.7069   -1.5048   -1.0551    5.1981   11.4967   15.0701   15.8738];
                mean_epoch_dfr_v = [0.3987    2.0608    2.3714    3.7228    5.3132    4.2332    5.4226    6.1551    0.2048  0.6405];
            end
                %% N*-1, N*-2, N*-3 
            % dfr_v = [-1.9331 -2.3789 -0.7867 -0.4078 3.7139 7.6448 11.1526 11.4609];
                % using 440ms(=200+250) after the shape onset.
            if 0
                dfr_v = [-2.2937   -4.1131   -2.1419   -1.7470    4.8592   10.1247   14.8966   15.3848];
                mean_epoch_dfr_v = [1.4215    1.5817    1.9503    3.7944    4.9660    4.1672    5.2233    5.4376    1.0164   -0.0795];
            end
            % dfr_v = [-2.7433   -3.3367   -1.6834   -1.2381    4.4063    9.7445   13.4855   13.8751];
            
                %% N*, N*-1, N*-2 
            % dfr_v = [-3.7664   -5.0942   -2.7689   -2.0950    5.0231   10.3736   14.2488   15.4976];
                        
            
                % 1st to N*-1
            % mean_epoch_dfr_v = [0.8, 2.2, 1.9, 2.7, 4.4, 3.9, 3.8, 4.7, 1.2, 0.4];
            % mean_epoch_dfr_v = [0.9, 1.6, 2.2, 3.3, 4.9, 4.0, 5.1, 5.8, 0.7, 0.2];
            % mean_epoch_dfr_v = [1.4, 1.6, 2.0, 3.8, 5.0, 4.2, 5.2, 5.4, 1.0, -0.1];
                % N*-1, N*-2, N*-3 
            % mean_epoch_dfr_v = [1.6, 2.3, 1.7, 2.6, 4.2, 3.8, 3.6, 4.6, 1.2, 0.1];
                
                % neutral firing rate @ logLR = 0
                % All trials
            % fr_neutral = [20 21.4, 23.8, 25.5, 27.9, 31.1, 32.8, 33.9, 34.6, 34.7, 30.6];
                % Tin only
            % fr_neutral = [20 21.1, 23.4, 24.9, 26.9, 31.3, 34.9, 37.0, 40.0, 42.0, 43.2];
            
        case 2
            fr_offset = 30;
            t_nd_sd = 60;
            t_nd_min = 200;
            all_epoch_pos_gain = 33;
            all_epoch_neg_gain = 33;
            %% 1st to N*
            % dfr_v = [-6.2310 -7.8367 -8.2396 -5.0321 10.0297 14.0943 14.0162 14.2624];
            
            %% 1st to N*-1
            % dfr_v = [-17.8019 -13.7730 -11.1763 -6.3642 12.4904 18.2569 18.1824 19.7452];
            
            %% N*-1 only
            % dfr_v = [-15.7411 -12.5903 -9.7030 -5.9058 12.7594 18.2565 18.9426 17.7148];
            
            %% N*-1, N*-2, N*-3 
                % using 440ms(=200+250) after the shape onset.
            dfr_v = [-22.4881  -16.1828  -12.6341   -7.4541   13.1486   21.5478   21.1620   18.5215];
            mean_epoch_dfr_v = [4.0311    0.9178   -0.6037   -3.3214    0.8682   -2.1809   -2.8114    2.4595, 0,0];
            
            % 1st to N*-1
            % mean_epoch_dfr_v = [3.3, 2.6, -1.5, -0.5, -1.0, -2.3, 2.2, 0, 0, 0];
            % mean_epoch_dfr_v = [3.3, 2.6, 0, 0, 0, 0, 0, 0, 0, 0];
    end
    % dfr_v = [info.dfr_mean_all(1:4),0,0,info.dfr_mean_all(5:8)];
    dfr_v = [dfr_v(1:4),0,0,dfr_v(5:8)];
    
    
    % figure; hold on
    % plot(WOE,[dfr_v(1:4),nan,nan,dfr_v(7:10)],'ko','MarkerSize',10,'MarkerFaceColor','k')
    
    mean_dfr_v = sum(dfr_v.*(shape_prob + fliplr(shape_prob))/2);
    dfr_v = dfr_v - mean_dfr_v;
    
    if exist('dfr_v2','var')
        dfr_v2 = [dfr_v2(1:4),0,0,dfr_v2(5:8)];
        mean_dfr_v2 = sum(dfr_v2.*(shape_prob + fliplr(shape_prob))/2);
        dfr_v2 = dfr_v2 - mean_dfr_v2;
    end
    % plot(WOE,[dfr_v(1:4),nan,nan,dfr_v(7:10)],'ro','MarkerSize',10,'MarkerFaceColor','r')
    % ylim([-12 12])
else
    fr_offset = 0;
    all_epoch_pos_gain = 10;
    all_epoch_neg_gain = 10;
end
    
    
%%

% Dealing with gain for the gain-change model.
% excluding N* shape
% pos_dfr_gain = 10.5+1*(1:10); % Tin only
pos_dfr_gain = 9.4+1.5*(1:10); % Tin & Tout
neg_dfr_gain = 2.9+0.0*(1:10);

gamma = pos_dfr_gain/all_epoch_pos_gain;
if strcmp(info.model,'gain')
    relative_gain = gamma;
else
%     gamma_norm = gamma./gamma(1);
%     all_epoch_pos_gain = all_epoch_pos_gain.*gamma(1);
%     all_epoch_neg_gain = all_epoch_neg_gain.*gamma(1);
%     urgency = [0 diff(B.*(gamma_norm - 1)./gamma_norm)];
    relative_gain = ones(1,num_step);
end

pos_dfr_offset = 0*ones(1,num_step); % shift in the baseline FR in each epoch (i.e. urgency)
neg_dfr_offset = 0*ones(1,num_step);

% Choose the noise type:
% either 'anti' (anti-correlated) or 'ind' (independent)
noise_type = info.noise_correlation_type;
noise_type = 'ind';

fig = 1;

if 0
    figure;hold on;
    dfr_v(5:6) = nan;
    h = ploterr(WOE,dfr_v,[],[],1,'ko','abshhy',0);
    set(h(1),'MarkerFaceColor','k','MarkerSize',6);
    dfr_v = dfr_v - nanmean(dfr_v);
    h = ploterr(WOE,dfr_v,[],dfr_sd,1,'ro','abshhy',0);
    set(h(1),'MarkerFaceColor','r','MarkerSize',12);
    out = [];
    return
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Numerical method (Monte-Carlo simulation)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preallocation
num_trial = 1e4;
x1_hfr = nan(num_trial,num_step+1);
x2_hfr = nan(num_trial,num_step+1);
x1_fr = nan(num_trial,num_step+1);
x2_fr = nan(num_trial,num_step+1);

cum_woe = nan(num_trial,num_step+1);
x1_correct_fr = nan(num_trial,1);
x2_correct_fr = nan(num_trial,1);
x1_error_fr = nan(num_trial,1);
x2_error_fr = nan(num_trial,1);
result = nan(num_trial,1);
choice = nan(num_trial,1);
num_accum = nan(num_trial,1);
cum_woe_correct = nan(num_trial,1);
cum_woe_error = nan(num_trial,1);
cum_woe_end = nan(num_trial,1);

%%
for k = 1:num_trial
    
    if 0
        u = B.*(1-exp(-a.*((0:num_step)-d)));
        u = a.*((0:num_step)-d);
        u(u<0) = 0;
        urgency = diff(u);
    elseif 0
        % urgency signal with trial-to-trial variation in slope
        % slope value is drawn from a gamma distribution
        a_var = gamrnd(a,b,1);
        % urgency signal is added after a time-offset: d
        u = a_var.*((0:num_step)-d);
        u((0:num_step)<d) = 0;
        % u = zeros(1,11); % without urgency
        urgency = diff(u);
    else
        % urgency = diff(fr_neutral) + randn(1);
        % urgency = 2.7*ones(1,10) + randn(1);
        urgency = mean_epoch_dfr_v;
    end
    
    % if mod(k,1e3)==0
    %     printf('%d\n',k)
    % end

    x1_hfr(k,1) = fr_offset + init_dfr_sd*randn(1);
    x2_hfr(k,1) = fr_offset + init_dfr_sd*randn(1);
    x1_fr(k,1) = x1_hfr(k,1);
    x2_fr(k,1) = x2_hfr(k,1);

    cum_woe(k,1) = 0;
    rev = 20;
    for ei = 1:num_step
            % WOE
            
        if symmetric_flag
            shape_ind = find(rand(1)<=cum_shape_prob,1,'first');
        else
            if k<=num_trial/2   
                shape_ind = find(rand(1)<=cum_shape_prob,1,'first');
            else
                shape_ind = find(rand(1)<=cum_shape_prob_flip,1,'first');
            end
        end
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
        
        if 0
            switch info.model
                case 'simple' % {'simple','urgency'}
                    if pos_ind
                        sig1 = shape_woe * 0.1 * pos_dfr_gain(ei) + noise1;
                        sig2 = -shape_woe * 0.1 * neg_dfr_gain(ei) + noise2;
                    else
                        sig1 = shape_woe * 0.1 * neg_dfr_gain(ei) + noise1;
                        sig2 = -shape_woe * 0.1 * pos_dfr_gain(ei) + noise2;
                    end
                otherwise
                    if pos_ind
                        sig1 = shape_woe * 0.1 * all_epoch_pos_gain + noise1;
                        sig2 = -shape_woe * 0.1 * all_epoch_neg_gain + noise2;
                    else
                        sig1 = shape_woe * 0.1 * all_epoch_neg_gain + noise1;
                        sig2 = -shape_woe * 0.1 * all_epoch_pos_gain + noise2;
                    end
            end
        else
            if 1 % x1_fr(k,ei)<35
                sig1 = dfr_v(shape_ind)  + noise1;
            else
                sig1 = dfr_v2(shape_ind)  + noise1;
            end
            if 1 % x2_fr(k,ei)<35
                if symmetric_flag
                    sig2 = dfr_v(11-shape_ind) + noise2;
                else
                    sig2 = -dfr_v(shape_ind) + noise2;
                end
            else
                if symmetric_flag
                    sig2 = dfr_v2(11-shape_ind) + noise2;
                else
                    sig2 = -dfr_v2(shape_ind) + noise2;
                end
            end
        end

        % x1_hfr(k,ei+1) = alpha * x1_hfr(k,ei) + sig1 + urgency(ei) - beta * x2_hfr(k,ei);
        % x2_hfr(k,ei+1) = alpha * x2_hfr(k,ei) + sig2 + urgency(ei) - beta * x1_hfr(k,ei);

        x1_hfr(k,ei+1) = x1_hfr(k,ei) + sig1 + urgency(ei);
        x2_hfr(k,ei+1) = x2_hfr(k,ei) + sig2 + urgency(ei);

        x1_fr(k,ei+1) = relative_gain(ei)*x1_hfr(k,ei+1);
        x2_fr(k,ei+1) = relative_gain(ei)*x2_hfr(k,ei+1);
        
        
        if x1_fr(k,ei)>30
            burst = 10;
            x1_fr(k,ei+1) = x1_fr(k,ei+1) + burst;
        end
        if x1_fr(k,ei)>30
            burst = 10;
            x2_fr(k,ei+1) = x2_fr(k,ei+1) + burst;
        end
        
        % Assuming the pupulation FR cannot go below min_fr
        % if x1_hfr(k,ei+1)<min_fr
        %     x1_hfr(k,ei+1)=min_fr;
        % end
        % if x2_hfr(k,ei+1)<min_fr
        %     x2_hfr(k,ei+1)=min_fr;
        % end
        
        % Assuming the pupulation FR cannot go below min_fr
        if 1
            if x1_fr(k,ei+1)<min_fr
                x1_fr(k,ei+1)=min_fr;
            end
            if x2_fr(k,ei+1)<min_fr
                x2_fr(k,ei+1)=min_fr;
            end
        end

        % if FR of either neuron reaches the bound...
        if (x1_fr(k,ei+1)>B) || (x2_fr(k,ei+1)>B)
            num_accum(k) = ei;
            cum_woe_end(k) = cum_woe(k,ei+1);
            rt(k) = 250*(ei-1) + 0 + (t_nd_min+t_nd_sd*2) + t_nd_sd*randn(1);
            if 0 % plot an example trial
                figure(1);
                set(gcf,'position',[400 200 200 400])
                subplot(2,1,1);hold on;
                plot(0:10,x1_hfr(k,:)-[0,cumsum(urgency)],'-ko','MarkerFaceColor','k');
                xlim([0 10])
                ylim([0 70])
                subplot(2,1,2);hold on;
                plot(0:10,x2_hfr(k,:)-[0,cumsum(urgency)],'-ko','MarkerFaceColor','k');
                if 1
                    subplot(2,1,1)
                    plot(0:10,x1_fr(k,:),'-ro','MarkerFaceColor','r','LineWidth',2);
                    plot([0 10],[B B],'k--')
                    subplot(2,1,2)
                    plot(0:6,[20*ones(1,6),x2_fr(k,7)],'-o','LineWidth',2,'color',0.7*[1,0,1],'MarkerFaceColor',0.7*[1,0,1]);
                    plot(0:10,x2_fr(k,:),'-ro','MarkerFaceColor','r');
                    plot(6:10,x2_fr(k,7:end),'-ro','LineWidth',2);
                    plot([0 10],[B B],'k--')
                end
                xlim([0 10])
                ylim([0 70])
                out = [];
                return
            end

            if x1_fr(k,ei+1)>x2_fr(k,ei+1)
                % simulate a saccadic burst
                x1_fr(k,ei+1) = x1_fr(k,ei+1) + burst;
                % x1_fr(k,ei+1) = B;
                choice(k) = 1;
                
                if symmetric_flag
                    x1_correct_fr(k) = x1_fr(k,ei+1);
                    x2_correct_fr(k) = x2_fr(k,ei+1);
                    cum_woe_correct(k) = cum_woe(k,ei+1);
                    result(k) = 1;
                else
                    if k <= num_trial/2
                        x1_correct_fr(k) = x1_fr(k,ei+1);
                        x2_correct_fr(k) = x2_fr(k,ei+1);
                        cum_woe_correct(k) = cum_woe(k,ei+1);
                        result(k) = 1;
                    else
                        x1_error_fr(k) = x1_fr(k,ei+1);
                        x2_error_fr(k) = x2_fr(k,ei+1);
                        cum_woe_error(k) = cum_woe(k,ei+1);
                        result(k) = 0;
                    end
                end
                
            elseif x1_fr(k,ei+1)<x2_fr(k,ei+1)
                % simulate a saccadic burst
                x2_fr(k,ei+1) = x2_fr(k,ei+1) + burst;
                % x2_fr(k,ei+1) = B;
                choice(k) = 2;
                
                if symmetric_flag
                    x1_error_fr(k) = x1_fr(k,ei+1);
                    x2_error_fr(k) = x2_fr(k,ei+1);
                    cum_woe_error(k) = cum_woe(k,ei+1);
                    result(k) = 0;
                else
                    if k <= num_trial/2
                        x1_error_fr(k) = x1_fr(k,ei+1);
                        x2_error_fr(k) = x2_fr(k,ei+1);
                        cum_woe_error(k) = cum_woe(k,ei+1);
                        result(k) = 0;
                    else
                        x1_correct_fr(k) = x1_fr(k,ei+1);
                        x2_correct_fr(k) = x2_fr(k,ei+1);
                        cum_woe_correct(k) = cum_woe(k,ei+1);
                        result(k) = 1;
                    end
                end
                    
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
                
                if rand(1)>0.5
                    choice(k) = 1;
                else
                    choice(k) = 2;
                end
            end
            if 0
                if x1_fr(k,ei+1)>B
                    x1_fr(k,ei+1) = B;
                end
                if x2_fr(k,ei+1)>B
                    x2_fr(k,ei+1) = B;
                end
            end
            break
        end
    end

end

bins = 1:10;
num_accum_hist = histc(num_accum,bins);
num_accum_hist_correct = histc(num_accum(result==1),bins)';
num_accum_hist_error = histc(num_accum(result==0),bins)';

n_accum_pdf = info.n_accum./sum(info.n_accum);
n_accum_correct_pdf = info.n_accum_correct./sum(info.n_accum);
n_accum_error_pdf = info.n_accum_wrong./sum(info.n_accum);

% add eps to avoid infinite error
out.RT_dist_correct = num_accum_hist_correct/num_trial + eps;
out.RT_dist_wrong = num_accum_hist_error/num_trial + eps;


% Psychometric curve
unique_cum_woe = unique([cum_woe_correct;cum_woe_error;-cum_woe_correct;-cum_woe_error]);
unique_cum_woe = unique_cum_woe(isfinite(unique_cum_woe));
for wi = 1:length(unique_cum_woe)
    if symmetric_flag
        R_num1 = sum(cum_woe_correct==unique_cum_woe(wi));
        L_num1 = sum(cum_woe_error==unique_cum_woe(wi));

        R_num2 = sum(cum_woe_error==-unique_cum_woe(wi));
        L_num2 = sum(cum_woe_correct==-unique_cum_woe(wi));

        R_num(wi) = R_num1 + R_num2;
        L_num(wi) = L_num1 + L_num2;
    else
        switch info.id
            case 1
                L_num(wi) = sum((cum_woe_end==unique_cum_woe(wi)) & choice==1);
                R_num(wi) = sum((cum_woe_end==unique_cum_woe(wi)) & choice==2);
            case 2
                L_num(wi) = sum((cum_woe_end==unique_cum_woe(wi)) & choice==2);
                R_num(wi) = sum((cum_woe_end==unique_cum_woe(wi)) & choice==1);
        end
    end
    PR(wi) = R_num(wi)/(R_num(wi) + L_num(wi));
end

[betaFit,errFit,xflg,oput] = fminsearch(@(beta) logistFitSK(beta,[R_num',(R_num + L_num)'],unique_cum_woe'/10),1);
out.betaFit = betaFit;

        
    %% Compute descriptive statistics of numerical solution
if fig_switch
    figure; hold on
    for ei = 5
        pick = logical(num_accum==ei & choice==1)
        % plot(x1_fr(1:100,:)')
        % plot(mean(x1_fr(pick,:))')
        plot(x1_fr(pick,:)')
    end
    set(gcf,'position',[200 200 1200 600])
    
    %%
    % RT histogram
    figure(fig);clf
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
    bar((1:10)+shift,num_accum_hist/num_trial,bar_width,'r')
    ylabel('Prob','FontSize',24,'FontWeight','bold')
    set(gca,'XTickMode','manual','XTick',1:10);
    set(gca,'FontSize',18,'FontWeight','bold','Box','OFF','TickDir','out');
    y = get(gca,'Ylim');
    xlim([0 11])

    subplot(3,1,2);hold on;
    bar((1:10)-shift,n_accum_correct_pdf(1:10),bar_width,'k')
    bar((1:10)+shift,num_accum_hist_correct/num_trial,bar_width,'r')
    ylabel('Prob','FontSize',24,'FontWeight','bold')
    set(gca,'XTickMode','manual','XTick',1:10);
    set(gca,'FontSize',18,'FontWeight','bold','Box','OFF','TickDir','out');
    xlim([0 11])
    ylim(y);

    subplot(3,1,3);hold on;
    bar((1:10)-shift,n_accum_error_pdf(1:10),bar_width,'k')
    bar((1:10)+shift,num_accum_hist_error/num_trial,bar_width,'r')
    xlabel('Number of shapes used for decision','FontSize',24,'FontWeight','bold')
    ylabel('Prob','FontSize',24,'FontWeight','bold')
    set(gca,'XTickMode','manual','XTick',1:10);
    set(gca,'FontSize',18,'FontWeight','bold','Box','OFF','TickDir','out');
    xlim([0 11])
    ylim(y);
    
    drawnow;
    % return
    fig = fig+1;

    %% Format figures
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
        
        % Survivor function
    figure(fig+3)
    set(gcf,'position',[500 550,500,400],'color','w')
    res_trial = num_trial-cumsum(num_accum_hist);
    bar(1:10,res_trial)
%%
        % mean and s.d. of bound reaching WOE
    figure(fig+4);hold on
    set(gcf,'position',[1000,550,800,500],'color','w')
    
    if symmetric_flag
        for ei = 1:10
            mean_cum_woe_correct(ei) = nanmean(cum_woe_correct(num_accum==ei));
            sd_cum_woe_correct(ei) = nanstd(cum_woe_correct(num_accum==ei));
            mean_cum_woe_error(ei) = nanmean(-cum_woe_error(num_accum==ei));
            sd_cum_woe_error(ei) = nanstd(-cum_woe_error(num_accum==ei));
            mean_cum_woe_all(ei) = nanmean([cum_woe_correct(num_accum==ei);-cum_woe_error(num_accum==ei)]);
            sd_cum_woe_all(ei) = nanstd([cum_woe_correct(num_accum==ei);-cum_woe_error(num_accum==ei)]);
        end
    else
        for ei = 1:10
            mean_cum_woe_correct(ei) = nanmean(cum_woe_end(num_accum==ei & result==1 & choice==1));
            sd_cum_woe_correct(ei) = nanstd(cum_woe_end(num_accum==ei & result==1 & choice==1));
            mean_cum_woe_error(ei) = nanmean(cum_woe_end(num_accum==ei & result==0 & choice==1));
            sd_cum_woe_error(ei) = nanstd(cum_woe_end(num_accum==ei & result==0 & choice==1));
            mean_cum_woe_all(ei) = nanmean(cum_woe_end(num_accum==ei & choice==1));
            sd_cum_woe_all(ei) = nanstd(cum_woe_end(num_accum==ei & choice==1));
        end
    end

    ploterr([1:num_step]-0.3,mean_cum_woe_correct/10,[],{mean_cum_woe_correct/10,mean_cum_woe_correct/10+sd_cum_woe_correct/10},2,'-r','abshhy',0)
    ploterr([1:num_step]+0.1,mean_cum_woe_error/10,[],{mean_cum_woe_error/10-sd_cum_woe_error/10,mean_cum_woe_error/10},2,'-b','abshhy',0)
    ploterr([1:num_step]-0.1,mean_cum_woe_all/10,[],{mean_cum_woe_all/10-sd_cum_woe_all/10,mean_cum_woe_all/10},2,'-k','abshhy',0);
    
    % ploterr([1:10]-0.2,mean_cum_woe_correct/10,[],sd_cum_woe_correct/10,1,'m-.','abshhy',0.1);
    % ploterr([1:10]-0.1,mean_cum_woe_error/10,[],sd_cum_woe_error/10,1,'c-.','abshhy',0.1);
    
    % errorbar([1:10]-0.2,mean_cum_woe_correct/10,sd_cum_woe_correct/10,'m-','LineWidth',2);
    % errorbar([1:10]-0.1,mean_cum_woe_error/10,sd_cum_woe_error/10,'c-','LineWidth',2);
    % errorbar([1:10]-0.1,mean_cum_woe_all/10,sd_cum_woe_all/10,'k-','LineWidth',2);

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
    
    ploterr([1:10]-0.2,meanCumLLR_correct,[],{meanCumLLR_correct+stdCumLLR_correct,meanCumLLR_correct},2,'--m','abshhy',0)
    ploterr([1:10]+0.2,meanCumLLR_error,[],{meanCumLLR_error-stdCumLLR_error,meanCumLLR_error},2,'--c','abshhy',0)
    h=ploterr([1:10],meanCumLLR_all,[],{meanCumLLR_all-stdCumLLR_all,meanCumLLR_all},2,'--k','abshhy',0);
    set(h(1),'color',[0.8,0.8,0.8])
    set(h(2),'color',[0.8,0.8,0.8])
    
    % errorbar([1:10],meanCumLLR_all,stdCumLLR_all,'--k',linewidth,2)
    % errorbar([1:10],meanCumLLR_correct,stdCumLLR_correct,'--r',linewidth,2)
    % errorbar([1:10],meanCumLLR_error,stdCumLLR_error,'--b',linewidth,2)

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
    set(gcf,'position',[500 50,800,400],'color','w')
    T_in_out = 1; % 1: Tin, 0: Tout
    group = 8;
    
    % color_map = colormap(jet(group));
    interval = floor((64-1)/(group-1));
    map = colormap;
    color_map = map(1:interval:(1+interval*(group-1)),:);
    
    if 0
        FM = [x1_fr;x2_fr];
        WM = [cum_woe;-cum_woe];
    else
        FM = x1_fr;
        WM = cum_woe;
    end
    
    for ei = 2:11
        pick = isfinite(FM(:,ei));
        num_per_group = floor(sum(pick)/group);
        M = [FM(pick,ei),WM(pick,ei)];
        M = sortrows(M,2);
        for j = 1:group
            pick_group = (num_per_group*(j-1)+1):(num_per_group*j);
            fr_mean(j) = mean(M(pick_group,1));
            fr_sd(j) = std(M(pick_group,1));
            fr_se(j) = fr_sd(j)/sqrt(num_per_group);
            mean_woe(j) = mean(M(pick_group,2));
            
            % model prediction
            subplot(2,5,ei-1);hold on;
            h = ploterr(mean_woe(j)/10,fr_mean(j),[],fr_se(j),1,'ko');
            set(h(1),'MarkerSize',8,'color','r','MarkerFaceColor','r');
            % set(h(1),'MarkerSize',8,'color',color_map(j,:),'MarkerFaceColor',color_map(j,:));
            
            % observed firing rate
            plot(info.fr_woe(ei-1,j),info.fr_mean(ei-1,j),'ko','MarkerSize',8,'MarkerFaceColor','k')
            % plot(info.fr_woe(ei-1,j),info.fr_mean(ei-1,j),'ko','color',color_map(j,:)*0.5,'MarkerFaceColor',color_map(j,:)*0.5)
            
        end
        
        [beta bint r rint stats] = regress([x1_fr(:,ei);x2_fr(:,ei)],[ones(size(cum_woe,1)*2,1),[cum_woe(:,ei);-cum_woe(:,ei)]],0.05);
        offset(ei) = beta(1);
        offset_ci(ei) = (bint(1,2)-bint(1,1))/2;
        slope(ei) = beta(2)*10;
        slope_ci(ei) = (bint(2,2)-bint(2,1))/2*10;
        p_values(ei) = stats(3);
        
        % model prediction
        plot(mean_woe/10,slope(ei).*mean_woe./10+offset(ei),'r-');
        
        % observed firing rate
        plot(mean_woe/10,info.fr_slope(ei-1).*mean_woe./10+info.fr_offset(ei-1),'k-');
        
        
        axis([-2.5 2.5 0 65]);
        maxY = 65;
        text(-2,maxY*1,sprintf('Slope'),'FontSize',12,'FontWeight','bold','Color','r');
        text(-2,maxY*0.9,sprintf('%.1f\\pm%.1f',slope(ei),slope_ci(ei)),'FontSize',12,'FontWeight','bold','Color','r');
        text(-2,maxY*0.8,sprintf('Offset'),'FontSize',12,'FontWeight','bold','Color','r');
        text(-2,maxY*0.7,sprintf('%.1f\\pm%.1f',offset(ei),offset_ci(ei)),'FontSize',12,'FontWeight','bold','Color','r');
        
        text(0.5,maxY*0.4,sprintf('Slope'),'FontSize',12,'FontWeight','bold');
        text(0.5,maxY*0.3,sprintf('%.1f\\pm%.1f',info.fr_slope(ei-1),info.fr_slope_ci(ei-1)),'FontSize',12,'FontWeight','bold');
        text(0.5,maxY*0.2,sprintf('Offset'),'FontSize',12,'FontWeight','bold');
        text(0.5,maxY*0.1,sprintf('%.1f\\pm%.1f',info.fr_offset(ei-1),info.fr_offset_ci(ei-1)),'FontSize',12,'FontWeight','bold');
    end
    
    if 0    
        for ei = 2:11
            unique_cum_woe = unique([cum_woe(:,ei);-cum_woe(:,ei)]);
            unique_cum_woe = unique_cum_woe(isfinite(unique_cum_woe));
            for wi = 1:length(unique_cum_woe)
                pick1 = logical(cum_woe(:,ei)==unique_cum_woe(wi));
                pick2 = logical(-cum_woe(:,ei)==unique_cum_woe(wi));
                % pick1 = logical(cum_woe(:,ei)==unique_cum_woe(wi) & result==T_in_out);
                % pick2 = logical(-cum_woe(:,ei)==unique_cum_woe(wi) & result==(1-T_in_out));
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
            axis([-2.5 2.5 0 65]);
            maxY = 65;
            text(-2,maxY*1,sprintf('Slope'),'FontSize',12,'FontWeight','bold');
            text(-2,maxY*0.85,sprintf('%.1f\\pm%.1f',slope(ei),slope_ci(ei)),'FontSize',16,'FontWeight','bold');
            text(0,maxY*0.25,sprintf('Offset'),'FontSize',12,'FontWeight','bold');
            text(0,maxY*0.1,sprintf('%.1f\\pm%.1f',offset(ei),offset_ci(ei)),'FontSize',16,'FontWeight','bold');
        end
    end
%%      
    % delta_fr vs delta_logLR

    figure(fig+7)
    set(gcf,'position',[1000 50,800,400],'color','w')


    T_in_out = 0; % 1: Tin, 0: Tout
    
    switch info.id
        case 1
            max_ei = 10;
        case 2
            max_ei = 8;
    end
    
    for ei = 2:max_ei
        all_delta_fr = [];
        all_delta_woe = [];
        for wi = 1:length(deci_WOE)
            if symmetric_flag
                delta_woe = cum_woe(:,ei) - cum_woe(:,ei-1);
                pick1 = logical(delta_woe==deci_WOE(wi));
                % pick1 = logical(delta_woe==deci_WOE(wi) & result==T_in_out);
                delta_x1_fr = x1_fr(pick1,ei)-x1_fr(pick1,ei-1);
                pick2 = logical(delta_woe==-deci_WOE(wi));
                % pick2 = logical(delta_woe==-deci_WOE(wi) & result==(1-T_in_out));
                delta_x2_fr = x2_fr(pick2,ei)-x2_fr(pick2,ei-1);
                mean_delta_fr(ei-1,wi) = nanmean([delta_x1_fr;delta_x2_fr]);
                sd_delta_fr(ei-1,wi) = nanstd([delta_x1_fr;delta_x2_fr]);
                se_delta_fr(ei-1,wi) = sd_delta_fr(ei-1,wi)/sqrt(sum(pick1)+sum(pick2));
                all_delta_fr = [all_delta_fr;delta_x1_fr;delta_x2_fr];
                all_delta_woe = [all_delta_woe;deci_WOE(wi)*ones(sum(pick1)+sum(pick2),1)];
            else
                delta_woe = cum_woe(:,ei) - cum_woe(:,ei-1);
                pick1 = logical(delta_woe==deci_WOE(wi));
                % pick1 = logical(delta_woe==deci_WOE(wi) & result==T_in_out);
                delta_x1_fr = x1_fr(pick1,ei)-x1_fr(pick1,ei-1);
                
                mean_delta_fr(ei-1,wi) = nanmean(delta_x1_fr);
                sd_delta_fr(ei-1,wi) = nanstd(delta_x1_fr);
                se_delta_fr(ei-1,wi) = sd_delta_fr(ei-1,wi)/sqrt(sum(pick1));
                all_delta_fr = [all_delta_fr;delta_x1_fr];
                all_delta_woe = [all_delta_woe;deci_WOE(wi)*ones(sum(pick1),1)];
            end
        end
        if sum(isfinite(all_delta_fr))<10
            continue
        end
        [beta bint r rint stats] = regress(all_delta_fr,[ones(size(all_delta_fr)),all_delta_woe],0.05);
        offset(ei) = beta(1);
        offset_ci(ei) = (bint(1,2)-bint(1,1))/2;
        slope(ei) = beta(2)*10;
        slope_ci(ei) = (bint(2,2)-bint(2,1))/2*10;
        p_values(ei) = stats(3);

        subplot(2,5,ei-1);hold on;
        % model prediction
        h = ploterr(WOE,mean_delta_fr(ei-1,:),[],se_delta_fr(ei-1,:),1,'or','abshhy',0);
        set(h(1),'MarkerFaceColor','r','MarkerSize',8);
        plot(deci_WOE/10,slope(ei).*deci_WOE./10+offset(ei),'r-');
        
        % observed firing rate
        h = ploterr(info.dfr_woe(ei-1,:),info.dfr_mean(ei-1,:),[],info.dfr_se(ei-1,:),1,'ko','abshhy',0);
        set(h(1),'MarkerFaceColor','k','MarkerSize',8);
        plot(deci_WOE/10,info.dfr_slope(ei-1).*deci_WOE./10+info.dfr_offset(ei-1),'k-');
        
        axis([-1 1 -20 30])
        y = get(gca,'Ylim');
        text(-0.9,y(1)+diff(y)*1,sprintf('Slope'),'FontSize',12,'FontWeight','bold','Color','r');
        text(-0.9,y(1)+diff(y)*0.9,sprintf('%.1f\\pm%.1f',slope(ei),slope_ci(ei)),'FontSize',12,'FontWeight','bold','Color','r');
        text(-0.9,y(1)+diff(y)*0.8,sprintf('Offset'),'FontSize',12,'FontWeight','bold','Color','r');
        text(-0.9,y(1)+diff(y)*0.7,sprintf('%.1f\\pm%.1f',offset(ei),offset_ci(ei)),'FontSize',12,'FontWeight','bold','Color','r');
        
        text(0,y(1)+diff(y)*0.4,sprintf('Slope'),'FontSize',12,'FontWeight','bold');
        text(0,y(1)+diff(y)*0.3,sprintf('%.1f\\pm%.1f',info.dfr_slope(ei-1),info.dfr_slope_ci(ei-1)),'FontSize',12,'FontWeight','bold');
        text(0,y(1)+diff(y)*0.2,sprintf('Offset'),'FontSize',12,'FontWeight','bold');
        text(0,y(1)+diff(y)*0.1,sprintf('%.1f\\pm%.1f',info.dfr_offset(ei-1),info.dfr_offset_ci(ei-1)),'FontSize',12,'FontWeight','bold');
        
    end
%%
    if 0
        % delta fr as a function of cumLLR
        
        figure(fig+8)
        set(gcf,'position',[1500 50,500,400],'color','w')
        T_in_out = 1;
        for ei = 2:11
            all_delta_fr = [];
            all_delta_woe = [];
            wi = find(deci_WOE==-3);
            cum_woe_list = (-9*(ei-2)):2:(9*(ei-2));
            for cwi = 1:length(cum_woe_list)
                delta_woe = cum_woe(:,ei) - cum_woe(:,ei-1);
                % pick1 = logical(delta_woe==deci_WOE(wi) & cum_woe(:,ei-1)==cum_woe_list(cwi));
                pick1 = logical(delta_woe==deci_WOE(wi) & cum_woe(:,ei-1)==cum_woe_list(cwi) & result==T_in_out);
                delta_x1_fr = x1_fr(pick1,ei)-x1_fr(pick1,ei-1);
                % pick2 = logical(delta_woe==-deci_WOE(wi) & cum_woe(:,ei-1)==-cum_woe_list(cwi));
                % delta_x2_fr = x2_fr(pick2,ei)-x2_fr(pick2,ei-1);
                pick2 = 0;
                delta_x2_fr = [];
                
                mean_delta_fr(ei-1,cwi) = nanmean([delta_x1_fr;delta_x2_fr]);
                sd_delta_fr(ei-1,cwi) = nanstd([delta_x1_fr;delta_x2_fr]);
                se_delta_fr(ei-1,cwi) = sd_delta_fr(ei-1,cwi)/sqrt(sum(pick1)+sum(pick2));
                all_delta_fr = [all_delta_fr;delta_x1_fr;delta_x2_fr];
                all_delta_woe = [all_delta_woe;cum_woe_list(cwi)*ones(sum(pick1)+sum(pick2),1)];
            end
            
            if sum(isfinite(all_delta_fr))<10
                continue
            end
            [beta bint r rint stats] = regress(all_delta_fr,[ones(size(all_delta_fr)),all_delta_woe],0.05);
            offset(ei) = beta(1);
            offset_ci(ei) = (bint(1,2)-bint(1,1))/2;
            slope(ei) = beta(2)*10;
            slope_ci(ei) = (bint(2,2)-bint(2,1))/2*10;
            p_values(ei) = stats(3);
            
            subplot(2,5,ei-1);hold on;
            ploterr(cum_woe_list/10,mean_delta_fr(ei-1,:),[],se_delta_fr(ei-1,:),1,'ok','abshhy',0);
            plot(cum_woe_list/10,slope(ei).*cum_woe_list./10+offset(ei),'k--');
            axis([-2 2 -20 20])
            y = get(gca,'Ylim');
            text(0,y(1)+diff(y)*0.9,sprintf('Slope'),'FontSize',16,'FontWeight','bold');
            text(0,y(1)+diff(y)*0.8,sprintf('%.1f\\pm%.1f',slope(ei),slope_ci(ei)),'FontSize',16,'FontWeight','bold');
            text(0,y(1)+diff(y)*0.6,sprintf('Offset'),'FontSize',16,'FontWeight','bold');
            text(0,y(1)+diff(y)*0.5,sprintf('%.1f\\pm%.1f',offset(ei),offset_ci(ei)),'FontSize',16,'FontWeight','bold');
        end
        
    end
        
%%        
    % Psychometric curve
    figure(fig+9);hold on;
    set(gcf,'position',[2000 50,500,400],'color','w')
    
    % These lines are already computed above
    % 
    % unique_cum_woe = unique([cum_woe_correct;cum_woe_error;-cum_woe_correct;-cum_woe_error]);
    % unique_cum_woe = unique_cum_woe(isfinite(unique_cum_woe));
    % for wi = 1:length(unique_cum_woe)
    %     R_num1 = sum(cum_woe_correct==unique_cum_woe(wi));
    %     L_num1 = sum(cum_woe_error==unique_cum_woe(wi));
    % 
    %     R_num2 = sum(cum_woe_error==-unique_cum_woe(wi));
    %     L_num2 = sum(cum_woe_correct==-unique_cum_woe(wi));
    % 
    %     R_num(wi) = R_num1 + R_num2;
    %     L_num(wi) = L_num1 + L_num2;
    %     PR(wi) = R_num(wi)/(R_num(wi) + L_num(wi));
    % end
    % 
    % [betaFit,errFit,xflg,oput] = fminsearch(@(beta) logistFitSK(beta,[R_num',(R_num + L_num)'],unique_cum_woe'/10),1);
    
    maxLLR = abs(info.uniqueLLR);
    cum_woe_axis = -maxLLR:maxLLR;
    p_pred = 1./(1+10.^(-betaFit.*(cum_woe_axis/10)))';
    
    plot(info.uniqueLLR/10,info.PR,'ko','MarkerFaceColor','k');
    
    if symmetric_flag
        plot(unique_cum_woe/10,PR,'ro','MarkerFaceColor','r');
        plot(cum_woe_axis/10,p_pred,'r--');
    else
        switch info.id
            case 1
                plot(-unique_cum_woe/10,PR,'ro','MarkerFaceColor','r');
                plot(-cum_woe_axis/10,p_pred,'r--');
            case 2
                plot(unique_cum_woe/10,PR,'ro','MarkerFaceColor','r');
                plot(cum_woe_axis/10,p_pred,'r--');
        end
    end
    xlabel('total logLR','FontSize',24,'FontWeight','bold')
    ylabel('Proportion of R choice','FontSize',24,'FontWeight','bold')
    set(gca,'FontSize',18,'FontWeight','bold','Box','OFF','TickDir','out');
    axis([-3 3 0 1]);
    drawnow
    num_fig = 10;
    
    fig = fig+num_fig;
    
%% RT distribution
    figure(fig+10);clf;hold on;
    set(gcf,'position',[0,0,1000,400],'color','w')
    rt_obs_hist = histc(info.RT,0:3000);
    rt_pred_hist = histc(rt,0:3000);
    
    p_rt_obs_hist = rt_obs_hist./sum(rt_obs_hist);
    p_rt_pred_hist = rt_pred_hist./sum(rt_pred_hist);
    p_rt_pred_hist(1) = 0; % force to zero

    % fill(0:3000,rt_obs_hist,'k');
    % fill(0:3000,rt_pred_hist,'r');
    subplot(2,1,1)
    fill(0:3000,p_rt_obs_hist,'k');
    subplot(2,1,2); hold on
    fill(0:3000,p_rt_obs_hist,'k');
    fill(0:3000,p_rt_pred_hist,'r','EdgeColor','r');
    % set(h,'EdgeColor','r')
    
    figure(fig+11);clf;hold on;
    plot(0:3000,cumsum(p_rt_obs_hist),'k-')
    plot(0:3000,cumsum(p_rt_pred_hist),'r-')
    

end

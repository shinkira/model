function out = frWoe2dDynSimCalc(theta,info,varargin)

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

num_trial = 5e3;
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
init_dfr_sd = theta(6);
min_fr = theta(7); % reflecting lower bound
alpha = theta(8); % self excitation
beta = theta(9); % mutual inhibition

dt = 10;
t = 0:dt:2500;
tau_r = 100;
burst = 0; % saccadic burst

WOE = -0.9:0.2:0.9;
deci_WOE = -9:2:9;

shape_prob = make790StimProb(1);
cum_shape_prob = cumsum(shape_prob);
cum_shape_prob_flip = cumsum(fliplr(shape_prob));

randX = rand(num_trial,num_step);
shape_woe = nan(size(randX));
shape_ind = nan(size(randX));
for i = 1:length(shape_prob);
    shape_woe(randX>0 & randX<=shape_prob(i)) = deci_WOE(i);
    shape_ind(randX>0 & randX<=shape_prob(i)) = i;
    randX = randX - shape_prob(i);
end

if 1
    % Shapes are sampled from the other distribution for the latter half of trials.
    shape_woe(num_trial/2+1:end,:) = -shape_woe(num_trial/2+1:end,:);
    shape_ind(num_trial/2+1:end,:) = 11-shape_ind(num_trial/2+1:end,:);
end


switch info.id
    case 1
        fr_offset = 20;
        t_nd_sd = 60;
        t_nd_min = 270;
        tnd_pre = 200; % lag between the shape onset and when the FR starts reflecting logLR
        tnd_post = 200; % lag between when the FR reaches to threshold and time of saccade

            %% N*-1 only (200-250)
        if 0
            dfr_v = [-2.8782   -3.4534   -0.9936   -0.2898    6.6742   13.3659   18.3843   19.3156];
            mean_epoch_dfr_v = [3.9, 6.3, 3.3, 4.4, 6.4, 5.1, 6.4, 6.9, -0.3, 5.7];
        end
            %% N*-1, N*-2
        if 1
            dfr_v = [-2.9878   -3.7069   -1.5048   -1.0551    5.1981   11.4967   15.0701   15.8738];
            mean_epoch_dfr_v = [0.3987    2.0608    2.3714    3.7228    5.3132    4.2332    5.4226    6.1551    0.2048  0.6405];
        end
            %% N*-1, N*-2, N*-3 
            % using 440ms(=200+250) after the shape onset.
        if 0
            dfr_v = [-2.2937   -4.1131   -2.1419   -1.7470    4.8592   10.1247   14.8966   15.3848];
            mean_epoch_dfr_v = [1.4215    1.5817    1.9503    3.7944    4.9660    4.1672    5.2233    5.4376    1.0164   -0.0795];
        end

    case 2
        fr_offset = 30;
        t_nd_sd = 60;
        t_nd_min = 200;
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

end

dfr_v = [dfr_v(1:4),0,0,dfr_v(5:8)];

mean_dfr_v = sum(dfr_v.*(shape_prob + fliplr(shape_prob))/2);
dfr_v = dfr_v - mean_dfr_v;

% Choose the noise type:
% either 'anti' (anti-correlated) or 'ind' (independent)
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

x1_fr = nan(num_trial,num_step+1);
x2_fr = nan(num_trial,num_step+1);

cum_woe = nan(num_trial,num_step+1);
x1_correct_fr = nan(num_trial,1);
x2_correct_fr = nan(num_trial,1);
x1_error_fr = nan(num_trial,1);
x2_error_fr = nan(num_trial,1);
result = nan(num_trial,1);
choice = nan(num_trial,1);
GoRT = nan(num_trial,1);
rt = nan(num_trial,1);
num_accum = nan(num_trial,1);
cum_woe_correct = nan(num_trial,1);
cum_woe_error = nan(num_trial,1);
cum_woe_end = nan(num_trial,1);

r1 = nan(num_trial,length(t));
r2 = nan(num_trial,length(t));

%%
urgency = cumsum([0,mean_epoch_dfr_v]);
urgency_dt = [];
for i = 1:10
    temp = linspace(urgency(i),urgency(i+1),26);
    urgency_dt = [urgency_dt,temp(1:25)];
end

shift_num = tnd_pre/10;
urgency_dt = circshift([urgency_dt,zeros(1,100)],[0,shift_num]);
urgency_dt = urgency_dt(1:251);

for k = 1:num_trial
    if mod(k,1e3)==0
        printf('%d\n',k)
    end

    x1_fr(k,1) = fr_offset + init_dfr_sd*randn(1);
    x2_fr(k,1) = fr_offset + init_dfr_sd*randn(1);
    
    plot_flag = 1;
    
    % hidden rate
    
    hr1 = nan(1,length(t));
    hr2 = nan(1,length(t));
    hr1(1) = x1_fr(k,1);
    hr2(1) = x2_fr(k,1);
    r1(k,1) = hr1(1);
    r2(k,1) = hr2(1);
    
    % observable rate (urgency added)

    % r1 = nan(1,length(t));
    % r2 = nan(1,length(t));

    cum_woe(k,1) = 0;
    N = 0;
    ei = 0;
    
    for ti = 1:length(t)-1
        
        if mod((t(ti) - tnd_pre),250)==0
            tic
            ei = ei+1;
                % FR
            noise1 = dfr_sd * randn(1);
            switch noise_type
                case 'anti'
                    noise2 = -noise1;
                case 'ind'
                    noise2 = dfr_sd * randn(1);
            end
            % toc
            pos_ind = heaviside(shape_woe);
            sig1 = dfr_v(shape_ind(k,ei))  + noise1;
            if symmetric_flag
                sig2 = dfr_v(11-shape_ind(k,ei)) + noise2;
            else
                sig2 = -dfr_v(shape_ind(k,ei)) + noise2;
            end
            x1_fr(k,ei+1) = x1_fr(k,ei) + sig1;
            x2_fr(k,ei+1) = x2_fr(k,ei) + sig2;
            % toc
            % Assuming the pupulation FR cannot go below min_fr
            if 1
                if x1_fr(k,ei+1)<min_fr
                    x1_fr(k,ei+1)=min_fr;
                end
                if x2_fr(k,ei+1)<min_fr
                    x2_fr(k,ei+1)=min_fr;
                end
            end
        end
        % toc
        hr1(ti+1) = hr1(ti) + dt/tau_r*(x1_fr(k,ei+1) - hr1(ti));
        hr2(ti+1) = hr2(ti) + dt/tau_r*(x2_fr(k,ei+1) - hr2(ti));
        
        r1(k,ti+1) = hr1(ti+1) + urgency_dt(ti);
        r2(k,ti+1) = hr2(ti+1) + urgency_dt(ti);
        % toc
        % if FR of either neuron reaches the bound...
        if (r1(k,ti+1)>B) || (r2(k,ti+1)>B)
            rt(k,1) = t(ti+1) + tnd_post;

            if rt(k,1)>=2500
                continue
            end
            
            N = floor(rt(k,1)/250)+1;
            GoRT(k,1) = mod(rt(k,1),250);
            num_accum(k) = ei;
            
            cum_woe(k,1:N) = cumsum(shape_woe(k,1:N));
            cum_woe_end(k) = cum_woe(k,ei);

            if (r1(k,ti+1)>B) > (r2(k,ti+1)>B)
                % simulate a saccadic burst

                choice(k) = 1;
                
                if k <= num_trial/2
                    r1_correct_fr(k) = r1(k,ti+1);
                    r2_correct_fr(k) = r2(k,ti+1);
                    cum_woe_correct(k) = cum_woe(k,ei+1);
                    result(k) = 1; % correct
                else
                    r1_error_fr(k) = r1(k,ti+1);
                    r2_error_fr(k) = r2(k,ti+1);
                    cum_woe_error(k) = cum_woe(k,ei+1);
                    result(k) = 0; % error
                end

            elseif (r1(k,ti+1)>B) < (r2(k,ti+1)>B)
                % simulate a saccadic burst
                
                choice(k) = 2;
                
                if k <= num_trial/2
                    r1_error_fr(k) = r1(k,ti+1);
                    r2_error_fr(k) = r2(k,ti+1);
                    cum_woe_error(k) = cum_woe(k,ei+1);
                    result(k) = 0; % error
                else
                    r1_correct_fr(k) = r1(k,ti+1);
                    r2_correct_fr(k) = r2(k,ti+1);
                    cum_woe_correct(k) = cum_woe(k,ei+1);
                    result(k) = 1; % correct
                end

            else % Tied case: this heppens very rarely, if any.
                if rand(1)>0.5
                    r1_correct_fr(k) = r1(k,ti+1);
                    r2_correct_fr(k) = r2(k,ti+1);
                    cum_woe_correct(k) = cum_woe(k,ei+1);
                    result(k) = 1; % correct
                else
                    r1_error_fr(k) = r1(k,ti+1);
                    r2_error_fr(k) = r2(k,ti+1);
                    cum_woe_error(k) = cum_woe(k,ei+1);
                    result(k) = 0; % error
                end

                if rand(1)>0.5
                    choice(k) = 1;
                else
                    choice(k) = 2;
                end
            end
            
            if 0
                if r1(k,ti+1)>B
                    r1(k,ti+1) = B;
                end
                if r2(k,ti+1)>B
                    r2(k,ti+1) = B;
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

%% RT distribution
    figure(fig+1);clf;hold on;
    set(gcf,'position',[0,0,1000,400],'color','w')
    rt_obs_hist = histc(info.RT,0:3000);
    rt_pred_hist = histc(rt,0:3000);
    
    p_rt_obs_hist = rt_obs_hist./sum(rt_obs_hist);
    cum_p_rt_obs_hist = cumsum(p_rt_obs_hist);
    p_rt_obs_hist_10ms = diff([0;cum_p_rt_obs_hist(t+1)]);
    
    p_rt_pred_hist = rt_pred_hist./sum(rt_pred_hist);
    p_rt_pred_hist(1) = 0; % force to zero

    subplot(2,1,1)
    % fill(0:3000,p_rt_obs_hist,'k');
    fill(t,p_rt_obs_hist_10ms','k-')
    subplot(2,1,2); hold on
    % fill(0:3000,p_rt_obs_hist,'k');
    fill(t,p_rt_obs_hist_10ms','k-')
    fill(0:3000,p_rt_pred_hist,'r','EdgeColor','r');
    % set(h,'EdgeColor','r')
    
    figure(fig+11);clf;hold on;
    plot(0:3000,cumsum(p_rt_obs_hist),'k-')
    plot(0:3000,cumsum(p_rt_pred_hist),'r-')
    
%%  
    % load sim0107.mat
    LLR{1} = diff(cum_woe,1,2);
    LLR{2} = nan(size(LLR{1}));
    for ei = 1:10
        pick = isfinite(LLR{1}(:,ei)); 
        num_shift = 10-ei;
        LLR{2}(pick,:) = fliplr(circshift(LLR{1}(pick,:),[0,num_shift]));
    end

    % *********** Backward Analysis ***********
    lastStim2Sac = GoRT;
    LLR_back = LLR{2}/10;
    shift_size = 10;
    bin_size = 40;
    num_bin = 60;
    trial_id = (1:num_trial)';
    
    for ti = 1:num_bin
        
        pick_before = logical(lastStim2Sac>=0);
        pick_after = logical(lastStim2Sac-bin_size>=0);
        pickDiff = logical((pick_before-pick_after)==1);
        
        % remove shapes for the next iteration
        lastStim2Sac = lastStim2Sac - shift_size;
        pick_res = logical(lastStim2Sac - shift_size>=0);
        pick_remove = logical((pick_before-pick_res)==1);
        lastStim2Sac(pick_remove) = lastStim2Sac(pick_remove) + 250;
        % LLR_picked = LLR_back(pickDiff,1);
        for ri = 1;
            
            pick = logical(pickDiff & isfinite(LLR_back(:,1)) & trial_id<=num_trial/2);
            
            LLR_picked = LLR_back(pick,1);
            resLLR_picked = nansum(LLR_back(pick,2:end),2); % sum of logLR for the rest of shapes
            choice_picked = choice(pick);
            [logitCoef,dev,stats] = glmfit([LLR_picked, resLLR_picked],2-choice_picked,'binomial','link','logit');
            % Convert the log base from exp to 10
            logitCoef = logitCoef/log(10);
            stats.se = stats.se/log(10);
            
            glm_w(ri,ti,:) = logitCoef;
            glm_w_ci(ri,ti,:) = stats.se * norminv(0.975,0,1);
            glm_w_lo(ri,ti,:) = logitCoef - stats.se * norminv(0.975,0,1);
            glm_w_hi(ri,ti,:) = logitCoef + stats.se * norminv(0.975,0,1);
            glm_p(ri,ti,:) = stats.p;
            
        end
        
        % shift the matrix for trials that were picked and replace with NaN.
        LLR_back(pick_remove,:) = circshift(LLR_back(pick_remove,:),[0 -1]);
        LLR_back(pick_remove,end) = nan;
        
    end
    figure(fig);shg
    subplot(1,2,2);hold on;
    t_bin_center = bin_size/2 + shift_size.*((1:num_bin)-1);
    
    if 0
        h1 = ploterr(-t_bin_center,squeeze(glm_w(1,:,2)),[],{squeeze(glm_w_lo(1,:,2)),squeeze(glm_w_hi(1,:,2))},1,'ro-','abshhy',0);
        h2 = ploterr(-t_bin_center,squeeze(glm_w(2,:,2)),[],{squeeze(glm_w_lo(2,:,2)),squeeze(glm_w_hi(2,:,2))},1,'bo-','abshhy',0);
        set(h1(1),'MarkerFaceColor','r');
        set(h2(1),'MarkerFaceColor','b');
    else
        fillTrace(-t_bin_center,squeeze(glm_w(1,:,2)),squeeze(glm_w_ci(1,:,2)),0.2*[1,1,1]);
        % h3 = ploterr(-t_bin_center,squeeze(glm_w(1,:,2)),[],{squeeze(glm_w_lo(1,:,2)),squeeze(glm_w_hi(1,:,2))},1,'ko-','abshhy',0);
        % set(h3(1),'MarkerFaceColor','k');
    end
    hold on
    plot(-t_bin_center,0*ones(1,num_bin),'k--');

    xlim([-600 0]);
    ylim([-0.5 3])
    set(gcf,'position',[500 100 400 300],'color','w');
    set(gca,'FontSize',24,'FontWeight','bold','Box','off','TickDir','out','YTick',[],'YAxisLocation','right')

end

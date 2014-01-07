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
num_trial = 1e4;
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
    
    urgency = mean_epoch_dfr_v;
    
    if mod(k,1e3)==0
        printf('%d\n',k)
    end

    x1_fr(k,1) = fr_offset + init_dfr_sd*randn(1);
    x2_fr(k,1) = fr_offset + init_dfr_sd*randn(1);
    
    r1 = nan(1,length(t));
    r2 = nan(1,length(t));

    plot_flag = 1;
    r1(1) = x1_fr(k,1);
    r2(1) = x2_fr(k,1);

    cum_woe(k,1) = 0;
    N = 0;
    ei = 0;
    
    for ti = 1:length(t)-1
        
        if mod((t(ti) - tnd_pre),250)==0
            N = N+1;
        end
        
        if mod((t(ti) - tnd_pre),250)==0
            
            ei = ei+1;
            
            % A new shape affects the firing rate 200ms after the onset.        

                % WOE
            if k<=num_trial/2   
                shape_ind = find(rand(1)<=cum_shape_prob,1,'first');
            else
                shape_ind = find(rand(1)<=cum_shape_prob_flip,1,'first');
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

            sig1 = dfr_v(shape_ind)  + noise1;
            if symmetric_flag
                sig2 = dfr_v(11-shape_ind) + noise2;
            else
                sig2 = -dfr_v(shape_ind) + noise2;
            end

            x1_fr(k,ei+1) = x1_fr(k,ei) + sig1 + urgency(ei);
            x2_fr(k,ei+1) = x2_fr(k,ei) + sig2 + urgency(ei);

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
        
        r1(ti+1) = r1(ti) + dt/tau_r*(x1_fr(k,ei+1) - r1(ti));
        r2(ti+1) = r2(ti) + dt/tau_r*(x2_fr(k,ei+1) - r2(ti));
        
        % if FR of either neuron reaches the bound...
        if (r1(ti+1)>B) || (r2(ti+1)>B)
            
            rt(k) = t(ti+1) + tnd_post;
            num_accum(k) = ei;
            cum_woe_end(k) = cum_woe(k,ei+1);

            if (r1(ti+1)>B) > (r2(ti+1)>B)
                % simulate a saccadic burst

                choice(k) = 1;
                
                if k <= num_trial/2
                    r1_correct_fr(k) = r1(ti+1);
                    r2_correct_fr(k) = r2(ti+1);
                    cum_woe_correct(k) = cum_woe(k,ei+1);
                    result(k) = 1; % correct
                else
                    r1_error_fr(k) = r1(ti+1);
                    r2_error_fr(k) = r2(ti+1);
                    cum_woe_error(k) = cum_woe(k,ei+1);
                    result(k) = 0; % error
                end

            elseif (r1(ti+1)>B) < (r2(ti+1)>B)
                % simulate a saccadic burst
                
                choice(k) = 2;
                
                if k <= num_trial/2
                    r1_error_fr(k) = r1(ti+1);
                    r2_error_fr(k) = r2(ti+1);
                    cum_woe_error(k) = cum_woe(k,ei+1);
                    result(k) = 0; % error
                else
                    r1_correct_fr(k) = r1(ti+1);
                    r2_correct_fr(k) = r2(ti+1);
                    cum_woe_correct(k) = cum_woe(k,ei+1);
                    result(k) = 1; % correct
                end

            else % Tied case: this heppens very rarely, if any.
                if rand(1)>0.5
                    r1_correct_fr(k) = r1(ti+1);
                    r2_correct_fr(k) = r2(ti+1);
                    cum_woe_correct(k) = cum_woe(k,ei+1);
                    result(k) = 1; % correct
                else
                    r1_error_fr(k) = r1(ti+1);
                    r2_error_fr(k) = r2(ti+1);
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
                if r1(ti+1)>B
                    r1(ti+1) = B;
                end
                if r2(ti+1)>B
                    r2(ti+1) = B;
                end
            end
            break
        end
        
        
        if 0
            % if FR of either neuron reaches the bound...
            if (x1_fr(k,ei+1)>B) || (x2_fr(k,ei+1)>B)
                num_accum(k) = ei;
                cum_woe_end(k) = cum_woe(k,ei+1);
                
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
    figure(fig+10);clf;hold on;
    set(gcf,'position',[0,0,1000,400],'color','w')
    rt_obs_hist = histc(info.RT,0:3000);
    rt_pred_hist = histc(rt,0:3000);
    
    p_rt_obs_hist = rt_obs_hist./sum(rt_obs_hist);
    cum_p_rt_obs_hist = cumsum(p_rt_obs_hist);
    p_rt_obs_hist_10ms = diff([0;cum_p_rt_obs_hist(t+1)]);
    
    p_rt_pred_hist = rt_pred_hist./sum(rt_pred_hist);
    p_rt_pred_hist(1) = 0; % force to zero

    % fill(0:3000,rt_obs_hist,'k');
    % fill(0:3000,rt_pred_hist,'r');
    subplot(2,1,1)
    % fill(0:3000,p_rt_obs_hist,'k');
    fill(t,p_rt_obs_hist_10ms','k-')
    subplot(2,1,2); hold on
    % fill(0:3000,p_rt_obs_hist,'k');
    fill(t,p_rt_obs_hist_10ms','k-')
    fill(0:3000,p_rt_pred_hist,'r','EdgeColor','r');
    % set(h,'EdgeColor','r')
    
    figure(fig+11);clf;hold on;
    % plot(0:3000,cumsum(p_rt_obs_hist),'k-')
    plot(t,p_rt_obs_hist_10ms','k-')
    plot(0:3000,cumsum(p_rt_pred_hist),'r-')
    
end



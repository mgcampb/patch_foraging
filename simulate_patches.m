% simulation parameters
opt.num_trials = 10000;
opt.tau = 80;
opt.iti = 150;
opt.max_patch_time = 500;
opt.patch_time = 0:opt.max_patch_time;

% struct for results
rez = struct;
rez.N0 = unifrnd(0,1,opt.num_trials,1); % draw N0 from uniform distribution
% rez.N0 = unifrnd(0,0.2,opt.num_trials,1) + binornd(1,0.5,opt.num_trials,1)*0.8;
% rez.N0 = binornd(1,0.5,opt.num_trials,1);
% N0 = [0.1 0.5 0.9]';
% rez.N0 = N0(randsample(3,opt.num_trials,true));
rez.rew_prob = rez.N0 * exp(-(0:opt.max_patch_time)/opt.tau);
rez.rew = binornd(1,rez.rew_prob);
rez.cum_rew = cumsum(rez.rew,2);

%% omniscient mouse
rez.omnisc = struct;
rew_rate_omnisc = rez.cum_rew./(opt.iti+repmat(opt.patch_time,opt.num_trials,1));
[~,max_idx_omnisc] = max(rew_rate_omnisc,[],2);
rez.omnisc.leave_time = opt.patch_time(max_idx_omnisc)';
rez.omnisc.rew = rez.cum_rew(sub2ind(size(rez.cum_rew),(1:opt.num_trials)',max_idx_omnisc));
rez.omnisc.total_rew = sum(rez.omnisc.rew);
rez.omnisc.total_time = sum(rez.omnisc.leave_time)+opt.num_trials*opt.iti;
rez.omnisc.rew_rate = rez.omnisc.total_rew./rez.omnisc.total_time;

%% strategy 1: constant wait time
rez.strat1 = struct;
rez.strat1.wait_time = opt.patch_time; % test all possible wait times
rez.strat1.rew = rez.cum_rew(:,1:numel(rez.strat1.wait_time));
rez.strat1.total_rew = sum(rez.strat1.rew);
rez.strat1.total_time = (rez.strat1.wait_time+opt.iti)*opt.num_trials;
rez.strat1.rew_rate = rez.strat1.total_rew./rez.strat1.total_time;

%% strategy 2: constant wait time after each reward
% step through each trial to decide when to leave
rez.strat2 = struct;
rez.strat2.wait_time = opt.patch_time;
rez.strat2.leave_time = opt.max_patch_time * ones(opt.num_trials,numel(rez.strat2.wait_time));
rez.strat2.rew = repmat(sum(rez.rew,2),1,numel(rez.strat2.wait_time));
for i = 1:opt.num_trials
    time_since_last_reward = 0;
    left_already = zeros(size(rez.strat2.wait_time));
    for t = 1:size(rez.rew,2)
        if rez.rew(i,t)==0
            time_since_last_reward = time_since_last_reward+1;
        else
            time_since_last_reward = 0;
        end
        leave = rez.strat2.wait_time <= time_since_last_reward;
        rez.strat2.leave_time(i,leave & ~left_already) = t-1;
        rez.strat2.rew(i,leave & ~left_already) = rez.cum_rew(i,t);
        left_already(leave) = 1;
    end
end
rez.strat2.total_rew = sum(rez.strat2.rew);
rez.strat2.total_time = sum(rez.strat2.leave_time)+opt.iti*opt.num_trials;
rez.strat2.rew_rate = rez.strat2.total_rew./rez.strat2.total_time;

%% strategy 3: integrate to bound
rez.strat3 = struct;
rez.strat3.bound = 1;
rez.strat3.slope = linspace(0,1,20);
rez.strat3.step = linspace(0,-1,20);
rez.strat3.leave_time = nan(opt.num_trials,numel(rez.strat3.slope),numel(rez.strat3.step));
rez.strat3.rew = nan(opt.num_trials,numel(rez.strat3.slope),numel(rez.strat3.step));

for j = 1:numel(rez.strat3.slope)
    for k = 1:numel(rez.strat3.step)
        time_ramp = rez.strat3.slope(j) * repmat(0:opt.max_patch_time,opt.num_trials,1);
        rew_ramp = rez.strat3.step(k) * rez.cum_rew;
        dec_var = time_ramp+rew_ramp;
        for i = 1:opt.num_trials
            leave_idx = find(dec_var(i,:)>=rez.strat3.bound,1);
            if ~isempty(leave_idx)
                rez.strat3.leave_time(i,j,k) = leave_idx-1;
                rez.strat3.rew(i,j,k) = rez.cum_rew(i,leave_idx);
            else
                rez.strat3.leave_time(i,j,k) = opt.max_patch_time;
                rez.strat3.rew(i,j,k) = rez.cum_rew(i,end);
            end
        end
    end
end

rez.strat3.total_rew = squeeze(sum(rez.strat3.rew));
rez.strat3.total_time = squeeze(sum(rez.strat3.leave_time)) + opt.num_trials*opt.iti;
rez.strat3.rew_rate = rez.strat3.total_rew./rez.strat3.total_time;

%% fig: rew, time, and rew rate for different strategies
hfig(1) = figure('Position',[200 200 300 500]);
hfig(1).Name = 'rew, time and rate vs wait time';

subplot(3,1,1); hold on;
plot(rez.strat1.wait_time,rez.strat1.total_rew,'-o')
plot(rez.strat2.wait_time,rez.strat2.total_rew,'-o')
legend({'strat1','strat2'});
ylabel('total rew');

subplot(3,1,2); hold on;
plot(rez.strat1.wait_time,rez.strat1.total_time,'-o')
plot(rez.strat2.wait_time,rez.strat2.total_time,'-o')
legend({'strat1','strat2'});
ylabel('total time');

subplot(3,1,3); hold on;
plot(rez.strat1.wait_time,rez.strat1.rew_rate,'-o')
plot(rez.strat2.wait_time,rez.strat2.rew_rate,'-o')
legend({'strat1','strat2'});
xlabel('wait sec');
ylabel('rew rate');

%% fig: dist of PRT
hfig(2) = figure; hold on;
hfig(2).Name = 'PRT dist for best model 2';
[~,max_idx2] = max(rez.strat2.rew_rate);
histogram(rez.strat2.leave_time(:,max_idx2));
[~,max_idx1] = max(rez.strat1.rew_rate);
plot([rez.strat1.wait_time(max_idx1) rez.strat1.wait_time(max_idx1)],ylim(),'r--');
xlabel('patch residency time (sec)');
ylabel('count');
legend({'best model 2','best model 1'});

%% fig: cum reward over time colored by N0
hfig(3) = figure; hold on;
hfig(3).Name = 'cumul reward over time colored by N0';
for i = 1:500
    plot(opt.patch_time,rez.cum_rew(i,:),'-','Color',[rez.N0(i) 0 0]);
end
map = [(0:0.01:1)' zeros(101,1) zeros(101,1)];
colormap(map);
cb = colorbar;
ylabel(cb,'N_0');
xlabel('sec');
ylabel('cum. rewards');

%% fig: PRT vs N0 for best strat2 model
hfig(4) = figure;
hfig(4).Name = 'PRT vs N0 for best model 2';
[~,max_idx2] = max(rez.strat2.rew_rate);
scatter(rez.N0,rez.strat2.leave_time(:,max_idx2),20,'k','MarkerEdgeAlpha',0.1);
xlabel('N_0');
ylabel('patch residency time (sec)');

%% fig: heatmap of rew rate vs slope/step params for strat3
hfig(5) = figure;
hfig(5).Name = 'heatmap of rew rate vs slope and step for strat3';
imagesc(rez.strat3.slope,rez.strat3.step,rez.strat3.rew_rate);
colorbar;

%% fig: reward matrix
hfig(6) = figure('Position',[400 400 400 600]);
hfig(6).Name = 'reward matrix';
imagesc(opt.patch_time,1:opt.num_trials,rez.rew);
xlabel('time in patch');
ylabel('patch number');
title('reward matrix');

%% fig: best model of each strategy
hfig(7) = figure('Position',[300 300 500 300]);
hfig(7).Name = 'best possible performance for each strategy';
rez.best_rew_rate = [max(rez.strat1.rew_rate),max(rez.strat2.rew_rate),max(max(rez.strat3.rew_rate)),rez.omnisc.rew_rate];
bar(rez.best_rew_rate)
xticklabels({'strat1','strat2','strat3','omnisc'});
ylabel('reward rate');
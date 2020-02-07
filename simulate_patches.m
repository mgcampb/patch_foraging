% simulation parameters
opt.num_trials = 10000;
opt.tau = 8;
opt.iti = 15;
opt.max_patch_time = 50;
opt.patch_time = 0:opt.max_patch_time;

% struct for results
rez = struct;
% rez.N0 = unifrnd(0,1,opt.num_trials,1); % draw N0 from uniform distribution
% rez.N0 = unifrnd(0,0.2,opt.num_trials,1) + binornd(1,0.5,opt.num_trials,1)*0.8;
% rez.N0 = binornd(1,0.5,opt.num_trials,1);
% N0 = [0.1 0.5 0.9]';
N0 = [0.125 0.25 0.5]';
% N0 = [0.1 0.9]';
rez.N0 = N0(randsample(numel(N0),opt.num_trials,true));
rez.rew_prob = [ones(opt.num_trials,1) rez.N0 * exp(-(0:opt.max_patch_time-1)/opt.tau)];
rez.rew = binornd(1,rez.rew_prob);
rez.cum_rew = cumsum(rez.rew,2);

rew_size = [1 2 4];

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
rez.strat3.slope = linspace(0.2,0.3,30);
rez.strat3.step = linspace(-0.2,-0.3,30);
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

%% strategy 4: omnicient (allow different wait time for each N0)
rez.omnisc = struct;
rez.omnisc.wait_time = opt.patch_time; % test all possible wait times
uniqN0 = unique(rez.N0);
rez.omnisc.total_rew = nan(numel(uniqN0),numel(rez.omnisc.wait_time));
for i = 1:numel(uniqN0)
    rez.omnisc.total_rew(i,:) = sum(rez.cum_rew(rez.N0==uniqN0(i),:));
    rez.omnisc.total_time(i,:) = (rez.omnisc.wait_time+opt.iti)*sum(rez.N0==uniqN0(i));
end
allcombs = 1:numel(rez.omnisc.wait_time);
for i = 1:numel(uniqN0)-1
    allcombs = combvec(allcombs,1:numel(rez.omnisc.wait_time));
end
rez.omnisc.rew_rate = nan(size(allcombs,2),1);
for i = 1:numel(rez.omnisc.rew_rate)    
    total_rew_this = 0;
    total_time_this = 0;
    for j = 1:numel(uniqN0)
        total_rew_this = total_rew_this+rez.omnisc.total_rew(j,allcombs(j,i));
        total_time_this = total_time_this+rez.omnisc.total_time(j,allcombs(j,i));
    end
    rez.omnisc.rew_rate(i) = total_rew_this/total_time_this;
end

%% random strategy (shuffle leave times from strat2)
rez.rand = struct();
[~,max_idx2] = max(rez.strat2.rew_rate);
optimal_leave_times = rez.strat2.leave_time(:,max_idx2);
rez.rand.leave_time = optimal_leave_times(randperm(numel(optimal_leave_times)));
rez.rand.rew = nan(size(rez.rand.leave_time));
for i = 1:numel(rez.rand.rew)
    rez.rand.rew(i) = rez.cum_rew(rez.rand.leave_time(i)+1);
end
rez.rand.rew_rate = sum(rez.rand.rew)/(sum(rez.rand.leave_time)+opt.num_trials*opt.iti);

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
scatter(rez.N0+randn(numel(rez.N0),1)*0.01,rez.strat2.leave_time(:,max_idx2),20,'k','MarkerEdgeAlpha',0.1);
xticks(uniqN0);
xlabel('N_0');
ylabel('patch residency time (sec)');

%% fig: heatmap of rew rate vs slope/step params for strat3
hfig(5) = figure;
hfig(5).Name = 'heatmap of rew rate vs slope and step for strat3';
imagesc(rez.strat3.slope,rez.strat3.step,rez.strat3.rew_rate);
xlabel('strat 3 slope');
ylabel('strat 3 step');
hcb = colorbar;
ylabel(hcb,'rew. rate');

%% fig: reward matrix
hfig(6) = figure('Position',[400 400 400 600]);
hfig(6).Name = 'reward matrix';
imagesc(opt.patch_time,1:opt.num_trials,rez.rew);
xlabel('time in patch');
ylabel('patch number');
title('reward matrix');

%% fig: reward matrix sorted by PRT (strat2)
hfig(7) = figure('Position',[400 400 400 600]);
hfig(7).Name = 'reward matrix (sorted by strat2 leave times)';
[~,max_idx2] = max(rez.strat2.rew_rate);
[sorted_leave_time,sort_idx] = sort(rez.strat2.leave_time(:,max_idx2));
imagesc(opt.patch_time,1:opt.num_trials,rez.rew(sort_idx,:));
hold on;
plot(sorted_leave_time,1:opt.num_trials,'r-');
xlabel('time in patch');
ylabel('sorted patch number');
title(sprintf('reward matrix\n(sorted by strat2 leave times)'));

%% fig: best model of each strategy
hfig(8) = figure('Position',[300 300 500 300]);
hfig(8).Name = 'best possible performance for each strategy';
rez.best_rew_rate = [rez.rand.rew_rate,max(rez.strat1.rew_rate),max(rez.strat2.rew_rate),max(max(rez.strat3.rew_rate)),max(rez.omnisc.rew_rate)];
bar(rez.best_rew_rate)
xticklabels({'random','strat1','strat2','strat3','omnisc'});
ylabel('reward rate');
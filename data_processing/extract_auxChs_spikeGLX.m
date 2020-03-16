spikeDir = pwd;
cd ..


%% location of data
aux_dir = pwd;

[~,main_name]=fileparts(aux_dir);
NIDAQ_file = fullfile(aux_dir,strcat(main_name,'_t0.nidq.bin'));
NIDAQ_config = fullfile(aux_dir,strcat(main_name,'_t0.nidq.meta'));

%% get neuropixels sync pulse times
fpNIDAQ=fopen(NIDAQ_file);
datNIDAQ=fread(fpNIDAQ,[9,Inf],'*int16');
fclose(fpNIDAQ);

%% get the nidaq sample rate
dat=textscan(fopen(NIDAQ_config),'%s %s','Delimiter','=');
names=dat{1};
vals=dat{2};
loc=contains(names,'niSampRate');
sync_sampling_rate=str2double(vals{loc});

cd(spikeDir)


ch_vel = 2; % signal used to compute velocity from rotary encoder
ch_licks = 3; % https://media1.tenor.com/images/15a6acd8b109696db10f8cf459d6ce6d/tenor.gif?itemid=16400378
ch_rews = 4; % reward valve opening
ch_events = 5; % 'events' signalled from within DAQ for patch cue appearing, etc

e = datNIDAQ(ch_events,:)/6000;
r = datNIDAQ(ch_rews,:)/6000;

e(1) = 0; % hack to deal w first trial, we fake detect it as the start of openephys and fix it below
t.patchCue = find(e(1:end-9) < 0.5 & e(2:end-8) > 0.5 & e(10:end) > 0.5) + 1;
t.rews = find(r(1:end-9) < 1.0 & r(2:end-8) > 1.0 & r(10:end) > 1.0) + 1;
t.leaves = find(e(1:end-9) > -0.5 & e(2:end-8) < -0.5 & e(10:end) < -0.5) + 1;

% analog signal for leave actually occurs 0.5 seconds after patch is left,
% correct for that here

%t.leaves = t.leaves - 15000;

t.sessionStart = t.rews(1); % session starts w reward delivery

if t.patchCue(1) < t.sessionStart
   t.patchCue(1) = []; 
   display('false start detected before session began - removed')
end

if t.leaves(1) < t.sessionStart
   t.leaves(1) = []; 
   display('false leave detected before session began - removed')
end

t.rews(1) = []; % remove it; it's fake

assert(all(t.patchCue(2:end) > t.patchCue(1:end-1)));
assert(all(t.rews(2:end) > t.rews(1:end-1)));
assert(all(t.leaves(2:end) > t.leaves(1:end-1)));

if t.patchCue(2) < t.leaves(1)
    t.patchCue = t.patchCue(2:end);
end

if length(t.patchCue) > length(t.leaves)
    t.patchCue(end) = [];
end

assert(all(length(t.patchCue) == length(t.leaves)));

didstop = [];
patchStop_indx = [];

maxVoltsID = [];

j = 1; k = 1;

for i = 1:length(t.patchCue)

    %fprintf('i (starts) = %d\n', i);
    volts = max(e(t.patchCue(i):t.patchCue(i)+1000)); % 1000 is arbitrary
    nextRew = min(t.rews(t.rews > t.patchCue(i)));
    
    if ~isempty(nextRew) && nextRew < t.leaves(i)
        
        % mouse stopped on patch
        patchCue_didstop_indx(j) = t.patchCue(i);
        patchLeave_didstop_indx(j) = t.leaves(i);
        patchStop_indx(j) = nextRew;
        eventVolts(j) = volts;
        j = j + 1;
        didstop(i) = 1;
        
    else

        % patch was skipped
        didstop(i) = 0;

    end
    t.leaves(i) = t.leaves(i);
    maxVoltsID(i) = volts;
end

assert(length(t.patchCue) == length(t.leaves));
assert(length(t.patchCue) == length(didstop));
assert(length(t.patchCue) == length(maxVoltsID));


patchCue_ts = t.patchCue / sync_sampling_rate * 1000;
patchLeave_ts = t.leaves / sync_sampling_rate * 1000;
patchStop_ts = patchStop_indx / sync_sampling_rate * 1000;

patchCue_didstop_ts = patchCue_didstop_indx / sync_sampling_rate * 1000;
patchLeave_didstop_ts = patchLeave_didstop_indx / sync_sampling_rate * 1000;

assert(all(patchCue_ts < patchLeave_ts));

trialID = round(maxVoltsID*20)';
trialID = floor(trialID / 10);
trialID(trialID==4)=3;

trialID_didstop = round(eventVolts*20)';
trialID_didstop = floor(trialID_didstop / 10);
trialID_didstop(trialID_didstop==4) = 3;

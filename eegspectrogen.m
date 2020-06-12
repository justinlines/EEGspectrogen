%Justin Lines 2016. this program analyzes individual EEG files for spectral
%analysis and spiking. For cortex change frequencies

%% first

clear
winlength=4096; %used to be 2048 in samples

overlap=0;
welch=5;
seconds=60;


itime=21.137920; %input for injection time in minutes


injection_time=itime*(seconds/60);
[ephysname pathway]=uigetfile(['C:\Users\justi\Desktop\Im workin here\*.abf']);
efile=[pathway ephysname];



%% second

[D,si,h]=abfload(efile);

% file=load('012816.mat');
% D=file.A;
% % D(6300000:6600000,1)=zeros(300001,1);
% % D(30400000:30600000,1)=zeros(200001,1);
% si=file.si;

D(:,2)=D(:,1);
% D(:,3)=D(:,1);




conversion=1e-6;
Ttime=si*conversion;    %convert sampling interval into seconds
Fs=1/Ttime;     %determine sampling frequency

% a=h.tags.timeSinceRecStart;

%% Filter time

LFPfilt=designfilt('lowpassfir','PassbandFrequency',250,'StopbandFrequency',300,'PassbandRipple',1,'StopbandAttenuation',30,'SampleRate',Fs);

y=filtfilt(LFPfilt,D(:,1));
D(:,1)=y;

% d=fdesign.highpass('N,F3dB,Ap',10,300,0.5,Fs);    %Highpass for spike analysis
% designmethods(d);
% Hd=design(d,'cheby1');
% y=filter(Hd,D(:,2));
% D(:,2)=y;

% Keep in mind this overwrites spike analysis
spikefilt=designfilt('highpassfir','PassbandFrequency',300,'StopbandFrequency',250,'PassbandRipple',1,'StopbandAttenuation',30,'SampleRate',Fs);

y=filtfilt(spikefilt,(D(:,2)));
D(:,2)=y;

clear y

time=[0:Ttime:si*conversion*size(D(:,1),1)];
windw=hamming(winlength);

segment=find(time>=seconds);
rounds=floor(length(time)/segment(1));

%%UP STATES ANALYSIS%%
threshold=mean(movmean(D(:,1).^2,10000))+1*std(movmean(D(:,1).^2,10000));   %USED TO BE 2*std!!!
totup=movmean(D(:,1).^2,10000)>threshold;   %check when upstates are up in total

check=0;
start=[];
done=[];
numup=zeros(size(D(:,1)));  %checks when upstates begin
for i=1:length(D)
    if D(i,1)>threshold && check == 0   %Check when signal is above threshold
        start(end+1)=i;
        check=1;
        numup(i)=1;
    elseif check==1 && D(i,1)<0         %Don't count again until the signal comes back down
        done(end+1)=i;
        check=0;
    end
end



%% Spectral Content
wave=zeros(1,rounds); %parameters for spike analysis
bigspec=[];
bigT=[];
numstd=4;
threshold=mean(D(:,2))-numstd*std(D(:,2));
wait=waitbar(0,'Doing it'); %perform welch's method
for ro=1:rounds  %make spectrograms in segments
    waitbar(ro/rounds,wait)
    ephystime=time((ro-1)*segment(1)+1:ro*segment(1));    %select time of 
    ephystart=(ro-1)*segment(1)+1;
    ephysend=ro*segment(1);
    ephysignal=D(ephystart+1:ephysend,1); %pick out signal to do spectrogram of
    
    %% The good stuff
    [S,F,T]=spectrogram(ephysignal,windw,overlap*winlength,16384,Fs);
    
    I=find(F>200);
    flimit=I(1);
    
    Sdb=20*log10(abs(S));
    
    I=find(F>30);
    J=find(F>80);
    gamma(ro)=mean(mean(Sdb(I(1):J(1),:)));
    I=find(F>1);
    J=find(F>4);
    delta(ro)=mean(mean(Sdb(I(1):J(1),:)));         
    I=find(F>4);
    J=find(F>12);
    theta(ro)=mean(mean(Sdb(I(1):J(1),:)));
    I=find(F>100);
    J=find(F>250);
    SWR(ro)=mean(mean(Sdb(I(1):J(1),:)));
    
    %%
    
    activity(ro)=mean(sqrt(D(ephystart+1:ephysend,1).^2));
    
    spectro(:,ro)=mean(Sdb(1:J(1),:)');
    J=find(F>300);
    spectrofreq=F(1:J(1));
    spectrotime(ro)=ro;
    
    bigspec(:,end+1:end+size(Sdb,2))=Sdb(1:J(1),:);
    bigT(end+1:end+length(T))=T;
    
    clear S
    
    numupmin(ro)=sum(numup(ephystart:ephysend));
    totupmin(ro)=sum(totup(ephystart:ephysend));
    
    %%Welches Method%%
%     for i=1:length(T)
%         
%         if i<=welch
%             welchSdb(:,i)=mean(Sdb(:,1:i+welch),2);
%         elseif i>=size(T,2)-welch
%             welchSdb(:,i)=mean(Sdb(:,i-welch:size(T,2)),2);
%         else
%             welchSdb(:,i)=mean(Sdb(:,i-welch:i+welch),2);
%         end
%     end
%     clear Sdb
%     
%     if ro==1
%         welchSmax=max(max(welchSdb));
%         clipvals=[welchSmax-100,welchSmax]; %make coloring good in figure
%     end
    
        %%Spike detection%%
    spikes=D(ephystart+1:ephysend,2);
    boundary=0;
    for i=9:length(spikes)-24
        if spikes(i)<threshold && boundary==0
            boundary=1;
            wave(ro)=wave(ro)+1;
            wv(1:32,wave(ro))=spikes(i-8:i+23);
%             t(wave)=i;
%             raster(i)=1;
        else
            boundary=0;
        end
    end

end
close(wait)





inject=round(injection_time);

%FFT of the 10 minutes prior to inject and the 70-80 minutes after
% 
% injectD=itime*segment(1);
% base_start=round(injectD-segment(1)*10);
% base_end=round(injectD);
% 
% stim_start=round(injectD+70*segment(1));
% stim_end=round(injectD+80*segment(1));
% 
% prebase=abs(fft(D(base_start:base_end,1)));
% w = [0:Fs/(length(prebase)):Fs-(Fs/length(prebase))];
% prestim=abs(fft(D(stim_start:stim_end,1)));
% 
% smooth_base=sgolayfilt(prebase,5,501); %Smooth the FFT to be presentable
% smooth_stim=sgolayfilt(prestim,5,501);
% 
% base=smooth_base(1:60002);
% stim=smooth_stim(1:60002);
% w=w(1:60002);

% save([pathway ephysname(1:length(ephysname)-4) 'power.mat'],'inject','gamma','SWR','delta','theta','spectrofreq','spectro','activity','numupmin','totupmin','start','done','base','stim','w')
save([pathway ephysname(1:length(ephysname)-4) 'power.mat'],'inject','gamma','SWR','delta','theta','spectrofreq','spectro','activity','numupmin','totupmin')

% gammaZ=(gamma-min(gamma))/(max(gamma)-min(gamma)); %try and make all values positive.... might be a terrible way to do this
% 
% % gammaZ=gammaZ(1:89);
% 
% gamma_baseline=mean(gammaZ(1:inject));
% 
% normalized_gamma=gammaZ/gamma_baseline;
% save([pathway ephysname(1:length(ephysname)-4) '_frequencies.mat'],'inject','gamma','delta','theta','SWR','ephysignal')

inside=0;
up=[];
uper=[];
wide=[];
for i=1:length(D(:,1))
    if inside==0 && totup(i)==1
        inside=1;
    end
    if inside==1 && totup(i)==0
        inside=0;
        uper(end+1)=max(up);
        wide(end+1)=length(up);
        up=[];
    end
    if inside==1 && totup(i)==1
        up(end+1)=D(i,1);
    end
end

mean(uper)
mean(wide)/Fs
length(uper)/(length(time)/Fs)


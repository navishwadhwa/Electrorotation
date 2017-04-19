% This is a ".m" file to determine the frequency among time domain for avi
% files. Use the "simpleConvertTDMS.m" or "ConvertTDMS.m" first to convert
% the TDMS files to mat format. To determine the steps, pottslab package is
% used here so JavaPath should be set first by run the "installPottslab.m"
% in the "Pottslab0.5" folder.This version include analysis of dissociation
% and recovery periods.

clc;
close all;
clear all;
%% Load Data by setting path and selecting data folder
warning('off', 'MATLAB:colon:nonIntegerIndex')

% Ask the user to select folder containing the .mat data files to be
% analyzed for freq. Use current working directory as a starting point.

start_path=pwd;
folder_name = uigetdir(start_path,'Select folder containing data files to be analyzed:');
cd(folder_name)
%----------------------------------------------------------------

%Collect a list of all the files in the current directory
Video_Files=dir('*.avi');
Tdms_Files=dir('*s.mat');
%Measure how long the list of files is
lengthFiles=length(Video_Files);
% for i=1: length(Video_Files)
% Video_Files(i).num=str2num(Video_Files(i).name(1:end-4));
% end
%% Run for each data file separately

   for kk = 1:lengthFiles
    %clearvars -except kk Files lengthFiles folder_name
    clc;
%% load information of the video.
%h=waitbar(0,'Loading video, please wait...');
video_name = Video_Files(kk).name;
video_path = strcat(folder_name,'\',video_name);
obj = VideoReader(video_path); %read the objects
% get the height and width for the video
vidHeight = obj.Height;
vidWidth = obj.Width;
% set a struct for store the data
s = struct('cdata',zeros(vidHeight,vidWidth,3,'uint8'),...
    'graydata',zeros(vidHeight,vidWidth,1,'uint8'),...
    'denoi_graydata',zeros(vidHeight,vidWidth,1,'uint8'),...
    'dvalue',zeros(vidHeight,vidWidth,1,'double'));
   numFrame = 1; 

while hasFrame(obj)
    %waitbar(numFram/)
    s(numFrame).cdata = readFrame(obj);
    s(numFrame).graydata = rgb2gray(s(numFrame).cdata); % convert rgb to gray image
    % use the averaging filter mode to denoise the image
    s(numFrame).denoi_graydata = filter2(fspecial('average',7),s(numFrame).graydata)/255; 
   %s(numFrame).bwdata = imbinarize(s(numFrame).denoi_graydata,thre);
   
    % convert the unit8 format to double for the gray image
    s(numFrame).dvalue = im2double(s(numFrame).denoi_graydata);
    numFrame = numFrame+1;
end
numFrame = numFrame-1; % minus one since count one more for the number of Frame

    %% Load related tdms for the video and get the information of the fields state.    
    %waitbar(1,'Loading the tdms files, please wait...');
    related_tdms= Tdms_Files(kk).name;
    tdms_path = strcat(folder_name,'\',related_tdms);
    load(tdms_path);
    related_tdms(find(related_tdms==char('_')))=[];
    fs = UntitledVoltage_2.Property.wf_samples; % sampling frequency of tdms
    F_state = UntitledVoltage_2.Data; % Shows whether the field was ON or OFF
    F_state = (F_state > 2); % convert to binary (0/1)
    Switch = diff(F_state); % find if the field was switched
    OnOff = find(Switch); % Indices of the switching data points
    %Define pairs in sample numbers describing the times when the field was on
    %and off. The first column is when the event startes and the second when it
    %ended. Frequency calculated between these two.    
    ON_pairs = [OnOff(1:2:end) OnOff(2:2:end)];
    OFF_pairs = [OnOff(2:2:end-2) OnOff(3:2:end)];
    ON_cnt = (ON_pairs(3,2)-ON_pairs(3,1)); % # samples in each ON event

    %!!!Notice!!! : some data need (OFF_pairs(1,2)-OFF_pairs(1,1))
    OFF_cnt = (OFF_pairs(3,2)-OFF_pairs(3,1)); % # samples in each OFF event 
    
    % The matrix above missed the last OFF event immediately following the last
    % ON event and before the recovery starts (two seconds after the last ON).
    % Adding that bit to the Off_pairs here.
    rec_flag = 0;
    if size(F_state,1)-ON_pairs(end,2)>OFF_cnt && F_state(end)==0
        %some data need (OFF_pairs(1,2)-OFF_pairs(1,1))
        OFF_pairs = vertcat([1 OnOff(1)],OFF_pairs,[OnOff(end) OnOff(end)+(OFF_pairs(3,2)-OFF_pairs(3,1))]);
        rec_flag = 1;
    elseif size(F_state,1)-ON_pairs(end,2)<OFF_cnt && F_state(end)==0
        OFF_pairs = vertcat([1 OnOff(1)],OFF_pairs);
    else
        OFF_pairs = vertcat([1 OnOff(1)],OFF_pairs);
    end
     % get REC pairs
     REC_array = OFF_pairs(end,2):OFF_cnt:size(F_state,1);
     REC_pairs = [REC_array(1:end-1)' REC_array(2:end)'];

% Since the sampling frequency between two files is different, get the real second for each pairs     
OFF_pairs=((OFF_pairs-1)./fs);
REC_pairs = ((REC_pairs-1)./fs);


 %% Tracking the center of the cell in the video
waitbar(3,'Loading the tdms files, please wait...');
pos=zeros(numFrame,2); %two columns, the first is x position, the second is y position.
thre=0.2; % the threshhold for center calculating 
for i=1:numFrame
    pic=s(i).dvalue;
    pos(i,:)=center(pic,thre); % center is a customed function to find the center of the cell in the image.
end

%% Apply wavelet transform 
value_min=2;
T = obj.Duration/numFrame;             % Sampling period
Fs = 1/T;            % Sampling frequency of video
L = numFrame;             % Length of signal
t = (0:L-1)*T;        % Time vector
wlet = 'morse';              % name of the mother wavelet
Scale_min=0;  % Set the minimum of frequency to 0 first, and search a proper one in the following steps
Scale_max=1e7; % Set the maximum of frequency to a large value first, and search a proper one in the following steps
    
% Got abs of coefficients for both X signals and Y signals. Plus them.
[coefsx, freqx] = cwt(pos(:,1),wlet,Fs,'TimeBandwidth',60,'VoicesPerOctave',48); 
[coefsy, freqy] = cwt(pos(:,2),wlet,Fs,'TimeBandwidth',60,'VoicesPerOctave',48);
Scale_min=max(Scale_min,find(freqx>20, 1, 'last' ));
Scale_max=min(Scale_max,length(freqx));
coepow=abs(coefsy)+abs(coefsx);
coepow(Scale_max+1:end,:)=[];
coepow(1:Scale_min-1,:)=[];
freq =freqx(Scale_min:1:Scale_max);
%% Dealing with dissociation period
OFF_period=zeros(size(OFF_pairs,1),2); % OFF pairs in the video 
coepow_off=[]; % store the coefficients of WT for dissociation 
real_time_off=[]; % correspond the real time for each period
inter_off=[]; % dashed lines inserted to denote each period

for i=1: size(OFF_pairs,1)
OFF_period(i,1) = find((t>=OFF_pairs(i,1)),1,'first'); % find the start point for OFF period i
OFF_period(i,2) = find((t<=OFF_pairs(i,2)),1,'last'); % find the end poingt for OFF period i
coepow_off=[coepow_off coepow(:,OFF_period(i,1):OFF_period(i,2))]; % get the coefficients for the OFF period i
tt=linspace(OFF_pairs(i,1),OFF_pairs(i,2),OFF_period(i,2)-OFF_period(i,1)+1); % get the real time
real_time_off=[real_time_off tt]; % collect real time for fields off
inter_off=[inter_off length(real_time_off)]; % collect lines position 
end;
 inter_off=[0 inter_off]; % insert the start time.
 inter_off(end)=[];
 %% peak finder for OFF
        value_min=3; % threshold of peak
        peak_sampling_off= 1; % sampling frequency for peak finder
        real_peak_off=[];
        % a customed function for peak finder
        peak_off=wltpeakfinder(coepow_off(:,1:peak_sampling_off:end),3,value_min);
       
        % correspond the real value of peaks for the WT scale
        tt = isnan(peak_off);
        for i=1: length(peak_off)
        if tt(i) == 0
            scale= peak_off(i);     
            real_peak_off(i)= freq(scale);     
        else  real_peak_off(i)= 0;
        end
        end
 
%% Dealing with recovery period
REC_period=zeros(size(REC_pairs,1),2); % REC pairs in the video 
coepow_rec=[]; % store the coefficients of WT for recovery 
real_time_rec=[]; % correspond the real time for each period
inter_rec=[];
for i=1: size(REC_pairs,1)
REC_period(i,1) = find((t>=REC_pairs(i,1)),1,'first'); % find the start point for REC period i
REC_period(i,2) = find((t<=REC_pairs(i,2)),1,'last'); % find the end poingt for REC period i
coepow_rec=[coepow_rec coepow(:,REC_period(i,1):REC_period(i,2))]; % get the coefficients for the REC period i
tt=linspace(REC_pairs(i,1),REC_pairs(i,2),REC_period(i,2)-REC_period(i,1)+1); % get the real time
real_time_rec=[real_time_rec tt]; % collect real time for fields recovery
inter_rec=[inter_rec length(real_time_rec)]; % collect lines position 
end;



 %% peak finder for REC
        value_min=3; % threshold of peak
        peak_sampling_rec= 20; % sampling frequency for peak finder
        real_peak_rec=[];
        % a customed function for peak finder
        peak_rec=wltpeakfinder(coepow_rec(:,1:peak_sampling_rec:end),3,value_min);
       
        % correspond the real value of peaks for the WT scale
        tt = isnan(peak_rec);
        for i=1: length(peak_rec)
        if tt(i) == 0
            scale= peak_rec(i);     
            real_peak_rec(i)= freq(scale);     
        else  real_peak_rec(i)= 0;
        end
        end
                
        
%% plot the figure 
 h = figure;
    A1 = axes('position',[0.08 0.2 0.4 0.6]);
    A2 = axes('position',[0.56 0.2 0.4 0.6]);

%% Dissociation 
% set the grid for the wavelet spectrum
spec_sampling=1; % sampling frequency for the wavelet spectrum

set(gcf,'CurrentAxes',A1)
x=(linspace(1,size(coepow_off(:,1:spec_sampling:end),2),size(coepow_off(:,1:spec_sampling:end),2)));
y=freq(1:end);
[X Y]=meshgrid(x,y);

   % plot wavelet spectrum
   surf(X, Y,coepow_off(:,1:spec_sampling:end),'LineStyle','none');
   c=colorbar;
   c.Label.String = 'Intensity';
   set(gca,'Clim',[value_min 6]);
   %plot period lines
  xlabel('Time (s)');
  ylabel('Frequency (Hz)');
  name=['Dissociation' ',' related_tdms];
  title(name);
  axis tight;
 % set x label for each pr
  label_step=size(OFF_pairs,1)/8; % the interbal of x axis label for off periods.
   % plot peak  
  inter_off=inter_off./spec_sampling;
  z=40*ones(1,length(peak_off));
  hold on;
  peaknum = linspace(1,length(peak_off),length(peak_off));
  plot3(peaknum(1:end),real_peak_off(1:end),z(1:end),'.','color','k','MarkerSize',8);
  legend('Intensity of WT','Peak','Location','best');
  set(gca,'xtick',[]);
  set(gca,'xtick',inter_off(1:label_step:end));
  set(gca,'xticklabel',{floor(OFF_pairs(1:label_step:end,1))});  
    yline=[min(freq) max(freq)];
   zline=[20 20];
   xline=[0 0];
   hold on; 
for i=1:length(inter_off)
    plot3(xline+inter_off(i),yline,zline,'--w');
end
  view([0 90]);
%% Recovery
% set the grid for the wavelet spectrum
set(gcf,'CurrentAxes',A2)
spec_sampling=20; % sampling frequency for the wavelet spectrum
x=(linspace(1,size(coepow_rec(:,1:spec_sampling:end),2),size(coepow_rec(:,1:spec_sampling:end),2)));
y=freq(1:end);
[X Y]=meshgrid(x,y);


   % plot wavelet spectrum
   surf(X, Y,coepow_rec(:,1:spec_sampling:end),'LineStyle','none');
   c=colorbar;
   c.Label.String = 'Intensity';
   set(gca,'Clim',[value_min 6]);
   
  xlabel('Time (s)');
  ylabel('Frequency (Hz)');
  name=['Recovery' ',' related_tdms];
  title(name);
  axis tight;
 % set x label for each pr

 % plot peak  
  z=40*ones(1,length(peak_rec));
  inter_rec=inter_rec./spec_sampling;
  hold on;
  peaknum = linspace(1,length(peak_rec),length(peak_rec));
  plot3(peaknum(1:end),real_peak_rec(1:end),z(1:end),'.','color','k','MarkerSize',8);
  view([0 90]);
  legend('Intensity of WT','Peak','Location','best');
  label_step=size(REC_pairs,1)/4;% the interbal of x axis label for off periods.
  set(gca,'xtick',[]);
  set(gca,'xtick',inter_rec(1:label_step:end));
  set(gca,'xticklabel',{floor(REC_pairs(1:label_step:end,1))});  

% save png file  
set(gcf,'PaperType','A4')
    set(gcf,'PaperUnits','centimeters')
    xSize = 28; ySize = 15;
    xLeft = (30-xSize)/2; yTop = (21-ySize)/2;
    set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
    set(gcf,'Position',[100 100 xSize*50 ySize*50])
    set(gcf,'PaperOrientation','landscape');
    plot_name = strcat(related_tdms(1:end-4),'video.png');
    saveas(gcf,plot_name)
   end  

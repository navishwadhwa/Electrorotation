% This is a ".m" file to determine the frequency among time domain for tdms
% files. Use the "simpleConvertTDMS.m" or "ConvertTDMS.m" first to convert
% the TDMS files to mat format. To determine the steps, pottslab package is
% used here so JavaPath should be set first by run the "installPottslab.m"
% in the "Pottslab0.5" folder.This version include analysis of dissociation
% and recovery periods.

clc;
close all;
clear all;
%% Load Data by setting path
warning('off', 'MATLAB:colon:nonIntegerIndex')

% Ask the user to select folder containing the .mat data files to be
% analyzed for freq. Use current working directory as a starting point.

start_path=pwd;
folder_name = uigetdir(start_path,'Select folder containing .mat files to be analyzed:');
cd(folder_name)
%----------------------------------------------------------------

%Collect a list of all the files in the current directory
Files=dir('*s.mat');

%Measure how long the list of files is
lengthFiles=length(Files);


%% Run for each data file separately

   for kk = 1:lengthFiles
    clc;
    %clearvars -except kk Files lengthFiles folder_name
    
    file_name = Files(kk).name;
    path_name = strcat(folder_name,'\',file_name);
    
    %Load the mat file containing the data. Eventually this should be in a loop
    %going through all files in a cerain folder(s).
    load(path_name)
    
    %% Start the main codes
    %% Load measured parameters and arrange them to pairs.
    
    X_signal = UntitledVoltage_0.Data; % Contains X-Photomultiplier signal
    Y_signal = UntitledVoltage_1.Data; % Contains Y-Photomultiplier signal
    fs = UntitledVoltage_2.Property.wf_samples; % number of samples per sec
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
    OFF_cnt = (OFF_pairs(3,2)-OFF_pairs(3,1)); % # samples in each OFF event
   

    % The matrix above missed the last OFF event immediately following the last
    % ON event and before the recovery starts (two seconds after the last ON).
    % Adding that bit to the Off_pairs here.
    rec_flag = 0;
    if size(F_state,1)-ON_pairs(end,2)>OFF_cnt && F_state(end)==0
        OFF_pairs = vertcat([1 OnOff(1)],OFF_pairs,[OnOff(end) OnOff(end)+(OFF_pairs(3,2)-OFF_pairs(3,1))]);
        rec_flag = 1;
    elseif size(F_state,1)-ON_pairs(end,2)<OFF_cnt && F_state(end)==0
        OFF_pairs = vertcat([1 OnOff(1)],OFF_pairs);
    else
        OFF_pairs = vertcat([1 OnOff(1)],OFF_pairs);
    end
    
    %% Dissociation Periods
    
    %% Recovery Periods        
        %Recovery starts when the final OFF event ends, until the end of sampling.
        %Divided that section of the data into chunks of the same size as OFF
        %events for frequency calculation. Then the same parameters as OFF can be
        %used for REC.
        if rec_flag == 1
           % Set parameters
           d_sample = 1e3; % The resolution selected of peak finding step.
           real_rec_time =[];
           all_rec_coepow=[];
           peak=[];
           real_peak=[];
           pottsL2=[]; %Step value.         
           h=waitbar(0,'Dealing with Recovery periods. Please wait...') %show a time bar
           Inter=0; % Insterted lines for different periods.
           REC_array = OFF_pairs(end,2):ON_cnt:size(X_signal,1);
           REC_pairs = [REC_array(1:end-1)' REC_array(2:end)'];
           Npairs = size(REC_pairs,1);
           Scale_min=0;
           Scale_max=1e7;
    for ii = 1:Npairs
         waitbar(ii/Npairs);   %Define the process of the bar
         %Get the "chunk" of data for this computation
         X_chunk = X_signal(REC_pairs(ii,1):REC_pairs(ii,2));    
         Y_chunk = Y_signal(REC_pairs(ii,1):REC_pairs(ii,2));     
         wlet = 'morse';              % name of the mother wavelet
         % Got abs of coefficients for both X signals and Y signals. Plus them.
        [coefsx, freqx] = cwt(X_chunk,wlet,fs,'TimeBandwidth',60,'VoicesPerOctave',48);
        [coefsy, freqy] = cwt(Y_chunk,wlet,fs,'TimeBandwidth',60,'VoicesPerOctave',48);
        % Got abs of coefficients for both X signals and Y signals. Plus them.
        coepow=abs(coefsx)+abs(coefsy);
        %Select the interested periods and pick out one data with a
        %d_sample distance.
        %Use wltpeakfinder function to find the step.
        % a modified peakfinder function, which can find the frequency with
        % maximum intensity at each instant in the wavelet transform graph
        % Parameter 1: wavelet transform result;
        % Parameter 2: threshold for the relative height of the peak: small
        % value means large threshold
        % Parameter 3: threshold for the peakfinder(): the peak has to be higher than it    
        all_rec_coepow = [all_rec_coepow coepow];%Collect all the coefficients of wavelet transform.
        %List all the real time.
        x=linspace(ceil(REC_pairs(ii,1)/fs),floor(REC_pairs(ii,2)/fs),(floor(REC_pairs(ii,2)/fs)-ceil(REC_pairs(ii,1)/fs)+1));
        real_rec_time=[real_rec_time x];
        Interi=size(coepow,2);
        Inter=[Inter Inter(end)+Interi];
        Scale_min=max(Scale_min,find(freqx>16, 1, 'last' ));
        Scale_max=min(Scale_max,length(freqx));
    end
   
   peak=wltpeakfinder(all_rec_coepow,3,2)+Scale_min-1;
   real_peak=[];
        tt = isnan(peak);
        for i=1: length(peak)
        if (tt(i) == 0) && (scale<Scale_max)
            scale= peak(i);
            real_peak(i)= freqx(scale);     
        else  real_peak(i)= nan;
        end
        end
       
        f=(real_peak)';
        nanpos=isnan(f);
        nanpos=find(nanpos==1);
        f(nanpos)=0;
        freq_int =freqx(Scale_min:1:Scale_max);
        all_rec_coepow(1:Scale_min-1,:)=[];
        all_rec_coepow(Scale_max+1:end,:)=[];
        % Select the range of frequency which we are interested.
        freq_int=[freq_int;0];
        sss=zeros(1,size(all_rec_coepow,2));

        all_rec_coepow=[all_rec_coepow;sss];  
        % A function uses potts model. Parameter 1: data
        % Parameter 2: level of peak. The bigger number will be more coarse.
        
      %% plots part
        pottsL2 = minL2Potts(f, 5);
        zeropos=find(pottsL2==0);
        pottsL2(zeropos)=nan;
        
   
      %Exclude NaN first

   all_rec_coepow_int=all_rec_coepow(Scale_min:Scale_max,:); 
   %plot wavelet transform power spectrum.
   wavelet_spectum(all_rec_coepow_int,freq_int,'Recovery',fileName);
   
   %Modify the steps. The second parameter is the thresh which considered as
   %the same steps;
        nanjud = isnan(pottsL2);
        nanpo = find(nanjud==1);
% for i=1:length(nanpo)
%     tt=1;
%     while isnan(pottsL2(nanpo(i)-tt)) 
%         tt=tt+1;
%     end
% pottsL2(nanpo(i)) = pottsL2 (nanpo(i)-tt);
% end
   step_modified = stepmodify(pottsL2,0.1,1);
 
   %Plot peaks and steps.
   ps_plot(peak,real_peak,pottsL2,step_modified,freq_int,Inter,real_rec_time,100,0,14);
   
   end;
   end;



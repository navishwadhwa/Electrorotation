% This is a ".m" file to determine the frequency among time domain for tdms
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
folder_name = uigetdir(start_path,'Select folder containing .mat files to be analyzed:');
cd(folder_name)
%----------------------------------------------------------------

%Collect a list of all the files in the current directory
Files=dir('*s.mat');

%Measure how long the list of files is
lengthFiles=length(Files);


%% Run for each data file separately

   for kk = 1:lengthFiles
    %clearvars -except kk Files lengthFiles folder_name
    clc;
    
    file_name = Files(kk).name;
    path_name = strcat(folder_name,'\',file_name);
    
    %delete "_" in the file name.
    
    
    %Load the mat file containing the data. Eventually this should be in a loop
    %going through all files in a cerain folder(s).
    load(path_name)
    fileName = file_name;
    fileName(find(fileName==char('_')))=[];
    
   %% --------------%% Start the main codes %%---------------------------------------
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
    
    %% ------------- %% Dissociation Periods %% --------------------------
     %% prepare the parameter
     d_sample_off = 1e2; % The resolution of peak finding step.
     Scale_min=0;  % Set the minimum of frequency to 0 first, and search a proper one in the following steps
     Scale_max=1e7; % Set the maximum of frequency to a large value first, and search a proper one in the following steps
     all_off_coepow=[];% Collect all coefficients of wavelet tranform for the dissociation periods.
     real_off_time=[]; % Arrange the real time.
     Npairs = size(OFF_pairs,1); % Numbers of off pairs.
     peak_off=[]; %peak found by peak finder in the coefficients of wavelet transform.
     real_peak_off=[];%corresponding real frequency value for peak 
     pottsL2_off=[]; %Step finder.         
     h=waitbar(0,'Dealing with dissocation periods. Please wait...') %show a time bar
     Insert_off=0; % Inserted lines for different periods.
     %% WT and peak\step finding for dissocation
     for ii=1: Npairs
         waitbar(ii/Npairs);   %Define the process of the bar
         
         %Get the "chunk" of data for this computation
         X_chunk = X_signal(OFF_pairs(ii,1):OFF_pairs(ii,2));  
         Y_chunk = Y_signal(OFF_pairs(ii,1):OFF_pairs(ii,2)); 
         wlet = 'morse';   % name of the mother wavelet
         % Got abs of coefficients for both X signals and Y signals. Plus them.
        [coefsx, freqx] = cwt(X_chunk,wlet,fs,'TimeBandwidth',60,'VoicesPerOctave',48); 
        [coefsy, freqy] = cwt(Y_chunk,wlet,fs,'TimeBandwidth',60,'VoicesPerOctave',48); 
        coepow=abs(coefsx)+abs(coefsy);
        
        %since the first cycle for dissociation is the original periods
        %, the scale is defferent with the following
        if (ii>1) &&(length(freqx)<size(all_off_coepow,1))
            all_off_coepow(length(freqx)+1:end,:)=[];
        else if (ii>1) && (length(freqx)>size(all_off_coepow,1))
            coepow(size(all_off_coepow,1)+1:end,:)=[];
            end
        end      
            
        all_off_coepow = [all_off_coepow coepow];
        
        %find the real time 
        x=linspace(ceil(OFF_pairs(ii,1)/fs),floor(OFF_pairs(ii,2)/fs),(floor(OFF_pairs(ii,2)/fs)-ceil(OFF_pairs(ii,1)/fs)+1));
        real_off_time=[real_off_time x];
        %calculate the inter line
        Inserti=size(coepow,2);
        Insert_off=[Insert_off Insert_off(end)+Inserti];
        %find the minimum and maximum of scale (set the maxiumum of the frequency is 16Hz,or another proper one)
        Scale_min=max(Scale_min,find(freqx>16, 1, 'last' ));
        Scale_max=min(Scale_max,length(freqx));
     end
       % delete useless value of coefficients 
       all_off_coepow(1:Scale_min-1,:)=[];
       all_off_coepow(Scale_max+1:end,:)=[];
       % select the frequency used 
       freq_off =freqx(Scale_min:1:Scale_max);
      
        %Use wltpeakfinder function to find the step.
        % a modified peakfinder function, which can find the frequency with
        % maximum intensity at each instant in the wavelet transform graph
        % Parameter 1: wavelet transform result;
        % Parameter 2: threshold for the relative height of the peak: small
        % value means large threshold
        % Parameter 3: threshold for the peakfinder(): the peak has to be higher than it    
       peak_off=wltpeakfinder(all_off_coepow(:,1:d_sample_off:end),3,2.3)+Scale_min-1;
       
       % correspond to real value of peak
       peaknan = isnan(peak_off);
        for i=1: length(peak_off)
        if peaknan(i) == 0  
            real_peak_off(i)= freqx(peak_off(i));     
        else  real_peak_off(i)= 0; %if the peak is not found, in the mostly situation the cell stops
        end
        end
        
        % Set the zero of frequency and wavelet spectrum.
        freq_off=[freq_off;0];
        zeroline=zeros(1,size(all_off_coepow,2));
        all_off_coepow=[all_off_coepow;zeroline];  
        % minL2Potts: A function uses potts model.
        % Parameter 1: data
        % Parameter 2: level of peak. The bigger number will be more coarse.
        pottsL2_off = minL2Potts(real_peak_off, 5);
        zeropos=find(pottsL2_off==0);
        pottsL2_off(zeropos)=nan;
        % for i=1:length(nanpo)
        %     tt=1;
        %     while isnan(pottsL2(nanpo(i)-tt)) 
        %         tt=tt+1;
        %     end
        % pottsL2(nanpo(i)) = pottsL2 (nanpo(i)-tt);
        % end
        
        % modify the step by "stepmodify" function. The second parameter is the thresh which considered as
        %the same steps;
        step_modified_off = stepmodify(pottsL2_off,2,10);
       %step_modified_off = pottsL2_off;
%         nanjud = isnan(pottsL2_off);
%         nanpo = find(nanjud==1);
    %% ------------- %% Recovery Periods  %% ----------------------------     
        %Recovery starts when the final OFF event ends, until the end of sampling.
        %Divided that section of the data into chunks of the same size as OFF
        %events for frequency calculation. Then the same parameters as OFF can be
        %used for REC.
        if rec_flag == 1
           % Set parameters
           d_sample_rec = 1e3; % The resolution selected of peak finding step.
           Scale_min=0;  % Set the minimum of frequency to 0 first, and search a proper one in the following steps
           Scale_max=1e7; % Set the maximum of frequency to a large value first, and search a proper one in the following steps
           real_rec_time =[]; % Arrange the real time.
           all_rec_coepow=[]; % Collect all coefficients of wavelet tranform for the recovery periods.
           peak_rec=[]; %peak found by peak finder in the coefficients of wavelet transform.
           real_peak_rec=[]; %corresponding real frequency value for peak
           pottsL2_rec=[]; %Step value.         
           h=waitbar(0,'Dealing with Recovery periods. Please wait...') %show a time bar
           Insert_rec=0; % Insterted lines for different periods.
           REC_array = OFF_pairs(end,2):OFF_cnt:size(X_signal,1);
           REC_pairs = [REC_array(1:end-1)' REC_array(2:end)'];
           Npairs = size(REC_pairs,1);
     %% WT and peak\step finding for recovery    
    for ii = 1:Npairs
         waitbar(ii/Npairs);   %Define the process of the bar
         %Get the "chunk" of data for this computation
         X_chunk = X_signal(REC_pairs(ii,1):REC_pairs(ii,2));    
         Y_chunk = Y_signal(REC_pairs(ii,1):REC_pairs(ii,2));     
         wlet = 'morse';   % name of the mother wavelet
         % Got abs of coefficients for both X signals and Y signals. Plus them.
        [coefsx, freqx] = cwt(X_chunk,wlet,fs,'TimeBandwidth',60,'VoicesPerOctave',48);
        [coefsy, freqy] = cwt(Y_chunk,wlet,fs,'TimeBandwidth',60,'VoicesPerOctave',48);
        % Got abs of coefficients for both X signals and Y signals. Plus them.
        coepow=abs(coefsx)+abs(coefsy);
        all_rec_coepow = [all_rec_coepow coepow];%Collect all the coefficients of wavelet transform.
 
        %List all the real time.
        x=linspace(ceil(REC_pairs(ii,1)/fs),floor(REC_pairs(ii,2)/fs),(floor(REC_pairs(ii,2)/fs)-ceil(REC_pairs(ii,1)/fs)+1));
        real_rec_time=[real_rec_time x];
        
        %calculate the inter line
        Inserti=size(coepow,2);
        Insert_rec=[Insert_rec Insert_rec(end)+Inserti];
        %find the minimum and maximum of scale (set the maxiumum of the frequency is 16Hz,or another proper one)
        Scale_min=max(Scale_min,find(freqx>16, 1, 'last' ));
        Scale_max=min(Scale_max,length(freqx));
    end
       % delete useless value of coefficients
       all_rec_coepow(1:Scale_min-1,:)=[];
       all_rec_coepow(Scale_max+1:end,:)=[];
       real_rec_time = real_rec_time-1;
       % select the frequency used
       freq_rec =freqx(Scale_min:1:Scale_max);

       % find the peak by wltpeakfinder function
       peak_rec=wltpeakfinder(all_rec_coepow(:,1:d_sample_rec:end),3,2)+Scale_min-1;
       
       % correspond to real value of peak
       peaknan = isnan(peak_rec);
        for i=1: length(peak_rec)
        if peaknan(i) == 0 
            real_peak_rec(i)= freqx(peak_rec(i));    
        else  real_peak_rec(i)= 0; %if the peak is not found, in the mostly situation the cell stops
        end
        end
        
%         f=(real_peak_rec)';
%         nanpos=isnan(f);
%         nanpos=find(nanpos==1);
%         f(nanpos)=0;

        % Set the zero of frequency and wavelet spectrum.
        freq_rec=[freq_rec;0];
        zeroline=zeros(1,size(all_rec_coepow,2));
        all_rec_coepow=[all_rec_coepow;zeroline];  
        
        % A function uses potts model. 
        % Parameter 1: data
        % Parameter 2: level of peak. The bigger number will be more coarse.
        pottsL2_rec = minL2Potts(real_peak_rec, 5);
        zeropos=find(pottsL2_rec==0);
        pottsL2_rec(zeropos)=nan;
        
   %Exclude NaN first
   %Modify the steps. The second parameter is the thresh which considered as
   %the same steps;
%         nanjud = isnan(pottsL2_rec);
%         nanpo = find(nanjud==1);
% for i=1:length(nanpo)
%     tt=1;
%     while isnan(pottsL2(nanpo(i)-tt)) 
%         tt=tt+1;
%     end
% pottsL2(nanpo(i)) = pottsL2 (nanpo(i)-tt);
% end
   step_modified_rec = stepmodify(pottsL2_rec,0.5);
   end;
   
   %% Plot the dissociation and recovery parts together
    h = figure;
    A1 = axes('position',[0.08 0.2 0.4 0.6]);
    A2 = axes('position',[0.56 0.2 0.4 0.6]);
    set(gcf,'CurrentAxes',A1)
    wavelet_spectum(all_off_coepow(:,1:d_sample_off:end),freq_off,'Dissociation',fileName);
    ps_plot(peak_off,real_peak_off,pottsL2_off,step_modified_off,freq_off,Insert_off,real_off_time,1,0,14,fileName,0);
    view([0 90]);
    
    set(gcf,'CurrentAxes',A2)
    wavelet_spectum(all_rec_coepow(:,1:d_sample_rec:end),freq_rec,'Recovery',fileName);
    ps_plot(peak_rec,real_peak_rec,pottsL2_rec,step_modified_rec,freq_rec,Insert_rec,real_rec_time,15,0,14,fileName,1);
    view([0 90]);
    
    set(gcf,'PaperType','A4')
    set(gcf,'PaperUnits','centimeters')
    xSize = 28; ySize = 15;
    xLeft = (30-xSize)/2; yTop = (21-ySize)/2;
    set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
    set(gcf,'Position',[100 100 xSize*50 ySize*50])
    set(gcf,'PaperOrientation','landscape');
    plot_name = strcat(file_name(1:end-5),'.png');
    saveas(gcf,plot_name)
  
   end;



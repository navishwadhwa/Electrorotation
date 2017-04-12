clc
close all
clear all

%% read frames for the whole video.
FileName = '9.avi';
obj = VideoReader(FileName); %read the objects
vidHeight = obj.Height;
vidWidth = obj.Width;
s = struct('cdata',zeros(vidHeight,vidWidth,3,'uint8'),...
    'graydata',zeros(vidHeight,vidWidth,1,'uint8'),...
    'denoi_graydata',zeros(vidHeight,vidWidth,1,'uint8'),...
    'bwdata',zeros(vidHeight,vidWidth,1,'logical'),...
    'dvalue',zeros(vidHeight,vidWidth,1,'double'));
numFrame = 1; thre = 0.2; 
while hasFrame(obj)
    s(numFrame).cdata = readFrame(obj);
    s(numFrame).graydata = rgb2gray(s(numFrame).cdata);
    s(numFrame).denoi_graydata = filter2(fspecial('average',7),s(numFrame).graydata)/255;
    s(numFrame).bwdata = imbinarize(s(numFrame).denoi_graydata,thre);
    s(numFrame).dvalue = im2double(s(numFrame).denoi_graydata);
    numFrame = numFrame+1;
end
numFrame = numFrame-1;
value_min=2;
 %% find center
pos=zeros(numFrame,2); %two columns, the first is x position, the second is y position.
thre=0.2; 
for i=1:numFrame
    pic=s(i).dvalue;
    pos(i,:)=center(pic,thre);
end

%%
T = obj.Duration/numFrame;             % Sampling period
Fs = 1/T;            % Sampling frequency
L = numFrame;             % Length of signal
t = (0:L-1)*T;        % Time vector
wlet = 'morse';              % name of the mother wavelet


% Got abs of coefficients for both X signals and Y signals. Plus them.
[coefsx, freqx] = cwt(pos(:,1),wlet,Fs,'TimeBandwidth',60,'VoicesPerOctave',48); 
[coefsy, freqy] = cwt(pos(:,2),wlet,Fs,'TimeBandwidth',60,'VoicesPerOctave',48); 
coepow=abs(coefsy)+abs(coefsx);
dsample=1;
x=(linspace(0,size(coepow,2)/Fs,size(coepow,2)/dsample));
y=freqx(1:end);
[X Y]=meshgrid(x,y);

   figure;
   surf(X, Y,coepow(:,1:dsample:end),'LineStyle','none');
   c=colorbar;
   c.Label.String = 'Intensity';
  % c.Limits=[1 6];
  
  set(gca,'Clim',[value_min 6]);
  xlabel('Time (s)');
  ylabel('Frequency (Hz)');
 % name=[my_title ',' filename];
 % title(name);
  axis tight;

  %% peak finder
        value_min=3;
        peak=wltpeakfinder(coepow(:,1:dsample:end),3,value_min);
        real_peak=[];
        tt = isnan(peak);
        for i=1: length(peak)
        if tt(i) == 0
            scale= peak(i);
            
            real_peak(i)= freqx(scale);     
        else  real_peak(i)= 0;
        end
        end
        
        f=(real_peak)';
        nanpos=isnan(f);
        nanpos=find(nanpos==1);
        f(nanpos)=0;
        z=40*ones(1,length(peak));
        hold on;
       % timescale = linspace(1,length(peak),length(peak));
        plot3(x(1:end),real_peak(1:end),z(1:end),'.','color','k','MarkerSize',8);


function peakseries=wltpeakfinder(wlt,peaksel,peakthre)
siz=size(wlt);
sel=(max(wlt)-min(wlt))/peaksel;
peakseries=zeros(1,siz(2));
for i=1:siz(2)
    peak=peakfinder(wlt(:,i),sel(i),peakthre);
    lenpeak=length(peak);
    if lenpeak==1
        peakseries(i)=peak;
    elseif lenpeak==0
        peakseries(i)=nan;
    else
        [x,j]=max(wlt(peak,i));
        peakseries(i)=peak(j);
    end
    
end
end

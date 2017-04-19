function pos=center(imag,threshold)
% pos=[x y]
% threshold=65;
region=(imag>threshold);
imag=imag.*region;
xy=size(imag);
x=xy(1);
y=xy(2);
pos=[0 0];
for i=1:x
    for j=1:y
       pos=pos+[j i]*imag(i,j);
    end
end
pos=pos/sum(sum(imag));
end

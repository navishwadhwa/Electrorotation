function pos=center(imag)
% pos=[x y]
% threshold=65;
region=imbinarize(imag);
imag(~region)=0;
imag = im2double(imag);
% imshow(imag)
xy=size(imag);
x=xy(1);
y=xy(2);
[X,Y]=meshgrid(1:x,1:y);
pos = [sum(sum(X'.*imag))/sum(sum(imag)) sum(sum(Y'.*imag))/sum(sum(imag))];
% 
% for i=1:x
%     for j=1:y
%        pos=pos+[j i]*imag(i,j);
%     end
% end
% pos=pos/sum(sum(imag));
end

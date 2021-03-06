function out=stepmodify(step,thresh,dmin)
step(find(isnan(step)==1))=0;
%% thresh max
% dvalue=diff(step);
% spike= find(abs(dvalue)>thrmax);
% for i=1:length(spike)
%     tt=1;
%     while (step(spike(i)-tt)>=thrmax)
%         tt=tt+1;
%     end
% step(spike(i))=step(spike(i)-tt);
% end
%% thresh min
% spike= find(step<thrmin);
% for i=1:length(spike)
%     tt=1;
%     while (step(spike(i)-tt)<thrmin)
%         tt=tt+1;
%     end
% step(spike(i))=step(spike(i)-tt);
% end
%% 
d=diff(step); % differentiate distance with each two adjacent steps;
%step = smooth1(step,jumpmax);
jump=find(d~=0)'; % Find positions are not the same.
dis=diff(jump);
pos=find((dis<=dmin));
for i =1 : length(pos)
step(jump(pos(i))+1:jump(pos(i)+1))=step(jump(pos(i)));
end
stepl=[1; jump(1:end-1)+1]; % each step starting point
stepr=[jump; length(step)]; % each step ending point
len=length(stepr); % The number of steps
for n=1:10
for i=1:len-n
    a=step(stepl(i):stepr(i+n));
    di=max(a)-min(a);
    if di<thresh
    step(stepl(i):stepr(i+n))=a*0+mean(a);
    end
%     if di>thrmax
%      step(stepl(i):stepr(i+n))=a*0+mean(a);
%     end
end
end

out=step;
end

function y=smooth1(x,threshold)
dx=diff(x);
jumpl=find(dx~=0)';
jump=[0 jumpl length(x)];
inter=diff(jump);
spikepos=find(inter<=threshold);
for i=1:length(spikepos)
x((jump(spikepos(i))+1):jump(spikepos(i)+1))=x(jump(spikepos(i)));
end
y=x;
end

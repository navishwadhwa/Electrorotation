function y = wavelet_spectum(all_coepow,freq_int,my_title,filename)

x=linspace(1,size(all_coepow,2),size(all_coepow,2));
y=freq_int;
[X Y]=meshgrid(x,y);

   surf(X, Y,all_coepow,'LineStyle','none');
   c=colorbar;
   c.Label.String = 'Intensity';
  % c.Limits=[3 8];
  set(gca,'Clim',[3 6]);
  xlabel('Time (s)');
  ylabel('Frequency (Hz)');
  name=[my_title ',' filename];
  title(name);
  axis tight;
 
end
      



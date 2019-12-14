function [hmout] = set_hm_color(hmin)
%akesheatmap object and spits out new one with modified colormap

color1=[0,158,254]';
color2=[255,255,255]';
color3 = [254,155,0]';

hmin.Colormap = [linspacen(color1,color2),linspacen(color2,color3)]'/255;

try
hmin.ColorLimits=repmat(max(abs((hmin.ColorLimits))),size(hmin.ColorLimits)).*sign(hmin.ColorLimits);
catch
   if all(hmin.ColorLimits >0)
       hmin.ColorLimits = [ 0, max(hmin.ColorLimits)];
       hmin.Colormap = linspacen(color2,color3)'/255;
   elseif all(hmin.ColorLimits <0)
    hmin.ColorLimits   = [ min(hmin.ColorLimits),0];
      hmin.Colormap = linspacen(color1,color2)'/255;
       
   else
   end
    
end


hmout = hmin;


end


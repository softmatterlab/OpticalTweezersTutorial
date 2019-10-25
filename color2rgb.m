function rgb=color2rgb(colorst)
 if strcmp(colorst,'wine')==1
    %wine for non-linear fits
     rgb=[73/255,4/255,10/255];
 elseif strcmp(colorst,'yellow')==1
      %yellow for linear fit
     
     rgb=[241/255,185/255,14/255];
 elseif strcmp(colorst, 'blue_gray')==1
     %blue-gray for experimental data 
     rgb=[7/255, 79/255, 129/255];
 else 
     error('select the color')
     rgb=[0,0,0];
 end

end
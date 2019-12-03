function [valr, dvalr, nsig]= round_significance(val, dval)
%function to round the values to the significative figures-digits given by
%the error
nsig=abs(floor(log10(dval)));
nsigval=abs(ceil(log10(val)));
%nsig2=abs(ceil(log10(dval)));
nsig=nsig+nsigval;
    
valr = num2str(round(val,nsig,'significant'));

dvalr=num2str(round(dval,1,'significant'));


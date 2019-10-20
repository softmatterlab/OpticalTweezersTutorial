function value=chiqu (p, t, y, R)

 value=(y - fexp(t, p)')' * R * (y - fexp(t, p)')
 if isnan(value)
     value=1;
 end
end
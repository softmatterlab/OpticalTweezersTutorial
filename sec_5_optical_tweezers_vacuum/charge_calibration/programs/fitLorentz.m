function [estimates, sse, model] = fitLorentz(xdata, ydata, start_point)
% Call fminsearch with a starting point.

% start_point = rand(1, 2);
% start_point= [5E4 1E2 500 0];
model = @Lorentz;
[estimates,sse] = fminsearch(model, start_point,optimset('MaxFunEvals',5000,'MaxIter',10000));
% expfun accepts curve parameters as inputs, and outputs sse,
% the Chi squared for A*exp(-lambda*xdata)-ydata,
% and the FittedCurve. FMINSEARCH only needs sse, but we want
% to plot the FittedCurve at the end.
    function [sse, FittedCurve] = Lorentz(params)
        a1 = params(1);
        a2 = params(2);
        a3 = params(3);
        offs = params(4);
        %lambda = params(2);
        %FittedCurve = 0.5*(A0+A1*cos(2*xdata)+A2*sin(2*xdata)+A3*cos(4*xdata)+A4*sin(4*xdata));
        %FittedCurve = A*Gamma./((xdata.^2-w0^2).^2+Gamma^2*w0^2)+offs;
        FittedCurve = abs(a1)./((xdata.^2-a3^2).^2+a2^2*xdata.^2)+abs(offs);
        % FittedCurve = A + B^2*cos(2*pi*f*xdata+phi);
        % sse = sum( FittedCurve - ydata.*log(FittedCurve) +log(factorial(ydata)) );     
        %sse = sum( FittedCurve - ydata.*log(FittedCurve)  );  % max
        %liklihood fit
        sse = sum( (FittedCurve-ydata).^2  );     % least squared
    end
end
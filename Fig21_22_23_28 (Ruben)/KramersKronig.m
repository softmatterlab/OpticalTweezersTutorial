function rechi=KramersKronig(omega,imchi)
%The program inputs are the vector of the frequency
%components and the vector of the imaginary
%part of the response function under examination.
%The two vectors must have the same length 
%and the frequency vector omega must be equispaced. 
%If not, apply MATLAB functions such as interp.
%The output is the estimate of the real part as obtained
%with K-K relations.


if size(omega,1)>size(omega,2);
omega=omega';
end; if size(imchi,1)>size(imchi,2);
imchi=imchi';
end;
%Here the program rearranges the two vectors so that,
%whichever their initial shape, they become row vectors.
g=size(omega,2);
%Size of the vectors.%
rechi=zeros(size(imchi));
%The output is initialized.
a=zeros(size(imchi));
b=zeros(size(imchi));
%Two vectors for intermediate calculations are initialized
deltaomega=omega(2)-omega(1);

for j=1:g;
    imch = imchi;
    om = omega;   
    int1 = imch.*om./(om.^2-omega(j)^2);   
    int1(j) = 0;   
    alphaj = sum( int1 );
    rechi(j)=2/pi*deltaomega*alphaj;   
end;

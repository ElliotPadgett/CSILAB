function [si,en] = loadEELS(fnm)
%loadEELS loads an EELS spectrum image in dm3 format
% inputs:
%   fnm -- the name of the dm3 file to be loaded
% outputs:
%   si -- the 3D spectrum image, ordered [x,y,en]
%   en -- a vector of energy values corresponding to the energy axis.
%       (Warning!  this is often wrong and you may need to calibrate
%       yourself)
%
%This function is part of the CSILAB Package written by the Muller Group 
%at Cornell University
%Contributors include: Elliot Padgett, Megan Holtz, Paul Cueva, Julia
%   Mundy, Huolin Xin, Peter Ercius, David Muller

data = dm3Reader(fnm);
si=data.SI;
en=data.en;

N=size(si); N(1)=N(2); N(2)=size(si,1);
tempsi=zeros(N);
for i=1:length(en)
    tempsi(:,:,i)=si(:,:,i)';
end
si=tempsi;

end
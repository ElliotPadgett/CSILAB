function [ a_siL, a_siH] = alignDualEELSZLP(enL,enH,siL,siH)
% alignDualEELSZLP aligns the energy for dual EELS by finding the ZLP
% and shifting the SIs appropriately.
%
% inputs:
%   enL -- Energy values for the lower energy range
%   enH -- Energy values for the lower energy range
%   siL -- Low loss spectrum image, assumed to contain the ZLP
%   siH -- High loss spectrum image
% outputs:
%   a_siL -- Aligned low loss spectrum image
%   a_siH -- Aligned high loss spectrum image
%
%This function is part of the CSILAB Package written by the Muller Group 
%at Cornell University
%Contributors include: Elliot Padgett, Megan Holtz, Paul Cueva, Julia
%   Mundy, Huolin Xin, Peter Ercius, David Muller

p = size(siL);

siL = reshape(siL,p(1)*p(2),p(3));
siH = reshape(siH,p(1)*p(2),p(3));

[~,maxInd] = max(siL,[], 2);
maxInd = enL(maxInd);
%offset = maxInd - round(mean(maxInd));
offset = -1*maxInd;

specL = 0*siL;
specH = 0*siH;

for i=1:p(1)*p(2)
    
    osL = [enL; siL(i,:)];
    osH = [enH; siH(i,:)];
    
    specL(i,:) = ShiftSpec(osL,enL,offset(i));
    specH(i,:) = ShiftSpec(osH,enH,offset(i));
    
end

a_siL = reshape(specL,p(1),p(2),p(3));
a_siH = reshape(specH,p(1),p(2),p(3));

function spec = ShiftSpec(oS,E,offset)
% function ShiftSpec(oS,E,offset)
% oS = oringal spectrum 2 by N matrix
% E = energy axis
% offset = energy shift (+ shift the oS right, - shift the oS left)
% by Huolin Xin

if size(oS,1) ~= 2
    error('Original spectrum must be in 2xN format: first row: energy, second row: counts');
end

if nargin == 2
    offset = 0;
elseif nargin == 1
    E = oS(:,1);
    offset = 0;
end

nE = E - offset;
oEmax = max(oS(1,:));
oEmin = min(oS(1,:));
nEmax = max(nE);
nEmin = min(nE);
%spec = interp1(oS(1,:),oS(2,:),nE,'linear','extrap');
if nEmin>oEmin && nEmax<oEmax
    spec = interp1(oS(1,:),oS(2,:),nE);
elseif (nEmin>oEmax)||(nEmax<oEmin)
    error('Energy out of range!');
else
    indexmid = find(nE>=oEmin & nE<=oEmax);
    indexbig = find(nE>oEmax);
    indexsml = find(nE<oEmin);
    spec = zeros(size(nE));
    spec(indexmid) = interp1(oS(1,:),oS(2,:),nE(indexmid));
    NoS = size(oS,2);
    Span = 5;
    spec(indexsml) = mean(oS(2,1:Span));
    spec(indexbig) = mean(oS(2,NoS-Span+1:NoS));
end
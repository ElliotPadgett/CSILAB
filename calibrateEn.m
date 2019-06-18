function [ en ] = calibrateEn( si, dispersion, en_idx, en_val )
%calibrateEn Return energy vector for EELS spectra
%   Calculates the energy vector either from function arguments or from
%   user input based on selecting a spectrum feature
% inputs:
%   si -- the 3D spectrum image, ordered [x,y,en]
%   dispersion -- energy interval per pixel
%   en_idx -- (enter if known) vector index in en of known feature
%   en_val -- (enter if known) energy value of known feature
% outputs:
%   en -- a vector of energy values corresponding to the energy axis.
%
%This function is part of the CSILAB Package written by the Muller Group 
%at Cornell University
%Contributors include: Elliot Padgett, Megan Holtz, Paul Cueva, Julia
%   Mundy, Huolin Xin, Peter Ercius, David Muller


%assume last non-singleton dimension is energy
NN=size(si); NN=NN(NN~=1); Ne=NN(end); 

if nargin<3
    %User has not input specific values.  Show them the spectrum for them
    %to pick a known point
    
    %make a 1D spectrum to plot
    if length(size(si)) == 3
        meanspec = squeeze(mean(mean(si,1),2));
    elseif length(NN) ==1
        meanspec = si;
    else 
        meanspec = mean(si,1);
    end
    
    %plot spectrum
    f = figure;
    plot(meanspec)
    title('Click on known energy position');
    [en_idx,~] = ginput(1);
    en_idx = round(en_idx);
    close(f)
    
    %ask for known energy value
    prompt = {'What is the known energy?'};
    name = 'Energy';
    numlines = 1;
    n = inputdlg(prompt,name,numlines,{''},'on');
    en_val = str2num(n{1});
end

if nargin<2
    %user did not enter a dispersion
    prompt = {'What is the dispersion?'};
    name = 'Dispersion';
    numlines = 1;
    n = inputdlg(prompt,name,numlines,{''},'on');
    dispersion = str2num(n{1});
end

%Calculate energy vector
en=calculateEn(dispersion,en_idx,en_val,Ne);

end

function [ en ] = calculateEn( dispersion, en_idx, en_val, N )

en = 1:N;
en = (en-en(en_idx))*dispersion+en_val;

end

function bg = fitBG(si,en,fitwin,type,oversample)
%fitbg Fits background for EELS spectrum image
% inputs:
%   si -- the 3D spectrum image, ordered [x,y,en]
%   en -- a vector of energy values corresponding to the energy axis.
%   fitwin -- energy values specifying energy range for background fit, in
%             format [enmin, enmax]
%   type -- type of function for background fit.  Options are 'linear',
%           'exponential', 'powerlaw', 'lcpl'.
%   oversample -- (optional) degree of oversampling in pixels.  Default is
%                 no oversampling
% outputs:
%   bg -- Calculated background from background fit, in same form as si.
%
%This function is part of the CSILAB Package written by the Muller Group 
%at Cornell University
%Contributors include: Elliot Padgett, Megan Holtz, Paul Cueva, Julia
%   Mundy, Huolin Xin, Peter Ercius, David Muller

if nargin<5
    oversample=0;
end

%set up pre-processing and parameterization for linear fit
switch type
    case 'linear'
        %no parameterization needed
        X = @(x) x;
        Y = @(y) y;
        bg_func = @(En,A,B) A + B.*En;
    case 'exponential'
        %fit to log of si
        X = @(x) x;
        Y = @(y) log(checkNegs(y));
        bg_func = @(En,A,B) exp(A).*exp(B.*En);
    case 'powerlaw'
        %fit to log of en and si
        X = @(x) log(x);
        Y = @(y) log(checkNegs(y));
        bg_func = @(En,A,B) exp(A).*(En.^B);
    case 'lcpl'
        %fit to log of en and si
        X = @(x) log(x);
        Y = @(y) log(checkNegs(y));
        %then take two representative powers for linear combination
        bg_func = @(En,A1,A2,B) A1.*(En.^B(1)) + A2.*(En.^B(2));
    otherwise
        error('Unknown fit type. Options are linear, exponential, powerlaw, or lcpl.')
end

% reshape the 3D si into 2D si
p = size(si);
if length(p)==3
    si = reshape(si, p(1)*p(2),p(3));
end

% select appropriate energy range
logen = en>=min(fitwin) & en<=max(fitwin);
si_bg = si(:,logen);
en_bg = en(logen);

% change dimension of energy axis to make 2nd dimension energy dimension
if size(en_bg,1)>size(en_bg,2)
    en_bg = en_bg';  en = en';
end

%Do linear fit on parametrized data
if oversample
    nE=length(en_bg);
    %get slope/shape parameter b from fit to smoothed data
    si_bg_smooth = oversample_smoothing(reshape(si_bg, p(1),p(2),nE),oversample); %get back to 3D space and smooth
    si_bg_smooth = reshape(si_bg_smooth, p(1)*p(2),nE); %return to 2D space for fit
    [~, b] = linearFit(X(en_bg'), Y(si_bg_smooth')); 
    %get shift/scale parameter from b and original data
    a = slopeConstrainedLinearFit(X(en_bg'),Y(si_bg'),b);
else
    [a, b] = linearFit(X(en_bg'), Y(si_bg')); % primed to get into Huolin format, backwards from our own
end

if strcmp(type,'lcpl')
    %calculate representative powers for lcpl
    perc_vals = [20,80]; 
    b = prctile(b, perc_vals );
    comp = [en_bg.^b(1); en_bg.^b(2)];
    A1=a; A2=a;
    
    for i=1:size(si_bg,1) % for each spectrum
        if oversample
            coef = linsolve(comp',si_bg_smooth(i,:)'); % do the fit for shape
            coef=coef*mean(si_bg(i,:)./(coef'*comp)); % do the fit for scale
        else
            coef = linsolve(comp',si_bg(i,:)'); % do the fit
        end
        A1(i) = coef(1); A2(i) = coef(2); 
    end
    

    %Calculate background from fit parameters
    en_grid = repmat(en, [size(si,1), 1]);
    A1 = repmat(A1', [1, size(si,2)]);  A2 = repmat(A2', [1, size(si,2)]);
    bg = bg_func(en_grid,A1,A2,b);
    
else
    %positive sloping background is unphysical.  Replace with flat
    b(b>0)=0;
    a=slopeConstrainedLinearFit(X(en_bg'),Y(si_bg'),b);
    
    %Calculate background from fit parameters
    en_grid = repmat(en, [size(si,1), 1]);
    a = repmat(a', [1, size(si,2)]);  b = repmat(b', [1, size(si,2)]);
    bg = bg_func(en_grid,a,b);
end

% if originally 3d, reshape back to 3d
if length(p)==3
    bg = reshape(bg, p(1),p(2),p(3));
end

if isinf(max(bg(:)))
    error('bg is inf!')
end

end

function [a,b] = linearFit(x,y)
% Does a linear regression least squares fit to a+bx
%x is [N by 1]
%y is [N by M]
% a and b are [1 by M]

[N,M] = size(y);

%calculate sums
sumx = sum(x).*ones(1,M);
sumy = sum(y,1);
sumx2 = sum(x.^2).*ones(1,M);
sumxy = sum((x*ones(1,M)).*y,1);

%calculate fit parameters
b = (N*sumxy - sumx.*sumy)./(N*sumx2 - sumx.^2);
a = (sumy-b.*sumx)/N;
end

function [a] = slopeConstrainedLinearFit(x,y,b)
%finds the intercept a of a best-fit line a+bx given b

[N,M] = size(y);

sumx = sum(x).*ones(1,M);
sumy = sum(y,1);

a = (sumy-b.*sumx)/N;
end

function [ si_smooth ] = oversample_smoothing( si, oversample )
% does oversampling smoothing to SI (oversampling fwhm in all directions except
% energy direction, with oversample in pixels)

si_smooth=si;
sigma = oversample / (2*sqrt(2*log(2))) ;
p = size(si);
if length(p)==3 % 3D SI
    
    si_smooth = imgaussfilt( si , sigma);
    
elseif length(p)==2 && p(2)==length(en) % energy in second dimension
    for i=1:p(2)
        si_smooth(:,i) = imgaussfilt( si(:,i) , sigma);
    end
end

end

function [ x ] = checkNegs( x )
%checkNegs replaces negatives and zeros with 1 and warns the user
%   Prevents crashes or NaNs during exponential or power law background 
%   fitting

if min(x(:))<=0
    warning('Negatives or zeros in background before exponential or power law fit');
    x(x<=0)=1;
end

end
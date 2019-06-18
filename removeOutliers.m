function imf = removeOutliers(im, box, level)
%removeOutliers removes outlier "hot pixels" from image or image stack
% inputs:
%   im -- an image or image stack to have outlier pixels removed
%   box -- radius of neighborhood to examine to determine if pixel is an
%          outlier
%   level -- significance level for identifying outlier, in number of
%            "sigmas"
% outputs:
%   imf -- filtered image with outliers removed
%
%This function is part of the CSILAB Package written by the Muller Group 
%at Cornell University
%Contributors include: Elliot Padgett, Megan Holtz, Paul Cueva, Julia
%   Mundy, Huolin Xin, Peter Ercius, David Muller

im = double(im);

% Make neighborhood box
nhood = ones(2*box+1);
nhood(box+1,box+1)=0;

mf=zeros(size(im)); st=zeros(size(im));

for i=1:size(im,3)  
mf(:,:,i) = medfilt2(im(:,:,i),size(nhood),'symmetric');
st(:,:,i) = stdfilt(im(:,:,i),nhood);
end

boundHigh = mf + level*st;
boundLow = mf - level*st;

outlier = im<boundLow | im>boundHigh;

imf = im;
imf(outlier) = mf(outlier);



end
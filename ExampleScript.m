%This is an example script showing standard usage of the functions in the
%CSILAB package for analysis of Electron Energy Loss Spectroscopy (EELS)
%data
%
%The CSILAB Package was written by the Muller Group at Cornell University
%Contributors include: Elliot Padgett, Megan Holtz, Paul Cueva, Julia
%   Mundy, Huolin Xin, Peter Ercius, David Muller

% read in data
fnm = 'ExampleData/01_SI_SI_430eV_4ms_HQ_0,25eVpx.dm3';
[si, en] = loadEELS(fnm);

% calibrate energy
dispersion = 0.25; % eV/px
%en_idx = 612; en_val = 532; %oxygen onset
%en=calibrateEn(si,dispersion,en_idx, en_val);
en = calibrateEn(si,dispersion);

% preprocess data
si = removeOutliers(si,3,8);

% fit background
type = 'lcpl';
fitwin = [780,830];
ovs = 3;
bg = fitBG(si,en,fitwin,type,ovs);

% subtract background
bgsubsi = si - bg;

% integrate signal
intwin = [837,874];
int_im = integrateSI(bgsubsi,en,intwin);

% Show results
 meanSpec = @(si3D) squeeze(mean(mean(si3D,2),1));

figure;
subplot(2,1,1), plot(en,meanSpec(si),...
                    en(en>fitwin(1)),meanSpec(bg(:,:,en>fitwin(1))),...
                    en(en>fitwin(1)),meanSpec(bgsubsi(:,:,en>fitwin(1))))
                xlabel('Energy loss (eV)'),ylabel('counts')
                xlim([min(en),max(en)])
subplot(2,1,2), imagesc(int_im); colormap(gray); axis image;

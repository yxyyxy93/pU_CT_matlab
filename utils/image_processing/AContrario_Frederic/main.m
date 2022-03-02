% AC-ARPENOS example script
% (download first the images) 
% code for the paper 
% "An a-contrario approach to quasi-periodic noise removal"
% IEEE International Conference on Image Processing (ICIP)
% Frederic Sur
% 2015


% Mandrill experiment
%%%%%%%%%%%%%%%%%%%%%
clc;
clear; 
close all;
im_gt = double(imread('mandril_gray.tif'));
[x,y] = meshgrid(1:size(im_gt,2),1:size(im_gt,1));
im = im_gt+(50*sin(2*pi*200/size(im_gt,1)*x).*sin(2*pi*2/size(im_gt,2)*y))+0*randn(size(im_gt));
im_noise_remov = ACARPENOS(im,128,0);

%%
% spectrum interpolation with TV minimization:
im_noise_remov=ACARPENOS(im,128,0,1);

%%
% It can be seen (in the "corrected power spectrum")
% that among the meaningful spikes, 
% some of them are not caused by periodic noise.
% They are due to the texture of the fur which is
% present in most patches, and invalidates the independence
% assumption. 
% The corresponding spectrum parts have however a 
% large NFA, well separated from the noise spikes 
% having a very low NFA. 
% It is possible to get rid of them with:

close all
im_noise_remov=ACARPENOS(im,128,-5);




% Boat experiment
%%%%%%%%%%%%%%%%%

clear all, close all
im_gt=double(imread('boat.png'));
[x,y]=meshgrid(1:size(im_gt,2),1:size(im_gt,1));
im=im_gt+(50*sin(2*pi*100/size(im_gt,1)*x).*sin(2*pi*40/size(im_gt,2)*y))+0*randn(size(im_gt));
im_noise_remov=ACARPENOS(im,128,0);
% see the discussion of the artefacts in 
% "Automated Removal of quasi-PEriodic NOise using
% frequency domain Statistics"
% IS&T / SPIE Journal of Electronic Imaging
% Frederic Sur, Michel Grediac
% 2015

% with TV minimization:
im_noise_remov=ACARPENOS(im,128,0,1);
% note that ringing artefacts have disappeared


% Mariner 4
%%%%%%%%%%%

clear all, close all
im=double((imread('mariner04_07b.gif')));
ACARPENOS(im,128,0);
ACARPENOS(im,128,0,1);
% note that some ringing can been deleted inside the black borders

% try also:
ACARPENOS(im,128,-1);


% Lunar Orbiter 
%%%%%%%%%%%%%%%

clear all, close all
im=double((imread('Striping_Noise_Sample.jpg')));
ACARPENOS(im,256,0);



% Mariner 6
%%%%%%%%%%%

clear all, close all
im=double(imread('mariner6_close.jpg'));
ACARPENOS(im,128,0);

ACARPENOS(im,128,0,1);
% note that ringing is removed in the lower part of the image

% also try:
ACARPENOS(im,256,0);


% Vicking
%%%%%%%%%

clear all, close all
im=double(imread('vicking1_f603a14a_crop.pgm'));
ACARPENOS(im,128,0);
ACARPENOS(im,128,-5);



% Scanline
%%%%%%%%%%

clear all, close all
im=double(imread('scanline.tif'));
im=im(:,:,1);
ACARPENOS(im,64,0);



% Halftone images
%%%%%%%%%%%%%%%%%

clear all, close all
im=double(imread('halftone.pgm'));
ACARPENOS(im,128,0);

clear all, close all
im=double(imread('halftone-eye.pgm'));
ACARPENOS(im,128,0);



% images from  Gonzalez & Woods textbook
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all, close all
im=double(imread('Fig0505(a)(applo17_boulder_noisy).tif'));
ACARPENOS(im,128,0);
ACARPENOS(im,128,0,1);

clear all, close all
im=double(imread('Fig0519(a)(florida_satellite_original).tif'));
ACARPENOS(im,256,0);
ACARPENOS(im,256,0,1);
% note that artefacts are removed from saturated, white parts

clear all, close all
im=double(imread('Fig0520(a)(NASA_Mariner6_Mars).tif'));
im=im(:,:,1);
ACARPENOS(im,128,0);

clear all, close all
im=double(imread('Fig0421(car_newsprint_sampled_at_75DPI).tif'));
ACARPENOS(im,100,0);

clear all, close all
im=double(imread('Fig0422(newspaper_shot_woman).tif'));
ACARPENOS(im,200,0);


% An image without any periodic noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

im=double(imread('lena512.tif'));
ACARPENOS(im,128,0);

% the spikes detected in the spectrum are actually false alarms
% with a NFA near the threshold 1 (log_10(NFA) = -0.719)
% Moreover, the amplitude of the noise component is below 1 gray level.

% With a secured threshold:
ACARPENOS(im,128,-1);

% no periodic noise is detected

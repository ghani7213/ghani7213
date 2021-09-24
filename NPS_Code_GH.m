%% This code is intended for 2D noise power spectra calculations which is
%an important image quality metric in medical imaging. Here I have also
%provided the test images. The code is reading those images( no object) and
%applying the NPS formula.
%%
close all; clear all; clc;
%Copy and paste the directory. My images are in the following dir
cd('C:\Users\usman\Desktop\Micro-CT vs Phase Contrast\Modified_NPS\NPS_GH');
%% Reading Multiple images to perform NPS calculations
% The loop reads the .tif files in the directory
files = dir('*.tif');
num_files = numel(files);
images = cell(1, num_files);
for k = 1:num_files
    images{k} =imread(files(k).name);
    images{k}=double(images{k});
end
%% Defining Pixel, voxel size and the corresponding f accd to Nyquist Crit
M=2; % Geometric Magnification. In most cases its 1, but here it is 2
PixelPitch=0.050; % specify the pixel pitch in mm. 
DeltaX=PixelPitch/2; %The effective pixel size in magnification
f=(1:256/2)/256/DeltaX; % Here we are defining the Nyquist Frequency 
%in 128 steps.
%% Perform Subraction
DI=images{1}-images{2}; %The idea is that imageA-imageB = noise residual
%% Extact 16 Non Overlapped ROIs
% You can take the NPS from the whole image. This cropping of the ROI makes
% it efficient and reliable to discard any residual dc value. 
for row=50:50:200;
   for col=50:50:200;
    rect=[col,row,255,255]; % The ROI size is 256x256
        DX{row,col}=imcrop(DI, rect);
        NPS{row,col}=fft2(DX{row,col}-mean(mean(DX{row,col})));
        abs_NPS{row,col}=(abs(NPS{row,col}).^2);
        % shift the negative part to the left of the zero frequency
        shifted_NPS{row,col}=fftshift(abs_NPS{row,col});
        % formula for the NPs i-e Zhou etal, Yang et
        TwoDNps{row,col}=DeltaX^2*shifted_NPS{row,col}/256^2;
        figure(2),imshow(TwoDNps{row,col},[]);
        a{row,col}=imcrop(TwoDNps{row,col},[128 121 127 15]);
        %Here we are taking 14 rows from the center for the average. This is
        %an IEC recommendation (7 rows below and 7 rows above the center
        a{row,col}=mean(a{row,col});
    end
end
a_new=(a{50,50}+a{50,100}+a{50,150}+a{50,200}+a{100,50}+a{100,100}+a{100,150}+a{100,200}+a{150,50}+a{150,100}+a{150,150}+a{150,200}+a{200,50}+a{200,100}+a{200,150}+a{200,200});
a_new=a_new/16;
a_new1=(TwoDNps{50,50}+TwoDNps{50,100}+TwoDNps{50,150}+TwoDNps{50,200}+TwoDNps{100,50}+TwoDNps{100,100}+TwoDNps{100,150}+TwoDNps{100,200}+TwoDNps{150,50}+TwoDNps{150,100}+TwoDNps{150,150}+TwoDNps{150,200}+TwoDNps{200,50}+TwoDNps{200,100}+TwoDNps{200,150}+TwoDNps{200,200});
a_new2=a_new1/16;
figure(3), plot (f,a_new); % This is the plot of the horizontal profile
%% Radial NPS averaging; as the NPS is 2D circular symmetric, we can do 
%radial averaging
b=rscan(a_new2);
[c,i]=max(b)
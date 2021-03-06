close all; clear all; clc;
%% This is a simple matlab code for flat field image correction
% Specify the Directory of Raw Images
cd('C:\Users\usman\Desktop\Journal\Bee');
files = dir('*.tif');  % Here the search is for tiff images 
num_files = numel(files);
Image = cell(1, num_files);
for k = 1:num_files
    Images{k} =imread(files(k).name);
    Images{k}=double(Images{k});
end
%% Specify the Directory for Dark Images  
cd('C:\Users\usman\Desktop\Journal\Dark');
files = dir('*.tif');
num_files = numel(files);
dark = cell(1, num_files);
for k = 1:num_files
    dark{k} =imread(files(k).name);
    dark{k}=double(dark{k});
end
% Here we are averaging the dark images
% (2316, 2316) is the array size of the image. Please change it according
% to your image
image=zeros(2316, 2316);
for i=1:1:k
    image = image + dark{i};
end
% Averaging is performed
image=image/i;


%% Specify the Directory Flat field images
cd('C:\Users\usman\Desktop\Journal\Flat');
files = dir('*.tif');
num_files = numel(files);
FF = cell(1, num_files);
for k = 1:num_files
    FF{k} =imread(files(k).name);
    FF{k}=double(FF{k});
end
% Here we are averaging the dark images
% (2316, 2316) is the array size of the image. Please change it according
% to your image
image1=zeros(2316, 2316);
for i=1:1:k
    image1 = image1 + FF{i};
end
%Averaging
image1=image1/i;

%% Flat Fielding process C=((Raw-Dark)./(Flat-Dark)).*(mean2(Flat-Dark));
%Return to the Raw image directory
cd('C:\Users\usman\Desktop\Journal\Bee');
files = dir('*.tif');
num_files = numel(files);
for k=1:num_files
    Differ=Images{k}-image;
    Differ1=image1-image;
    m=mean(mean(Differ1));
    C{k}=(Differ./Differ1).*m;
    % Here we are doing complementary flaipping to save 
    C{k}=imrotate(C{k},-90);
    C{k}=fliplr(C{k});
    figure; imshow(C{k},[]);
    %Writng the Processed Corrected Data
    a=strcat('',num2str(k),'.raw');
    fid=fopen(a,'w');
    fwrite(fid,C{k},'uint16');
    fclose(fid);
end




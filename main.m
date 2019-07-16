clc; 
clear; 
close all;

%% read in images and make them into a video
% read in all jpg images in the file
file_path = '/Users/xiaochengguo/Desktop/RedChair/';  
img_path_list = dir(strcat(file_path,'*.jpg'));      
img_num = length(img_path_list);
image = zeros(240,320,img_num);            

% make the images as vedio frames
% videoName = 'redchair.avi';
% fps = 25;
% aviobj=VideoWriter(videoName);  
% aviobj.FrameRate=fps;
% open(aviobj);

if img_num >0 
    for j=1:img_num
        image_name = img_path_list(j).name;
        im = imread(strcat(file_path,image_name));
        im = rgb2gray(im);
        image(:,:,j) = im;
        %writeVideo(aviobj,im);
        fprintf('%d %s\n',j,strcat(file_path,image_name));
    end
end
%close(aviobj);

%% method 1：use next frame minus last frame to generate differential matrix
%dif = image(:,:,1);
Th = 30;                    %Threshold

mask = zeros(240,320,img_num);
for k=2:img_num
   dif = image(:,:,k) - image(:,:,k-1);
   mask(:,:,k) = abs(dif);
   for i=1:240
       for j = 1:320
           if mask(i,j,k) <= Th
              mask(i,j,k) = 0;
           else
              mask(i,j,k) = 1;
           end
       end
   end
% writeVideo(aviobj,mask);
 end
% close(aviobj);

%% method 2: use three frames to generate mask
% mask1 = zeros(240,320,img_num);
% mask2 = zeros(240,320,img_num);
box_filter = 1/9 * ones(3,3);      %box filter
% Gaussian filter
sig = 2.2;
m = 5 * sig;
Gaussian_filter = zeros(m,m);
for x = -(m - 1)/2:1:(m - 1)/2
    for y = -(m - 1)/2:1:(m - 1)/2
        Gaussian_filter(x+(m+1)/2,y+(m+1)/2) = exp(-(x*x + y*y)/(2*sig*sig));
    end
end
% making members to integer
Gaussian_filter = round(Gaussian_filter/Gaussian_filter(1,1));
K = sum(Gaussian_filter(:));
Gaussian_filter = (1/K) * Gaussian_filter;
%use filter to image
for q = 1:img_num
    image(:,:,q) = imfilter(image(:,:,q),Gaussian_filter,'corr');
end
videoName = 'mask.avi';     %Output video
fps = 25;                    %rate of frames
aviobj=VideoWriter(videoName);  
aviobj.FrameRate=fps;
open(aviobj);
for k=2:img_num-1
   dif1 = image(:,:,k) - image(:,:,k-1);
   dif2 = image(:,:,k+1) - image(:,:,k);
%    A = dif1(dif1<35 & dif1>-35);           %find the pixels belongs to stationary area in dif1
%    B = dif2(dif2<35 & dif2>-35);           %find the pixels belongs to stationary area in dif2
%    s1 = std2(A);                           %calculate sigma in stationary area in dif1
%    s2 = std2(B);                           %calculate sigma in stationary area in dif2                          
%    s = (s1+s2)/2;                          %average two sigmas
   Th = 20;                               %calculate threshold for current frame
   mask1 = abs(dif1);
   mask2 = abs(dif2);
   for i=1:240
       for j = 1:320
           if mask1(i,j) <= Th
              mask1(i,j) = 0;
           else
              mask1(i,j) = 1;
           end
           if mask2(i,j) <= Th
              mask2(i,j) = 0;
           else
              mask2(i,j) = 1;
           end
       end
   end
 mask = mask1 & mask2;
 mask = double(mask);
 writeVideo(aviobj,mask);
end

% figure(1)
% subplot(2,3,3)
% imshow(mask(:,:,220));
% title(['Threshold = ',num2str(Th)],'FontName','Times New Roman','FontSize',30,'Color','r');
close(aviobj);

% figure(1)
% subplot(1,2,1)
% imshow(image(:,:,220),[0,255]);
% title('Moment frame','FontName','Times New Roman','FontSize',45,'Color','r');
% subplot(1,2,2)
% imshow(mask(:,:,220));
% title(['Best Threshold = ',num2str(Th)],'FontName','Times New Roman','FontSize',45,'Color','r');



%% 0.5[-1,0,1] filter
mask = zeros(240,320,img_num);
for k=2:img_num-1
   dif =0.5 *(image(:,:,k+1) - image(:,:,k-1));
   mask(:,:,k) = abs(dif);
   for i=1:240
       for j = 1:320
           if mask(i,j,k) <= Th
              mask(i,j,k) = 0;
           else
              mask(i,j,k) = 1;
           end
       end
   end
% writeVideo(aviobj,mask); 
end

%% 1D Gaussian filter
% generate gaussian filter
background = image(:,:,(225:225+10));
background = mean(background,3);
sig = 5;
m = 5 * sig;
Gaussian_filter = zeros(1,m);
for x = -(m - 1)/2:1:(m - 1)/2
        Gaussian_filter(x+(m+1)/2) = exp(-(x*x)/(2*sig*sig));
end
% making members to integer
Gaussian_filter = round(Gaussian_filter/Gaussian_filter(1));
K = sum(Gaussian_filter(:));
% oprate on images
mask = zeros(240,320,img_num);
dif = zeros(240,320);
for k=(m+1)/2:img_num-(m+3)/2
    for n = 1:m
    dif = dif + Gaussian_filter(n)*image(:,:,k+n-((m+1)/2));
    end
    dif = dif/K;
    mask(:,:,k) = abs(dif);
    for i=1:240
       for j = 1:320
           if mask(i,j,k) >= Th + 120
              mask(i,j,k) = 1;
           else
              mask(i,j,k) = 0;
           end
       end
   end
% writeVideo(aviobj,mask); 
end

%% plot part

for n = 215:225
figure(n-214)
subplot(1,2,1)
imshow(image(:,:,n),[0,255]);
title('Moment frame','FontName','Times New Roman','FontSize',45,'Color','r');
subplot(1,2,2)
imshow(mask(:,:,n));
title('Foreground','FontName','Times New Roman','FontSize',45,'Color','r');
end

    
    

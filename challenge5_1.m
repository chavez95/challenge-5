%{
    **********************  README  ************************
    
    INPUT: 
    filename: holds file name of video to be processed.
    
    OUTPUT:
    cartoon.gif: processed video is saved to this file 

    ********************************************************
%}

clc
clear all
filename = 'originalVideos/golf8.mp4';
v = VideoReader(filename); %read in video to an object

% vWrite = VideoWriter('cartoon.avi','Motion JPEG AVI');
% vWrite.FrameRate = 10;
% open(vWrite);

M = v.Height; %extract parameters
N = v.Width;
D = floor(v.Duration * v.FrameRate); %num of frames
final = zeros(M, N, 1, D); %outputs
final2 = zeros(M, N, 1, D);
vertSobel = [1, 2, 1; 0, 0, 0; -1, -2, -1]; %sobel filters
horzSobel = [1, 0, -1; 2, 0, -2; 1, 0, -1];
% im = imread('lake.jpg');
% M = size(im,1);
% N = size(im,2);
% final = zeros(M, N);
% im = rgb2gray(im);
% im = double(im);
% padIm = zeros(M+2, N+2);
% padIm(2:M+1, 2:N+1) = im;
sigmaD = 3;
sigmaR = 10;
initial = -1 : 1;
X = zeros(3); %meshgrid for distance gaussian filter
Y = zeros(3);
for i = 1:3
    X(i, :) = initial;
    Y(:, i) = initial;
end
GaussFilt = exp(-(X.^2+Y.^2)/(2*sigmaD^2)); %gaussian filter
i = 1;
% for m = 2:M+1
%     for n = 2:N+1
%         intMat = padIm(m-1:m+1,n-1:n+1);
%         intWeights = exp(-(intMat-padIm(m,n)).^2/(2*sigmaR^2));
%         biFilt = intWeights.*GaussFilt;
%         final(m-1,n-1) = sum((sum(padIm(m-1:m+1, n-1:n+1).*biFilt))');
%     end
% end
while hasFrame(v)
    frame = rgb2gray(readFrame(v));     %extract one frame
    padFrame = zeros(M+2, N+2);         %zero padding each frame to make convolution easier
    padFrame(2:M+1,2:N+1) = frame;      %inserting frame into zero padded frame
    horzFrame = zeros(M+2, N+2);        %output frames for each filter
    vertFrame = zeros(M+2, N+2);
    biFrame = zeros(M+2, N+2);          % creating space for bilateral filter
    for m = 2:M+1
        for n = 2:N+1
            horzFrame(m,n) = sum((sum(padFrame(m-1:m+1, n-1:n+1).*horzSobel))');    % applying horizontal sobel filter. detects horizontal edges
            vertFrame(m,n) = sum((sum(padFrame(m-1:m+1, n-1:n+1).*vertSobel))');    % applying vertical sobel filter. detects vertical edges
            
            intMat = padFrame(m-1:m+1,n-1:n+1);                         %intensity matrix
            intWeights = exp(-(abs(intMat-padFrame(m,n))).^2/(2*sigmaR^2));  %intensity filter
            biFilt = intWeights.*GaussFilt;                             %populating bilateral filter
            biFilt = biFilt./sum(biFilt(:));
            biFrame(m,n) = sum((sum(padFrame(m-1:m+1, n-1:n+1).*biFilt))'); % applying bilateral filter four pixels at a time
        end
    end
    final2(:, :, :, i) = (biFrame(2:M+1,2:N+1));                  % inserting each bilateral filtered frame into final output
    magFrame = (horzFrame.^2) + (vertFrame.^2);                 % magnitude sobel filter
    magFrame = sqrt(magFrame);
    final(:, :, :, i) = final2(:, :, :, i) - (magFrame(2:M+1,2:N+1) );                  % inseting each sobel filter frame into final output
    i = i + 1;
    %break;
end
% writeVideo(vWrite, uint8(final));
% close(vWrite);
%imshow(uint8(horzFrame));
imwrite(uint8(final), 'cartoon.gif', 'DelayTime', 1/10, 'LoopCount', 100);
% imwrite(uint8(final2), 'cartoonbilateralfilter.gif', 'DelayTime', 1/10);
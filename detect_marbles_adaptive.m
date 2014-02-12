
function [ Images_marbles,success ] = detect_marbles_adaptive( directory )
%Function detect_marbles_adaptive
%AV1 Practical 1
% This function detects marbles in images located in the provided directory
% and outputs their detection in an array, it returns the number of
% marbles detected and a flag indicating success or failure.
% Parts of code based detectp.m
% Input :  directory - The directory where images files are located in
% format %d.jpg
% Output : Images_marbles - An array of images detected
%          sucess - Number of detected marbles, 0 if none which corresponds
%          to failure, non-zero is sucess
%Robert Evans, Francisco Aguirre

Images_marbles=true;
success=true;

%alpha=0.8;
%beta=1.2;
alpha=0.97;
beta=1.03;
tau=0.9;
iNumBackgroundImages=funcConfig('nImagesAverageBackground');

%If we can't calculate the background, return error
if ~funcCalculateAdaptiveBackground(directory)
    success=false;
end;


%Counting how many data images we have
num_Images=length(dir(strcat(directory,'*jpg')));

% loop over all images
%num_Images=8;
for i = 1 : num_Images
    imageCurrentRGB = imread([strcat(directory,int2str(i)), '.jpg'],'jpg');
    if i==1
        [rows,cols,~] = size(imageCurrentRGB);
        bgsigma=zeros(rows,cols);
        bghistory=zeros(rows,cols,3,funcConfig('nImagesAverageBackground'));
        %Load the background statistics
        load(funcConfig('fileAdaptiveBackground'),'imageSigma','imagesBackGround');
        bgsigma=imageSigma;
        bghistory=imagesBackGround;
    end
    
    
    for r = 1:rows %200:250%1:rows/10 %14:14 %1:rows/10
        %[i,r];
        for c = 1 : cols%200:250%1 : cols/10 %25:25 %1 : cols/10
            % compute cromaticity coordinates
            R=double(imageCurrentRGB(r,c,1));
            G=double(imageCurrentRGB(r,c,2));
            B=double(imageCurrentRGB(r,c,3));
            rgbsum=(R+G+B);
            sat = rgbsum/3;  % and saturation
            if rgbsum == 0
                rgbsum=1;
            end
            red = R/rgbsum;
            green = G/rgbsum;
            sigma = bgsigma(r,c);
            if sigma == 0
                sigma=0.05;
            end
            
            % search history at pixel for samples that satisfy saturation ratio test
            prob = 0;
            count=0;
            for k = 1 : iNumBackgroundImages
                if bghistory(r,c,3,k) > 0
                    ratio = sat/bghistory(r,c,3,k);
                else
                    ratio = sat/10;
                end
                if alpha < ratio && ratio < beta
                    count = count+1;
                    prob = prob + kernel(sigma,(red-bghistory(r,c,1,k)))* ...
                        kernel(sigma,(green-bghistory(r,c,2,k)));
                end
            end
            if count > 0
                p_x_b = (prob/count);
                p_b = 0.99;     % rough estimate of % a pixel is BG.
                p_x_f= 0.01;   % rough estimate of prob that this foreground pixel value
                % is chosen. Assumes all FG values possible uniformly and
                % there are 1000 FG values.
                p_b_x = (p_x_b * p_b) / (p_x_f * (1-p_b) + p_x_b * p_b);
                
                fgprob(r,c)=1-p_b_x; % prob of foreground
                if fgprob(r,c) > (1-tau)
                    % colour different enough than foreground
                    bgthresh(r,c)=1;
                else
                    % found background pixel, so update background model
                    bgthresh(r,c)=0;
                    bghistory(r,c,1,mod(i,50)+1) = red;
                    bghistory(r,c,2,mod(i,50)+1) = green;
                    bghistory(r,c,3,mod(i,50)+1) = sat;
                end
            else
                % saturation different enough than foreground
                fgprob(r,c)=1;
                bgthresh(r,c)=1;
            end
        end
    end
    
    
    figure(1)
    image(fgprob*50)
    colormap(gray)
    figure(2)
    image(bgthresh*100)
    colormap(gray)
    
    % Output the detections binary image for each frame.
    imwrite(bgthresh,['RESULTS/BGTHRESH/',int2str(i),'.jpg'],'jpg');
    
end


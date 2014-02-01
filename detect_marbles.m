
function [ Images_marbles,success ] = detect_marbles( directory )
%Function detect_marbles
%AV1 Practical 1
% This function detects marbles in images located in the provided directory
% and outputs its detection in an array, it returns too the number of
% marbles detected
% Parts of code based detectp.m
% Input :  directory - The directory where images files are located in
% format %d.jpg
% Output : Images_marbles - An array of images detected
%          sucess - Number of detected marbles, 0 if none which corresponds
%          to failure, non-zero is sucess
%Robert Evans, Francisco Aguirre
%Temporary
Images_marbles=true;
success=true;

alpha=0.8;
beta=1.2;
tau=0.05;
SIGMAFUDGE=6;

rows=480;
cols=640;

%If we can't calculate the background, return error
if ~funcCalculateAdaptiveBackground(directory,1)
    success=false;
end;

%Load the background statistics
load(funcConfig('fileAdaptiveBackground'),'bgsigma','bghistory');

%Counting how many data images we got
num_Images=length(dir(strcat(directory,'*jpg')));

% loop over all images
for i = 1 : num_Images
    imageCurrentRGB = imread([strcat(directory,int2str(i)), '.jpg'],'jpg');
    for r = 1:rows %200:250%1:rows/10 %14:14 %1:rows/10
[i,r]
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
            count = 0;
            for k = 1 : 50
            if bghistory(r,c,3,k) > 0
                ratio = sat/bghistory(r,c,3,k);
            else
                ratio = sat/10; 
                Im = imread([strcat(directory,int2str(i)), '.jpg'],'jpg');
            end
            if alpha < ratio && ratio < beta
            count = count+1;
            prob = prob + kernel(sigma,red-bghistory(r,c,1,k))* ...
                   kernel(sigma,green-bghistory(r,c,2,k));
            end
            end
            if count > 0 
                p_x_b = (prob/count);
                p_b = 0.99;     % rough estimate of % a pixel is BG. Should compute.
                p_x_f= 0.001;   % rough estimate of prob that this foreground pixel value
                                % is chosen. Assumes all FG values possible uniformly and
                                % there are 1000 FG values. Ought to compute this given
                                % some (r,g) samples
                p_b_x = (p_x_b * p_b) / (p_x_f * (1-p_b) + p_x_b * p_b);

                fgprob(r,c)=1-p_b_x; % prob of foreground
            if fgprob(r,c) > (1-tau)
            % colour different enough that foreground
                bgthresh(r,c)=1;
            else
            % found background pixel, so update background model
                bgthresh(r,c)=0;
                bghistory(r,c,1,mod(i,50)+1) = red;
                bghistory(r,c,2,mod(i,50)+1) = green;
                bghistory(r,c,3,mod(i,50)+1) = sat;
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
            for k = 1 : 50
                if bghistory(r,c,3,k) > 0
                ratio = sat/bghistory(r,c,3,k);
                else
                    ratio = sat/10;
                end
                if alpha < ratio && ratio < beta
                    count = count+1;
                    prob = prob + kernel(sigma,red-bghistory(r,c,1,k))* ...
                           kernel(sigma,green-bghistory(r,c,2,k));
                end
            end
            if count > 0 
                p_x_b = (prob/count);
                p_b = 0.99;     % rough estimate of % a pixel is BG. Should compute.
                p_x_f= 0.001;   % rough estimate of prob that this foreground pixel value
                              % is chosen. Assumes all FG values possible uniformly and
                              % there are 1000 FG values. Ought to compute this given
                              % some (r,g) samples
                p_b_x = (p_x_b * p_b) / (p_x_f * (1-p_b) + p_x_b * p_b);

                fgprob(r,c)=1-p_b_x; % prob of foreground
                if fgprob(r,c) > (1-tau)
                % colour different enough that foreground
                bgthresh(r,c)=1;
            else
                % found background pixel, so update background model
                bgthresh(r,c)=0;
                bghistory(r,c,1,mod(i,50)+1) = red;
                bghistory(r,c,2,mod(i,50)+1) = green;
                bghistory(r,c,3,mod(i,50)+1) = sat;
            end
            else
                % saturation different enough that foreground
                fgprob(r,c)=1;
                bgthresh(r,c)=1;
            end
        end
    end
end

    figure(1)  
    image(fgprob*50)
    colormap(gray)
    figure(2)  
    image(bgthresh*100)
    colormap(gray)

    imwrite(bgthresh,['RESULTS/BGTHRESH/',int2str(i),'.jpg'],'jpg');
    else
        % saturation different enough that foreground
        fgprob(r,c)=1;
        bgthresh(r,c)=1;
    end
    end


    figure(1)  
    image(fgprob*50)
    colormap(gray)
    figure(2)  
    image(bgthresh*100)
    colormap(gray)

    imwrite(bgthresh,['RESULTS/BGTHRESH/',int2str(i),'.jpg'],'jpg');    
end    
    % End of detectp
    
  %{
    [MR,MC,Dim] = size(imgBackGround);
    
    % load image
  Im = imread([strcat(directory,int2str(i)), '.jpg'],'jpg');
  if fig1 > 0
    figure(fig1)
    clf
    imshow(Im)
  end
  Imwork = double(Im);

  %Extract marbles from image
  %marble_n is the n detected marble, and it is not yet a tracking id
  %cc contains the centres of the marbles (<marble_n>,<frame>)
  %cr contains the radii of the marbles (<marble_n>,<frame>)
  %radii is an array containing the radius of each marble as is indexed by
  %(<marble_n)
  %n_detected is the number of detected marbles
  
  [cc(:,i),cr(:,i),radii,n_detected]=extract_marbles(Imwork,imgBackGround,fig1,fig2,fig3,fig15,i);
  if flag==0
    continue
  end
%{
  if fig1 > 0
    figure(fig1)
    hold on
    for c = -0.97*radius: radius/20 : 0.97*radius
      r = sqrt(radius^2-c^2);
      plot(cc(i)+c,cr(i)+r,'g.')
      plot(cc(i)+c,cr(i)-r,'g.')
    end
    %eval(['saveas(gcf,''TRACK/trk',int2str(i-1),'.jpg'',''jpg'')']);  
  end
%}
      pause(0.3)
end

% show positions
if fig4 > 0
  figure(fig4)
  hold on
  clf
  plot(cc,'r*')
  plot(cr,'g*')
end

%}
end


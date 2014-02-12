
function [ matMarbles, n_detected ] = detect_marbles( directory, bShowImages )
%Function detect_marbles
%AV1 Practical 1
% This function detects marbles in images located in the provided directory
% and outputs their detections in a cell array, it also returns the number of
% marbles detected in each frame.
% Parts of code based detectp.m
% Input :  directory - The directory where images files are located in
% format %d.jpg
% Output : matMarbles - cell array of detected marbles, indexed by frame,
%                       containing x, y and radius
%          n_detected - vector containing the number of detected marbles per frame
%Robert Evans, Francisco Aguirre

imgBackGround = double(imread(strcat(directory,'1.','jpg')));
[MR,MC,Dim] = size(imgBackGround);

%Count how many data images we have
num_Images=length(dir('SEQ1/*jpg'));

% Configure which figures to show.
fig1=1;
fig2=0; %show image after erosion
fig15=0; %debug, show foreground minus background
fig3=3; %debug, show detected objects
fig4=0;

% Foreground detection for all images.
for i = 1 : num_Images
    % load image
  Im = imread([strcat(directory,int2str(i)), '.jpg'],'jpg');
  if fig1 > 0
    figure(fig1)
    clf
    imshow(Im)
  end
  Imwork = double(Im);

  %Extract marbles from image
  %matMarbles is an matrix, column 1 is x, column 2 is y
  %columns 3 is radius for each marbles
  %columns 4 to (256+4) are the histogram of the marble
  %each rows has a detected object (hopefully a marble)
  %n_detected is the number of detected marbles
  [matMarbles{i}, n_detected(i), extracted_background]=extract_marbles(Imwork,imgBackGround,fig1,fig2,fig3,bShowImages);

  % Add the extracted background to the existing background, so that it
  % adapts to new conditions.
  imgBackGround = averageInExtractedBackground(imgBackGround, extracted_background);

  if fig1 > 0
    figure(fig1)
    hold on
    for i2=1:n_detected(i)
        plot(matMarbles{i}(i2,1),matMarbles{i}(i2,2),'g*');
    end
    
    %eval(['saveas(gcf,''TRACK/trk',int2str(i-1),'.jpg'',''jpg'')']);  
    pause(0.3)
  end


end

% show positions
if fig4 > 0
  figure(fig4)
  hold on
  clf
%  plot(matMarbles(i),'r*')
%  plot(cr,'g*')
end

end

function averagedBackground = averageInExtractedBackground(currentBackground, newBackground)
mask_idx = newBackground ~= 0;
currentBackground(mask_idx) = (currentBackground(mask_idx) + newBackground(mask_idx)) / 2;
averagedBackground = currentBackground;
end
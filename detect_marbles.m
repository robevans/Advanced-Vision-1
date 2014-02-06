
function [ matMarbles,success ] = detect_marbles( directory )
%Function detect_marbles
%AV1 Practical 1
% This function detects marbles in images located in the provided directory
% and outputs its detection in an array, it returns too the number of
% marbles detected
% Parts of code based detectp.m
% Input :  directory - The directory where images files are located in
% format %d.jpg
% Output : matMarbles - cell array of detected marbles, indexed by frame,
%                       containing x, y and radius
%          sucess - Number of detected marbles, 0 if none which corresponds
%          to failure, non-zero is sucess
%Robert Evans, Francisco Aguirre
success=true;

imgBackGround = double(imread(strcat(directory,'1.','jpg')));
[MR,MC,Dim] = size(imgBackGround);

%Counting how many data images we got
num_Images=length(dir('SEQ1/*jpg'));

 % loop over all images
fig1=1;
fig2=0; %show image after erosion
fig15=0; %debug, show foreground minus background
fig3=3; %debug, show detected objects
fig4=0;

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
  %each rows has a detected object (hopefully a marble)
  %n_detected is the number of detected marbles
  

  [matMarbles{i},n_detected(i)]=extract_marbles(Imwork,imgBackGround,fig1,fig2,fig3);


  if fig1 > 0
    figure(fig1)
    hold on
    for i2=1:n_detected(i)
        plot(matMarbles{i}(i2,1),matMarbles{i}(i2,2),'g*');
    end
%{    
    for c = -0.97*radius: radius/20 : 0.97*radius
      r = sqrt(radius^2-c^2);
      plot(cc(i)+c,cr(i)+r,'g.')
      plot(cc(i)+c,cr(i)-r,'g.')
    end
%}
    
    %eval(['saveas(gcf,''TRACK/trk',int2str(i-1),'.jpg'',''jpg'')']);  
  end

  
      pause(0.3)
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


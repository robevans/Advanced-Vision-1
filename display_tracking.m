function [ success ] = display_tracking( directory,nTrackedMarbles,matMarblesPosition)
%Function display_tracking
%This function display trajectories of detected marbles
%AV Practical 1 20140210
%Authors: Robert Evans, Francisco Aguirre
%Input:    directory - Directory where images are located
%           matMarblesPosition - matrix of marbles position
%           (frame,marble_id,x,y)
%Output:   success - false if there was an error
success=true;
iTrackFigure=10;
sColors=['y','m','c','r','g','b','w','k'];

%Counting how many data images we got
%nFrames=length(dir(strcat(directory,'*jpg')));

imgFrame = imread(strcat(directory,'1.jpg'),'jpg');
%[iRows,iColumns,~] = size (imgFrame);

figure(iTrackFigure)
clf
imshow(imgFrame);
hold on;

for iMarble=1:nTrackedMarbles
    
    sColor=sColors(mod(iMarble,length(sColors)-1)+1);
    idxValidPositions=find(matMarblesPosition(:,iMarble,1)>0);
    xDetected=matMarblesPosition(idxValidPositions,iMarble,1);
    yDetected=matMarblesPosition(idxValidPositions,iMarble,2);
    
    plot(xDetected,yDetected,'Color',sColor,'LineStyle','-','Marker','*');
end


end


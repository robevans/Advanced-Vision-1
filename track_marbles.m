function [ strTrackedMarbles,matMarblesPosition,nTrackedMarbles,success ] = track_marbles( matMarbles, directory, vecnDetected )
%function track_marbles
%AV Practical 1
%This function receives centers and radius of detected marbles and
%returns an struct of the tracked marbles, indexed by id.
%Portions of code taken from condense.m
%Input:  matMarbles, a cell array containing, per frame
%         a 3-column matrix with marbles in rows
%        and xcenter, ycenter and radius in the columns
%        directory - Image directory
%        vecnDetected - vector with the number of detected marbles by frame
%Output:  strTrackedMarbles - A cell array containing by frame an array of
%x and y positions of tracked marbles plus the radius and the histogram,
%256 values fot the later.
%         success - non-zero if there were problems, 0 otherwise
strTrackedMarbles=struct('frame',{},'xcenter',{},'ycenter',{});
nTrackedMarbles=0;
success=true;

%Real number of max marbles plus a slack of wrongly detected objects
%which need to be cleaned later
nMaxMarbles=funcConfig('nMaxMarbles')+10;
nNumHypothesis=funcConfig('nNumHypothesis');
sColors=['y','m','c','r','g','b','w','k'];

radius=10;
% Percentage of velocity we keep after colliding with something
% might be slower
fCollisionVelocityLoss = 0.9;
pstop=0.01;      % probability of stopping vertical motion %moved from 0.05 to 0.01
pCollision=0.05;    % probability of bouncing at current state (overestimated), changed from 0.3 to 0.15
pOutBounds=0.4;  %probability that the marbles goes off the frame
%weight used when we detected something but the hypothesis said it was
%out of bounds
dWeightImprobable=0.001;


%Columns of matHypothesis, not used currently
% 1 - frame
% 2 - hypothesis number
% 3 - marble_id
% 4 - x
% 5 - y
% 6 - vx
% 7 - vy
% 8 - state
%        1 - moving
%        2 - collision
%        3 - stop
%        4 - off-limits
% 9 - weight
%matHypothesis=zeros(nMaxMarbles,nNumHypothesis,9);

%Count how many data images we got
nFrames=length(dir(strcat(directory,'*jpg')));

%Debugging variables
debugFigure=1; %Figure to draw on
bdebugAllHypothesis=0; %Show all hypothesis
bdebugSelectedHypothesis=1; %Show selected hypothesis
bdebugCleanTracking=1; %Show debug info for cleaning operations

matMarblesPosition=zeros(nFrames,nMaxMarbles,2);

%matState matrix columns
%1 - frame
%2 - marble_id
%3 - hypothesis
%4 - x
%5 - y
%6 - vx
%7 - vy
matState=zeros(nFrames,nMaxMarbles,nNumHypothesis,4);
% est. probability of state, we keep two weights, one based on x and y
% distance and another on bhattacharryya distance,
% the first column is for weighted average of both
% the second for the vector based
% the third for the bhattacharryya one
weights=zeros(nFrames,nMaxMarbles,nNumHypothesis,3);
trackstate=zeros(nFrames,nMaxMarbles,nNumHypothesis); % state in each marble, hyp and frame

nSamplesHypothesis=10;
%The minimum number of frames we need to see something to consider it an
%object
nMinimumFramesVisible=4;
%idcount, max index of the matIdent array
idcount=zeros(nMaxMarbles,1);
%matIdent matrix keeps probabilities for each hypothesis
%the more rows it has, the more probable it is to pick a given
%hypothesis
matIdent=zeros(nNumHypothesis*nSamplesHypothesis,2);

for iFrame = 1 : nFrames
    %  fprintf('track_marbles. Current Frame: %d\n',iFrame);
    % load image
    imgFrame = imread([strcat(directory,int2str(iFrame)), '.jpg'],'jpg');
    if iFrame==1
        [iRows,iColumns,~] = size (imgFrame);
        %pre-initialize state and weights
        %Fill hypothesis, if it is first frame we just guess
        for iMarble=1:nMaxMarbles
            for iHyp=1:nNumHypothesis
                for i2Frame=1:nFrames
                    matState(i2Frame,iMarble,iHyp,:) = [floor(iColumns*rand(1)),floor(iRows*rand(1)),0,0]';
                    weights(i2Frame,iMarble,iHyp,1)=1/nNumHypothesis;
                    weights(i2Frame,iMarble,iHyp,2)=1/nNumHypothesis;
                    weights(i2Frame,iMarble,iHyp,3)=1/nNumHypothesis;
                end
            end
        end
        
    end
    
    if bdebugAllHypothesis||bdebugSelectedHypothesis > 0
        figure(debugFigure)
        clf
        imshow(imgFrame);
        %      pause(0.5);
    end
    
    
    for iMarble=1:nMaxMarbles
        
        % condensation tracking
        % generate nNumHypothesis new hypotheses from current nNumHypothesis hypotheses
        % first create an auxiliary array matIdent() containing state vector j
        % SAMPLE*p_k times, where p is the estimated probability of the
        % sample j
        if iFrame ~= 1
            %matIdent if an arrary with column 1 is the hypothesis, column 2 is
            %marble_id
            idcount(iMarble)=0;
            for iHyp = 1 : nNumHypothesis    % generate sampling distribution
                num=floor(nSamplesHypothesis*nNumHypothesis*weights(iFrame-1,iMarble,iHyp,1));  % number of samples to generate
                if num > 0
                    matIdent(idcount(iMarble)+1:idcount(iMarble)+num,iMarble) = iHyp*ones(1,num);
                    idcount(iMarble)=idcount(iMarble)+num;
                end
            end
        end
    end
    
    
    %Find suitable marbles which were detected
    matDetectedMarbles = matMarbles{iFrame};
    
    for iMarble=1:max([1,nTrackedMarbles,vecnDetected(iFrame)])
        if iMarble<=vecnDetected(iFrame)
            xDetected = matDetectedMarbles(iMarble,1);
            yDetected = matDetectedMarbles(iMarble,2);
            histDetected = (matDetectedMarbles(iMarble,4:(3+256)))';
        else
            xDetected = 0;
            yDetected = 0;
            histDetected = 0;
        end
        
        for iHyp = 1 : nNumHypothesis
            %      fprintf('Current marble %d. Current hyp %d\n',iMarble,iHyp);
            
            % sample randomly from the auxiliary array matIdent()
            if iFrame==1 % make a random vector
                % x and y are random within the image rows and columns
                % velocity and visibility are zero
                vecCurrent = [floor(iColumns*rand(1)),floor(iRows*rand(1)),0,0]';
            else
                iOldSample = matIdent(ceil(idcount(iMarble)*rand(1)),iMarble); % select which old sample
                vecCurrent(1) = matState(iFrame-1,iMarble,iOldSample,1);  % get its state vector
                vecCurrent(2) = matState(iFrame-1,iMarble,iOldSample,2);
                vecCurrent(3) = matState(iFrame-1,iMarble,iOldSample,3);
                vecCurrent(4) = matState(iFrame-1,iMarble,iOldSample,4);
                
                % sample about this vector randomly
                for n = 1 : 4
                    vecCurrent(n) = vecCurrent(n) + 1*randn(1);
                end
            end
            
            % hypothesize next state
            vecPredicted = vecCurrent;   % no process at start
            if iFrame == 1    % initial time - assume marble is moving
                trackstate(iFrame,iMarble,iHyp)=1;
            else
                if trackstate(iFrame-1,iMarble,iHyp)==3 || trackstate(iFrame-1,iMarble,iHyp)==4 % if already stopped or off-limits
                    vecPredicted(3) = 0;
                    vecPredicted(4) = 0;
                    if trackstate(iFrame-1,iMarble,iHyp)==4
                        vecPredicted(1)=0;
                        vecPredicted(2)=0;
                    end
                    % stay in the same state it was before
                    trackstate(iFrame,iMarble,iHyp)=trackstate(iFrame-1,iMarble,iHyp);
                else
                    r=rand(1);   % random number for state selection
                    if r < pstop
                        %if we are to stop then velocity becomes zero
                        vecPredicted(4) = 0;
                        vecPredicted(3) = 0;
                        trackstate(iFrame,iMarble,iHyp)=3;
                    elseif r < (pCollision + pstop)
                        %We are in the collision state
                        %One of many types of collision can happen
                        %We could have hit an stationary object
                        %In that case, we bounce from our current direction
                        %We could have been hitten from either side and the
                        %velocities of the incoming object sum vectorially to ours
                        rCollisionAngle=2*pi*rand(1);
                        
                        vecPredicted(3) = vecCurrent(3) * fCollisionVelocityLoss * cos(rCollisionAngle);
                        vecPredicted(4) = vecCurrent(4) * fCollisionVelocityLoss * sin(rCollisionAngle);
                        
                        vecPredicted(1) = vecCurrent(1) + 3*abs(vecCurrent(3))*(rand(1)-0.5);
                        vecPredicted(2) = vecCurrent(2) + 3*abs(vecCurrent(4))*(rand(1)-0.5);
                        
                        trackstate(iFrame,iMarble,iHyp)=2;  % set into collision state
                    elseif r < (pCollision + pstop + pOutBounds)
                        %We're out of bounds
                        vecPredicted(1) = 0;
                        vecPredicted(2) = 0;
                        vecPredicted(3) = 0;
                        vecPredicted(4) = 0;
                        trackstate(iFrame,iMarble,iHyp) = 4;
                    else % normal motion
                        vecPredicted(1) = vecCurrent(1) + vecCurrent(3);
                        vecPredicted(2) = vecCurrent(2) + vecCurrent(4);
                        trackstate(iFrame,iMarble,iHyp)=1;
                    end
                end
            end
            
            matState(iFrame,iMarble,iHyp,:) = vecPredicted;
            
            % weight hypothesis by distance from observed data
            %only if we detected something and if the hypothesis says it is
            %not out of bounds.
            if xDetected~=0 && yDetected~=0 && trackstate(iFrame,iMarble,iHyp)==4
                weights(iFrame,iMarble,iHyp,2) = dWeightImprobable;
                weights(iFrame,iMarble,iHyp,3) = dWeightImprobable;
            elseif xDetected~=0 && yDetected~=0
                %Get the histogram for the hypothesis
                histHypothesis=histogramOfCircleAroundPoint(matState(iFrame,iMarble,iHyp,1),...
                    matState(iFrame,iMarble,iHyp,2),radius,imgFrame);
                dDistanceHistograms=bhattacharyya_distance(histDetected,histHypothesis);
                weights(iFrame,iMarble,iHyp,3) = (1-dDistanceHistograms)+dWeightImprobable;
                
                %This distance is proportional to the distance between the detection and
                %the hypothesis location, the greater the distance, the weight goes much
                %lower
                dDistanceCentroids = [xDetected,yDetected] - [matState(iFrame,iMarble,iHyp,1),matState(iFrame,iMarble,iHyp,2)];
                weights(iFrame,iMarble,iHyp,2) = 1/(dDistanceCentroids*dDistanceCentroids');
                
                % print hypohtesis details
                if bdebugAllHypothesis > 0
                    fprintf('Frame %d, Marble %d Hyp %d x %f y %f state %d\n avgHistDetected %f avgHistogramHyp %f BhaWeight %.8f VecWeight %.8f\n',iFrame,iMarble,iHyp,...
                        matState(iFrame,iMarble,iHyp,1),matState(iFrame,iMarble,iHyp,2),...
                        trackstate(iFrame,iMarble,iHyp),...
                        mean(histHypothesis),mean(histDetected),...
                        weights(iFrame,iMarble,iHyp,3),weights(iFrame,iMarble,iHyp,2));
                    
                    %                    %{                %    figure(debugFigure)
                    %                     hold on
                    %
                    %                     sColor=sColors(trackstate(iFrame,iMarble,iHyp));
                    %                     rectangle('Position',[matState(iFrame,iMarble,iHyp,1),matState(iFrame,iMarble,iHyp,2),radius,radius],'Curvature',[1,1],...
                    %                             'EdgeColor',sColor);
                    %                         }%
                end
                %Correct for observed data
                matState(iFrame,iMarble,iHyp,1)=xDetected;
                matState(iFrame,iMarble,iHyp,2)=yDetected;
            elseif xDetected==0 && yDetected==0
                %We are tracking something that is not on the screen
                %so it if out of bounds
                if trackstate(iFrame,iMarble,iHyp)==4
                    weights(iFrame,iMarble,iHyp,2)=1;
                    weights(iFrame,iMarble,iHyp,3)=dWeightImprobable;
                elseif iFrame~=1
                    weights(iFrame,iMarble,iHyp,2)=dWeightImprobable;
                    weights(iFrame,iMarble,iHyp,3)=dWeightImprobable;
                end
            elseif iFrame~=1
                weights(iFrame,iMarble,iHyp,2) = dWeightImprobable;
                weights(iFrame,iMarble,iHyp,3) = dWeightImprobable;
            end
            %Check that the weights were set
            if weights(iFrame,iMarble,iHyp,2)==0 || weights(iFrame,iMarble,iHyp,3)==0
                error('track_marbles error. weights not set. Frame %d iMarble %d iHyp %d\n',iFrame,iMarble,iHyp);
            end
            if bdebugAllHypothesis
                fprintf('Frame %d iMarble %d iHyp %d weight(1) %.7f weight(2 %.7f weight(3) %.7f\n',...
                    iFrame,iMarble,iHyp,weights(iFrame,iMarble,iHyp,1),weights(iFrame,iMarble,iHyp,2),...
                    weights(iFrame,iMarble,iHyp,3));
            end
        end
    end
    
    
    for iMarble=1:max(nTrackedMarbles,vecnDetected(iFrame))
        % rescale new hypothesis weights
        
        totalw=sum(weights(iFrame,iMarble,:,2));
        %         %If all hypothesis failed, the marble probably is out of bounds
        %         if totalw==0
        %             totalw=0.0000001;
        %         end
        if isnan(totalw)
            error('track_marbles error. total vector distance weight is nan.\n');
        end
        weights(iFrame,iMarble,:,2)=weights(iFrame,iMarble,:,2)/totalw;
        totalw=sum(weights(iFrame,iMarble,:,3));
        if isnan(totalw)
            error('track_marbles error. Total bhattacharrya weight is nand.\n');
        end
        weights(iFrame,iMarble,:,3)=weights(iFrame,iMarble,:,3)/totalw;
        weights(iFrame,iMarble,:,1)=(0.6*weights(iFrame,iMarble,:,2)+0.4*weights(iFrame,iMarble,:,3))/2;
        
        % select top hypothesis to draw
        subset=weights(iFrame,iMarble,:,1);
        top = find(subset == max(subset));
        if length(top)>1
            top=top(1);
        end
        
        %Add to our tracking
        %If the winning hypothesis says that the marble is out of bounds
        %then we don't add coordinates, we don't really track out of bounds
        %marbles.
        if trackstate(iFrame,iMarble,top)==4 %|| weights(iFrame,iMarble,top,1)<(1/nNumHypothesis)
            %            if weights(iFrame,iMarble,top,1)<(1/nNumHypothesis)
            %                fprintf('Weight is about chance. Not adding to tracking. Weight %.7f\n',weights(iFrame,iMarble,top,1));
            %            end
            continue;
        end
        
        matMarblesPosition(iFrame,iMarble,1)=matState(iFrame,iMarble,top,1);
        matMarblesPosition(iFrame,iMarble,2)=matState(iFrame,iMarble,top,2);
        if nTrackedMarbles<vecnDetected(iFrame)
            nTrackedMarbles=nTrackedMarbles+1;
        end
        
        %        trackstate(iFrame,iMarble,top);
        % display final top hypothesis
        if bdebugSelectedHypothesis > 0
            %            fprintf('Added to tracking. Frame: %d. iMarble %d. x:%d y:%d state:%d\n',iFrame,iMarble,...
            %                matMarblesPosition(iFrame,iMarble,1),matMarblesPosition(iFrame,iMarble,2),...
            %                trackstate(iFrame,iMarble,top));
            figure(debugFigure)
            hold on
            sColor=sColors(trackstate(iFrame,iMarble,top));
            rectangle('Position',[matMarblesPosition(iFrame,iMarble,1),...
                matMarblesPosition(iFrame,iMarble,2),radius,radius],...
                'Curvature',[1,1],'Edgecolor',sColor);
        end
        
    end
    
    %Generate cell array for ground truth comparison
    %if length(strTrackedMarbles)>=iFrame
    %    xcenter=[strTrackedMarbles(iFrame).xcenter matMarblesPosition(iFrame,:,2)];
    %     ycenter=[strTrackedMarbles(iFrame).ycenter matMarblesPosition(iFrame,:,3)];
    %  else
    if matMarblesPosition(iFrame)>0
        xcenter=[matMarblesPosition(iFrame,:,1)];
        ycenter=[matMarblesPosition(iFrame,:,2)];
    else
        xcenter=[0];
        ycenter=[0];
    end
    
    strTrackedMarbles(iFrame)=struct('frame',iFrame,...
        'xcenter',xcenter, ...
        'ycenter',ycenter);
    %strTrackedMarbles{iFrame}.xcenter=matMarblesPosition(iFrame,2);
    %strTrackedMarbles{iFrame}.ycenter=matMarblesPosition(iFrame,3);
    %    pause(0.5) % wait a bit for the display to catch up
    
    figure(3)
    clf
    imshow(imgFrame);
    hold on;
    
    for i2Marble=1:nTrackedMarbles
        
        sColor=sColors(mod(i2Marble,length(sColors)-1)+1);
        idxValidPositions=find(matMarblesPosition(:,i2Marble,1)>0);
        x2Detected=matMarblesPosition(idxValidPositions,i2Marble,1);
        y2Detected=matMarblesPosition(idxValidPositions,i2Marble,2);
        
        plot(x2Detected,y2Detected,'Color',sColor,'LineStyle','-','Marker','*');
        
    end
    pause(1);
end

%Filter out objects detected for less than a number of frames
for iMarble=1:nMaxMarbles
    if length(matMarblesPosition(:,iMarble,1))<nMinimumFramesVisible
        if bdebugCleanTracking
            fprintf('Tracking cleaning. Removing marble id %d as it doesn''t have the minimum frames.\n',iMarble);
        end
        %delete from detected objects
        matMarblesPosition(:,iMarble,1)=0;
        matMarblesPosition(:,iMarble,2)=0;
        for iFrame=1:nFrames
            strTrackedMarbles(iFrame).xcenter(iMarble)=0;
            strTrackedMarbles(iFrame).ycenter(iMarble)=0;
        end
        
    end
end

end


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
%x and y positions of tracked marbles
%         success - non-zero if there were problems, 0 otherwise
strTrackedMarbles=struct('frame',{},'xcenter',{},'ycenter',{});
nTrackedMarbles=0;
success=true;

%Real number of max marbles plus a slack of wrongly detected objects
%which need to be cleaned later
nMaxMarbles=funcConfig('nMaxMarbles')+10;
nNumHypothesis=funcConfig('nNumHypothesis');
sColors=['c','y','b','k'];

radius=10;
% Percentage of velocity we keep after colliding with something
% might be slower
fCollisionVelocityLoss = 0.9;
pstop=0.01;      % probability of stopping vertical motion %moved from 0.05 to 0.01
pCollision=0.15;    % probability of bouncing at current state (overestimated), changed from 0.3 to 0.15
pOutBounds=0.4;  %probability that the marbles goes off the frame


%Columns of matHypothesis
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
matHypothesis=zeros(nMaxMarbles,nNumHypothesis,9);

%Counting how many data images we got
nFrames=length(dir(strcat(directory,'*jpg')));

%Debugging variables
debugFigure=1; %Figure to draw on
bdebugAllHypothesis=0; %Show all hypothesis
bdebugSelectedHypothesis=1; %Show selected hypothesis

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
weights=zeros(nFrames,nNumHypothesis,nMaxMarbles);    % est. probability of state
trackstate=zeros(nFrames,nMaxMarbles,nNumHypothesis); % state in each marble, hyp and frame

nSamplesHypothesis=10;
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
                weights(i2Frame,iHyp,iMarble)=1/nNumHypothesis;
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
                num=floor(nSamplesHypothesis*nNumHypothesis*weights(iFrame-1,iHyp,iMarble));  % number of samples to generate
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
        if nTrackedMarbles<=vecnDetected(iFrame)
            xDetected = matDetectedMarbles(iMarble,1);
            yDetected = matDetectedMarbles(iMarble,2);
        else
            xDetected = 0;
            yDetected = 0;
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
            vecCurrent(n) = vecCurrent(n) + 5*randn(1);
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
                vecPredicted(3) = 0;
                vecPredicted(4) = 0;
                trackstate(iFrame,iMarble,iHyp) = 4;
            else % normal motion
                %          A = A3;
                %          Bu = Bu3;
                vecPredicted(1) = vecCurrent(1) + vecCurrent(3);
                vecPredicted(2) = vecCurrent(2) + vecCurrent(4);
                trackstate(iFrame,iMarble,iHyp)=1;
            end
        end
        %      xp=A*xc + Bu;      % predict next state vector
        %xp = A*xc;
    end
            
            matState(iFrame,iMarble,iHyp,:) = vecPredicted;

            % weight hypothesis by distance from observed data
            %only if we detected something
            if xDetected~=0 && yDetected~=0
                dvec = [xDetected,yDetected] - [matState(iFrame,iMarble,iHyp,1),matState(iFrame,iMarble,iHyp,2)];
                weights(iFrame,iHyp,iMarble) = 1/(dvec*dvec');
            
                % draw some samples over one image
                if bdebugAllHypothesis > 0
                    figure(debugFigure)
                    hold on

                    sColor=sColors(trackstate(iFrame,iMarble,iHyp));                
                    rectangle('Position',[matState(iFrame,iMarble,iHyp,1),matState(iFrame,iMarble,iHyp,2),radius,radius],'Curvature',[1,1],...
                            'EdgeColor',sColor);
                end   
            elseif iFrame~=1
                weights(iFrame,iHyp,iMarble) = weights(iFrame-1,iHyp,iMarble);
            end
            
  end
  end
  
        
  for iMarble=1:max(nTrackedMarbles,vecnDetected(iFrame))
        % rescale new hypothesis weights
        
        totalw=sum(weights(iFrame,:,iMarble)');
        weights(iFrame,:,iMarble)=weights(iFrame,:,iMarble)/totalw;
        
        % select top hypothesis to draw
        subset=weights(iFrame,:,iMarble);
        top = find(subset == max(subset));
        if length(top)>1
            top=top(1);
        end
        
        %Add to our tracking
        matMarblesPosition(iFrame,iMarble,1)=matState(iFrame,iMarble,top,1);
        matMarblesPosition(iFrame,iMarble,2)=matState(iFrame,iMarble,top,2);
        if nTrackedMarbles<vecnDetected(iFrame)
            nTrackedMarbles=nTrackedMarbles+1;
        end
        
%        trackstate(iFrame,iMarble,top);
        % display final top hypothesis
        if bdebugSelectedHypothesis > 0
            figure(debugFigure)
            hold on
            sColor=sColors(trackstate(iFrame,iMarble,iHyp));  
            rectangle('Position',[matState(iFrame,iMarble,top,1),matState(iFrame,iMarble,top,2),radius,radius],'Curvature',[1,1],...
                        'Edgecolor',sColor);
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
  
%  end
  
  strTrackedMarbles(iFrame)=struct('frame',iFrame,...
            'xcenter',xcenter, ...
            'ycenter',ycenter);
  %strTrackedMarbles{iFrame}.xcenter=matMarblesPosition(iFrame,2);
  %strTrackedMarbles{iFrame}.ycenter=matMarblesPosition(iFrame,3);  
%    pause(0.5) % wait a bit for the display to catch up   
                
end

end


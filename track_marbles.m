function [ strTrackedMarbles,success ] = track_marbles( matMarbles, directory, vecnDetected )
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
%Output:  strTrackedMarbles - A cell array containing by frame
%         success - non-zero if there were problems, 0 otherwise
strTrackedMarbles=true;
success=true;

% Kalman filter initialization
R=[[0.2845,0.0045]',[0.0045,0.0455]'];
H=[[1,0]',[0,1]',[0,0]',[0,0]'];
Q=0.01*eye(4);
P = 100*eye(4);
dt=1;
A=[[1,0,0,0]',[0,1,0,0]',[dt,0,1,0]',[0,dt,0,1]'];
%g = 6; % pixels^2/time step
%Bu = [0,0,0,g]';
%TODO eliminate unused variables
kfinit=0;

% multiple condensation states
nCondensationSamples=100;          % number of condensation samples
nFrames=71;        % number of time frames
radius=10;         % radius of a marble
fCollisionVelocityLoss = 0.9; % Percentage of velocity we keep after colliding with something
% might be slower

vecState=zeros(nCondensationSamples,nFrames,4); % State vectors
%[ x, y, velocity_x, velocity_y]
weights=zeros(nCondensationSamples,nFrames);    % est. probability of state
trackstate=zeros(nCondensationSamples,nFrames);         % state=1,2,3,4;
% 1 - moving
% 2 - collision
% 3 - stop
% 4 - off-limits
P=zeros(nCondensationSamples,nFrames,4,4);    % est. covariance of state vec.
for i = 1 : nCondensationSamples              % initialize estimated covariance
    for j = 1 : nFrames
        P(i,j,1,1) = 100;
        P(i,j,2,2) = 100;
        P(i,j,3,3) = 100;
        P(i,j,4,4) = 100;
    end
end
pstop=0.01;      % probability of stopping vertical motion %moved from 0.05 to 0.01
pCollision=0.15;    % probability of bouncing at current state (overestimated), changed from 0.3 to 0.15
pOutBounds=0.4;  %probability that the marbles goes off the frame
xc=zeros(4,1);   % selected state
TP=zeros(4,4);   % predicted covariance


% loop over all images
fig1=1; % Show green circle around detected marble
fig2=0;
fig15=0;
fig3=0;

imgBackGround = double(imread(strcat(directory,'1.','jpg')));
[MR,MC,~] = size(imgBackGround);

%Cycle through images
for i = 1 : nFrames
    
  if fig1 > 0
      Im = imread(strcat(directory,int2str(i),'.jpg'),'jpg');
      figure(fig1)
      clf
      imshow(Im)
  end
  
    %{
  % load image
  if i < 11
    Im = (imread(['DATA/ball0000010',int2str(i-1), '.jpg'],'jpg'));
  else
    Im = (imread(['DATA/ball000001',int2str(i-1), '.jpg'],'jpg'));
  end
  if fig1 > 0
    figure(fig1)
    clf
    imshow(Im)
  end
  Imwork = double(Im);
    %}
    
    
    % extract ball
    %  [cc(i),cr(i),radius,flag]=extractball(Imwork,Imback,fig1,fig2,fig3,fig15,i);
    
    if vecnDetected(i)==0
        
        for k = 1 : nCondensationSamples
            vecState(k,i,:) = [floor(MC*rand(1)),floor(MR*rand(1)),0,0]';
            weights(k,i)=1/nCondensationSamples;
        end
        continue
    end
    
    %Find suitable marbles which were detected
    
    matDetectedMarbles = matMarbles{i};

    
    % display green estimated ball circle over original image
    %{
    if fig1 > 0
        figure(fig1)
        hold on
        for c = -0.99*radius: radius/10 : 0.99*radius           
            r = sqrt(radius^2-c^2);
            
            plot(cc(i)+c,cr(i)+r,'g.')
            plot(cc(i)+c,cr(i)-r,'g.')
        end
    end
    %}
    
    % condensation tracking
    % generate NCON new hypotheses from current NCON hypotheses
    % first create an auxiliary array ident() containing state vector j
    % SAMPLE*p_k times, where p is the estimated probability of j
    if i ~= 1
        SAMPLE=10;
        ident=zeros(100*SAMPLE,1);
        idcount=0;
        for j = 1 : nCondensationSamples    % generate sampling distribution
            num=floor(SAMPLE*100*weights(j,i-1));  % number of samples to generate
            if num > 0
                ident(idcount+1:idcount+num) = j*ones(1,num);
                idcount=idcount+num;
            end
        end
    end
    
    % generate NCON new state vectors
    for j = 1 : nCondensationSamples
        
        % sample randomly from the auxiliary array ident()
        if i==1 % make a random vector
            % x and y are random within the image rows and columns
            % velocity and visibility are zero
            xc = [floor(MC*rand(1)),floor(MR*rand(1)),0,0]';
        else
            k = ident(ceil(idcount*rand(1))); % select which old sample
            xc(1) = vecState(k,i-1,1);  % get its state vector
            xc(2) = vecState(k,i-1,2);
            xc(3) = vecState(k,i-1,3);
            xc(4) = vecState(k,i-1,4);
            
            % sample about this vector from the distribution (assume no covariance)
            for n = 1 : 4
                xc(n) = xc(n) + 5*sqrt(P(j,i-1,n,n))*randn(1);
            end
        end
        
        % hypothesize if it is going into a bounce or tabletop state
        if i == 1    % initial time - assume marble is moving
            xp = xc;   % no process at start
            %      A = A3;
            %      Bu = Bu3;
            trackstate(j,i)=1;
        else
            if trackstate(k,i-1)==3 || trackstate(k,i-1)==4 % if already stopped or off-limits
                %        A = A1;
                %        Bu = Bu1;
                xc(3) = 0;
                xc(4) = 0;
                % stay in the same state it was before
                trackstate(j,i)=trackstate(k,i-1);
            else
                r=rand(1);   % random number for state selection
                if r < pstop
                    %          A = A1;
                    %          Bu = Bu1;
                    %if we are to stop then velocity becomes zero
                    xc(4) = 0;
                    xc(3) = 0;
                    trackstate(j,i)=3;
                elseif r < (pCollision + pstop)
                    %          A = A2;
                    %          Bu = Bu2;
                    %We are in the collision state
                    %One of many types of collision can happen
                    %We could have hit an stationary object
                    %In that case, we bounce from our current direction
                    %We could have been hitten from either side and the
                    %velocities of the incoming object sum vectorially to ours
                    rCollisionAngle=2*pi*rand(1);
                    
                    xc(3) = xc(3) * fCollisionVelocityLoss * cos(rCollisionAngle);
                    xc(4) = xc(4) * fCollisionVelocityLoss * sin(rCollisionAngle);
                    
                    xc(1) = xc(1) + 3*abs(xc(3))*(rand(1)-0.5);
                    xc(2) = xc(2) + 3*abs(xc(4))*(rand(1)-0.5);
                    
                    trackstate(j,i)=2;  % set into bounce state
                elseif r < (pCollision + pstop + pOutBounds)
                    %We're out of bounds
                    xc(3) = 0;
                    xc(4) = 0;
                    trackstate(j,i) = 4;
                else % normal motion
                    %          A = A3;
                    %          Bu = Bu3;
                    %xc(1) = xc(1) + xc(3);
                    %xc(2) = xc(2) + xc(4);
                    trackstate(j,i)=1;
                end
            end
            %      xp=A*xc + Bu;      % predict next state vector
            xp = A*xc;
        end
        
        % update & evaluate new hypotheses via Kalman filter
        % predictions
        for u = 1 : 4 % extract old P()
            for v = 1 : 4
                TP(u,v)=P(k,i-1,u,v);
            end
        end
        PP = A*TP*A' + Q;    % predicted error
        % corrections
        K = PP*H'*inv(H*PP*H'+R);      % gain
                
        for iMarble=1:length(matDetectedMarbles)
            xDetected = matDetectedMarbles(iMarble,1);
            yDetected = matDetectedMarbles(iMarble,2);
            
            vecState(j,i,:) = (xp + K*([xDetected,yDetected]' - H*xp))';    % corrected state
            P(j,i,:,:) = (eye(4)-K*H)*PP;                    % corrected error
            
            % weight hypothesis by distance from observed data
            dvec = [xDetected,yDetected] - [vecState(j,i,1),vecState(j,i,2)];
            weights(j,i) = 1/(dvec*dvec');
            
            % draw some samples over one image
            if fig1 > 0
                figure(fig1)
                hold on
                for c = -0.99*radius: radius/10 : 0.99*radius
                    r = sqrt(radius^2-c^2);
                    if trackstate(j,i)==1                 % moving
                        plot(vecState(j,i,1)+c,vecState(j,i,2)+r,'c.')
                        plot(vecState(j,i,1)+c,vecState(j,i,2)-r,'c.')
                    elseif trackstate(j,i)==2             % collision
                        plot(vecState(j,i,1)+c,vecState(j,i,2)+r,'y.')
                        plot(vecState(j,i,1)+c,vecState(j,i,2)-r,'y.')
                    elseif trackstate(j,i)==3             % stop
                        plot(vecState(j,i,1)+c,vecState(j,i,2)+r,'b.');
                        plot(vecState(j,i,1)+c,vecState(j,i,2)-r,'b.');
                    else                                  % out of bounds
                        plot(vecState(j,i,1)+c,vecState(j,i,2)+r,'k.')
                        plot(vecState(j,i,1)+c,vecState(j,i,2)-r,'k.')
                    end
                end
            end
        end
        
        % rescale new hypothesis weights
        totalw=sum(weights(:,i)');
        weights(:,i)=weights(:,i)/totalw;
        
        % select top hypothesis to draw
        subset=weights(:,i);
        top = find(subset == max(subset));
        trackstate(top,i);
        
 
    end
        % display final top hypothesis
        if fig1 > 0
            figure(fig1)
            hold on
            for c = -0.99*radius: radius/10 : 0.99*radius
                r = sqrt(radius^2-c^2);
                %      plot(x(top,i,1)+c,x(top,i,2)+r+1,'b.')
                %      plot(x(top,i,1)+c,x(top,i,2)+r,'y.')
                plot(vecState(top,i,1)+c,vecState(top,i,2)+r,'r.')
                plot(vecState(top,i,1)+c,vecState(top,i,2)-r,'r.')
                %      plot(x(top,i,1)+c,x(top,i,2)-r,'y.')
                %      plot(x(top,i,1)+c,x(top,i,2)-r-1,'b.')
            end
            %    eval(['saveas(gcf,''COND/cond',int2str(i-1),'.jpg'',''jpg'')']);
        end
        
        pause(0.5) % wait a bit for the display to catch up   
    
    
end


end


function [ strTrackedMarbles,success ] = track_marbles( matMarbles )
%function track_marbles
%AV Practical 1
%This function receives centers and radius of detected marbles and 
%returns an struct of the tracked marbles, indexed by id.
%Portions of code taken from condense.m
%Input:  matMarbles, a cell array containing, per frame
%         a 3-column matrix with marbles in rows
%        and xcenter, ycenter and radius in the columns
%Output:  strTrackedMarbles - A cell array containing by frame
%         success - non-zero if there were problems, 0 otherwise
Images_marbles=true;
Tracking_marbles=true;
success=true;


% multiple condensation states
NCON=100;          % number of condensation samples
MAXTIME=60;        % number of time frames
x=zeros(NCON,MAXTIME,4);      % state vectors
weights=zeros(NCON,MAXTIME);  % est. probability of state
trackstate=zeros(NCON,MAXTIME);         % state=1,2,3;
P=zeros(NCON,MAXTIME,4,4);    % est. covariance of state vec.
for i = 1 : NCON              % initialize estimated covariance
for j = 1 : MAXTIME
  P(i,j,1,1) = 100;
  P(i,j,2,2) = 100;
  P(i,j,3,3) = 100;
  P(i,j,4,4) = 100;
end
end
pstop=0.05;      % probability of stopping vertical motion
pbounce=0.30;    % probability of bouncing at current state (overestimated)
xc=zeros(4,1);   % selected state
TP=zeros(4,4);   % predicted covariance

% loop over all images
fig1=1;
fig2=0;
fig15=0;
fig3=0;
for i = 1 : MAXTIME

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

  % extract ball
  [cc(i),cr(i),radius,flag]=extractball(Imwork,Imback,fig1,fig2,fig3,fig15,i);
  if flag==0
    for k = 1 : NCON
      x(k,i,:) = [floor(MC*rand(1)),floor(MR*rand(1)),0,0]';
      weights(k,i)=1/NCON;
    end
    continue
  end

  % display green estimated ball circle over original image
  if fig1 > 0
    figure(fig1)
    hold on
    for c = -0.99*radius: radius/10 : 0.99*radius
      r = sqrt(radius^2-c^2);
      plot(cc(i)+c,cr(i)+r,'g.')
      plot(cc(i)+c,cr(i)-r,'g.')
    end
  end

  % condensation tracking
  % generate NCON new hypotheses from current NCON hypotheses
  % first create an auxiliary array ident() containing state vector j
  % SAMPLE*p_k times, where p is the estimated probability of j
  if i ~= 1
    SAMPLE=10;
    ident=zeros(100*SAMPLE,1);
    idcount=0;
    for j = 1 : NCON    % generate sampling distribution
      num=floor(SAMPLE*100*weights(j,i-1));  % number of samples to generate
      if num > 0
        ident(idcount+1:idcount+num) = j*ones(1,num);
        idcount=idcount+num;
      end
    end
  end

  % generate NCON new state vectors
  for j = 1 : NCON

    % sample randomly from the auxiliary array ident()
    if i==1 % make a random vector
      xc = [floor(MC*rand(1)),floor(MR*rand(1)),0,0]';
    else
      k = ident(ceil(idcount*rand(1))); % select which old sample
      xc(1) = x(k,i-1,1);  % get its state vector
      xc(2) = x(k,i-1,2);
      xc(3) = x(k,i-1,3);
      xc(4) = x(k,i-1,4);

      % sample about this vector from the distribution (assume no covariance)
      for n = 1 : 4
        xc(n) = xc(n) + 5*sqrt(P(j,i-1,n,n))*randn(1);
      end
    end

    % hypothesize if it is going into a bounce or tabletop state
    if i == 1    % initial time - assume falling
      xp = xc;   % no process at start
      A = A3;
      Bu = Bu3;
      trackstate(j,i)=3;
    else
      if trackstate(k,i-1)==1  % if already stopped bouncing
        A = A1;
        Bu = Bu1;
        xc(4) = 0;
        trackstate(j,i)=1;     % stay stopped bouncing
      else
        r=rand(1);   % random number for state selection
        if r < pstop
          A = A1;
          Bu = Bu1;
          xc(4) = 0;
          trackstate(j,i)=1;
        elseif r < (pbounce + pstop)
          A = A2;
          Bu = Bu2;
          % add some random vertical motion due to imprecision
          % about time of bounce
          xc(2) = xc(2) + 3*abs(xc(4))*(rand(1)-0.5); 
          xc(4) = -xc(4)*loss;  % invert vertical velocity (lossy)
          trackstate(j,i)=2;  % set into bounce state
        else % normal motion
          A = A3;
          Bu = Bu3;
          trackstate(j,i)=3;
        end
      end
      xp=A*xc + Bu;      % predict next state vector
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
    x(j,i,:) = (xp + K*([cc(i),cr(i)]' - H*xp))';    % corrected state
    P(j,i,:,:) = (eye(4)-K*H)*PP;                    % corrected error

    % weight hypothesis by distance from observed data
    dvec = [cc(i),cr(i)] - [x(j,i,1),x(j,i,2)];
    weights(j,i) = 1/(dvec*dvec');

    % draw some samples over one image
    if i == 15 & fig1 > 0
      figure(fig1)
      hold on
      for c = -0.99*radius: radius/10 : 0.99*radius
        r = sqrt(radius^2-c^2);
        if trackstate(j,i)==1                 % stop
          plot(x(j,i,1)+c,x(j,i,2)+r,'c.')
          plot(x(j,i,1)+c,x(j,i,2)-r,'c.')
        elseif trackstate(j,i)==2             % bounce
          plot(x(j,i,1)+c,x(j,i,2)+r,'y.')
          plot(x(j,i,1)+c,x(j,i,2)-r,'y.')
        else                                  % normal
          plot(x(j,i,1)+c,x(j,i,2)+r,'k.')
          plot(x(j,i,1)+c,x(j,i,2)-r,'k.')
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

  % display final top hypothesis
  if fig1 > 0
    figure(fig1)
    hold on
    for c = -0.99*radius: radius/10 : 0.99*radius
      r = sqrt(radius^2-c^2);
%      plot(x(top,i,1)+c,x(top,i,2)+r+1,'b.')
%      plot(x(top,i,1)+c,x(top,i,2)+r,'y.')
      plot(x(top,i,1)+c,x(top,i,2)+r,'r.')
      plot(x(top,i,1)+c,x(top,i,2)-r,'r.')
%      plot(x(top,i,1)+c,x(top,i,2)-r,'y.')
%      plot(x(top,i,1)+c,x(top,i,2)-r-1,'b.')
    end
%    eval(['saveas(gcf,''COND/cond',int2str(i-1),'.jpg'',''jpg'')']);  
  end

      pause(0.5) % wait a bit for the display to catch up
end


end


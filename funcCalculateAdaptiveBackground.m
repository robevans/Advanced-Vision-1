function [success] = funcCalculateAdaptiveBackground (sDirectory, bForceRecalculation, bShowPlots)
%function funcCalculateAdaptiveBackground
%AV 1 Practical 1 20140201
%This function calculates a background image using a probabilistic method
%using chromaticity
%It takes the number of images to use from the funcConfig method parameter
%'nImagesAverageBackground'
%Input:  sDirectory - Directory where images are located
%        bForceRecalculation - force recalculation, otherwise use existing file is one exists
%        bShowPlots - display plots about the background distribution
%Output:  imageBackGround - An image of the same size of the initial image
%                           with the background calculation
%Authors: Robert Evans, Francisco Aguirre


success=false;

if nargin<2
    bForceRecalculation=false;
end;

if nargin<3
    bShowPlots=false;
end;

if ~bForceRecalculation
    if exist(funcConfig('fileAdaptiveBackground'),'file')==2
        success=true;
        return;
    end;
end;

iNumBackgroundImages=funcConfig('nImagesAverageBackground');

% use N images to estimate initial statistics of background model.
for i = 1 : iNumBackgroundImages
    
    % load image
    imageCurrentRGB = imread([strcat(sDirectory,int2str(i)), '.jpg'],'jpg');
    if i==1
        [iRows,iColumns,iDimension] = size(imageCurrentRGB);
        imagesBackGround = zeros(iRows,iColumns,3,iNumBackgroundImages);
        imageSigma = zeros (iRows, iColumns);
        uint8BackGround=uint8(zeros(iRows,iColumns,3,iNumBackgroundImages));
    end
    uint8BackGround(:,:,:,i+1)=imageCurrentRGB;
    
    % initialize background model for (normalised red, normalised green, saturation)
    if i==(iNumBackgroundImages-1)
        
        % Convert to chromaticity
        rgbsum = sum(uint8BackGround,3);
        rgbsum(rgbsum(:,:,1,:)==0) = 1; % Avoid divide by zero errors.
        imagesBackGround(:,:,1,:) = double(uint8BackGround(:,:,1,:)) ./ rgbsum; % Red Chromaticity
        imagesBackGround(:,:,2,:) = double(uint8BackGround(:,:,2,:)) ./ rgbsum; % Green Chromaticity
        imagesBackGround(:,:,3,:) = rgbsum / 3; % Saturation
        
        % estimate sigma using robust estimator to avoid outliers
        bgdiff = abs(imagesBackGround(:,:,1,1:(iNumBackgroundImages-1)) - imagesBackGround(:,:,1,2:iNumBackgroundImages));
        % compute sigma
        imageSigma = median(bgdiff,4) ./ (0.68*sqrt(2));
        
        % hypothesise initial history with median value in case foreground present
        medians = median(imagesBackGround,4);
        imagesBackGround = repmat(medians,1,1,1,iNumBackgroundImages)+repmat(imageSigma,1,1,3,iNumBackgroundImages).*randn(iRows,iColumns,3,iNumBackgroundImages);
    end
end

if bShowPlots
    % show std deviation across image
    figure(1)
    image(imageSigma*2000)
    colormap(gray)
    figure(2)
    
    % plot reconstructed distribution
    for k=1:iNumBackgroundImages
        x(k) = imagesBackGround(100,100,1,k);
        xc(k) = k/iNumBackgroundImages;
    end
    xh = hist(x,xc);
    plot(xc,xh)
end

%save background model for detection
save(funcConfig('fileAdaptiveBackground'),'imageSigma','imagesBackGround');

success=true;


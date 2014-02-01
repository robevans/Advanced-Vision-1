function [success] = funcCalculateAdaptiveBackground (sDirectory, iCurrentImage, bForceRecalculation)
%function funcCalculateAdaptiveBackground
%AV 1 Practical 1 20140201
%This function calculates a background image using an adaptive method
%using chromaticity
%It takes the number of images to use from the funcConfig method parameter
%'nImagesAverageBackground'
%Input:  sDirectory - Directory where images are located
%        iCurrentImage - index number of the current image to process
%        bForceRecalculation - force recalculation, otherwise use existing file is one exists
%Output:  imageBackGround - An image of the same size of the initial image
%                           with the background calculation
%Authors: Robert Evans, Francisco Aguirre


success=false;

if nargin<3
    bForceRecalculation=false;
end;    

if ~bForceRecalculation
    if exist(funcConfig('fileAdaptiveBackground'),'file')==2
       success=true;
       return;
    end;
end;    

iTotalBackgroundImages=funcConfig('nImagesAverageBackground');
bgdiff=zeros((iTotalBackgroundImages-1),1);
iCurrentRow=0;

% use first 50 images to estimate statistics
for i = 1 : iTotalBackgroundImages
% load image
    imageCurrentRGB = imread([strcat(sDirectory,int2str(i)), '.jpg'],'jpg');
    if i==1 
        [iRows,iColumns,iDimension] = size(imageCurrentRGB);       
        imagesBackGround = zeros(iRows,iColumns,3,iTotalBackgroundImages);
        imageSigma = zeros (iRows, iColumns);
        uint8BackGround=uint8(zeros(iRows,iColumns,3,iTotalBackgroundImages));       
    end
    uint8BackGround(:,:,:,i+1)=imageCurrentRGB;

% initialize background model for (normalised red, normalised green, saturation)
    if i==(iTotalBackgroundImages-1)
      for r = 1:iRows
          iCurrentRow = r
        for c = 1 : iColumns
          % convert to chromaticity
	  for k = 1 : iTotalBackgroundImages
	    R=double(uint8BackGround(r,c,1,k)); 
	    G=double(uint8BackGround(r,c,2,k)); 
	    B=double(uint8BackGround(r,c,3,k)); 
            rgbsum=(R+G+B);
            sat = rgbsum/3;
            if rgbsum == 0
                rgbsum=1;
            end
            red = R/rgbsum;
            green = G/rgbsum;
            imagesBackGround(r,c,1,k)= red;
            imagesBackGround(r,c,2,k)= green;
            imagesBackGround(r,c,3,k)= sat;
          end

          % estimate sigma using robust estimator to avoid outliers
	  for k = 1 : (iTotalBackgroundImages-1)
	    bgdiff(k) = abs(imagesBackGround(r,c,1,k) - imagesBackGround(r,c,1,k+1));
          end

          % compute sigma
          imageSigma(r,c)=median(bgdiff)/(0.68*sqrt(2));
          
          % hypothesise initial history with median value in case foreground present
          rmed = median(imagesBackGround(r,c,1,:));
          imagesBackGround(r,c,1,:) = rmed + imageSigma(r,c)*randn(50,1);
          gmed = median(imagesBackGround(r,c,2,:));
          imagesBackGround(r,c,2,:) = gmed + imageSigma(r,c)*randn(50,1);
          smed = median(imagesBackGround(r,c,3,:));
          imagesBackGround(r,c,3,:) = smed + imageSigma(r,c)*randn(50,1);
        end
      end

      % show std deviation across image
      figure(1)
      image(imageSigma*2000)
      colormap(gray)
      figure(2)

      % plot reconstructed distribution
      for k=1:iTotalBackgroundImages
  	    x(k) = imagesBackGround(100,100,1,k);
        xc(k) = k/50;
      end
      xh = hist(x,xc);
      plot(xc,xh)

      %save background model for detection
      save(funcConfig('fileAdaptiveBackground'),'imageSigma','imagesBackGround');
    end
success=true;
end

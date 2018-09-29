function outputFactor = ppbkc_outputFactorCorrection(outputFactor,primaryFluence,kernelExtension,kernelResolution,FWHM)

sigma = FWHM/(sqrt(8*log(2))); % mm
sigmaVox = sigma/kernelResolution; % vox
iCenter = kernelExtension/2;

[X,Y] = meshgrid(-kernelExtension/2+1:kernelExtension/2);
Z = sqrt(X.^2 + Y.^2)*0.5;
primflu = interp1(primaryFluence(:,1),primaryFluence(:,2),Z,'linear',0);

for i = 1:numel(outputFactor(:,1))
        
    if outputFactor(i,1) < 5.4 * sqrt(2)*sigma
        
        lowerLimit = round(iCenter - outputFactor(i,1)/2/kernelResolution + 1);
        upperLimit = floor(iCenter + outputFactor(i,1)/2/kernelResolution);

        fieldShape = 0*primflu;
        fieldShape(lowerLimit:upperLimit,lowerLimit:upperLimit) = 1;

         gaussFilter = 1/(sqrt(2*pi)*sigmaVox) .* exp( -X.^2 / (2*sigmaVox^2) ) ...
                    .* 1/(sqrt(2*pi)*sigmaVox) .* exp( -Y.^2 / (2*sigmaVox^2) );
               

        convRes = fftshift( ifft2(    fft2(fieldShape .*primflu,kernelExtension,kernelExtension) ...
                                   .* fft2(gaussFilter,kernelExtension,kernelExtension) ) );
                               
        outputFactor(i,2) = outputFactor(i,2) / convRes(iCenter-1,iCenter-1);
                     
    end
    
end

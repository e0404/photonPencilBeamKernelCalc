function fRNorm = ppbkc_calcKernelNorm(kernelExtension,kernelResolution,primaryFluence)

mode = 'radialInt'; % can be either 'pdc' or 'radialInt'

if strcmp(mode,'radialInt')
    
    dr    = kernelResolution;
    r     = dr*[0:(kernelExtension/2-1)];
    r_mid = (r(1:end-1)+r(2:end))/2;

    % compute integral
    int = dr*2*pi*r_mid.*interp1(primaryFluence(:,1), ...
                                 primaryFluence(:,2), ...
                                 r_mid);

    % interpolate data on original grid
    fRNorm = 4*interp1(r_mid,int,r,'spline')';

    % wo kommt der bekackte faktor 4 her? wieso brauch ich den, damit es
    % wieder passt? --> vermutlich fehlt ein faktor 1/4 beim anderen kern!
    
elseif strcmp(mode,'pdc')
    
    iMaxirradtiation_mm  = kernelExtension/2*kernelResolution;
    iNumofRadii          = iMaxirradtiation_mm/kernelResolution;

    [X,Y] = meshgrid(-kernelExtension/2+1:kernelExtension/2);
    Z = sqrt(X.^2 + Y.^2)*0.5;
    primflu = interp1(primaryFluence(:,1),primaryFluence(:,2),Z,'linear',0);

    fResolutionFactor = 1/kernelResolution;

    fRNorm = zeros(iNumofRadii,1);

    % this is a radial integration over the primary fluence!
    % this may be an option to improve the numerical stability of the code here
    % by replacing the integration with a sound approach!
    for i = 1:kernelExtension
        for j=1:kernelExtension
            fRadius_mm(i,j) = sqrt(((kernelExtension/2-i+1)*kernelResolution)^2 + ((kernelExtension/2-j+1)*kernelResolution)^2);
            fRadius_RF(i,j) = fRadius_mm(i,j).*fResolutionFactor;
            iRadius_RF = round(fRadius_RF(i,j),0);
            viRadius_RF(i,j) = iRadius_RF;
            if iRadius_RF < iNumofRadii
                % it it not ok to add a the primary fluence in full. actually
                % you need to scale this by the area of the pixel i.2.
                % fKernelResolution^2!
                fRNorm(iRadius_RF+1) = fRNorm(iRadius_RF+1)+primflu(i,j);
            end
        end
    end
    
end
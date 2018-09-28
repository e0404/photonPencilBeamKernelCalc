%% run import script to read in the data
ppbkc_readData

%% set parameters

SAD = params('SAD');
SSD = 870;

%% scale tpr

fieldSizeScaleFactorTpr = (SSD+tprDepths)/SAD;

scaledTpr = interp2(tprFieldSizes,tprDepths,tpr, ...
                fieldSizeScaleFactorTpr*tprFieldSizes, ...
                tprDepths*ones(1,numel(tprFieldSizes)), ...
                'spline');
            
%% normalize TPR_0, determine maxima and mu

[tprMax,tprMaxIx] = max(tpr);
meanMaxPos_mm     = ceil(mean(tprDepths(tprMaxIx)));

tpr_0 = tpr(:,1)/tprMax(1);

% compute mu for TPR0 with exponential fit, only use data points behind max
% (neglects build-up) 
[~,ix] = min(abs(tprDepths-meanMaxPos_mm));

fSx  = sum(tprDepths(ix+1:end));
fSxx = sum(tprDepths(ix+1:end).^2);
ftmp = -log(tpr_0(ix+1:end));
fSy  = sum(ftmp);
fSxy = sum(ftmp.*tprDepths(ix+1:end));

% mu = 0.005066; % reference value from literature
mu = (fSxy-((fSx*fSy)/length(tpr_0(ix+1:end)))) / ...
            (fSxx-(fSx^2/length(tpr_0(ix+1:end))));
        

%% compute betas

maxPos_fun = @(x) (log(mu)-log(x))/(mu-x);

beta(1) = fmincon(@(x) (maxPos_fun(x) - meanMaxPos_mm).^2,1,[],[],[],[],0,1000);
beta(2) = fmincon(@(x) (maxPos_fun(x) - (meanMaxPos_mm+1/mu)/2).^2,1,[],[],[],[],0,1000);
beta(3) = fmincon(@(x) (maxPos_fun(x) - 1/mu).^2,1,[],[],[],[],0,1000);

%% get weights for indivdual compoents of all TPRs
D_1 = (beta(1)/(beta(1)-mu))*(exp(-mu*tprDepths(ix+1:end))-exp(-beta(1)*tprDepths(ix+1:end)));
D_2 = (beta(2)/(beta(2)-mu))*(exp(-mu*tprDepths(ix+1:end))-exp(-beta(2)*tprDepths(ix+1:end)));
D_3 = (beta(3)/(beta(3)-mu))*(exp(-mu*tprDepths(ix+1:end))-exp(-beta(3)*tprDepths(ix+1:end)));

mx1 = [D_1 D_2 D_3]'*[D_1 D_2 D_3];
mx2 = [D_1 D_2 D_3]'*scaledTpr(ix+1:end,:);

W_ri = (mx1\mx2)';

%% compute normalization for kernel
kernelExtension = 512; % pixel
kernelResolution = 0.5; % mm
fRNorm = ppbkc_calcKernelNorm(kernelExtension,kernelResolution,primaryFluence);

%% computation of correction factors for output factor

FWHM  = params('fwhm_gauss'); % mm
sigma = FWHM/(sqrt(8*log(2))); % mm
sigmaVox = sigma/kernelResolution; % vox
iCenter = kernelExtension/2;

correctionFactors = ones(size(outputFactor,1),1);

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
        correctionFactors(i) = 1/convRes(iCenter-1,iCenter-1);
                     
    end
    
end

%% compute equivalent field size for circular fields!
fEquivalentFieldSize = [1:kernelExtension/2].*kernelResolution*sqrt(pi);

%% calculate corrected output factors at equivalent field size
correctedOutputFactor = interp1(outputFactor(:,1),outputFactor(:,2).*correctionFactors,fEquivalentFieldSize, 'linear', 'extrap');

%% compute kernels
D_1_spline = interp1(tprFieldSizes,W_ri(:,1),fEquivalentFieldSize, 'spline');
D_2_spline = interp1(tprFieldSizes,W_ri(:,2),fEquivalentFieldSize, 'spline');
D_3_spline = interp1(tprFieldSizes,W_ri(:,3),fEquivalentFieldSize, 'spline');


fGradFitPar_1 = [correctedOutputFactor(1)*D_1_spline(1),diff(correctedOutputFactor.*D_1_spline)];
fGradFitPar_2 = [correctedOutputFactor(1)*D_2_spline(1),diff(correctedOutputFactor.*D_2_spline)];
fGradFitPar_3 = [correctedOutputFactor(1)*D_3_spline(1),diff(correctedOutputFactor.*D_3_spline)];

kernel1 = fGradFitPar_1'./fRNorm;
kernel2 = fGradFitPar_2'./fRNorm;
kernel3 = fGradFitPar_3'./fRNorm;

%% plot
matRadRoot = 'C:\Home\Bangertm\Git\matRad\';
load([matRadRoot 'photons_Generic.mat']);

matRadSsdIx = find([machine.data.kernel.SSD]==SSD);

figure
semilogx(machine.data.kernelPos,machine.data.kernel(matRadSsdIx).kernel1,'r')
hold on
plot(machine.data.kernelPos,machine.data.kernel(matRadSsdIx).kernel2,'b')
plot(machine.data.kernelPos,machine.data.kernel(matRadSsdIx).kernel3,'g')
plot(machine.data.kernelPos(1:256),kernel1,'r--')
plot(machine.data.kernelPos(1:256),kernel2,'b--')
plot(machine.data.kernelPos(1:256),kernel3,'g--')

axis([0 250 -1e-4 10e-4])

xlabel('mm')
ylabel('a.u.')
box on
grid minor
legend({'kernel 1 (matRad)','kernel 2 (matRad)','kernel 3 (matRad)' ...
       ,'kernel 1 (calc.)','kernel 2 (calc.)','kernel 3 (calc.)'});

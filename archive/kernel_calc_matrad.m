clear all
tic
load 'artiste_dkfz_data.mat'

% scale TPRs to current SSD
ssd = 870;
attSAD = 1000;
tpr_dataOld = tpr_data;

for i = 1:numel(fieldsize)
    for j = 1:numel(tpr_data(:,1))
        
        scaleFactor = (ssd+tpr_data(j,1))/attSAD;
        
        tpr_data(j,i+1) = interp2(fieldsize,tpr_data(:,1),tpr_dataOld(:,2:end),scaleFactor*fieldsize(i),tpr_data(j,1),'spline');
        
    end
end

% use new output factor
outputFactorOld = outputFactor;
outputFactor = [5.7433	 0.67650
10.1112	 0.73272
12.0144	 0.76148
13.9094	 0.78703
15.8798	 0.80749
17.8449	 0.8246
19.8299	 0.8352
29.7900	 0.8742
40.0000	 0.9032
50.0000	 0.9250
60.0000	 0.9447
80.0000	 0.9763
100.0000 1.0000
150.0000 1.0381
200.0000 1.0659
300.0000 1.0949
400.0000 1.1065];

% use new primary fluence
primaryFluence = [0.0	1.0
14.14	1.001290323
28.28	1.023125781
42.43	1.041875469
56.57	1.048037009
70.71	1.050857714
84.85	1.054963741
98.99	1.061915479
113.14	1.067296824
127.28	1.075788947
141.42	1.082460615
155.56	1.079379845
169.71	1.073728432
183.85	1.045476369
197.99	0.937324331
212.13	0.833028257
226.27	0.769582396
240.42	0.66324081
254.56	0.438434609
261.63	0.295868967
264.46	0.247421855
267.29	0.206576644
270.11	0.173953488
272.94	0.146446612
275.77	0.123295824
278.6	0.102315579
281.43	0.083020755
284.26	0.066886722
287.09	0.055458865
289.91	0.048987247
292.74	0.044596149
295.57	0.040970243
296.98	0.039689922
302.64	0.033933483
308.3	0.028872218
313.96	0.024456114
319.61	0.020525131
325.27	0.017494374
350.0   0.0];

% correct negative values

% TPR-Kurven normieren ,Maximas und mu bestimmen
tpr_size = size(tpr_data);
for i = 2:tpr_size(2)
    idx                  = find(tpr_data(:,i) == max(tpr_data(:,i)));
    Dmax_value(i-1)      = mean(tpr_data(idx,i));
    dmax_mm(i-1)         = mean(tpr_data(idx,1));
    tpr_data_renorm(:,i) = tpr_data(:,i)/Dmax_value(i-1);
    tpr_data_renorm(:,1) = tpr_data(:,1);
end

% compute mu for TPR0 with exponential fit, only use data points behind max
% (neglecgts build-up) 
d_max_average = ceil(mean(dmax_mm));
tpr_0         = tpr_data_renorm(:,1:2);
[c index]     = min(abs(tpr_0(:,1)-d_max_average));
fSx           = sum(tpr_0(index+1:end,1));
fSxx          = sum(tpr_0(index+1:end,1).^2);
ftmp          = -log(tpr_0(index+1:end,2));
fSy           = sum(ftmp);
fSxy          = sum(ftmp.*tpr_0(index+1:end,1));
fAttCoeff     = (fSxy-((fSx*fSy)/length(tpr_0(index+1:end,1))))/(fSxx-(fSx^2/length(tpr_0(index+1:end,1))));
% fAttCoeff = 0.005066; % reference value from literature

% get betas for the three components of individual tprs
kernel = [d_max_average; 0.5*(d_max_average+(1/fAttCoeff)); 1/fAttCoeff];
for j = 1:3
    s = 1/kernel(j);
    fDMax = 2;
    fDMaxOld = 0;
    i=0;
    while (abs((fDMaxOld-fDMax)/fDMax) > 0.000001)
        i=i+1;
        fDMaxOld = fDMax;
        fDMax = -1 / s * log(fAttCoeff / (fAttCoeff + s));
        s = s * fDMax / kernel(j);
    end
    beta(j) = s + fAttCoeff;
end

% get weights for indivdual compoents of all TPRs
D_1=(beta(1)/(beta(1)-fAttCoeff))*(exp(-fAttCoeff*tpr_data(:,1))-exp(-beta(1)*tpr_data(:,1)));
D_2=(beta(2)/(beta(2)-fAttCoeff))*(exp(-fAttCoeff*tpr_data(:,1))-exp(-beta(2)*tpr_data(:,1)));
D_3=(beta(3)/(beta(3)-fAttCoeff))*(exp(-fAttCoeff*tpr_data(:,1))-exp(-beta(3)*tpr_data(:,1)));

vec_D1_sig = D_1(index+1:end);
vec_D2_sig = D_2(index+1:end);
vec_D3_sig = D_3(index+1:end);
data_sig = tpr_data(index+1:end,2:end);

mx1 = [vec_D1_sig vec_D2_sig vec_D3_sig]'*[vec_D1_sig vec_D2_sig vec_D3_sig];
mx2 = [vec_D1_sig vec_D2_sig vec_D3_sig]'*data_sig;

W_ri = (mx1\mx2)';

%%


iKernelExtension     = 512; %512
distance             = nan(iKernelExtension,iKernelExtension);
iCenterX             = iKernelExtension/2;
iCenterY             = iCenterX;
kd.fKernelResolution = 0.5; % [mm]
iMaxirradtiation_mm  = iKernelExtension/2*kd.fKernelResolution;
iNumofRadii          = iMaxirradtiation_mm/kd.fKernelResolution;
maxDist              = primaryFluence(end,1);



x = -iKernelExtension/2+1:iKernelExtension/2;
y = x;
[X,Y] = meshgrid(x,y);
% sqrt((iCenterX-i+1)^2 + (iCenterY-j+1)^2)* kd.fKernelResolution;
Z = sqrt(X.^2 + Y.^2)*0.5;
primflu = interp1(primaryFluence(:,1),primaryFluence(:,2),Z);
primflu(Z > maxDist) = 0;



fResolutionFactor = 1/kd.fKernelResolution;
iNumOfRadii_step  = 1:iNumofRadii;


fRNorm = zeros(iNumofRadii,1);

% this is a radial integration over the primary fluence!
% this may be an option to improve the numerical stability of the code here
% by replacing the integration with a sound approach!
for i = 1:iKernelExtension
    for j=1:iKernelExtension
        fRadius_mm(i,j) = sqrt(((iKernelExtension/2-i+1)*kd.fKernelResolution)^2 + ((iKernelExtension/2-j+1)*kd.fKernelResolution)^2);
        fRadius_RF(i,j) = fRadius_mm(i,j).*fResolutionFactor;
    %     iRadius_RF = ceil(fRadius_RF+0.5);
        iRadius_RF = round(fRadius_RF(i,j),0);
        viRadius_RF(i,j)=iRadius_RF;
        if iRadius_RF < iNumofRadii
            fRNorm(iRadius_RF+1) = fRNorm(iRadius_RF+1)+primflu(i,j);
        end
    end
end

%% computation of correction factors for output factor

FWHM  = 5; % mm
sigma = FWHM/(sqrt(8*log(2))); % mm
sigmaVox = sigma/kd.fKernelResolution; % vox

for i = 1:numel(outputFactor(:,1))
        
    if outputFactor(i,1) < 5.4 * sqrt(2)*sigma
        
        lowerLimit = round(iCenterX - outputFactor(i,1)/2/kd.fKernelResolution + 1);
        upperLimit = floor(iCenterX+outputFactor(i,1)/2/kd.fKernelResolution);

        fieldShape = 0*primflu;
        fieldShape(lowerLimit:upperLimit,lowerLimit:upperLimit) = 1;

         gaussFilter = 1/(sqrt(2*pi)*sigmaVox) .* exp( -X.^2 / (2*sigmaVox^2) ) ...
                    .* 1/(sqrt(2*pi)*sigmaVox) .* exp( -Y.^2 / (2*sigmaVox^2) );
               

        convRes = fftshift( ifft2(    fft2(fieldShape .*primflu,iKernelExtension,iKernelExtension) ...
                                   .* fft2(gaussFilter,iKernelExtension,iKernelExtension) ) );
        correctionFactor = 1/convRes(iCenterX-1,iCenterY-1)            
                     
    end
    
end

%% compute equivalent field size for circular fields!
fEquivalentFieldSize = iNumOfRadii_step.*kd.fKernelResolution*sqrt(pi);

fScp = interp1(outputFactor(:,1),outputFactor(:,2),fEquivalentFieldSize, 'linear', 'extrap');

D_1_spline=interp1(fieldsize,W_ri(:,1),fEquivalentFieldSize, 'spline');
D_2_spline=interp1(fieldsize,W_ri(:,2),fEquivalentFieldSize, 'spline');
D_3_spline=interp1(fieldsize,W_ri(:,3),fEquivalentFieldSize, 'spline');


fGradFitPar_1=[fScp(1)*D_1_spline(1),diff(fScp.*D_1_spline)];
fGradFitPar_2=[fScp(1)*D_2_spline(1),diff(fScp.*D_2_spline)];
fGradFitPar_3=[fScp(1)*D_3_spline(1),diff(fScp.*D_3_spline)];

Kern1 = fGradFitPar_1'./fRNorm;
Kern2 = fGradFitPar_2'./fRNorm;
Kern3 = fGradFitPar_3'./fRNorm;

toc;

%%
load C:\Home\Bangertm\Git\matRad\photons_Generic.mat

if ssd == 950
    
    load pdc_kernels_ssd950.mat
    for i = 1:size(kernelnewssd950,1)

        refKernelPdc{kernelnewssd950(i,1)+1}(kernelnewssd950(i,2)+1,kernelnewssd950(i,3)+1) = kernelnewssd950(i,4);

    end

elseif ssd == 870
    
    load pdc_kernels_ssd870.mat
    for i = 1:size(kernelnewssd870,1)

        refKernelPdc{kernelnewssd870(i,1)+1}(kernelnewssd870(i,2)+1,kernelnewssd870(i,3)+1) = kernelnewssd870(i,4);

    end
end    
    
figure
subplot(2,1,1)
hold on
plot(machine.data.kernelPos(1:256),refKernelPdc{1}(1,1:256),'or')
plot(machine.data.kernelPos(1:256),Kern1,'r')
plot(machine.data.kernelPos(1:256),refKernelPdc{2}(1,1:256),'bo')
plot(machine.data.kernelPos(1:256),Kern2,'b')
plot(machine.data.kernelPos(1:256),refKernelPdc{3}(1,1:256),'go')
plot(machine.data.kernelPos(1:256),Kern3,'g')
axis([0 250 -1e-4 10e-4])

title(['ssd ' num2str(ssd) ' mm'])
xlabel('mm')
ylabel('a.u.')
box on
grid minor
legend({'pdc kernel 1','matRad kernel 1','pdc kernel 2','matRad kernel 2','pdc kernel 3','matRad kernel 3'})

subplot(2,1,2)
semilogy(machine.data.kernelPos(1:256),abs(refKernelPdc{1}(1,1:256)),'or')
hold on
plot(machine.data.kernelPos(1:256),abs(Kern1),'r')
plot(machine.data.kernelPos(1:256),abs(refKernelPdc{2}(1,1:256)),'bo')
plot(machine.data.kernelPos(1:256),abs(Kern2),'b')
plot(machine.data.kernelPos(1:256),abs(refKernelPdc{3}(1,1:256)),'go')
plot(machine.data.kernelPos(1:256),abs(Kern3),'g')
axis([0 250 -1e-4 10e-4])

title(['ssd ' num2str(ssd) ' mm'])
xlabel('mm')
ylabel('a.u.')
box on
grid minor
legend({'pdc kernel 1','matRad kernel 1','pdc kernel 2','matRad kernel 2','pdc kernel 3','matRad kernel 3'})
function machine = ppbkc_generateBaseData(name,inputDir)

%% read parameter file
fileHandle = fopen([inputDir filesep 'params.dat']);
tmp = textscan(fileHandle,'%s %f','CommentStyle',{'#'});
fclose(fileHandle);

params = containers.Map(tmp{1},tmp{2});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% write meta information
machine.meta.radiationMode = 'photons';
machine.meta.dataType      = '-';
machine.meta.created_on    = date;
machine.meta.created_by    = mfilename;
machine.meta.description   = ['photon pencil beam kernels calculated with ' mfilename ' for custom data'];
machine.meta.name          = name;
machine.meta.SAD           = params('SAD');
machine.meta.SCD           = params('source_collimator_distance');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% write data
machine.data.energy         = params('photon_energy');
machine.data.kernelPos      = [0:0.5:179.5];
machine.data.fwhm           = params('fwhm_gauss');
machine.data.iNumKernel     = params('i_num_kernel');
machine.data.electronRangeIntensity     = params('electron_range_intensity');

%% load primary fluence
fileHandle = fopen([inputDir filesep 'primflu.dat']);
machine.data.primaryFluence = cell2mat(textscan(fileHandle,'%f %f','CommentStyle',{'#'}));
fclose(fileHandle);

%% compute attenuation coefficent
% load tpr
tprTmp = [];
fileHandle = fopen([inputDir filesep 'tpr.dat']);
while ~feof(fileHandle)
    currLine = fgetl(fileHandle);
    if currLine(1) ~= '#'
        % fprintf([currLine '\n']);
        tprTmp = [tprTmp; str2num(currLine)];
    end
end
fclose(fileHandle);

% check if field size 0 included, if not interpolate
if tprTmp(1,2) ~= 0
    
    % interpolate
    tprZero = NaN*ones(size(tprTmp,1)-1,1);
    for i = 1:numel(tprZero)
        tprZero(i) = interp1(tprTmp(1,2:end),tprTmp(i+1,2:end),0,'linear','extrap');        
    end

    % insert
    tprTmp = [tprTmp(:,1) [0;tprZero] tprTmp(:,2:end)];
    
end

% refactor into field size, depths and tpr values
tprFieldSizes = tprTmp(1,2:end);
tprDepths     = tprTmp(2:end,1);
tpr           = tprTmp(2:end,2:end);

machine.data.surfaceDose = tpr(1,1);

% uncomment for debug plot
% semilogx(tprDepths,tpr)
% legend(arrayfun(@(n) ['field size = ' num2str(n)],tprFieldSizes,'Uni',0))
% xlabel('[mm]')
% ylabel('rel. dose')

% find max positions
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

% machine.data.m = 0.005066; % reference value for 6MV from literature
machine.data.m = ( fSxy - ( (fSx*fSy) / length(tpr_0(ix+1:end)) ) ) / ...
                 ( fSxx - ( (fSx^2)   / length(tpr_0(ix+1:end)) ) );

%% compute betas
maxPos_fun = @(x) (log(machine.data.m)-log(x))/(machine.data.m-x);

options = optimoptions('fmincon','Display','off');

machine.data.betas(1) = fmincon(@(x) (maxPos_fun(x) - meanMaxPos_mm).^2,1,[],[],[],[],0,1000,[],options);
machine.data.betas(2) = fmincon(@(x) (maxPos_fun(x) - (meanMaxPos_mm+1/machine.data.m)/2).^2,1,[],[],[],[],0,1000,[],options);
machine.data.betas(3) = fmincon(@(x) (maxPos_fun(x) - 1/machine.data.m).^2,1,[],[],[],[],0,1000,[],options);

%% compute normalization for kernel
kernelExtension  = 720; % pixel
kernelResolution = 0.5; % mm
kernelNorm       = ppbkc_calcKernelNorm(kernelExtension,kernelResolution,machine.data.primaryFluence);

%% output factor
% read data
fileHandle = fopen([inputDir filesep 'of.dat']);
outputFactor = cell2mat(textscan(fileHandle,'%f %f','CommentStyle',{'#'}));
fclose(fileHandle);

% make correction for small fields
outputFactor = ppbkc_outputFactorCorrection(outputFactor, ...
                                            machine.data.primaryFluence, ...
                                            kernelExtension, ...
                                            kernelResolution, ...
                                            params('fwhm_gauss'));
                                        
%% compute equivalent field size for circular fields!
equivalentFieldSizes = [1:kernelExtension/2].*kernelResolution*sqrt(pi);

%% calculate corrected output factors at equivalent field sizes
correctedOutputFactorAtEquiFieldSizes = interp1(outputFactor(:,1),outputFactor(:,2),equivalentFieldSizes, 'linear', 'extrap');

%% compute kernels
machine.data.kernelPos = kernelResolution * [0:(kernelExtension/2-1)];

for i = 1:501

    % log SSD
    machine.data.kernel(i).SSD = i + 499; % [mm]
    
    % scale tpr
    fieldSizeScaleFactorTpr = (machine.data.kernel(i).SSD+tprDepths)/machine.meta.SAD;

    scaledTpr = interp2(tprFieldSizes,tprDepths,tpr, ...
                  fieldSizeScaleFactorTpr*tprFieldSizes, ...
                  tprDepths*ones(1,numel(tprFieldSizes)), ...
                  'spline');
              
    % compute weights of depth dose components 
    D_1 = (machine.data.betas(1)/(machine.data.betas(1)-machine.data.m)) * ...
        (exp(-machine.data.m*tprDepths(ix+1:end))-exp(-machine.data.betas(1)*tprDepths(ix+1:end)));
    D_2 = (machine.data.betas(2)/(machine.data.betas(2)-machine.data.m)) * ...
        (exp(-machine.data.m*tprDepths(ix+1:end))-exp(-machine.data.betas(2)*tprDepths(ix+1:end)));
    D_3 = (machine.data.betas(3)/(machine.data.betas(3)-machine.data.m)) * ...
        (exp(-machine.data.m*tprDepths(ix+1:end))-exp(-machine.data.betas(3)*tprDepths(ix+1:end)));

    mx1 = [D_1 D_2 D_3]'*[D_1 D_2 D_3];
    mx2 = [D_1 D_2 D_3]'*scaledTpr(ix+1:end,:);

    W_ri = (mx1\mx2)';
    W_ri(:,4) = tpr(1,:)./tpr(1,1);
    
    D_1_spline = interp1(tprFieldSizes,W_ri(:,1),equivalentFieldSizes, 'spline');
    D_2_spline = interp1(tprFieldSizes,W_ri(:,2),equivalentFieldSizes, 'spline');
    D_3_spline = interp1(tprFieldSizes,W_ri(:,3),equivalentFieldSizes, 'spline');
    D_4_spline = interp1(tprFieldSizes,W_ri(:,4),equivalentFieldSizes, 'spline');
    
    fGradFitPar_1 = [correctedOutputFactorAtEquiFieldSizes(1)*D_1_spline(1), ...
                     diff(correctedOutputFactorAtEquiFieldSizes.*D_1_spline)];
    fGradFitPar_2 = [correctedOutputFactorAtEquiFieldSizes(1)*D_2_spline(1), ...
                     diff(correctedOutputFactorAtEquiFieldSizes.*D_2_spline)];
    fGradFitPar_3 = [correctedOutputFactorAtEquiFieldSizes(1)*D_3_spline(1), ...
                     diff(correctedOutputFactorAtEquiFieldSizes.*D_3_spline)];
    fGradFitPar_4 = [correctedOutputFactorAtEquiFieldSizes(1)*D_4_spline(1), ...
                         diff(correctedOutputFactorAtEquiFieldSizes.*D_4_spline)];
        
    machine.data.kernel(i).kernel1 = fGradFitPar_1' ./ kernelNorm;
    machine.data.kernel(i).kernel2 = fGradFitPar_2' ./ kernelNorm;
    machine.data.kernel(i).kernel3 = fGradFitPar_3' ./ kernelNorm;
    machine.data.kernel(i).kernel4 = fGradFitPar_4' ./ kernelNorm;
    
end

% save file
save(strcat('photons_',name),'machine');
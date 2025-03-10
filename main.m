
%% Creating a mask for background using z-score normalization
IntensityMaskGUI(profiles);
%% Saving the mask and threshold for low resolution mask
intense_mask = double(intense_mask);
intense_mask_nand = intense_mask;
intense_mask_nand(intense_mask<1) = nan;
mask_trs = threshold;

%% Defining scan parameters
tr = 3.4; te = tr/2;
pc_step = round(360/size(profiles,3));
pcnum = 360/pc_step;
fa = 35;
field_strength = 2.89;

%% Phase correction step
profiles = ...
        profiles.*exp(-1i*angle(mean(profiles,3))); % Rotating profiles to real axis first

%%%%%% Creating a dictionary for dictionary-based phase correction %%%%%%%%
rtrarr = [5]; 
t2arr = 80;
% Only a set of 1 relaxation parameters are enoug for phase correction

freq_step_num = round(360/pc_step);
[b0, rtr, t2] = ndgrid(linspace(0,1e3/tr,freq_step_num),rtrarr,t2arr);
D1 = zeros(round(360/pc_step), numel(rtr));


for i = 1:numel(rtr)
    p = ...
        bSSFPAnalytic(rtr(i)*t2(i), ...
                        t2(i), ...
                        te, ...
                        tr, ...
                        fa, ...
                        0:pc_step:359, ...
                        b0(i));
    % Analytic equation of bSSFP was used to simulate dictionary entries
    D1(:,i) = conj((squeeze((p))));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[profiles_phase_corr, offset_map1] = ...
    PhaseCorrection((profiles), conj(D1));

profiles = profiles_phase_corr; % Saving phase-corrected profiles


%% Selecting the spectral model for fat

%%%%% LIVER %%%%%
% freq = [-3.800 -3.400 -2.600 -1.940 -0.390 0.600];
% amp = [0.087 0.693 0.128 0.004 0.039 0.048];
%%%%%%%%%%%%%%%%%

%%%% Calimetrix phantom (commercially avaliable) %%%%% 
% freq = [-3.80, -3.40, -2.60, -1.94, -0.39, 0.60];   
% amp = [0.087 0.693 0.128 0.004 0.039 0.048];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% BUTTERS phantom (custom-made) %%%%%%%%%%%%%%%%%%%%
freq = [-3.950 -3.538 -3.266 -2.809 -2.595 -2.079 -0.751 -0.549 0.490];
amp = [0.09 0.607 0.047 0.058 0.036 0.025 0.031 0.038 0.069];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% KNEE %%%%%%%%%%%%%%%%%%%%%%%%
% freq = [-3.8 -3.4 -3.11 -2.7 -2.45 -1.93 -0.5 0.61];
% amp = [0.087 0.568 0.058 0.092 0.058 0.027 0.038 0.073];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Going to lower resolution for faster field mapping

scale_factor = 4;
profiles(isnan(profiles)) = 0;
clear profiles_lowres
for i = 1:size(profiles,3)
    profiles_lowres(:,:,i) = imresize(squeeze(profiles(:,:,i)), 1/scale_factor);
end


%% Calculating residuals 

freqhz = freq*42.6*field_strength; % Convert ppm to Hz

b0range = linspace(-400,400,101); % Defining the inhomogeneity range (Steps of 8Hz works well)

ResidualCalculation % This is the main function where residuals are calculated

%% Graph-cuts

%%% Here we initialize the graph-cut parameters
imDataParams.images = profiles_lowres;
Nx = size(imDataParams.images,1);
Ny = size(imDataParams.images,2);
imDataParams.TR = tr*1e-3;
algoParams.noise_bias_correction = 1;
algoParams.species(1).name = 'water';
algoParams.species(1).frequency = [0];
algoParams.species(1).relAmps = [1];
algoParams.species(2).name = 'fat';
algoParams.species(2).frequency = freq;
algoParams.species(2).relAmps = amp;
res(isnan(res)) = 0;
algoParams.residual = permute(res, [3 1 2]);
algoParams.bssfp_flag = 1;
algoParams.size_clique = 1; %Size of MRF neighborhood (1 uses an 8-neighborhood, common in 2D)
algoParams.NUM_R2STARS = 10; %Number of R2* values for quantization
algoParams.range_fm = [min(b0range) max(b0range)]; %Range of field map values
algoParams.NUM_FMS = size(res,3); %Number of field map values to discretize
algoParams.NUM_ITERS = 100; %Number of graph cut iterations
algoParams.SUBSAMPLE = 1; %Spatial subsampling for field map estimation (for speed)
algoParams.DO_OT = 0; %0,1 flag to enable optimization transfer descent (final stage of field map estimation)
algoParams.LMAP_POWER = 2; %Spatially-varying regularization (2 gives ~ uniformn resolution)
algoParams.lambda = 1e-2; %Regularization parameter
algoParams.LMAP_EXTRA = 1e-2; %More smoothing for low-signal regions
algoParams.TRY_PERIODIC_RESIDUAL = 0; %Take advantage of periodic residual if uniform TEs (will change range_fm)
imDataParams.PrecessionIsClockwise = 1;
imDataParams.FieldStrength = field_strength;

fm = GraphcutbSSFP( imDataParams, algoParams, zeros(Nx, Ny) ); % Applying graph-cuts
fm = imresize(fm,scale_factor); % Rescaling the field map to the full resolution again
imagesc((fm).*intense_mask_nand,[min(b0range), max(b0range)]), colormap hot, colorbar
axis image
axis off 
title("Estimated \Delta B_0 Map", "FontName","Arial","FontSize",16);

% Chemical shift map is estimated as well, it might be wrapped around 1/TR,
% nevertheless it is a good first sanity check to see if fat and water are
% well-seperated
cs = ((angle(mean(profiles,3))))*(1e3/tr)/(pi);
cs = mod(cs, 1e3/tr);
figure
imagesc((mod(cs-fm, 1e3/tr)).*intense_mask_nand), axis image, colormap hot, colorbar
title("Chemical shift map", "FontName","Arial","FontSize",16);

%% Initial 2-compartment fixed relaxation fat-water seperation 

% Fixing water and fat relaxation times
t1w = 16*80;
t2w = 80;
t1f = 6*80;
t2f = 80;

d = [];
k = 1;


%%% Fitting to 2-component dictionaries 
clear FF_GC

Nx = size(profiles,1);
Ny = size(profiles,2);
for x = 1:Nx
    for y = 1:Ny

        p = (squeeze(profiles(x, y, :)));
        p = ApplyB0Shift(p, -1*fm(x,y), tr); % Here the profiles are corrected for inhomogeneity field

        %%% Simulating fat profile
        f = zeros(pcnum,1);
        
        for i = 1:length(freq)
            t = (bSSFPAnalytic(t1f, t2f, te, tr, fa, (0:pc_step:359)', freqhz(i)));
             f = f + amp(i)*t;
        end

        %%% Simulating water profile
        w = (bSSFPAnalytic(t1w, t2w, te, tr, fa, (0:pc_step:359)', 0));

        %%% "da" is the small dictionary
        da = [w, f];

        if (cs(x,y)-fm(x,y))<0 || angle(mean(p))*(1e3/tr)/(pi)<0
            da = -1*da; 
        end
            %%% Initial weights and FF is found via direct least squares
            %%% fitting to the small dictionary
            weights = pinv(reimconcat(da))*reimconcat(p);
            if weights(1)==0 && weights(2)==0
                FF_GC(x,y) = 0;
            else
                FF_GC(x,y) = abs(weights(2))/abs(sum((weights)));
            end
    end
end


FF_GC(FF_GC>1) = 1;
FF_GC(FF_GC<0) = 0;

figure
imagesc(FF_GC.*intense_mask_nand), axis image, axis off, colormap parula, colorbar
title("Initial PDFF Estimation","FontName","Arial","FontSize",16)



%% DEPENDENCY FUNCTION DEFINITIONS %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [profiles, phase_offset] = PhaseCorrection(profiles, dictionary)

D = dictionary;

M = pinv(real(D'*D));

corrected_profiles = [];
phase_offset = [];

for x = 1:size(profiles,1)
    for y = 1:size(profiles,2)
        p = squeeze(profiles(x,y,:));
        h = D'*p;
        c = h.'*M*h;
        b = angle(c);
        a = b;
        phase_offset(x,y) = a;
    end
%     x
end

unwrapped = phase_offset;
unwrapped(isnan(unwrapped))=0;

for x = 1:size(profiles,1)
    for y = 1:size(profiles,2)
        p = squeeze(profiles(x,y,:));
        a = (unwrapped(x,y))/2;
        cp = p*exp(-1i*a);
        corrected_profiles(x,y,:) = cp;
        phase_offset(x,y) = -a;
    end
end
profiles = corrected_profiles;

end


function signal = bSSFPAnalytic(t1,t2,te,tr,alpha,pc,off_res)

    e1 = exp(-tr/t1);
    e2 = exp(-tr/t2);

    alpha = alpha*pi/180;
    d = (1-e1*cos(alpha)-(e2^2)*(e1-cos(alpha)));
    b = e2*(1-e1)*(1+cos(alpha))/d;
    a = e2;
    M = (1-e1)*sin(alpha)/d;
    delt = pc*pi/180 ;
    theta0 = 2*pi*off_res*(tr*1e-3);
    theta = theta0-delt;
    phi = theta0*te/tr;
    signal = (M*(1-a*exp(1i*theta))./(1-b*cos(theta)))*exp(-1i*phi);
    signal = conj(signal)*exp(-te/t2);

end


function shifted_profile = ApplyB0Shift(profile, b0, tr)

    if size(profile, 1)<size(profile, 2)
        profile = profile.';
    end
    
    bw = 1e3/tr;
    b0 = b0/bw;
    Npc = length(profile);
    pc_step = 360/Npc;
    pc = 0:pc_step:359;
    
    [F, modenums] = FTMat(Npc-1, pc);
    
    v = @(x) exp((1i)*2*pi*(2*modenums+1)*x/2);
    shifted_profile = F'*diag(v(b0))*F*profile;
    end
    
function ff= Weights2FFFast(weights, b0, water_span)

    [Nx,Ny,~] = size(weights);
    b0arr = unique(b0);
    for x = 1:Nx
        for y = 1:Ny
            w = squeeze(weights(x,y,:));
            w = sum(reshape(abs(w),length(b0arr),length(b0)/length(b0arr)),2);
            ff(x,y) = ((b0arr<(max(b0arr)-water_span)) + (b0arr>(water_span)) - 1)'*abs(w) / sum(abs(w));
    
        end
    end

end



function [weights, bestFitSignal, regMatrix] = NNLS_Laplace(profile, lambda, lib, parameters)

    arrayX = unique(parameters(:,1));
    arrayY = unique(parameters(:,2));
    
    acquiredSignalWithReg = [profile; zeros(size(lib,2),1)];
    
    sz = length(arrayX)*length(arrayY);
    lgY = length(arrayY);
    lgX = length(arrayX);
    
    
    
    projection_matrix = zeros(lgX, lgX*lgY);
    for i = 1:lgX
        projection_matrix(i, ((1:lgX:(lgX*lgY-lgX))+(i-1))) = 1;
    end
    
    dfVector = ones(sz-1,1);
    dfVector(lgX:lgX:end) = 0;
    dfVectorWrap = zeros(lgX*lgY-lgX+1,1);
    dfVectorWrap(1:lgX:end) = 1;
    rtrVector = ones(sz-lgX,1);
    diagVector = -4*ones(sz,1);
    diagVector(1:1:lgX) = -3;
    diagVector(end-lgX:1:end) = -3;
    
    dfLaplacian = diag(diagVector,0)+diag(dfVector,1)+diag(dfVector,-1)+diag(dfVectorWrap,lgX-1)+diag(dfVectorWrap,-lgX+1);
    rtrLaplacian = diag(rtrVector,lgX)+diag(rtrVector,-lgX);
    
    regMatrix = lambda*(dfLaplacian+rtrLaplacian);
    dictionaryWithReg = [lib; regMatrix];
    
    
    weights = lsqnonneg(dictionaryWithReg,acquiredSignalWithReg);
    bestFitSignal = lib*weights; 

end






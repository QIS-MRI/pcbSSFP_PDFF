%% Function: graphCutIterations
%% Description: graph cut iterations for field map estimation
%% Author: Diego Hernando, June 3, 2009 (modified for bSSFP)

function fm = graphCutIterations(imDataParams,algoParams,residual,lmap,cur_ind )

DEBUG = 0;
SMOOTH_NOSIGNAL = 1;
STARTBIG = 1;
dkg = 15;
DISPLAY_ITER = 0;

gyro = 42.58;
deltaF = [0 ; gyro*(algoParams.species(2).frequency(:))*(imDataParams.FieldStrength)];
lambda = algoParams.lambda;
try
    if algoParams.bssfp_flag == 1
        dt = imDataParams.TR;
    end
catch
    te = imDataParams.TE;
    dt = te(2)-te(1);
end
period = 1/dt;
[sx,sy,N,C,num_acqs] = size(imDataParams.images);
fms = linspace(algoParams.range_fm(1),algoParams.range_fm(2),algoParams.NUM_FMS);
dfm = fms(2)-fms(1);
resoffset = [0:(sx*sy-1)]'*algoParams.NUM_FMS;
[masksignal,resLocalMinima,numMinimaPerVoxel] = findLocalMinima( residual, 0.06 );
numLocalMin = size(resLocalMinima,1);
stepoffset = [0:(sx*sy-1)]'*numLocalMin;
clear cur_ind2

fm = zeros(sx,sy);
for kg=1:algoParams.NUM_ITERS

  fmiters(:,:,kg) = fm;

  if kg == 1 & STARTBIG==1
    lambdamap = lambda*lmap;
    prob_bigJump = 1;
  elseif (kg==dkg & SMOOTH_NOSIGNAL==1) | STARTBIG == 0
    lambdamap = lambda*lmap;
    prob_bigJump = 0.5;
  end

  cur_sign = (-1)^(kg);

  if rand < prob_bigJump
    cur_ind2(1,:,:) = cur_ind;
    repCurInd = repmat(cur_ind2,[numLocalMin,1,1]);
    if cur_sign>0
      stepLocator = (repCurInd+20/dfm>=resLocalMinima) & (resLocalMinima>0);
      stepLocator = squeeze(sum(stepLocator,1))+1 ;
      validStep = masksignal>0 & stepLocator<=numMinimaPerVoxel;
    else
      stepLocator = (repCurInd-20/dfm>resLocalMinima) & (resLocalMinima>0);
      stepLocator = squeeze(sum(stepLocator,1));
      validStep = masksignal>0 & stepLocator>=1;
    end
    nextValue = zeros(sx,sy);
    nextValue(validStep) = resLocalMinima(stepoffset(validStep) + stepLocator(validStep));
    cur_step = zeros(sx,sy);
    cur_step(validStep) = nextValue(validStep) - cur_ind(validStep);

    if rand < 0.5
      nosignal_jump = cur_sign*round(abs(deltaF(2))/dfm);
    else
      nosignal_jump = cur_sign*abs(round((period - abs(deltaF(2)))/dfm));
    end
    cur_step(~validStep) = nosignal_jump;

  else
    all_jump = cur_sign*ceil(abs(randn*3));
    cur_step = all_jump*ones(sx,sy);
    nextValue = cur_ind + cur_step;

    if cur_sign>0
      cur_step(nextValue(:)>length(fms)) = length(fms) - cur_ind(nextValue(:)>length(fms));
    else
      cur_step(nextValue(:)<1) = 1 - cur_ind(nextValue(:)<1);
    end
  end

  [A] = createExpansionGraphVARPRO_fast( residual, dfm, lambdamap, algoParams.size_clique,cur_ind, cur_step);
  A(A<0)=0;

  [flowvalTS,cut_TS,RTS,FTS] = max_flow(A',size(A,1),1);

  cut1 = (cut_TS==-1);
  cut1b = 0*cut1 + 1;
  cut1b(end) = 0;

  if sum(sum(A(cut1b==1, cut1b==0))) <= sum(sum(A(cut1==1, cut1==0)))
    cur_indST = cur_ind;
  else
    cut = reshape(cut1(2:end-1)==0,sx,sy);
    cur_indST = cur_ind + cur_step.*(cut);
  end

  cur_ind = cur_indST;
  cur_ind(cur_ind<1) = 1;
  cur_ind(cur_ind>length(fms)) = length(fms);
  fm = fms(cur_ind);

  if DISPLAY_ITER == 1
    imagesc(fm,[-600 600]);axis off equal tight,colormap gray;colorbar;
    title(['Iteration: ' num2str(kg)],'FontSize',24);
    drawnow;
  end

end

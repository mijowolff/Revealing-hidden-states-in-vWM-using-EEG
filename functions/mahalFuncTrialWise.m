function d = mahalFuncTrialWise(dat,lab)

% calculate relative Mahalanobis of each trial between two condition
% specific centroids 
% dat [trials, features, samples]
% e.g. dat = rand(150,20,500); % random data, 150 trials, 20 channels, 500 time points
% lab is a vector of binary condition values [1 or 2]
% output d is the distance to own centroid, relative to other condition
% centroid
% written by Mark Stokes, 2014
%%

nTrials = size(dat,1);
nSamples = size(dat,3);

d = zeros(nTrials,nSamples);% variable for saving simulated data
trl_ind = 1:nTrials;
for trl=1:nTrials
%     fprintf(['Doing ' num2str(trl) ' out of ' num2str(nTrials) ' trials' '\n'])
    tmpDat = dat(setdiff(trl_ind,trl),:,:);
    tmpLab = lab(setdiff(trl_ind,trl));
    for n=1:nSamples
        sigma = pinv(pCov(squeeze(tmpDat(tmpLab==1,:,n)),squeeze(tmpDat(tmpLab==2,:,n))));        
        P1 = mean(tmpDat(tmpLab==1,:,n),1);
        P2 = dat(trl,:,n);
        d1 = sqrt((P1-P2)*sigma*(P1-P2)');        
        P1 = mean(tmpDat(tmpLab==2,:,n),1);
        d2 = sqrt((P1-P2)*sigma*(P1-P2)');
        if lab(trl) == 1
            d(trl,n) = d2 - d1;
        elseif lab(trl) == 2
            d(trl,n) = d1 - d2;
        end
    end
end
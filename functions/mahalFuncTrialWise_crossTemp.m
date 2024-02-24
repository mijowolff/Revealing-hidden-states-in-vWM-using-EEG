function d = mahalFuncTrialWise_crossTemp(dat,lab)

% calculate relative Mahalanobis of each trial between two condition
% specific centroids for all pissible time-point time-point combinations
% dat [trials, features, samples]
% e.g. dat = rand(150,20,500); % random data, 150 trials, 20 channels, 500 time points
% lab is a vector of binary condition values [1 or 2]
% output d is the distance to own centroid, relative to other condition
% centroid for each time-point combination
% written by Mark Stokes, 2014, adapted for cross-temporal analysis by
% Michael Wolff, 2015
%%

nTrials = size(dat,1);
nSamples = size(dat,3);

d = zeros(nTrials,nSamples,nSamples);% variable for saving simulated data
trl_ind = 1:nTrials;
for trl=1:nTrials
%     fprintf(['Doing ' num2str(trl) ' out of ' num2str(nTrials) ' trials' '\n'])
    tmpDat = dat(setdiff(trl_ind,trl),:,:);
    tmpLab = lab(setdiff(trl_ind,trl));
    for n1=1:nSamples
        sigma = pinv(pCov(squeeze(tmpDat(tmpLab==1,:,n1)),squeeze(tmpDat(tmpLab==2,:,n1))));
        P1_1 = mean(tmpDat(tmpLab==1,:,n1),1);
        P1_2 = mean(tmpDat(tmpLab==2,:,n1),1);
        for n2=1:nSamples
            P2 = dat(trl,:,n2);                    
            d1 = sqrt((P1_1-P2)*sigma*(P1_1-P2)');           
            d2 = sqrt((P1_2-P2)*sigma*(P1_2-P2)');
            if lab(trl) == 1
                d(trl,n2,n1) = d2 - d1;
            elseif lab(trl) == 2
                d(trl,n2,n1) = d1 - d2;
            end
        end
    end
end


%cohort_wbn is 3D matrix that tabulates the weighted brain networks
% of the whole cohort
[no_subjs N N]=size(cohort_wbn);

threshold = 0.2;% could be adjusted


hubs_list=zeros(no_subjs,N);
hubs_list_score=zeros(no_subjs,N);

 for su=1:no_subjs
     wbn=squeeze(cohort_wbn(su,:,:));
    [hubs_list(su,:) hubs_list_score(su,:)]=detect_hubs_bn(wbn,threshold);
 end
  
  
%% CREATE A CO-OCCURENCE MATRIX OF SIZE ROIS X ROIS THAT COUNTS
%% HOW MANY TIMES EACH PAIR OF ROIS IS ENCOUNTERED SIMULTANEOUSLY 
%% AS HUB ACROSS THE COHORT

co_occurence=zeros(N,N);

for su=1:no_subjs
    rr=find(hubs_list(su,:)==1);
    
    for k=1:length(rr)
        for l=(k+1):length(rr)
            co_occurence(rr(k),rr(l))=co_occurence(rr(k),rr(l))+1;
            co_occurence(rr(l),rr(k))=co_occurence(rr(k),rr(l));
        end
    end
end
    
no_iter=100;
thres=0.1;
% winners will return consistent hubs across a cohort
[winners replicators_evolution]=replicator_dynamics(co_occurence,no_iter,thres);

%%get consistent hubs across the cohort
consistent_hubs=find(winners==1);




%%% add the pathway to your BCT software
addpath(genpath('C:\Users\mpnsd\Desktop\SOFTWARE\BCT\2019_03_03_BCT'))


%% run an example per subject

threshold = 0.1;% could be adjusted
no_subjs=160;% for example
N=90; % AAL

hubs_list=zeros(no_subjs,N);

 for su=1:no_subjs
    hubs_list(su,:)=detect_hubs_bn(wbn,threshold);
 end
 

%% CREATE A CO-OCCURENCE MATRIX OF SIZE ROIS X ROIS THAT COUNTS
%% HOW MANY TIMES EACH PAIR OF ROIS IS ENCOUNTERED AS HUB

co_occurence=zeros(N,N);

for k=1:N
    r1=find(hubs_list(:,k)=1);
    for l=1:N
        r2=find(hubs_list(:,l)=1);
        r3=intersect(r1,r2);
        co_occurence(k,l)=length(r3)/no_subjs;
        co_occurence(k,l)=co_occurence(l,k);
    end
end

no_iter=100;
thres=0.1;
% winners will return consistent hubs across a cohort

[winners replicators_evolution]=replicator_dynamics(co_occurence,no_iter,thres);








function hubs_list=detect_hubs_bn(wbn,threshold)


%% detect hubs from a brain network
% Below I summarize a simple pipeline:
% 
% 1) estimate nodal BC, nodal strength, nodal local efficiency, and nodal weighted clustering coefficient.
% 2) rank the 90s vectors of nodal BC and strength from the highest to the lowest
% 3) rank the 90s vectors of nodal local efficiency and weighted clustering coefficient from the lowest to the highest
% A node is a global hub when it is more global than local!!
% 
% 4) we can set a threshold of 20% in the 4 90s vectors refer to the highest or the lowest values
% 5) Scoring every node with values from 0 up to 4 determined by the total number of hub criteria fulfilled.
% 
% for example, if an ROI in DMN is in the 20% of the highest score in BC and strength and in the 20% of the lowest
%  values in LE and wCC than its score will be 4!
% 
% 6) Finally, we will get a score per ROI independently per participant 
%    We will prefer a ROI to be at least in one top K % (global:eithe BC or strength) and in at
%    least on lowest K % (local: either CC or LE)
%
% 7)We can further detect consistent hubs across cohort based on
%   replicated dynamics and graph theory

%% INPUT : wbn       = a weighed structural or functional brain network
%%       : threshold = set a threshold like 0.2 (20%)

%% OUPUT : hubs_list = a vector equals the size of ROIs with 1s (HUBS) 
%%            AND 0s (NON-HUBS)


%Dimitriadis Stavros 2020
%http://users.auth.gr/~stdimitr/index.html

%Please cited this m-file as:
%Dimitriadis SI, Laskaris NA, Tsirka V, Vourkas V, Micheloyannis S, Fotopoulos S. 
%Tracking brain dynamics via time-dependent network analysis. 
%Journal of Neuroscience Methods Volume 193, Issue 1, 30 October 2010,
%Pages 145-155

% % The  betweenness centrality and global efficiency is computed using an auxiliary connection-length
%   matrix L, defined as L_ij = 1/W_ij - 1 for all nonzero L_ij.
%This has an intuitive interpretation, as higher connection weights intuitively
%   correspond to shorter lengths.

%%NUMBER OF NODES SELECTED FROM THE ORDERED LISTS
[N,N]=size(wbn);
no=round(N*threshold);


%get the nonzeros values
A = wbn > 0;                                                  % 
% transformed weights to length
G=zeros(N,N);
G=wbn;
G(A)=1./G(A)-1;

%1.ESTIMATE NODAL BC
BC=betweenness_wei(G);

% NORMALIZED NODAL BC
nBC=zeros(1,N);
nBC=BC/((N-1)*(N-2));

%% RANKED nBC
[values,nBCsorted]=sort(nBC,'descend');

nBCsorted_sel=nBCsorted(1:no);

%2. ESTIMATE NODAL STRENGTH
nstr=zeros(1,N);
nstr=sum(wbn);

%% RANKED NODAL STRENGTH
[values,nstrsorted]=sort(nstr,'descend');

nSTRsorted_sel=nstrsorted(1:no);

% 3. ESTIMATE NODAL LOCAL EFFICIENCY
nle=zeros(1,N);
nle=efficiency_wei(wbn,2); % 2 is recommended instead of 1

%% RANKED NODAL STRENGTH
[values,nlesorted]=sort(nle);

nLEsorted_sel=nlesorted(1:no);


% 4. ESTIMATE NODAL CLUSTERING COEFFICIENT
ncc=zeros(1,N);
ncc=clustering_coef_wu(wbn);


%% RANKED NODAL STRENGTH
[values,nccsorted]=sort(ncc);

nCCsorted_sel=nccsorted(1:no);


%% FROM THE FOUR VECTORS, WE WILL 
%% nBCsorted_sel
%% nstrsorted_sel
%% nlesorted_sel
%% nccsorted_sel

vec1=zeros(1,N);
vec2=zeros(1,N);
vec3=zeros(1,N);
vec4=zeros(1,N);

vec1(nBCsorted_sel)=1;
vec2(nSTRsorted_sel)=1;
vec3(nLEsorted_sel)=1;
vec4(nCCsorted_sel)=1;

%% A NODE COULD BE A HUB IF AT LEAST IS WITHIN THE SORTED
%% LIST OF ONE GLOBAL AND ONE LOCAL FEATURE
vec12=sum([vec1; vec2]);
vec34=sum([vec3; vec4]);

%% THE UNION OF THE TWO LISTS DEFINE THE HUB LIST
r1=find(vec12==2);
r2=find(vec34==2);

% THE RETUREDN HUBS LIST WITH 1s (HUBS) AND 0s (NON-HUBS)
hubs_list=zeros(1,N);
hubs_list(union(r1,r2))=1;







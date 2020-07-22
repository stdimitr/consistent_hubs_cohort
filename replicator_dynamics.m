function [winners replicators_evolution]=replicator_dynamics(w,no_iter,threshold)

%implementation of replicator dynamics according to
%Neumann et al. 2005 "Meta-analysis of Functional Imaging data using replicator dynamics"

%INPUT: w = co-occurence matrix - (co-occurence matrix = z*z' where z is the binary representation - see table 1)
% no_iter = # of iterations
%threshold= define a threshold to distinquish winners for the rest (e.g. 0.1)
%OUTPUT : winners = binary vector (size = 1 x # of replicators) which defines the winner replicators 
%         replicators_evolution = show the evolution of the replicators across iteration 

%Dimitriadis Stavros 2009
%http://users.auth.gr/~stdimitr/index.html

%Please cited this m-file as:
%Dimitriadis SI, Laskaris NA, Tsirka V, Vourkas V, Micheloyannis S, Fotopoulos S. 
%Tracking brain dynamics via time-dependent network analysis. 
%Journal of Neuroscience Methods Volume 193, Issue 1, 30 October 2010,
%Pages 145-155

[x y]=size(w);

%initialize vector
X_new(1:no_iter+1,1:x)=0;
%initialization for the first iteration
X_new(1,:)=1/x;



for i=2:no_iter+1
    for j=1:x
        X_new(i,j)=X_new(i-1,j)*(sum(w(j,:).*X_new(i-1,:))/(X_new(i-1,:)*w*X_new(i-1,:)')); 
    end
end


winners=X_new(no_iter+1,:)> threshold;
replicators_evolution=X_new;





%co-occurence-matrix= z'*z where z is Table 1.
%Co-occurence matrix of TABLE 1 . P.168 
w=[ 0 6 4 1 1 0 ; 6 0 2 2 0 1; 4 2 0 2 0 0 ; 1 2 2 0 2 1;1 0 0 2 0 0 ; 0 1 0 1 0 0];

%calling replicator_dynamics - winners are the 2 first foci in accordance
%with the paper
[winners replicators_evolution]=replicator_dynamics(w,100,0.1);


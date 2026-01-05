
function [F,L]=normalize(F0,L0)

% F0 is the T by r matrix of initial factors
% L0 is the N by r matrix of initial factor loadings
% F is the T by r matrix of normalized factors 
% L is the N by r matrix of normalized factor loadings

[T,r]=size(F0);
N=size(L0,1);

[U,V,W]=svd(F0*L0');
F=U(:,1:r)*sqrt(T);
L=(V(1:r,:)*W')'/sqrt(T);

end



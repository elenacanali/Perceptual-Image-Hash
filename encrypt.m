function F_E=encrypt(F, K)
% input: secret seed K and intermediate hash F
%% Produce vector of key with logistic map
% size(F{3}, 1) is the cardinality of the feature point set extract by
% SIFT algorithm
m=size(F{3}, 1);
n_blocks=size(F{2},2);
L=m + n_blocks + 2*m + 2;
k=zeros(1, L);
%k(1)=K(1);
k(1) = K;
for i=2:L %logistic map
    k(i)=k(i-1)-k(i-1)^2;
end

%% Encryption: scalar multiplication 
tmp=diag(k(1:m));
F_E{1}=F{1}*tmp;
%tmp=diag(k(m+1: m + K(2)));
tmp=diag(k(m+1: m + n_blocks));
F_E{2}=F{2}*tmp;
F_E{3}=zeros(m,2);
%tmp=diag(k( m+K(2)+1 : 2 : L-3 ));
tmp=diag(k(m+n_blocks+1 : 2 : L-3 ));
F_E{3}(:,1)=tmp*F{3}(:, 1); % x coordinates
%tmp=diag(k( m+K(2)+2 : 2 : L-2 ));
tmp=diag(k( m+n_blocks+2 : 2 : L-2 ));
F_E{3}(:,2)=tmp*F{3}(:,2); % y coordinates
F_E{4}=F{4}*k(L-1);
F_E{5}=F{5}*k(L);


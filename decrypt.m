function F_D=decrypt(F_E, K)
% input: secret seed K and encrypted hash F_E
%% Produce vector of key with logistic map
% size(F{3}, 1) is the cardinality of the feature point sent extract by
% SIFT algorithm
n_blocks=size(F_E{2},2);
m=size(F_E{3}, 1);
L=m + n_blocks + 2*m + 2;
k=zeros(1, L);
%k(1)=K(1);
k(1) = K;
for i=2:L %logistic map
    k(i)=k(i-1)-k(i-1)^2;
end
k=1./k;

%% Decryption
tmp=diag(k(1:m));
F_D{1}=F_E{1}*tmp;
%tmp=diag(k(m+1: m + K(2)));
tmp=diag(k(m+1: m + n_blocks));
F_D{2}=F_E{2}*tmp;
F_D{3}=zeros(m,2);
%tmp=diag(k( m+K(2)+1 : 2 : L-3 ));
tmp=diag(k(m+n_blocks+1 : 2 : L-3 ));
F_D{3}(:,1)=tmp*F_E{3}(:, 1); % x coordinates
%tmp=diag(k( m+K(2)+2 : 2 : L-2 ));
tmp=diag(k( m+n_blocks+2 : 2 : L-2 ));
F_D{3}(:,2)=tmp*F_E{3}(:,2); % y coordinates
F_D{4}=F_E{4}*k(L-1);
F_D{5}=F_E{5}*k(L);










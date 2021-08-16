function [pcs,explained] = myPCA(X)

[V, E] = eig(cov(X)); %[eigvec,eigval]
[E, S] = sort(diag(E),'descend');
pcs = V(:,S);
explained = (E / sum(E))*100;

end
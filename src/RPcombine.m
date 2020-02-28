function f=RPcombine(G)%,fileName)%(B,W)%between(I,J), within(I,J)
%G=[B(2:end,:);W];
%lg=length(G);
lg=max(max(G));
Gsp=sparse(G(:,1),G(:,2),1,lg,lg);
Gsp=Gsp+Gsp';
Gsp(Gsp>1)=1;
Gsp=Gsp-diag(diag(Gsp));
%save(fileName,'G')
f=Gsp;
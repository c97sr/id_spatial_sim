function [f,g]=DHCMall(G,N)
n=length(G);
%
measureA=cell(N,1); measureB=measureA;
%[C1,C2,C]=clust_coeff(G);
measureA{1}=G; measureB{1}=G;
A=G; B=G;
X=zeros(n,N);
%Standard clustering:
for i=1:n
    vi=A(i,:);
    vfind=find(vi);
    lv=length(vfind);
    Wi=B(vfind,vfind);
    X(i,1)=sum(sum(Wi))/lv/(lv-1);%length(find(Wi))/lv/(lv-1);
end
%N-clustering:
for m=2:N
    nextA=A*G;
    nextA(speye(n)==1)=0;
    nextA(B==1)=0;
    nextA(nextA>1)=1;
    nextB=B+nextA;
    nextB(speye(n)==1)=0;
    nextB(nextB>1)=1;
    measureA{m}=nextA; measureB{m}=nextB;
    A=nextA; B=nextB;
    for i=1:n
        vi=B(i,:);%Ax if exactly N steps from i, Bx if up to N steps
        vfind=find(vi);
        lv=length(vfind);
        Wi=B(vfind,vfind);%B1 for neighbours, Bx for up to N steps between, BNm1 for N-1
        X(i,m)=sum(sum(Wi))/lv/(lv-1);%length(find(Wi))/lv/(lv-1);
    end
end
X=X';
f=[nanmean(X,2),nanvar(X,0,2),min(X,[],2),prctile(X,[5,25,50,75,95],2),max(X,[],2)];
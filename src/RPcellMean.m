function f=RPcellMean(c1,c2,c3)
%ncells=3;
%lc=length(c1);
[a,b]=size(c1{1});
x=zeros(a,b,ncells);
cout=cell(lc,1);
for i=1:lc
    x(:,:,1)=c1{i};
    x(:,:,2)=c2{i};
    x(:,:,3)=c3{i};
    cout{i}=nanmean(x,3);
end
f=cout;
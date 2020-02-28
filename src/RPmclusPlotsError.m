function f=RPmclusPlotsError(X)
%Input is RPcellMean(cCell)
%numClus=1;
clusmax=4;
la=length(X);
m1=zeros(clusmax,la);
M1=zeros(clusmax,la,2);
for i=1:la
    Xi=X{i};
    Xi(:,5)=Xi(:,6)-Xi(:,5);
    Xi(:,7)=Xi(:,7)-Xi(:,6);
    M1(:,i,:)=Xi(:,[5,7]);
    m1(:,i)=Xi(:,6);
end
x=repmat((0:6),clusmax,1);

cmap=lines(clusmax);
%cmap=[.2,.2,.2;1,.2,.2;.2,1,.2;.2,.2,1];
%cmap=.8*[0,0,0;1,0,0;0,1,0;0,0,1];
colors=cell(clusmax,1);
for i=1:clusmax
    colors{i}=cmap(i,:);
end
lineProps.col=colors;
lineProps.style='-';
%
fs=12; lw=2;
figure
%lineProps=[];
%lineProps.col=lines(clusmax);
hold on
%
x1=0:6;
%
mseb(x,m1,M1,lineProps,1);
h1=scatter(x1,m1(1,:),[],cmap(1,:),'o','linewidth',lw);
h2=scatter(x1,m1(2,:),[],cmap(2,:),'x','linewidth',lw);
h3=scatter(x1,m1(3,:),[],cmap(3,:),'^','linewidth',lw);
h4=scatter(x1,m1(4,:),[],cmap(4,:),'s','linewidth',lw);
hold off
set(gca,'fontsize',fs)
axis ([0,6,0,.7])%.35])
xlabel('Distance power \alpha')
ylabel(strcat('CC^m (nodes)'))
legend([h1,h2,h3,h4],'m=1','m=2','m=3','m=4','location','NW')
%legend([h1,h2],'m=1','m=2','location','NW')
grid on
grid minor
box on
%}
%{
figure
mseb(x2,h1',H1,lineProps,1);
set(gca,'fontsize',12)
axis ([0,6,0,.35])
xlabel('Distance power \alpha')
ylabel(strcat('CC^m (households)'))
legend('m=1','m=2','location','NW')
grid on
grid minor
box on
%}
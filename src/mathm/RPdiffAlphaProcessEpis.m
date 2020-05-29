function [f,g,h]=RPdiffAlphaProcessEpis
a=(0:1:6);%Alpha values
n=10000; lruns=10;
la=length(a);
P=zeros(la,lruns); Z=P;%5=max poss runs;
epi=cell(la,1);
%
lepi=300;%(Max) length of incidence vector
epiz=zeros(lepi,lruns);
[epi{:}]=deal(epiz);
fpath='..\..\..\Downloads\SRebola\id_spatial_sim-master\scenarios\ebola\output\epiDHh4w30pp2333r2p2_alpha';
for i=1:la
    ai=a(i);
    %events=importdata(strcat('DH_h4w120pp05_alpha',num2str(ai),'_pset_0_Events.out'));%,' ',1,0);
    events=importdata(fullfile(strcat(fpath,num2str(ai),'_pset_0_Events.out')));
    %DH1: DH4_30_p3_2_alpha
    %DH2: DH_h4w100pp08_alpha
    %DH_h4w120pp05_alpha
    %
    run=events.data(:,1);%events(2:end,1);
    day=events.data(:,2);%events(2:end,2);
    event=events.data(:,3);%events(3:10:end);%(2:end,3);
    %index=events(2:end,4);
    runs=unique(run);
    lruns=length(runs);
    for j=1:lruns
        runj=runs(j);
        d=day(run==runj)+1;
        ev=event(run==runj);
        lev=length(ev); zev=zeros(lev,1);
        zev(ev==0)=1;
        %zev=ev; zev(zev>0)=1;
        x=accumarray(d,zev);
        [peak,~]=max(x);
        z=sum(zev);
        P(i,j)=peak/n; Z(i,j)=z/n;
        %For full epidemics:
        epi{i}(1:length(x),j)=x;
    end
    zevAll=event; zevAll(zevAll>0)=1; dayAll=day+1;
    x=accumarray(dayAll,zevAll);
    epi{i}=x/n/lruns;
    clear events zevAll dayAll
end
f=P; g=Z; h=epi;
save('allEpis.mat','P','Z','epi')
function f=RPdiffAlphaProcessNets
a=(0:6);%Alpha values
%ai=0;
%w=[60,100,140,180,220];
%lw=length(w);
n=50000;%
la=length(a);
reps=1;%Num of nets for each alpha;
Gcell=cell(la,reps);%XX
fpath='C:\Users\dhaw\Documents\GitHub\id_spatial_sim\subExp\ebola\network\netDHflath3w200pp04A_alpha';
for i=1:la
    ai=a(i);%XX
    for j=1:reps
    %
        fname1=strcat(fpath,num2str(ai),'_arcs.out');
        fname2=strcat(fpath,num2str(ai),'_nodes.out');
        arcs=load(fullfile(fname1)); %Households
        arcs=arcs(:,1);%1st column
        nodes=load(fullfile(fname2));
        nodes=nodes(:,5);%5th column
        %DH1: net_h3w32pp25_alpha
        %DH2: netDH4_100_p08_2_alpha
        %DH3: netDH4_120_p05_2_alpha
        %DH4: netDH4_20_p35_2_alpha (h=3)
        Gcell{i,j}=RPmakeNets(arcs,nodes,n);
    end
    %
    %{
    fname1a=strcat(fpath,num2str(ai),'w',num2str(wi),'a_arcs.out');
    fname2a=strcat(fpath,num2str(ai),'w',num2str(wi),'a_nodes.out');
    fname1b=strcat(fpath,num2str(ai),'w',num2str(wi),'b_arcs.out');
    fname2b=strcat(fpath,num2str(ai),'w',num2str(wi),'b_nodes.out');
    fname1c=strcat(fpath,num2str(ai),'w',num2str(wi),'c_arcs.out');
    fname2c=strcat(fpath,num2str(ai),'w',num2str(wi),'c_nodes.out');
    arcs=load(fullfile(fname1a)); arcs=arcs(:,1); nodes=load(fullfile(fname2a)); nodes=nodes(:,5);
    Gcell{i,1}=RPmakeNets(arcs,nodes,n);
    arcs=load(fullfile(fname1b)); arcs=arcs(:,1); nodes=load(fullfile(fname2b)); nodes=nodes(:,5);
    Gcell{i,2}=RPmakeNets(arcs,nodes,n);
    arcs=load(fullfile(fname1c)); arcs=arcs(:,1); nodes=load(fullfile(fname2c)); nodes=nodes(:,5);
    Gcell{i,3}=RPmakeNets(arcs,nodes,n);
    %}
end
f=Gcell;
%save('allNetsDH4.mat','Gcell','-v7.3')%h_w_pw_R0
save('DH2netsA.mat','Gcell','-v7.3')%h_w_pw_R0
function f=DHCMsamplingSpeed(G,mmax,sampleNumber)%Or input v
%G must be simple: undirected with no self loops
n=length(G);
v=randsample(n,sampleNumber);%v is indices
%%
lv=length(v);
clusterVec1=zeros(lv,mmax);

%Gnbrs=cell{n,1};
[aye,jay]=find(G);
GG=accumarray(aye,jay,[size(G,1),1],@(x) {sort(x).'});

%%
imvw1=v';%node
imvw2=ones(lv,1);%info order
Cv=cell(n,1);%{};%Cell array of data for each node
lengthC=lv;

%Standard clustering:
for i=1:lv
    %Find neighbours of node i:
    vx=v(i);
    v1=mNbr(GG,vx,1,vx);%Indices of neighbours
    %Update stored data:
    Cv{i}=v1;
    %
    Gred=G(v1,v1);%Adj mat of neighbourhood
    lv1=length(v1);%Number of neighbours
    C1=sum(sum(Gred))/lv1/(lv1-1);%How many links exist/how many could exist
    clusterVec1(i,1)=C1;
end
%m-clustering, m>1:
if mmax>1
for m=2:mmax
    for i=1:lv
        %m-neighbours of v(i):
        vi=v(i);
        vx=Cv{i};%m-1 OR m (if already xcalculated in loop) nbrs of v(i)
        imvw2i=imvw2(i);
        if imvw2i<m%i.e.==m-1
            vm=mNbr(GG,vx,m-imvw2i,vi);%m-neighbours of i
            %vm=vm(vm~=vi) here?
            imvw2(i)=m;
            Cv{i}=vm;
        else
            vm=vx;
        end
        vmx=vm;%(vm~=vi);%Exclude node i
        lvm=length(vmx);
        %m-clustering:
        %
        links1=0;
        %m-neighbours of each m-neighbour:
        for j=1:lvm%parfor
            vmj=vmx(j);%vmOrder(j);
            %if vmj~=vi
            %Faster to order or to find all links?
            if ismember(vmj,imvw1)%If have already calculated some nbhds
                index=find(imvw1==vmj);%How many
                thisFar=imvw2(index);
                vmjSoFar=Cv{index};
                if thisFar<m
                    mNj=mNbr(GG,vmjSoFar,m-thisFar,vmj);
                    imvw2(index)=m;
                    Cv{index}=mNj;
                else
                    mNj=Cv{index}; %mNjOnly=Cw{index};
                end
            else
                mNj=mNbr(GG,vmj,m,vmj);
                lengthC=lengthC+1;
                index=lengthC;
                imvw1(index)=vmj;
                imvw2(index)=m;
                Cv{index}=mNj;
            end
            %Other method: find subgraph
            %mNj=mNj(~vi); mNjOnly=mNjOnly(~vi);
            mNjx=mNj;%(mNj~=vmj);%m-neighbours of node j, not including j
            links1=links1+length(intersect(vmx,mNjx));
            %end
        end
        clusterVec1(i,m)=links1/lvm/(lvm-1);
    end
end
%****Save output****
%f=clusterVec1;
end
X=clusterVec1';
f=[nanmean(X,2),nanvar(X,0,2),min(X,[],2),prctile(X,[5,25,50,75,95],2),max(X,[],2)];
end

function f=mNbr(GG,vx1,mMore,vi)
    vx=vx1;
    G1cell=GG(vx);
    v1=unique([G1cell{:}]);%Indices of next neighbours
    
    v1=union(vx,v1);%up-to-m-neighbours of i
    
    v1=v1(v1~=vi);%Exclude original node
    vx=v1;%In case want wx etc.
    if mMore>1
        k=2;
        while k<=mMore
            Gmcell=GG(vx);
            vm=unique([Gmcell{:}]);
            vm=union(vx,vm);%up-to-m-neighbours of i
            vm=vm(vm~=vi);
            vx=vm;
            k=k+1;
        end
    end
    f=vx;
end
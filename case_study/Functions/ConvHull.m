function  SetConvHull=ConvHull(P)
para.n=size(P,2);
NumVertex_P=size(P,1);
SetConvHull=[];

yalmip('clear')
M=sdpvar(para.n,1);
Input=sdpvar(1,para.n);
z=sdpvar(1,1);
indicator=sdpvar(NumVertex_P,1);
cnsmin=[];
cnsmax=[];

cnsmin=[cnsmin,P*M-Input*M<=-1.0000e-03*ones(NumVertex_P,1)+indicator*z];
cnsmax=[cnsmax,P*M-Input*M>=1.0000e-03*ones(NumVertex_P,1)-indicator*z];
cnsmin=[cnsmin,z>=0];
cnsmax=[cnsmax,z>=0];
ops = sdpsettings('verbose',0);
OPTmin=optimizer(cnsmin,z,ops,{Input,indicator},M);
OPTmax=optimizer(cnsmax,z,ops,{Input,indicator},M);

parfor i=1:NumVertex_P
    indicator=zeros(NumVertex_P,1);
    indicator(i)=1;
    Mmin(i,:)=OPTmin({P(i,:),indicator});
    Mmax(i,:)=OPTmax({P(i,:),indicator});
end
evalmin=any(isnan(Mmin),2); 
evalmax=any(isnan(Mmax),2);
idxmin=find(evalmin<=0.5);
idxmax=find(evalmax<=0.5);
SetConvHull=[SetConvHull;P(idxmin,:);P(idxmax,:)];
SetConvHull=unique(SetConvHull,'rows');
end
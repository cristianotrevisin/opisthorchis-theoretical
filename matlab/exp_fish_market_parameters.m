clc
clearvars
close all


load ../dataOCN/subcatchments

COLS = ["#ff595e", "#ff924c", "#ffca3a", "#8ac926", "#1982c4", "#6a4c93"];

NX = ceil(max(X));
NY = ceil(max(Y));

CCR = randi(25,[1 max(CTC)]);

CTD = zeros(1,length(CTC));
for i = 1:length(CTD)
    CTD(i) = CCR(CTC(i));
end

CTM = zeros(NX,NY); MSC = zeros(NX,NY);
XI = ceil(X);
YI = ceil(Y);

for i = 1:length(X)
    CTM(XI(i),YI(i)) = CTD(i);
    MSC(XI(i),YI(i)) = CTC(i);
end


rng(2021)


load ../dataOCN/FD
load ../dataOCN/SC
FD.nNodes = length(CTD);
FD.outlet = 1921;
FD.A = FD_A; FD.X = FD_X; FD.Y = FD_Y; FD.downNode = FD_downNode;
% AvailableNodes = setdiff(1:OCN$FD$nNodes,OCN$FD$outlet)
thrA = 120; cellsize = 1;



CTA = zeros(NX,NY);

load ../dataOCN/AG
for i = 1:length(X)
    CTA(XI(i),YI(i)) = A(CTC(i));
end

% population is sorted randomly
ranks = randperm(length(A));
% Populations are Zipf distributed
P = zipf(ranks, 2, 250);

% Fish population scaling 0.3
Nf = (A*1000).^0.3;
hf = ones(SC.nNodes,1)*0.001;


%% HYDROLOGICAL CONNECTIVITY
W = zeros(SC.nNodes,SC.nNodes);
W2 = zeros(SC.nNodes,SC.nNodes);
for nn = 1:SC.nNodes
    temp = find(downNode==nn);
    W(nn,temp) = 1;
    W2(nn,temp) = 1;
    W2(temp,nn) = 1;
end

%% EPIMODEL
% define parameters
par.beta_FH = 4.898e-5; 
par.beta_HS = 9.160e-11;
par.beta_SF = 3.477e-5;

par.mu_H = 1/365/10;
par.mu_S = 1/365;
par.mu_F = 1/365/2.5;

par.mS = 0;
par.mF = 0;

par.phi = 1;

%define setup
setup.nNodes = SC.nNodes;
setup.Cf = hf;
setup.W = W;
setup.W2 = W2;
setup.H = P';
setup.Nf = Nf;
setup.Ns = 20000*ones(SC.nNodes,1);

%define equivalent implicit network
setup_I.H = sum(P);
setup_I.Nf = sum(Nf);
setup_I.Ns = sum(setup.Ns);
setup_I.M = 1;
setup_I.Cf = 1;
setup_I.W = 0;
setup_I.W2 = 0;
setup_I.nNodes = 1;

Time = 1:100*365;

%initial solution
y0=zeros(3,SC.nNodes);
%[a, b, c] = EE_OPI(par,setup,3);
y0(1:3,43) = [1/P(43); 0; 0];
%y0(1,253)=0; %initial infected in node 3

%initial solution for implicit
y0_I=zeros(3,1);
%[a, b, c] = EE_OPI(par,setup,3);
y0_I(1:3,1) = [1/sum(P); 0; 0];
%y0(1,253)=0; %initial infected in node 3



DST = zeros(SC.nNodes,SC.nNodes);
MD = zeros(SC.nNodes,SC.nNodes);
for sc1 = 1:SC.nNodes
    for sc2 = 1:SC.nNodes
        DST(sc1,sc2) = sqrt((SCX(sc1)-SCX(sc2))^2+(SCY(sc1)-SCY(sc2))^2);
        if sc1~=sc2
            MD(sc1,sc2) = P(sc1)*P(sc2)/(DST(sc1,sc2)^par.phi);
        end
    end
end

chis = 0:10;

for loop_over_chi = 1:length(chis);

    LOC = chis(loop_over_chi)*0.1*ones(1,SC.nNodes);
    

    M = MD./(sum(MD,1)).*(1-LOC);
    for sc = 1:SC.nNodes; M(sc,sc) = LOC(sc); end
    
    setup.M = M;    

    y = model(Time,par,setup,y0);
    y_I = model(Time,par,setup_I,y0_I);

    WH(:,:,loop_over_chi) = y(:,1:3:end);
    SI(:,:,loop_over_chi) = y(:,2:3:end);
    FI(:,:,loop_over_chi) = y(:,3:3:end);

end
%%
figure()
%plot(Time/365,y_I(:,1))
hold on
for plt = 1:length(chis);
    TMP = squeeze(WH(:,:,plt));
plot(Time/365, TMP*P'/(sum(P)))
text(Time(end)/365, TMP(end,:)*P'/(sum(P)),num2str(0.1*chis(plt)))
end

xlabel('Year')
ylabel('W_H')





function M = assignSC(MSC,X)
    M = zeros(size(MSC));
    for iz = 1:length(X)
        M(MSC==iz) = X(iz);
    end
end


function P = zipf(rank, expn, minP)
%     ranks = 1:1:N;
%     pmf = (ranks.^(-expn))/sum(ranks.^(-expn));
%     samples = rand(1,M);
%     p = cumsum(pmf(:));
%     [~,x] = histc(samples,[0;p/p(end)])

    H = sum(1./(rank.^expn));
    
    
    P = 1./rank.^expn./H;
    P = P/min(P)*minP;
    
end
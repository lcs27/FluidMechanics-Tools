clear,clc,close all;
CArl_Dyn_assemble_matrix

st=5000*dtB;

outer=st/dtA;
inner=dtA/dtB;

UA=zeros(N,1);
UB=zeros(N+1,1);
UdA=zeros(N,1);
UdB=zeros(N+1,1);
UddA=zeros(N,1);
UddB=zeros(N+1,1);

prev_UA=zeros(N,1);
UA_inter=zeros(N,1);

UA_free=zeros(N,1);
UB_free=zeros(N+1,1);
UdA_free=zeros(N,1);
UdB_free=zeros(N+1,1);
UddA_free=zeros(N,1);
UddB_free=zeros(N+1,1);

UA_link=zeros(N,1);
UB_link=zeros(N+1,1);
UdA_link=zeros(N,1);
UdB_link=zeros(N+1,1);
UddA_link=zeros(N,1);
UddB_link=zeros(N+1,1);

pUA=zeros(N,1);
pUB=zeros(N+1,1);
pUdA=zeros(N,1);
pUdB=zeros(N+1,1);

RHSA_free=FA-KA*pUA;
RHSB_free=FB-KB*pUB;

uAs=zeros(N,outer);
uBs=zeros(N+1,inner*outer);
residual=zeros(outer,1);
for jj=1:outer
    RHSA_free=FA-KA*pUA;
    UddA_free=MA\RHSA_free;

    UdA_free=pUdA+gammaA*dtA*UddA_free;
    UA_free=pUA+betaA*dtA*dtA*UddA_free;

    for kk=1:inner

        UddB_free=MB\RHSB_free;

        UdB_free=pUdB+gammaB*dtB*UddB_free;
        UB_free=pUB+betaB*dtB*dtB*UddB_free;

        UA_inter=(1-kk/inner)*prev_UA+kk/inner*UA_free;
        RHS_coupling=CA*UA_inter-CB*UB_free;

        lambda=H\RHS_coupling;

        UddB_link=MB\(transpose(CB)*lambda);
        UdB_link=gammaB*dtB*UddB_link;
        UB_link=betaB*dtB*dtB*UddB_link;

        UB=UB_link+UB_free;
        UdB=UdB_link+UdB_free;
        UddB=UddB_link+UddB_free;

        pUB=UB+dtB*UdB+dtB*dtB*(1/2-betaB)*UddB;
        pUdB=UdB+(1-gammaB)*dtB*UddB;

        RHSB_free=FB-KB*pUB;

        uBs(:,(jj-1)*inner+(kk-1)+1)=UB;
    end

    UddA_link=MA\(-transpose(CA)*lambda);
    UdA_link=gammaA*dtA*UddA_link;
    UA_link=betaA*dtA*dtA*UddA_link;

    UA=UA_link+UA_free;
    UdA=UdA_link+UdA_free;
    UddA=UddA_link+UddA_free;

    pUA=UA+dtA*UdA+dtA*dtA*(1/2-betaA)*UddA;
    pUdA=UdA+(1-gammaA)*dtA*UddA;
    prev_UA=UA;

    RHSA_free=FA-KA*pUA;

    uAs(:,jj)=UA;

    residual(jj)=norm(CA*UA-CB*UB)/norm(CA*UA);

end
figure
hold on
plot(dtA*(1:outer),uAs(N/7,:),'r')
plot(dtA*(1:outer),uAs(3*N/7,:),'r')
plot(dtA*(1:outer),uAs(5*N/7,:),'r')
plot(dtA*(1:outer),uAs(6*N/7,:),'r--','LineWidth',3)
plot(dtA*(1:outer),uAs(6*N/7+4,:),'r*','LineWidth',3)
plot(dtA*(1:outer),uAs(N,:),'r','LineWidth',3)
plot(dtB*(1:(inner*outer)),uBs(1,:),'b--','LineWidth',3)
plot(dtB*(1:(inner*outer)),uBs(5,:),'b*','LineWidth',3)
plot(dtB*(1:(inner*outer)),uBs(N/7+1,:),'b','LineWidth',3)
plot(dtB*(1:(inner*outer)),uBs(3*N/7+1,:),'b')
plot(dtB*(1:(inner*outer)),uBs(5*N/7+1,:),'b')
plot(dtB*(1:(inner*outer)),uBs(N+1,:),'b')

figure
plot(xA,UA)
hold on
plot(xB,UB)
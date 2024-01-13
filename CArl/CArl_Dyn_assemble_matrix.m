N=70;
L=3.5;
e=0.5;
E=200000;
rho=1000;
Newmark_param
kappa=E;
es=e;

xA=(1:N)*L/N;
xB=(L-e)+(0:N)*L/N;
%% Assembly KA MA
KA=zeros(N,N);
MA=zeros(N,N);

for ii=1:N
    if ii<((L-e)/L*N)
        KA(ii,ii)=2*E*N/L;
        MA(ii,ii)=2*rho*L/(3*N);

        if (ii-1)>0
            KA(ii-1,ii)=-E*N/L;
            MA(ii-1,ii)=rho*L/(6*N);
        end

        KA(ii+1,ii)=-E*N/L;
        MA(ii+1,ii)=rho*L./(6*N);

    elseif ii==((L-e)/L*N)
        KA(ii-1,ii)=-E*N/L;
        MA(ii-1,ii)=rho*L/(6*N);

        KA(ii,ii)=3*E*N/(2*L);
        MA(ii,ii)=1*rho*L/(2*N);

        KA(ii+1,ii)=-E*N/(2*L);
        MA(ii+1,ii)=rho*L/(12*N);

    elseif ii<N

        KA(ii,ii)=E*N/(L);
        MA(ii,ii)=rho*L/(3*N);

        KA(ii-1,ii)=-E*N/(2*L);
        MA(ii-1,ii)=rho*L/(12*N);

        KA(ii+1,ii)=-E*N/(2*L);
        MA(ii+1,ii)=rho*L/(12*N);
    else
        KA(ii,ii)=E*N/(2*L);
        MA(ii,ii)=rho*L/(6*N);

        KA(ii-1,ii)=-E*N/(2*L);
        MA(ii-1,ii)=rho*L/(12*N);
    end
end

%% Assembly KB MB
KB=zeros(N+1,N+1);
MB=zeros(N+1,N+1);

for ii=1:(N+1)
    if ii==1
        KB(ii+1,ii)=-E*N/(2*L);
        MB(ii+1,ii)=rho*L/(12*N);
        KB(ii,ii)=E*N/(2*L);
        MB(ii,ii)=rho*L/(6*N);

    elseif ii<(e/L*N+1)

        KB(ii,ii)=E*N/(L);
        MB(ii,ii)=rho*L/(3*N);

        KB(ii-1,ii)=-E*N/(2*L);
        MB(ii-1,ii)=rho*L/(12*N);

        KB(ii+1,ii)=-E*N/(2*L);
        MB(ii+1,ii)=rho*L/(12*N);

    elseif ii==(e/L*N+1)
        KB(ii-1,ii)=-E*N/(2*L);
        MB(ii-1,ii)=rho*L/(12*N);
        KB(ii+1,ii)=-E*N/L;
        MB(ii+1,ii)=rho*L/(6*N);
        KB(ii,ii)=3*E*N/(2*L);
        MB(ii,ii)=rho*L/(2*N);

    elseif ii<(N+1)

        KB(ii,ii)=2*E*N/L;
        MB(ii,ii)=2*rho*L/(3*N);

        KB(ii-1,ii)=-E*N/L;
        MB(ii-1,ii)=rho*L/(6*N);

        KB(ii+1,ii)=-E*N/L;
        MB(ii+1,ii)=rho*L./(6*N);

    else
        KB(ii,ii)=E*N/L;
        MB(ii,ii)=rho*L/(3*N);

        KB(ii-1,ii)=-E*N/L;
        MB(ii-1,ii)=rho*L./(6*N);
    end
end

MA=MA+betaA*dtA*dtA*KA;
MB=MB+betaB*dtB*dtB*KB;
%% Assembly CA,CB
Ne=e/L*N;

CA=zeros(Ne+1,N);
CB=zeros(Ne+1,N+1);

for ii=1:(Ne+1)
    if ii==1
        CA(ii,N-Ne)=kappa*((N)/L+1/(es^2)*L/(3*N));
        CA(ii,N-Ne+1)=kappa*(-N/L+1/(es^2)*L/(6*N));
    elseif ii<=Ne
        CA(ii,ii+N-Ne-2)=kappa*(-N/L+1/(es^2)*L/(6*N));
        CA(ii,ii+N-Ne-1)=kappa*((2*N)/L+1/(es^2)*2*L/(3*N));
        CA(ii,ii+N-Ne)=kappa*(-N/L+1/(es^2)*L/(6*N));
    else
        CA(ii,N-1)=kappa*(-N/L+1/(es^2)*L/(6*N));
        CA(ii,N)=kappa*((N)/L+1/(es^2)*L/(3*N));
    end
end

for ii=1:(Ne+1)
    if ii==1
        CB(ii,1)=kappa*((N)/L+1/(es^2)*L/(3*N));
        CB(ii,2)=kappa*(-N/L+1/(es^2)*L/(6*N));
    else
        CB(ii,ii-1)=kappa*(-N/L+1/(es^2)*L/(6*N));
        CB(ii,ii)=kappa*((2*N)/L+1/(es^2)*2*L/(3*N));
        CB(ii,ii+1)=kappa*(-N/L+1/(es^2)*L/(6*N));
    end
end

%% Assembly FA,FB
FA=zeros(N,1);
FB=zeros(N+1,1);

FB(N+1)=10;

%% Assembly H
%H=CA*inv(MA)*transpose(CA)+CB*inv(MB)*transpose(CB);
H=betaA*dtA*dtA*CA*inv(MA)*transpose(CA)+betaB*dtB*dtB*CB*inv(MB)*transpose(CB);
%% LightCompetitionStrict
%
% This code computes the equilibrium density of a light competition model with different tradeoffs. 
% Supplement of 
%  " _Maintenance of high diversity in mechanistic forest dynamics models of
% competition for light_ "
%
% Matteo Detto, Jonathan M. Levine and Stephen W. Pacala
% Department of Ecology and Evolutionary Biology, Princeton University
%


%% Description
% Inputs
%
% * _t0_: fixed disturbace interval
% * _S_ : number of species in the pool
% * options: a structure contain model specific parameters
%
% Outputs
%
% * _k_ : a vector of parameter _k_
% * _n_ : a vector of equilibrium species density
% * _t_ : a vector of closing times
% * _ks_ : a vector of establishment conditions ( _k>ks_ )



function [k,n,ks,t,kk,nk,n0,k0,kmax,S0] = LightCompetitionStrict(t0,S,options)

S = 100;
gamma = 1.5;
n = zeros(S,1);
t = zeros(S,1);
ks = zeros(S,1);
k = zeros(S,1);
opts = optimset('TolX',1e-18,'TolFun',1e-18);

if strcmp(options.model,'out vs. up')
%% growth out vs. up
    F = options.F;
    
    c1 = gamma/(gamma+1);
    c2 = gamma+1;
    n0 = F/c2;
    k0 = c2/(F*t0^gamma);
    kmax = k0*options.kmax;
    lambda = kmax./(S+1);
    k(1) = k0+lambda;
    t1 = fminbnd(@(x) abs(F*c1*(t0-x)/t0 + F/c2 - (x^(-gamma))/k(1)),0,t0,opts);
    ks1 = c2/F/t1^c2*t0;
    S0 = (kmax-ks1)./(2*lambda)+1;
    

    n=zeros(S,1);
    t=zeros(S,1);
    ks = zeros(S,1);
    k = zeros(S,1);
    kr = rand(S,1)*kmax;
    k(1) = min(kr(kr>k0));
    t(1) = fminbnd(@(x) abs(F*c1*(t0-x)/t0 + F/c2 - (x^(-gamma))/k(1)),0,t0,opts);
    n(1) = F*c1*(t0 - t(1))/t0 + F/c2;
    Z = k(1)*n(1);
    ks(1) = c2/F/t(1)^c2*t0;
    
    i=1;
    while ks(i)<max(kr)
        i=i+1;
        k(i) = min(kr(kr>ks(i-1)));
        
        t(i)  = fminbnd(@(x) abs(F*c1*(t(i-1)-x)/t0 - (x^(-gamma) - Z)/k(i)),...
            0,(t0*c2./k(i)./F).^(1/c2),opts);
        n(i) = F*c1*(t(i-1) - t(i))/t0;
        Z = Z+k(i)*n(i);
        ks(i) = c2/F/t(i)^c2*t0;
    end
   
    n=n(1:i);
    k=k(1:i);
    ks=ks(1:i);
    t=t(1:i);
    
    kk=linspace(k0,kmax,1000);
    nk = c1*(n0/t0)^c1*kk.^(-(gamma+2)/(gamma+1));
    

elseif strcmp(options.model,'growth vs. fecundity')
%% growth vs. fecundity

    ns = options.ns;
    m = options.m;
    c = options.c;
    c1 = gamma/(gamma+1);
    c2 = gamma+1;
%     c  = 0.11;
%     m  = 100;
%     ns = 1000;

    
    fi0    = fzero(@(x) ns*x*m*(c-x)^gamma*t0^gamma - c2,[0 c/c2],opts);
    fimax = c/c2;
    fir = fi0 + rand(S,1)*(fimax - fi0);
            
    n=zeros(S,1);
    t=zeros(S,1);
    fi = zeros(S,1);
    fis = zeros(S,1);
    k = zeros(S,1);
    

    fi(1) = min(fir);
    k(1) = m*(c-fi(1)).^gamma;
    F = ns*fi(1);
    t(1) = fminbnd(@(x) abs(F*c1*(t0-x)/t0 + F(1)/c2 - (x^(-gamma))/k(1)),0,t0,opts);
    n(1) = F*c1*(t0 - t(1))/t0 + F/c2;
    Z = k(1)*n(1);
    
  
try
    fis(1) = fzero(@(x) ns*x*m*(c-x)^gamma*t(1)^c2 - c2*t0,[fi0 c/c2],opts);
catch ME
    if (strcmp(ME.identifier,'MATLAB:fzero:ValuesAtEndPtsSameSign'))
        fis(1) = fimax;
    end
end

i=1;
while fis(i) < max(fir) && i<S
    i=i+1;
    fi(i) = min(fir(fir>fis(i-1)));
    k(i) = m*(c-fi(i)).^gamma;
    F = ns*fi(i);
    t(i) = fminbnd(@(x) abs(F*c1*(t(i-1)-x)/t0 - (x^(-gamma) - Z)/k(i)),0,t(i-1),opts);
    n(i) = F*c1*(t(i-1) - t(i))/t0;
    Z = Z+k(i)*n(i);
    try
        fis(i) = fzero(@(x) ns*x*m*(c-x)^gamma*t(i)^c2 - c2*t0,[fis(i-1) c/c2],opts);
    catch ME
        if (strcmp(ME.identifier,'MATLAB:fzero:ValuesAtEndPtsSameSign'))
            fis(i) = fimax;
        end
    end
end

n=n(1:i-1);
fi=fi(1:i-1);
fis=fis(1:i-1);
t=t(1:i-1);
    
%% annaual plants (fecundity ~ crown area)
elseif strcmp(options.model,'annual') 

    F = options.F;
    k0 = 1/(F*t0^(1/gamma));
    kmax = k0*options.kmax;
    lambda = kmax./(S+1);
    k(1) = k0+lambda;


    t1 = fminbnd(@(x) abs(F*gamma*(log(t0)-log(x))+F - (x^(-gamma))/k(1)),0,t0,opts);
    ks1=t1.^-gamma/F;
    S0 = (kmax-ks1)./(2*lambda)+1;
    n0=F;
    

    kr = rand(S,1)*kmax;
    k(1) = min(kr(kr>k0));    
    t(1) = fminbnd(@(x) abs(F*gamma*(log(t0)-log(x))+F - (x^(-gamma))/k(1)),0,t0,opts);
    n(1) = F*gamma*log(t0/t(1))+F;
    Z = k(1)*n(1);
    ks(1)=Z/F;
    ti=t(1);
    i=1;
    while ks(i)<max(kr)
        i=i+1;
        k(i) = min(kr(kr>ks(i-1)));
        t(i) = fminbnd(@(x) abs(F*gamma*(log(t(i-1))-log(x)) - (x^(-gamma) - Z)/k(i)),0,ti,opts);
        n(i) = F*gamma*(log(ti)-log(t(i)));
        ti = t(i);
        Z = Z+k(i)*n(i);
        ks(i)=Z/F;
    end

    
    n=n(1:i);
    k=k(1:i);
    ks=ks(1:i);
    t=t(1:i);

    kk=linspace(k0,kmax,1000);
    nk = F./kk;
    

elseif strcmp(options.model,'understory growth')
%% understory growth  

    nu = options.nu;
    F = options.F;
    k0 = 1/(F*t0^gamma);
    kmax = k0*options.kmax;
    lambda = kmax./(S+1);
    k(1) = k0+lambda;
    
    H = @(t,ti) integral(@(x) gamma*((1-nu)*x + nu*t0).^gamma./x.^(gamma+1),t,ti);
       
    t1 = fminbnd(@(x) abs(F*H(x,t0) + F - x^(-gamma)/k(1)),0,t0,opts);
    ks1=((1-nu)*t1 + nu*t0)^(-gamma)/F;
    S0 = (kmax-k0-ks1)./(2*lambda)+1;
    n0=F;
    
    kr = rand(S,1)*kmax;
    k(1) = min(kr(kr>k0));
    t(1) = fminbnd(@(x) abs(F*H(x,t0) + F - x^(-gamma)/k(1)),0,t0,opts);
%     n(1) = F*(H(t0) - H(t(1)))+F;
    n(1) = F*H(t(1),t0)+F;
    ks(1)=((1-nu)*t(1) + nu*t0)^(-gamma)/F;
    Z = k(1)*n(1);
    i=1;
    while ks(i)<max(kr)
        i=i+1;
        k(i) = min(kr(kr>ks(i-1)));
        t(i) = fminbnd(@(x) abs(F*H(x,t(i-1)) - (x^(-gamma) - Z)/k(i)),0,t(i-1),opts);
%         n(i) = F*(H(t(i-1))-H(t(i)));
        n(i) = F*H(t(i),t(i-1));
        Z = Z+k(i)*n(i);
        ks(i)=((1-nu)*t(i)+nu*t0)^(-gamma)/F;
    end
    

    n=n(1:i);
    k=k(1:i);
    ks=ks(1:i);
    t=t(1:i);
    
    kk=linspace(k0,kmax,1000);
    nk = (1-nu)^gamma ./ (1 - nu*(kk/k0).^(1/gamma)).^(gamma+1).*(n0./kk);
    

elseif strcmp(options.model,'growth vs. maturation')
%%  growth vs. maturation
c = gamma+1;
n=zeros(S,1);
z=zeros(S,1);
t=zeros(S,1);
zs = zeros(S,1);
k = 20;
F = 1;
opts = optimset('TolX',1e-12,'TolFun',1e-6);
td = 1;

z0 = fzero(@(y) F*k*y^gamma*(td-y)-1,[td*gamma/c td]);
n0=F*(td-z0);

zmin = fzero(@(y) F*k*y^gamma*(td-y)-1,0.1);
lambda = (z0-zmin)/S;

    
% species #1
z(1)  = z0 + log(rand)*lambda;    
t(1)  = (F*k*(td-z(1))).^(-1/gamma);
n(1)  = F*(td-z(1));
zs(1) = fzero(@(y) F*k*y^gamma*(t(1)-y)-1,[t(1)*gamma/c t(1)]);

% species #2 
z(2) = zs(1) + log(rand)*lambda;
y0(1) = (1/F-k*n(1)*(1/(k*z(2)^gamma)+z(2))^gamma)/(k*z(2)^gamma);
y0(2) = (1/k - n(1)*z(2)^gamma)/z(2).^gamma;

n(2) = fzero(@Solution2,y0,opts);
t(2) = ((1/k - n(2).*z(2).^gamma)/n(1)).^(1/gamma);
zs(2) = fzero(@(y) F*k*y^gamma*(t(2)-y)-1,[t(2)*gamma/c t(2)]);
    
       
i=2;
while t(i)>z(2)
    i=i+1;
    z(i) = zs(i-1) + log(rand)*lambda;
    y0(1)   = n(1)*(t(i-1)/z(i))^gamma - n(1)*(1+1/(F*k*z(i)^c))^gamma;
    y0(2)   = (1/k - sum(n(1:i-1))*z(i)^gamma)/z(i).^gamma;   
    n(i) = fzero(@Solution,y0,opts);
    t(i) = ((1/k - sum(n(2:i).*z(2:i).^gamma))/n(1)).^(1/gamma);
    zs(i)= fzero(@(y) F*k*y^gamma*(t(i)-y)-1,[t(i)*gamma/c t(i)]);

end


n=n(1:i-1);
z=z(1:i-1);
zs=zs(1:i-1);
t=t(1:i-1);



elseif strcmp(options.model,'growth vs. survival')
%% growth vs. survival

%inputs
% t0: constant distrbance interval
% c = exponent of the survival function Fmax./(1+(k./kstar).^c)
% Fmax: max fecundity
% kstar: scaling parameter for growth
kmax = options.kmax;
c = options.c;
Fmax = options.Fmax;
omega0 = (gamma+1)./t0^gamma;
kmin = kmax*(gamma/(c + gamma))^gamma;

k0 = fzero(@(x) Fmax*x.*(1-(x./kmax).^(1/gamma)).^c  - omega0,[kmin kmax],opts);

W = @(t) t./t0.*gamma/(gamma+1);

kr = kmin + rand(S,1)*(k0 - kmin);

n  = zeros(S,1);
t  = zeros(S,1);
F  = zeros(S,1);
ks = zeros(S,1);
k  = zeros(S,1);


k(1) = max(kr);
F(1) = Fmax.*(1-(k(1)./kmax).^(1/gamma)).^c;

t(1) = fzero(@(x) 1 - W(x) - x^(-gamma)/k(1)/F(1),[0.001*t0 t0],opts);
n(1) = t(1)^(-gamma)/k(1);
Z = k(1)*n(1);
omega = omega0*(t0./t(1)).^(gamma+1);
ks(1) = fzero(@(x) Fmax*x.*(1-(x./kmax).^(1/gamma)).^c - omega,[kmin k(1)],opts);
    
i=1;
while ks(i)>min(kr)
    i=i+1;
    k(i) = max(kr(kr<ks(i-1)));
    F(i) = Fmax.*(1-(k(i)./kmax).^(1/gamma)).^c;
    W1 = W(t(i-1));
    t(i) = fminbnd(@(x) abs(W1 - W(x) - (x^(-gamma) - Z)/k(i)/F(i)),0,t(i-1),opts);
    n(i) = (t(i)^(-gamma)-Z)/k(i);
    Z = Z+k(i)*n(i);
    omega = omega0*(t0./t(i)).^(gamma+1);   
    try
    ks(i) = fzero(@(x) Fmax*x.*(1-(x./kmax).^(1/gamma)).^c - omega,[kmin k(i)],opts);
    catch ME
        if (strcmp(ME.identifier,'MATLAB:fzero:ValuesAtEndPtsSameSign'))
            ks(i) = kmin;
        end
    end
end
    
    
k=k(1:i);
n=n(1:i);
t=t(1:i);
F=F(1:i);
ks=ks(1:i);

end

%% plotting
if options.graph && (strcmp(options.model,'out vs. up') || strcmp(options.model,'understory growth') || strcmp(options.model,'annual'))
    
    
    subplot(options.sbpl(1))
    semilogy(k(1:end)/k0,n(1:end),'o','linewidth',2,'markersize',5)
    hold all;plot(k(1)/k0,n(1),'bo','linewidth',2,'markersize',5)
    ylabel('{\itN_i}','FontName','Cambria Math','interpreter','tex')
    title(options.model)
    
    subplot(options.sbpl(2))
    dk=diff(ks);
    plot(k(2:end)/k0,n(2:end)./dk,'o',kk/k0,nk,'r','linewidth',2,'markersize',5)
    hold all
    ylabel('{\itN_i} /  \Delta{\itk_i^*}','FontName','Cambria Math')
    xlabel('\it{k / k}_0','FontName','Cambria Math')
    
    
elseif options.graph && strcmp(options.model,'growth vs. fecundity')
    
    ff = linspace(min(fir),c/c2, 1000);
    nf = ns*gamma*(c2^(2*gamma+1)*t0^gamma.*ns.*ff.*m.*(c-ff).^(2*gamma+1)).^(-1/c2).*(c - ff*c2);
    
    subplot(options.sbpl(1))
    semilogy(fi(2:end)/c*100,n(2:end),'o','linewidth',2,'markersize',5)
    hold all
    plot(fi(1)/c*100,n(1),'bo','linewidth',2,'markersize',5)
    ylabel('{\itN_i}','FontName','Cambria Math','interpreter','tex')
    title('growth vs. fecundity')
    
    subplot(options.sbpl(2))
    df=diff(fis);
    plot(fi(2:end)/c*100,n(2:end)./df,'o','linewidth',2,'markersize',5);
    hold all
    plot(ff/c*100,nf,'r','linewidth',2)
    ylabel('{\itN_i} / \Delta ?_{\iti}^*','FontName','Cambria Math','interpreter','tex')
    xlabel('allocation to reproduction (%)','FontName','Cambria Math','interpreter','tex')
    
elseif options.graph && strcmp(options.model,'growth vs. survival')
    
    kk = linspace(k(end),k(2),1000);
    Fk = Fmax.*(1-(kk./kmax).^(1/gamma)).^c;
    muk = Fmax./Fk;
    
    dFdmu = -Fmax./muk.^2;
    tk = ((gamma+1)*t0./(Fk.*kk)).^(1./(gamma+1));
    dFdk = -Fmax*c*(1 - (kk/kmax).^(1/gamma)).^(c - 1).*(kk/kmax).^(1/gamma - 1)/(gamma*kmax);
    dtdk = -tk.^(-gamma).*t0./(Fk.*kk).^2.*(Fk+kk.*dFdk);
    nF = -gamma./(kk.*tk.^(gamma+1)).*dtdk./dFdk.*dFdmu;
    Fs = Fmax.*(1-(ks./kmax).^(1/gamma)).^c;
    mus = Fmax./Fs;
    
    subplot(options.sbpl(1))
    loglog(Fmax./F,n,'o','linewidth',2,'markersize',5)
    ylabel('{\itN_i}','FontName','Cambria Math','interpreter','tex')
    title('growth in light vs. survival in shade')
    
    subplot(options.sbpl(2))
    semilogx(Fmax./F(2:end),-n(2:end)./diff(log(mus)),'o','linewidth',2,'markersize',5)
    hold all
    semilogx(muk,-nF.*muk,'linewidth',2)
    xlabel('mortality in the shade (\it{m}_0 / \it{m})','FontName','Cambria Math','interpreter','tex')
    ylabel('{\itN_i} /  \Delta{log\it(m_i^*)}','FontName','Cambria Math')
    
elseif options.graph && strcmp(options.model,'growth vs. maturation')
    
    
subplot(options.sbpl(1))
semilogy(z(1:end)/z0,n(1:end),'o')
hold all;plot(z(1)/z0,n(1),'bo','linewidth',2)
ylabel('n_i')
xlabel('plant reproductive threshold (\it{z / z}_0)')
ylabel('\it{n}')
xlim([floor(z(end)/z0*100)/100 1])
title({'forest with species-specific','reproductive threshold'})

subplot(options.sbpl(2))
dz=diff(zs);
plot(z(2:end)/z0,-n(2:end)./dz,'o','linewidth',2)
hold all

zmin = fzero(@(y) y+1/(F*k*y^gamma)-z(2),z(2));
x=linspace(zmin,zs(1),1000);
mod = gamma*n(1)./x.^gamma.*(x+1./(F*k*x.^gamma)).^(gamma-1).*(1-gamma./(F*k*x.^c));

plot(x/z0,mod,'-','linewidth',2)
xlim([floor(z(end)/z0*100)/100 1])
ylabel('n(k)')
xlabel('plant reproductive threshold (\it{z / z}_0)')
ylabel('\it{n / dz}')
    
end

%% function for i=2
function y = Solution2(ni)
   
    ti = ((1/k - ni.*z(2).^gamma)/n(1)).^(1/gamma);
    h  = ni.*k.*z(2).^gamma.*(ti-z(2)) + t(1)-ti-n(1)*k*(t(1)^c-ti.^c)/c;
    y  = ni - F*h;
end

%% function for i>2
function y = Solution(ni)
   
    ti = ((1/k - sum(n(2:i-1).*z(2:i-1).^gamma) - ni*z(i).^gamma)/n(1)).^(1/gamma);
    h = ni*k*z(i)^gamma*(ti-z(i));
    h = h + (1-k*sum(n(2:i-1).*z(2:i-1).^gamma)).*(t(i-1) - ti) - n(1)*k*(t(i-1).^c - ti.^c)/c;
    y = ni - F*h;
end

end

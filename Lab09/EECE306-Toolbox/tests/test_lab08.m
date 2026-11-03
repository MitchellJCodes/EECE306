function test_08 = test_lab08()
%TEST_LAB08
e0=em.const.eps0();

%% Parallel plate capacitance
n=81; G=em.solve.grid2d([0 1],[0 1],n,n);
iy1=round(.25*(n-1))+1; iy2=round(.75*(n-1))+1;
mask=false(n,n); mask(iy1,:)=true; mask(iy2,:)=true;
Vfix=zeros(n,n); Vfix(iy2,:)=1;
C=em.derive.capacitance(G,1,mask,Vfix,struct('method','direct'));
Cref=e0/0.5;
assert(abs(C-Cref)/Cref<0.01,'Parallel-plate capacitance error exceeds 1%%.');

%% Two dielectric layers: series formula and refinement
% Even-sized grids put y=0.5 on a face, where the material interface belongs.
Ns=[42 82 162]; Ch=zeros(size(Ns)); Ca=Ch; Crefs=Ch;
for k=1:numel(Ns)
    nn=Ns(k); Gk=em.solve.grid2d([0 1],[0 1],nn,nn);
    [~,a]=min(abs(Gk.y-.25)); [~,b]=min(abs(Gk.y-.75));
    mk=false(nn,nn); mk(a,:)=true; mk(b,:)=true;
    vf=zeros(nn,nn); vf(b,:)=1;
    er=ones(nn,nn); er(Gk.Y<.5)=4;
    d1=.5-Gk.y(a); d2=Gk.y(b)-.5;
    Crefs(k)=e0/(d1/4+d2);
    Ch(k)=em.derive.capacitance(Gk,er,mk,vf,struct('method','direct','faceMean','harmonic'));
    Ca(k)=em.derive.capacitance(Gk,er,mk,vf,struct('method','direct','faceMean','arithmetic'));
end
assert(max(abs(Ch-Crefs)./Crefs)<0.01,'Two-layer harmonic capacitance misses series formula by >1%%.');
ea=abs(Ca-Crefs);
pa=log(ea(1:end-1)./ea(2:end))./log((Ns(2:end)-1)./(Ns(1:end-1)-1));
pA=mean(pa);
% For this aligned 1D layered case, harmonic face conductance reproduces
% the series resistance to roundoff, so its empirical order is not resolvable.
% The arithmetic-face error does refine; report that measured order.
assert(pA<1.7,sprintf('Arithmetic mean did not show degraded convergence (p=%.3f).',pA));
pH=NaN;

%% Coaxial capacitance on a Cartesian grid
% Refine the Cartesian grid to reduce staircasing error at circular conductors.
n=201; Gc=em.solve.grid2d([-1 1],[-1 1],n,n);
RR=hypot(Gc.X,Gc.Y); a=.3; b=.9;
mc=RR<=a | RR>=b; vfc=zeros(n,n); vfc(RR<=a)=1;
Cc=em.derive.capacitance(Gc,1,mc,vfc,struct('method','direct'));
CcRef=2*pi*e0/log(b/a);
assert(abs(Cc-CcRef)/CcRef<0.02,sprintf('Coax capacitance error %.2f%% exceeds 2%%.',100*abs(Cc-CcRef)/CcRef));

%% Capacitance independent of applied voltage
v100=vfc; v100(RR<=a)=100;
Cc100=em.derive.capacitance(Gc,1,mc,v100,struct('method','direct'));
assert(abs(Cc100-Cc)/Cc<1e-10,'Capacitance depends on applied voltage.');

%% Filling with a uniform dielectric multiplies capacitance by epsr
Cfill=em.derive.capacitance(Gc,3,mc,vfc,struct('method','direct'));
assert(abs(Cfill/Cc-3)<1e-10,'Uniform dielectric did not scale capacitance by epsr.');

%% Resistance of a uniform bar
Gr=em.solve.grid2d([0 1],[0 1],81,81);
mr=false(81,81); mr(:,1)=true; mr(:,end)=true;
vr=zeros(81,81); vr(:,1)=1; sigma=100;
Rb=em.derive.resistance(Gr,sigma,mr,vr,struct('method','direct'));
Rref=1/(sigma*1);
assert(abs(Rb-Rref)/Rref<0.01,sprintf('Bar resistance error %.3g.',abs(Rb-Rref)/Rref));

%% Tangential E is continuous across a vertical dielectric interface
Gt=em.solve.grid2d([-1 1],[0 1],81,81);
ert=ones(81,81); ert(:,Gt.x<0)=4;
mt=false(81,81); mt(1,:)=true; mt(end,:)=true;
vt=zeros(81,81); vt(end,:)=1;
Vt=em.solve.poisson2d(Gt,0,ert,mt,vt,struct('method','direct'));
% Tangential direction is y for the vertical interface x=0
iy=round(.4*(81-1))+1; ixL=find(Gt.x<0,1,'last'); ixR=find(Gt.x>0,1,'first');
EyL=-(Vt(iy+1,ixL)-Vt(iy-1,ixL))/(2*Gt.hy);
EyR=-(Vt(iy+1,ixR)-Vt(iy-1,ixR))/(2*Gt.hy);
assert(abs(EyL-EyR)<1e-8,'Tangential electric field is discontinuous across interface.');

fprintf('  Lab 8 parallel-plate rel err %.3g; two-layer harmonic order unresolved (aligned 1-D case); arithmetic order %.3f; coax rel err %.3g\n', ...
    abs(C-Cref)/Cref,pA,abs(Cc-CcRef)/CcRef);
disp('  test_lab08 checks complete');
end

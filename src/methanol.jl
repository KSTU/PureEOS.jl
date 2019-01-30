"""
Calculate methanol critical properties
T --- temperature in Kelvin
p --- pressure in bar
ro --- density in mol/l
"""
function methanolCritical()
	sub=newCritical()
	sub.T=512.60
	sub.p=8.1035*10	#MPa to bar
	sub.ro=275.56/32.04	#mol/l #275.56 kg/m-3
	return sub
end

"""
Calculate mathanol properties at triple point
T --- temperature in Kelvin
p --- pressure in bar
vRo --- vapor density in mol/l
lRo --- liquid density in mol/l
url = https://www.engineeringtoolbox.com/methanol-methyl-alcohol-properties-CH3OH-d_2031.html
"""
function methanolTriple()
	sub=newTriple()
	sub.T=175.50
	sub.p=1.86*10^(-8)	#MPa to bar
	sub.lRo=-1.0	#mol/l
	sub.vRo=-1.0	#mol/l
	return sub
end


"""
calculate mathanol properties at saturation line

T --- temperature in Kelvin
p --- pressure in bar
vRo --- vapor density in mol/l
lRo --- liquid density in mol/l
"""
function methanolSaturate(T)
A=[-8.8570247, 2.4072447, -2.6452501, -1.5044111]
pA=[1.0, 1.5, 2.0, 3.5]

B=[-1.131507, 7.463057, -16.52860, -16.93269, -77.01263]
pB=[0.2, 0.8, 2.5, 4.5, 8.0]

C=[-0.0770237, 1.686259, -0.5426160, -0.3579365]
pC=[1.0/6.0, 0.3, 1.0 , 2.4]

sub=newSaturate()
crit=methanolCritical()
if (T<crit.T)
	v=1.0-T/crit.T
	sub.p=crit.p*exp( 1.0/(1.0-v) * (A[1]*v^pA[1] + A[2]*v^pA[2] + A[3]*v^pA[3] + A[4]*v^pA[4]) )
	sub.vRo=crit.ro*exp(B[1]*v^pB[1] + B[2]*v^pB[2] + B[3]*v^pB[3] + B[4]*v^pB[4] + B[5]*v^pB[5])
	sub.lRo=crit.ro*exp(C[1]*v^pC[1] + C[2]*v^pC[2] + C[3]*v^pC[3] + C[4]*v^pC[4])
	sub.T=T
else
	sub.p=-1
	sub.lRo=-1
	sub.vRo=-1
	sub.T=-1
end
return sub
end

"""
calculate methanol properties at single phase


T --- temperature, Kelvin
p --- pressure, bar
ro --- density, mol/l
z ---
s --- entropy,
s0 --- ideal part of entropy,
ds --- nonideal part derived k_B, dimensionless
g --- Gibbs energy
g0 --- ideal part of Gibbs energy
dg --- nonideal part of Gibbs energy
f ---  Helmholtz free energy
f0 --- ideal part
df --- nonideal part, dimensionless
u --- internal energy
du --- nonideal part
h --- entalpy
dh --- nonideal part
w --- speed of sound,
cp ---
cv ---
"""
function methanolSinglephase(T,ro)

sub=newSingle()
crit=methanolCritical()
delta=ro/crit.ro
tau=crit.T/T
R=8.3144598
mm=0.03204

nOne=6
nTwo=10
nThree=0
nFour=0


c=[0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 1.000, 1.000, 1.000, 1.000,
1.000, 1.000, 1.000, 1.000, 1.000, 
2.000, 2.000, 2.000, 2.000, 2.000,
2.000, 2.000, 2.000, 2.000, 3.000,
3.000, 4.000]

d=[1.0, 1.0, 2.0, 2.0, 3.0,
3.0, 1.0, 2.0, 3.0, 3.0,
4.0, 5.0, 6.0, 7.0, 8.0,
1.0, 1.0, 2.0, 2.0, 3.0,
4.0, 5.0, 6.0, 7.0, 7.0,
6.0, 5.0]

t=[0.5, 0,750, 0,125, 1.5, 0,375,
1.750, 0.000, 0.500, 1.000, 3.750, 
2.000, 2.500, 0.000, 2.500, 4.500,
3.000, 4.000, 1.000, 3.000, 5.000,
6.000, 2.500, 5.000, 6.000, 12.000,
10.000, 15.000]

n=[0.12622395e02, -0.83224516e01,
-0.14647501e01, -0.12954522e1,
0.22417697, 0.22830533,
-0.46549039e01, -0.41099957e01,
0.70421007, 0.81617251e-01,
-0.37777607, 0.19627811,
0.45571723e-02, -0.11777859e-01,
-0.17277890e-04, 0.19096856e01,
-0.29551319e01, -0.28958480,
0.18169967e01, -0.96254996,
-0.11885503, -0.10730710,
0.15487654e-01, 0.66239025e-03,
0.15286750e-01, 0.31218155e-02,
0.14740469e-01]

al=[0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000, 
0.000, 1.075, 0.463, 0.876, 1.108, 
0.741, 4.032, 2.453, 2.300, 3.143]

bt=[0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 1.207, 0.090, 0.581, 0.947,
2.356, 27.010, 4.542, 1.287, 3.090]

gm=[0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 1.194, 1.986, 1.583, 0.756,
0.495, 1.002, 1.077, 1.493, 1.542]

ep=[0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.000, 0.000, 0.000, 0.000,
0.000, 0.779, 0.805, 1.869, 0.694,
1.312, 2.054, 0.441, 0.793, 0.313]

phidelta=0.0
phir=0.0
phitau=0.0
phidelta2=0.0
phitau2=0.0
dphideltadtau=0.0

i=0
for ni=1:nOne
	i=i+1
	phidelta=phidelta+n[i]*d[i]*delta^(d[i]-1.0)*tau^t[i]
	phir=phir+n[i]*delta^d[i]*tau^t[i]
	phitau=phitau+n[i]*t[i]*delta^d[i]*tau^(t[i]-1)
	phidelta2=phidelta2+n[i]*d[i]*(d[i]-1.0)*delta^(d[i]-2)*tau^t[i]
	phitau2=phitau2+n[i]*t[i]*(t[i]-1.0)*delta^d[i]*tau^(t[i]-2.0)
	dphideltadtau=dphideltadtau+n[i]*d[i]*t[i]*delta^(d[i]-1.0)*tau^(t[i]-1.0)
end
for ni=1:nTwo
	i+=1
	phidelta=phidelta+n[i]*exp(-delta^c[i])*(delta^(d[i]-1.0)*tau^t[i]*(d[i]-c[i]*delta^c[i]))
	phir=phir+n[i]*delta^(d[i])*tau^t[i]*exp(-delta^c[i])
	phitau=phitau+n[i]*t[i]*delta^d[i]*tau^(t[i]-1.0)*exp(-delta^c[i]);
	phidelta2=phidelta2+n[i]*exp(-delta^c[i])*(delta^(d[i]-2)*tau^(t[i])*((d[i]-c[i]*delta^c[i])*(d[i]-1.0-c[i]*delta^c[i])-c[i]^2*delta^c[i]))
	phitau2=phitau2+n[i]*t[i]*(t[i]-1.0)*delta^d[i]*tau^(t[i]-2)*exp(-delta^c[i])
	dphideltadtau=dphideltadtau+n[i]*t[i]*delta^(d[i]-1.0)*tau^(t[i]-1.0)*(d[i]-c[i]*delta^c[i])*exp(-delta^c[i])
end
for pn=1:nThree
	i+=1
	phidelta=phidelta+n[i]*delta^d[i]*tau^t[i]*exp(-al[i]*(delta-ep[i])^2-bt[i]*(tau-gm[i])^2)*(d[i]/delta-2.0*al[i]*(delta-ep[i]))
	phir=phir+n[i]*delta^d[i]*tau^t[i]*exp(-al[i]*(delta-ep[i])^2-bt[i]*(tau-gm[i])^2)
	phitau=phitau+n[i]*delta^d[i]*tau^t[i]*exp(-al[i]*(delta-ep[i])^2-bt[i]*(tau-gm[i])^2)*(t[i]/tau-2*bt[i]*(tau-gm[i]))
	phidelta2=phidelta2+n[i]*tau^t[i]*exp(-al[i]*(delta-ep[i])^2-bt[i]*(tau-gm[i])^2)*(-2*al[i]*delta^d[i]+4*al[i]^2*delta^d[i]*(delta-ep[i])^2-4*d[i]*al[i]*delta^(d[i]-1)*(delta-ep[i])+d[i]*(d[i]-1)*delta^(d[i]-2))
	phitau2=phitau2+n[i]*delta^d[i]*tau^t[i]*exp(-al[i]*(delta-ep[i])^2-bt[i]*(tau-gm[i])^2)*((t[i]/tau-2*bt[i]*(tau-gm[i]))^2-t[i]/tau^2-2*bt[i])
	dphideltadtau=dphideltadtau+n[i]*delta^d[i]*tau^t[i]*exp(-al[i]*(delta-ep[i])^2-bt[i]*(tau-gm[i])^2)*(d[i]/delta -2*al[i]*(delta-ep[i]))*(t[i]/tau-2*bt[i]*(tau-gm[i]))
end

sub.T=T
sub.ro=ro
sub.p=(1.0+phidelta*delta)*ro*R*T/100.0	#bar
sub.z=1.0+phidelta*delta	#
#ideal part
n0=[2.496674887, 2.900791185, -62.57135350, 10.99267739, 18.33682995, -16.36600476, -6.223234762, 2.803536282, 1.077809894, 0.969656970]
gamma0=[0.0, 0.0, 0.0, 4.119785, 3.264999, 3.769463, 2.931493, 8.225557, 10.31627, 0.5324892]

phi0 = n0[1] + n0[2]*log(tau) + n0[3]*tau + n0[4]*log(exp(gamma0[4]*tau)-1.0) + n0[5]*log(exp(gamma0[5]*tau)-1.0) + n0[6]*log(exp(gamma0[6]*tau)-1.0) + n0[7]*log(exp(gamma0[7]*tau)-1.0) + n0[8]*log(exp(gamma0[8]*tau)-1.0) + n0[9]*log(exp(gamma0[9]*tau)-1.0) + n0[10]*log(exp(gamma0*tau)-1.0)
phi0tau = n0[2]/tau + n0[3] + n0[4]*gamma0[4]*exp(gamma0[4]*tau)/(exp(gamma0[4]*tau)-1.0)+ n0[5]*gamma0[5]*exp(gamma0[5]*tau)/(exp(gamma0[5]*tau)-1.0) + n0[6]*gamma0[6]*exp(gamma0[6]*tau)/(exp(gamma0[6]*tau)-1.0) + n0[7]*gamma0[7]*exp(gamma0[7]*tau)/(exp(gamma0[7]*tau)-1.0) + n0[8]*gamma0[8]*exp(gamma0[8]*tau)/(exp(gamma0[8]*tau)-1.0) + n0[9]*gamma0[9]*exp(gamma0[9]*tau)/(exp(gamma0[9]*tau)-1.0) + n0[10]*gamma0[10]*exp(gamma0[10]*tau)/(exp(gamma0[10]*tau)-1.0)
phi0tau2 = -n0[2]/tau^2 + gamma0[4]*n0[4]*exp(gamma0[4]*tau)/(exp(gamma0[4]*tau)-1.0) - gamma0[4]^2*n0[4]*exp(2.0*gamma0[4]*tau)/(exp(gamma0[4]*tau)-1.0)^2 + gamma0[5]*n0[5]*exp(gamma0[5]*tau)/(exp(gamma0[5]*tau)-1.0) - gamma0[5]^2*n0[5]*exp(2.0*gamma0[5]*tau)/(exp(gamma0[5]*tau)-1.0)^2 + gamma0[6]*n0[6]*exp(gamma0[6]*tau)/(exp(gamma0[6]*tau)-1.0) - gamma0[6]^2*n0[6]*exp(2.0*gamma0[6]*tau)/(exp(gamma0[6]*tau)-1.0)^2 + gamma0[7]*n0[7]*exp(gamma0[7]*tau)/(exp(gamma0[7]*tau)-1.0) - gamma0[7]^2*n0[7]*exp(2.0*gamma0[7]*tau)/(exp(gamma0[7]*tau)-1.0)^2 + gamma0[8]*n0[8]*exp(gamma0[8]*tau)/(exp(gamma0[8]*tau)-1.0) - gamma0[8]^2*n0[8]*exp(2.0*gamma0[8]*tau)/(exp(gamma0[8]*tau)-1.0)^2 + gamma0[9]*n0[9]*exp(gamma0[9]*tau)/(exp(gamma0[9]*tau)-1.0) - gamma0[9]^2*n0[9]*exp(2.0*gamma0[9]*tau)/(exp(gamma0[9]*tau)-1.0)^2 + gamma0[10]*n0[10]*exp(gamma0[10]*tau)/(exp(gamma0[10]*tau)-1.0) - gamma0[10]^2*n0[10]*exp(2.0*gamma0[10]*tau)/(exp(gamma0[10]*tau)-1.0)^2
 
phi0delta = 0.0
phi0delta2 = 0.0
phi0taudelta=0.0
#entropy
#phi0=log(delta)+n0[1]+n0[2]*tau+n0[3]*log(tau)+n0[4]*log(1.0-exp(-gamma0[4]*tau))+n0[5]*log(1.0-exp(-gamma0[5]*tau))+n0[6]*log(1.0-exp(-gamma0[6]*tau))+n0[7]*log(1.0-exp(-gamma0[7]*tau))
#phi0tau=n0[2]+n0[3]/tau+n0[4]*gamma0[4]*(1.0/(1.0-exp(-gamma0[4]*tau))-1.0)+n0[5]*gamma0[5]*(1.0/(1.0-exp(-gamma0[5]*tau))-1.0)+n0[6]*gamma0[6]*(1.0/(1.0-exp(-gamma0[6]*tau))-1.0)+n0[7]*gamma0[7]*(1.0/(1.0-exp(-gamma0[7]*tau))-1.0)
#phi0delta=1/delta
#phi0delta2=-1/delta^2
#phi0tau2=-n0[3]/tau^2-n0[4]*gamma0[4]^2*exp(-gamma0[4]*tau)*(1.0-exp(-gamma0[4]*tau))^(-2)-n0[5]*gamma0[5]^2*exp(-gamma0[5]*tau)*(1.0-exp(-gamma0[5]*tau))^(-2)-n0[6]*gamma0[6]^2*exp(-gamma0[6]*tau)*(1.0-exp(-gamma0[6]*tau))^(-2)-n0[7]*gamma0[7]^2*exp(-gamma0[7]*tau)*(1.0-exp(-gamma0[7]*tau))^(-2)
#phi0taudelta=0.0

#Entropy
sub.s=(tau*(phi0tau+phitau)-phi0-phir)*R
sub.s0=(tau*(phi0tau)-phi0)*R
sub.ds=(tau*(phitau)-phir)	#devidet to R

#out.entrfrac=(tau*(phi0tau+phitau)-phi0-phir)/(tau*(phi0tau)-phi0)

dtaudtau=-tau/T
dphidt=dtaudtau*dphideltadtau

#sub.dpdt=(delta*cur_temp*dphidt+delta*phidelta+1)*R*T/100;
#phi0
#phir
#phidelta
#phidelta*delta
#(1+phi0+phir+phidelta*delta)
sub.g=(1+phi0+phir+phidelta*delta)*R*T	#R*cur_temp*
sub.g0=(1+phi0)*R*T
sub.dg=(phir+phidelta*delta)	#deriver RT
#
sub.f=(phi0+phir)*R*T
sub.f0=phi0*R*T
sub.df=phir

#Energy
sub.u=tau*(phi0tau+phitau)*R*T
sub.du=tau*(phitau)
#
sub.h=(1+tau*(phi0tau+phitau)+delta*phidelta)*R*T
sub.dh=1+tau*(phitau)+delta*phidelta;
#speed of sound
sub.w=sqrt(R/mm*T*(1+2*delta*phidelta+delta^2*phidelta2-(1+delta*phidelta-delta*tau*dphideltadtau)^2/(tau^2*(phi0tau2+phitau2))))	#0.018 перевод в килограмы
#
sub.cv=-tau^2*(phi0tau2+phitau2)*R
sub.cp=(-tau^2*(phi0tau2+phitau2)+(1+delta*phidelta-delta*tau*dphideltadtau)^2/(1+2*delta*phidelta+delta^2*phidelta2))*R
sub.dpdt=R/100.0*T*(1.0+2*delta*phidelta+delta^2*phidelta2)
#debug
# println("phi0 $(phi0)")
# println("phi0delta $(phi0delta)")
# println("phi0delta2 $(phi0delta2)")
# println("phi0tau $(phi0tau)")
# println("phi0tau2 $(phi0tau2)")
# println("phi0taudelta $(phi0taudelta)")
#
# println("phir $(phir)")
# println("phidelta $(phidelta)")
# println("phidelta2 $(phidelta2)")
# println("phitau $(phitau)")
# println("phitau2 $(phitau2)")
# println("dphideltadtau $(dphideltadtau)")
#enddebug
return sub
end 


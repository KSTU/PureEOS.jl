"""
Calculate ethane critical properties
T --- temperature in Kelvin
p --- pressure in bar
ro --- density in mol/l
"""
function LJ2Critical()
	sub=newCritical()
	sub.T=1.32
	sub.p=0.13006	#MPa to bar
	sub.ro=0.31	#mol/l #322 kg/m^3
	return sub
end

"""
Calculate ethane properties at triple point
T --- temperature in Kelvin
p --- pressure in bar
vRo --- vapor density in mol/l
lRo --- liquid density in mol/l
"""
function LJ2Triple()
	sub=newTriple()
	sub.T=90.368
	sub.p=1.14*10^-5	#Pa to bar
#	sub.lRo=999.793/18.015268	#kg/m^3 to mol/l
#	sub.vRo=0.00485458/18.015268	#kg/m^3s to mol/l
	return sub
end

"""
calculate water properties at saturation line
Equation of State for the Lennard-Jones Fluid
Article in Journal of Physical and Chemical Reference Data · June 2016
DOI: 10.1063/1.4945000
T --- temperature
p --- pressure
vRo --- vapor density in mol/l
lRo --- liquid density in mol/l
"""
function ethaneSaturate(T)
a=[-6.48647577, 1.47010078, -1.66261122, 3.57898378, -4.79105705]
b=[1.56138026, -0.381552776, 0.0785372040, 0.0370315089]
c=[-1.89879145, -3.65459262, 0.850562745, 0.363965487, -1.50005943, -2.26690389]
sub=newSaturate()
crit=ethaneCritical()
if (T<crit.T)
	v=1.0-T/crit.T
	sub.p=crit.p*exp(crit.T/T*(a[1]*v+a[2]*v^1.5+a[3]*v^2.5+a[4]*v^3.5+a[5]*v^4.0))
	sub.lRo=crit.ro*exp(b[1]*v^0.329+b[2]*v^(4.0/6.0)+b[3]*v^(8.0/6.0)+b[4]*v^(19.0/6.0))
	sub.vRo=crit.ro*exp(crit.T/T*(c[1]*v^(0.346)+c[2]*v^(5.0/6.0)+c[3]*v+c[4]*v^2.0+c[5]*v^3.0+c[6]*v^5))
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
function LJ2Singlephase(T,ro)

sub=newSingle()
crit=LJ2Critical()
delta=ro/crit.ro
tau=crit.T/T
R=1 #8.3144598
mm=0.03006904

nOne=6
nTwo=6
nThree=11


c=[0.00, 0.00, 0.00, 0.00, 0.00,
   0.00, 1.00, 2.00, 1.00, 2.00,	#10
   2.00, 1.00, 0.00, 0.00, 0.00,
   0.00, 0.00, 0.00, 0.00, 0.00,	#20
   0.00, 0.00, 0.00]

d=[4.00, 1.00, 1.00, 2.00, 2.00,
   3.00, 5.00, 2.00, 2.00, 3.00,	#10
   1.00, 1.00, 1.00, 1.00, 2.00,
   3.00, 3.00, 2.00, 1.00, 2.00,	#20
   3.00, 1.00, 1.00]

t=[ 1.000, 0.320, 0.505, 0.672, 0.843,
	0.898, 1.294, 2.590, 1.786, 2.770,	#10
	1.786, 1.205, 2.830, 2.548, 4.650,
	1.385, 1.460, 1.351, 0.660, 1.496,	#20
	1.830, 1.616, 4.970]

n=[0.52080730e-2, 0.21862520e1,
    -0.21610160e1, 0.14527000e1,
    -0.20417920e1, 0.18695286e0,
    -0.90988445e-1, -0.49745610e0,
    0.10901431e0, -0.80055922e0,
    -0.56883900e0, -0.62086250e0,
    -0.14667177e1, 0.18914690e1,
    -0.13837010e0, -0.38696450e0,
    0.12657020e0, 0.60578100e0,
    0.11791890e1, -0.47732679e0,
    -0.99218575e1, -0.57479320e0,
    0.37729230e-2]

al=[0.0, 0.0, 0.0, 0.0, 0.0,
	0.0, 0.0, 0.0, 0.0, 0.0,	#10
	0.0, 0.0, 2.067, 1.522, 8.820,
	1.722, 0.679, 1.883, 3.925, 2.461,	#20
	28.20, 0.753, 0.820]

bt=[0.0, 0.0, 0.0, 0.0, 0.0,
	0.0, 0.0, 0.0, 0.0, 0.0,	#10
	0.0, 0.0, 0.625, 0.638, 3.910,
	0.156, 0.157, 0.153, 1.160, 1.730,	#20
	383.0, 0.112, 0.119]

gm=[0.0, 0.0, 0.0, 0.0, 0.0,
	0.0, 0.0, 0.0, 0.0, 0.0,	#10
	0.0, 0.0, 0.710, 0.860, 1.940,
	1.480, 1.490, 1.945, 3.020, 1.110,	#20
	1.170, 1.330, 0.240]

ep=[0.0, 0.0, 0.0, 0.0, 0.0,
	0.0, 0.0, 0.0, 0.0, 0.0,	#10
	0.0, 0.0, 0.2053, 0.4090, 0.6000,
	1.2030, 1.8290, 1.3970, 1.3900, 0.5390,	#20
	0.9340, 2.3690, 2.4300]

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
	i=i+1
	phidelta=phidelta+n[i]*exp(-delta^c[i])*(delta^(d[i]-1.0)*tau^t[i]*(d[i]-c[i]*delta^c[i]))
	phir=phir+n[i]*delta^(d[i])*tau^t[i]*exp(-delta^c[i])
	phitau=phitau+n[i]*t[i]*delta^d[i]*tau^(t[i]-1.0)*exp(-delta^c[i]);
	phidelta2=phidelta2+n[i]*exp(-delta^c[i])*(delta^(d[i]-2)*tau^(t[i])*((d[i]-c[i]*delta^c[i])*(d[i]-1.0-c[i]*delta^c[i])-c[i]^2*delta^c[i]))
	phitau2=phitau2+n[i]*t[i]*(t[i]-1.0)*delta^d[i]*tau^(t[i]-2)*exp(-delta^c[i])
	dphideltadtau=dphideltadtau+n[i]*t[i]*delta^(d[i]-1.0)*tau^(t[i]-1.0)*(d[i]-c[i]*delta^c[i])*exp(-delta^c[i])
end
for pn=1:nThree
	i=i+1
	phidelta=phidelta+n[i]*delta^d[i]*tau^t[i]*exp(-al[i]*(delta-ep[i])^2-bt[i]*(tau-gm[i])^2)*(d[i]/delta-2.0*al[i]*(delta-ep[i]))
	phir=phir+n[i]*delta^d[i]*tau^t[i]*exp(-al[i]*(delta-ep[i])^2-bt[i]*(tau-gm[i])^2)
	phitau=phitau+n[i]*delta^d[i]*tau^t[i]*exp(-al[i]*(delta-ep[i])^2-bt[i]*(tau-gm[i])^2)*(t[i]/tau-2*bt[i]*(tau-gm[i]))
	phidelta2=phidelta2+n[i]*tau^t[i]*exp(-al[i]*(delta-ep[i])^2-bt[i]*(tau-gm[i])^2)*(-2*al[i]*delta^d[i]+4*al[i]^2*delta^d[i]*(delta-ep[i])^2-4*d[i]*al[i]*delta^(d[i]-1)*(delta-ep[i])+d[i]*(d[i]-1)*delta^(d[i]-2))
	phitau2=phitau2+n[i]*delta^d[i]*tau^t[i]*exp(-al[i]*(delta-ep[i])^2-bt[i]*(tau-gm[i])^2)*((t[i]/tau-2*bt[i]*(tau-gm[i]))^2-t[i]/tau^2-2*bt[i])
	dphideltadtau=dphideltadtau+n[i]*delta^d[i]*tau^t[i]*exp(-al[i]*(delta-ep[i])^2-bt[i]*(tau-gm[i])^2)*(d[i]/delta -2*al[i]*(delta-ep[i]))*(t[i]/tau-2*bt[i]*(tau-gm[i]))
end


sub.T=T
sub.ro=ro
sub.p=(1.0+phidelta*delta)*ro*T	#bar
sub.z=1.0+phidelta*delta	#
#ideal part
n0=[6.262265814, 
-1.515151515,
1.5,
0.0,
0.0,
0.0,
0.0,
0.0]
gamma0=[0.0, 0.0, 0.0, 1.4091052332, 4.0099170712, 6.5967098342, 13.9798102659,1.0]

#entropy
phi0=log(delta)+n0[1]+n0[2]*tau+n0[3]*log(tau)

phi0tau=n0[2]+n0[3]/tau
phi0delta=1/delta
phi0delta2=-1/delta^2
phi0tau2=-n0[3]/tau^2
phi0taudelta=0

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
sub.g=(1+phi0+phir+phidelta*delta)*T	#R*cur_temp*
sub.g0=(1+phi0)*T
sub.dg=(phir+phidelta*delta)	#deriver RT
#
sub.f=(phi0+phir)*T
sub.f0=phi0*T
sub.df=phir

#Energy
sub.u=tau*(phi0tau+phitau)*T
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



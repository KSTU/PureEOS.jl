"""
Calculate water critical properties
T --- temperature in Kelvin
p --- pressure in bar
ro --- density in mol/l
doi: 10.1063/1.555786
"""
function methanolCritical()
	sub=newCritical()
	sub.T=512.50    #K
	sub.p=80.9464   #bar
	sub.ro= 8.40    #mol/l #322 kg/m^3
	return sub
end

"""
Calculate water properties at triple point
T --- temperature in Kelvin
p --- pressure in bar
vRo --- vapor density in mol/l
lRo --- liquid density in mol/l
doi: 10.1063/1.555786
"""
function methanolTriple()
	sub=newTriple()
	sub.T=175.59
	sub.p=1.835*10.0^-6	#bar
	sub.lRo=28.226	#kg/m^3 to mol/l
	sub.vRo=1.271*10.0^(-7)	#kg/m^3s to mol/l
	return sub
end

"""
calculate methanol properties at saturation line
T --- temperature in Kelvin
p --- pressure in bar
vRo --- vapor density in mol/l
lRo --- liquid density in mol/l
doi: 10.1063/1.555786
"""
function waterSaturate(T)
a=[-10.75284879, 16.758206642, -3.603424623, 4.373231941, -2.381377449, 4.572198698]    #pressure
b=[2.51730906, -2.46669419, 3.06681768, -1.32507680]    #liquid density
c=[10.61966850, -2.55668202, 3.81845421, 4.79556752]    #vapor density
sub=newSaturate()
crit=methanolCritical()

if (T<crit.T)
	v=1.0-T/crit.T
	sub.p=exp(a[1]/v + a[2] + a[3]*v + a[4]*v^2 + a[5]*v^3 + a[6]*(1-v)^1.7)
	sub.lRo=crit.ro*(1+b[1]*(1-v)^0.35 + b[2]*(v-1) + b[3]*(v^2-1) + b[4]*(v^3-1))
	sub.vRo=crit.ro*exp(c[1]*(1 - 1/v) + c[2]*(1-v)^0.35 + c[3]*(1-v) + c[4]*(1-v)^2)
	sub.T=T
else
	sub.p=-1
	sub.lRo=-1
	sub.vRo=-1
	sub.T=-1
end
return sub
end




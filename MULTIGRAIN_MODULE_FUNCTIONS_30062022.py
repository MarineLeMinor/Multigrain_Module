# -*-coding:Latin-1 -*

# C:\Users\mleminor\AppData\Local\Programs\Python

###################################################################
### FUNCTIONS - MULTI-GRAIN/MULTI-MODE SEDIMENT TRANSPORT MODEL ###
###################################################################


# Import modules
import numpy
import math


# Median diameter
def calculateD50(dBins,stratTSLp):

	dBinsCenters=dBins[:-1]+0.5*(dBins[1:]-dBins[:-1])	
	
	dist=stratTSLp/numpy.sum(stratTSLp)	
	cumdist=numpy.cumsum(dist)
	diffdist=cumdist-0.5

	#d50=dBinsCenters[numpy.abs(diffdist)==numpy.amin(numpy.abs(diffdist))]

	##d50=numpy.sum(dBinsCenters*stratTSLp)/numpy.sum(stratTSLp)

	#if len(stratTSLp[stratTSLp!=0])==1:
	#	d50=dBinsCenters[stratTSLp!=0]

	#d50=d50[0]

	ind=numpy.arange(0,len(dBinsCenters),1)

	if len(stratTSLp[stratTSLp!=0])==1:
		d50=dBinsCenters[stratTSLp!=0]
		d50=d50[0]

	else:

		lim=diffdist[numpy.abs(diffdist)==numpy.amin(numpy.abs(diffdist))]

		if lim[0]<0:
			#indL=ind[numpy.abs(diffdist)==numpy.amin(numpy.abs(diffdist))]
			indL=ind[diffdist==lim[0]]
			indL=indL[-1]

			if indL!=ind[-1]:
				indR=indL+1
				d50=dBinsCenters[indL]+(0.5-cumdist[indL])/(cumdist[indR]-cumdist[indL])*(dBinsCenters[indR]-dBinsCenters[indL])

			else:
				print('problem',dBinsCenters[indL],(0.5-cumdist[indL]),(1-cumdist[indL]),(dBins[-1]-dBinsCenters[indL]))
				print(indL,lim[0],cumdist)
				d50=dBinsCenters[indL]+(0.5-cumdist[indL])/(1-cumdist[indL])*(dBins[-1]-dBinsCenters[indL])

		elif lim[0]>0:
			#indR=ind[numpy.abs(diffdist)==numpy.amin(numpy.abs(diffdist))]
			indR=ind[diffdist==lim[0]]
			indR=indR[0]

			if indR!=ind[0]:
				indL=indR-1
				d50=dBinsCenters[indL]+(0.5-cumdist[indL])/(cumdist[indR]-cumdist[indL])*(dBinsCenters[indR]-dBinsCenters[indL])

			else:
				d50=dBins[0]+(0.5-0)/(cumdist[indR]-0)*(dBinsCenters[indR]-dBins[0])

		elif lim[0]==0:
			d50=dBinsCenters[numpy.abs(diffdist)==numpy.amin(numpy.abs(diffdist))]


	return d50


# Diameter of 90-percentile
def calculateD90(dBinsCenters,stratTSLp):

	dist=stratTSLp/numpy.sum(stratTSLp)	
	cumdist=numpy.cumsum(dist)
	diffdist=cumdist-0.9

	d90=dBinsCenters[numpy.abs(diffdist)==numpy.amin(numpy.abs(diffdist))]

	if len(stratTSLp[stratTSLp!=0])==1:
		d90=dBinsCenters[stratTSLp!=0]

	d90=d90[0]

	return d90


# Sediment load in precipiton after abrasion
def calculatePrecipitonLoadAbrasion(dBins,dBinsCenters,hsp,dCoeff,Lp):

	d=dBinsCenters
	d_reduced=d*numpy.exp(-Lp/dCoeff)
	hsp_reduced=hsp*numpy.exp(-3*Lp/dCoeff)

	hsp_updated=numpy.zeros(len(dBinsCenters))

	for i in range(0,len(dBinsCenters)):

		ind=d_reduced[(d_reduced>dBins[i]) & (d_reduced<dBins[i+1])]

		if len(ind)!=0:

			hsp_updated[i]=hsp_updated[i]+numpy.sum(hsp_reduced[(d_reduced>dBins[i]) & (d_reduced<dBins[i+1])])

	return hsp_updated


# Shields number theta [-] (Shields, 1936)
def calculateShieldsNumber(rhos,rho,g,d,tau):

	theta=tau/((rhos-rho)*g*d)

	return theta


# Critical shear stress thetac [-] (Soulsby, 1997)
def calculateShieldsNumberSoulsby(rhos,rho,g,d,nu):

	dstar=((rhos-rho)*g*(d**3)/(rho*nu**2))**(1/3)

	thetac=0.3/(1+1.2*dstar)+0.055*(1-math.exp(-0.02*dstar))

	return thetac


# Shear stress tau [Pa]
def calculateShearStress(rhos,rho,g,d,theta):

	tau=theta*(rhos-rho)*g*d

	return tau


# Hiding-exposure factor zeta [-]
def calculateHEfactor(d,d50,gamma):

	zeta=(d/d50)**(1-gamma)

	return zeta


# Settling velocity Vs [m/s]
def calculateSettlingVelocity(rhos,rho,g,d,nu):

	# Specific density R [-]
	R=(rhos-rho)/rho

	# Settling velocity Vs [m/s] (Ferguson and Church, 2004)
	ws=R*g*(d**2)/(18*nu+(0.75*0.4*R*g*(d**3))**(1/2))
	
	return ws


# Saltation height hsalt [m]
def calculateHsalt(d,tauc,tau,h):

	# Transport stage Tstar [-]
	Tstar=tau/tauc-1

	# Saltation height hsalt [m] (Auel et al., 2017)
	hsalt=min(0.6*d+0.025*d*Tstar,h)

	return hsalt


# Rouse number P [-] (Rouse, 1937)
def calculateRouseNumber(ws,K,beta,ustar):

	P=ws/(K*beta*ustar)

	return P


# Vertical gradient of sediment distribution
def calculateGradient(h,hsalt,ustar,z0,P):

	# Reference height for the Rouse profile zref [m]
	zref=z0

	# Vertical increment
	z=numpy.linspace(z0,h,1000)
	dz=z[1]-z[0]
	zz=z[:-1]+0.5*dz

	# Logarithmic velocity profile
	u=numpy.log(zz/z0)

	# Zero velocity below reference height z0 [m] of velocity profile
	u[zz<=z0]=0

	# Integration of the velocity profile from the saltation height to the water depth: water discharge per unit width
	A=numpy.sum(u)*dz

	# Rouse profile
	c=numpy.power((h/zz-1)*zref/(h-zref),P)
	

	# Integration of the sediment load from the saltation height to the water depth: sediment discharge per unit width
	cu=numpy.multiply(c,u)
	B=numpy.sum(cu)*dz

	# Calculation of the gradient of sediment distribution
	r0=A/B

	return r0


# Characteristic transport height hstar [m]
def calculateHs(hsalt,h,r0):

	hs=hsalt+(h-hsalt)/r0

	return hs


# Characteristic transport velocity vstar [m/s]
def calculateVs(theta,thetac,ustar,K,hs,z0,u,P,rhos,rho,g,d):
	
	vlayer=ustar/K*(numpy.log(hs/z0)-1+z0/hs)

	if P<2.5:
		vs=vlayer

	else: 
		R=rhos/rho-1
		Tstar=theta/thetac-1
		vsalt=1.46*math.sqrt(R*g*d)*(Tstar**0.5)
		vs=min(vsalt,vlayer)

	vs=min(vs,u)

	return vs


# Characteristic transport length xi [m]
def calculateTransportLength(hs,vs,ws):

	xi=hs*vs/ws

	return xi


# Erosion coefficient ke [m2*s/kg]
def calculateErosionCoefficient(rhos,rho,g,d,ws):

	ce=18.82

	ke=math.pi/6*ce/(rhos*ws)

	return ke


# Erosion rate dote [m/s]
def calculateErosionRate(ke,tauc,tau,F):

	dote=F*ke*(tau-tauc)

	return dote


# Transport rate per unit width qs [m2/s]
def calculateQs(xi,E):

	qs=xi*E

	return qs


# Dimensionless transport capacity qsstar per unit width [-]
def calculateEinsteinNumber(rhos,rho,g,d,qs):

	qsstar=qs/(((rhos-rho)/rho*g*(d**3))**(1/2))

	return qsstar
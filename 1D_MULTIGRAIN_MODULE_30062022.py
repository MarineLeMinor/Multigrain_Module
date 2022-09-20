# -*-coding:Latin-1 -*

# C:\Users\mleminor\AppData\Local\Programs\Python

###############################################################################
###                   MULTI-GRAIN SEDIMENT TRANSPORT - 1D                   ###
###############################################################################

# Marine LE MINOR
# Geosciences Rennes, University of Rennes 1, CNRS, UMR 6118, Rennes, France
# 30 June 2022 


# Import modules
import numpy    
import math
import decimal
import csv
from matplotlib import pyplot
from matplotlib import cm
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy import interpolate, stats
import matplotlib
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

from MULTIGRAIN_MODULE_FUNCTIONS import *

##############################################################################
# CONSTANTS

# Gravitational constant [m/s2]
g=9.81

# Water density [kg/m3]
rho=1000
# Water viscosity [m2/s]
nu=1e-6

# Sediment density [kg/m3]
rhos=2650
# Sediment porosity [-]
porosity=0.6

# Van Karman constant [-]
K=0.41

# Exponent for Manning equation
alpha=2/3
# Factor for Manning equation
C=25
# Manning coefficient
n=1/C

# Threshold to merge two adjacent layers based on grain size distribution
RMSE_lim=0

# Size reduction coefficient in response to abrasion [m]
dCoeff=1e5


# Number of grain size classes, i.e., bins in the pdf
Nd=2#2
# Bins boundary in the pdf [m]
dBins=numpy.array([1e-4,1e-3,1e-2])#numpy.logspace(-5,-1,num=Nd+1)
#dBins=numpy.array([1e-4,1e-3])
#dBins=numpy.array([1e-3,1e-2])
# Bins centers [m]
dBinsCenters=dBins[:-1]+0.5*(dBins[1:]-dBins[:-1])

#print(calculateErosionCoefficient(rhos,rho,g,dBinsCenters[0],calculateSettlingVelocity(rhos,rho,g,dBinsCenters[0],nu)))
print(calculateErosionCoefficient(rhos,rho,g,dBinsCenters[1],calculateSettlingVelocity(rhos,rho,g,dBinsCenters[1],nu)))

print('Grain size classes',dBinsCenters)
# Hiding-exposure exponent
# Grain idependence gamma=0
# Equal mobility gamma=1
gamma=0

# Diffusivity ratio [-]
beta=1



##############################################################################
# GRID

# Channel length [m]
Lriver=10000
Llake=10000
L=Lriver+Llake # River length + Lake length
# Pixel length [m]
dx=250
# Number of pixels in longitudinal direction [-]
Nx=int(L/dx)
Nxriver=int(Lriver/dx)
print('Nxriver',Nxriver)
# Longitudinal pixel coordinates [m]
x=numpy.arange(dx/2,Nx*dx+dx/2,dx)

# Channel width [m]
W=50

# Number of pixels in transversal direction [-]
Ny=1
# Pixel width [-]
dy=int(W/Ny)
# Transversal pixel coordinates [m]
y=dy/2*numpy.ones(Nx)

# Pixel area [m2]
A=dx*dy

# Number of cells [-]
Ncells=Nx*Ny



alpha=0
minThick=4e-2

##############################################################################
# TIME

# Starting time [s]
t_start=0
# Time increment between two successive precipitons "launched" on the grid [s]
dt=1
# Ending time [s]
t_end=1e6

# Number of time steps ~ number of precipitons "launched" on the grid [-]
Nt=int((t_end-t_start)/dt)



##############################################################################
# INITIAL CONDITIONS

# Initial sediment load in the precipiton, i.e., upstream sediment discharge [m3/s]
QsInit=numpy.zeros(Nd)
QsInit[:]=1e-2
# Initial water discharge [m3/s]
Q=100

#6: 1e-3 100
#7: 1e-3 75
#8: 1e-1 75
#9: 1e-5 75
#10: 1e-3 75 s1=0.0003
#11: 1e-3 50 s1=0.0003
#12: 1e-3 50 s1=0.0003 save tau


#CaseM25: 1e-3 25 s1=0.0003
#CaseM50: 1e-3 50 s1=0.0003

#CaseM75: 1e-3 75 s1=0.0003
#CaseS75: 1e-3 75 s1=0.0003
#CaseS50: 1e-3 50 s1=0.0003
#CaseS25: 1e-3 25 s1=0.0003

# Initial bed slope [-]
Spmean=0.0001
# Initial bed topography [m]
z_down=0
z=z_down*numpy.ones(Nx)

for i in range(0,Nx):
	z[i]=z_down+Spmean*(x[-1]-x[i])

# Initial bed slope [-]
Spmean=0.001
# Initial bed topography [m]
z_down=z[int(Nxriver/2)]

for i in range(0,int(Nxriver/2)):
	z[i]=z_down+Spmean*(x[int(Nxriver/2)]-x[i])

for i in range(Nxriver,Nx):
	z[i]=z[Nxriver-1]*math.exp(-(x[i]-Lriver)/dx)

pyplot.figure()
pyplot.plot(x,z)
pyplot.show()
pyplot.close()

Zp=z.copy()

# Initial basin stratigraphy [m] - porosity not considered
strat=numpy.array([numpy.zeros(Nd)])
#strat[0,8:24]=0.04

#strat=numpy.zeros((100,Nd))
#for i in range(0,len(strat)):
#	strat[i,8:24]=0.0004

basin=[]
thickinit=numpy.zeros(Nx)

for i in range(0,Nx): 
	if i<Nxriver:
		strat=numpy.zeros((2,Nd))
		for j in range(0,len(strat)):
			strat[j,:]=0.8
			#strat[j,0]=1.6
			#strat[j,1]=0
	else:
		strat=numpy.array([numpy.zeros(Nd)])

	basin.append(strat)
	thickinit[i]=numpy.sum(basin[i])/(1-porosity)

#print('init',z,strat,thickinit)

# Age record

chrono=[]

for i in range(0,Nx): 
	#chrono.append(numpy.zeros(1))
	chrono.append(numpy.zeros(2))



# Grain-size specific parameters that will not change during the simulation
# Settling velocity [m/s]
Vs=numpy.zeros(Nd)
# Erosion coefficient [m2s/kg]
ke=numpy.zeros(Nd)

for i in range(0,Nd):
	Vs[i]=calculateSettlingVelocity(rhos,rho,g,dBinsCenters[i],nu)
	ke[i]=calculateErosionCoefficient(rhos,rho,g,dBinsCenters[i],Vs[i])


###############################################################################
# ITERATIONS

countDeletions=0

for t in range(1,Nt+1):

	#print(Nt-t)

	Zpmoinsun=Zp.copy()


	d50p=numpy.zeros(Nx)
	taup=numpy.zeros(Nx)
	taucCoarse=numpy.zeros(Nx)
	taucFine=numpy.zeros(Nx)
	dotdCoarse=numpy.zeros(Nx)
	dotdFine=numpy.zeros(Nx)
	doteCoarse=numpy.zeros(Nx)
	doteFine=numpy.zeros(Nx)
	FCoarse=numpy.zeros(Nx)
	xiCoarse=numpy.zeros(Nx)
	xiFine=numpy.zeros(Nx)

	#Q=50

	#if t<Nt/4:
	#	Q=25

	#elif t>=Nt/4 and t<2*Nt/4:
	#	Q=50

	#elif t>=Nt/2 and t<3*Nt/4:
	#	Q=75

	#if t>Nt/2:
	#	QsInit=numpy.zeros(Nd)
	#	QsInit[14:18]=0e-4

	# Initial heights of sediment on the precipiton [m]
	hsp=QsInit*dt/A

	for p in range(0,Nx):
		if p<Nxriver:

			# Calculate the slope.

			#Sp=-(Zp[p+1]-Zp[p])/dx

			if p!=Nxriver-1:
				Sp=-(Zp[p+1]-Zp[p])/dx

			else:
				Sp=-(Zp[p]-Zpmoinsun[p-1])/dx

			# Calculate the water depth.
			Hp=((n*Q/W)**(3/5))*(Sp**(-3/10))


			# Calculate the shear stress.
			tau=rho*g*Hp*Sp

			taup[p]=tau

			# Calculate shear velocity [m/s].
			ustar=(tau/rho)**0.5

			# Initiate transport length [m]
			xi=numpy.zeros(Nd)
			tauc=numpy.zeros(Nd)

			for i in range(0,Nd):
				tauc[i]=calculateShearStress(rhos,rho,g,dBinsCenters[i],calculateShieldsNumberSoulsby(rhos,rho,g,dBinsCenters[i],nu))

			# Calculate the excess of shear stress to identify grain sizes present in the coming precipiton and in the surficial stratigraphic layer that are in motion.
			excessTau=tau-tauc

			for i in range(0,Nd):
				# For mobile grain sizes:
				if excessTau[i]>0:

					# Calculate the bed characteristic diameter.
					if numpy.sum(numpy.sum(basin[p]))!=0:
						bed=numpy.copy(basin[p])
						d90=calculateD90(dBinsCenters,bed[0])

					else:
						d90=1e-3#dBinsCenters[-16]

					# Calculate the bed roughness.
					#z0=3*dBinsCenters[-16]/30
					z0=3*d90/30

					# Calculate the saltation height.
					hsalt=calculateHsalt(dBinsCenters[i],tauc[i],tau)

					# Calculate the Rouse number.
					P=calculateRouseNumber(Vs[i],K,beta,ustar)

					# Calculate the vertical gradient of distribution.
					r0=calculateGradient(hsalt,Hp,ustar,z0,P)

					# Calculate the characteristic transport height.
					hstar=calculateHstar(hsalt,Hp,r0)

					# Calculate the characteristic transport velocity.
					vstar=calculateVstar(hstar,ustar,z0,K)

					# Calculate the transport length.
					xi[i]=calculateTransportLength(hstar,vstar,Vs[i])

		# Select pixel stratigraphy
		stratp=numpy.copy(basin[p])

		agep=numpy.copy(chrono[p])

		dt_temp=0
		dhsp_tot=numpy.zeros(Nd)

		while dt_temp<dt:

			dhsp_eff=numpy.zeros(Nd)
			dhsp=numpy.zeros(Nd)

			# Initiate erosion flux [m/s]
			dote=numpy.zeros(Nd)

			if p<Nxriver:

				# Make sure that the topmost storage layer has a thickness of at least the diameter of the finest grains it contains.
				if len(stratp)>=2 and numpy.sum(stratp[0,:])<numpy.amin(dBinsCenters[stratp[0,:]!=0]):
					portion=numpy.amin(dBinsCenters[stratp[0,:]!=0])-numpy.sum(stratp[0,:])/numpy.sum(stratp[:,1])
					stratp[0,:]=stratp[0,:]+portion*stratp[1,:]
					stratp[1,:]=stratp[1,:]-portion*stratp[1,:]
					#print('MOD TSL')

				# Select topmost storage layer.
				stratTSLp=stratp[0,:]

				if numpy.sum(stratTSLp)>0:

					# Calculate the characteristic diameters, i.e., d50 and d90.
					d50=calculateD50(dBinsCenters,stratTSLp)
					d90=calculateD90(dBinsCenters,stratTSLp)

					# Calculate the threshold of motion of the median grain size.
					theta50=calculateShieldsNumberSoulsby(rhos,rho,g,d50,nu)
					tauc50=calculateShearStress(rhos,rho,g,d50,theta50)

					# Calculate the hiding-exposure factor.
					zeta=numpy.zeros(Nd)

					for i in range(Nd):
						gamma=1-0.67/(1+math.exp(1.5-dBinsCenters[i]/d50))
						zeta[i]=calculateHEfactor(dBinsCenters[i],d50,gamma)

					# Calculate the threshold of motion  for grain sizes present in the coming precipiton and in the surficial stratigraphic layer.
					tauc_corr=zeta*tauc50

					if dt_temp==0:
						taucCoarse[p]=tauc_corr[1]
						taucFine[p]=tauc_corr[0]

					# Calculate the excess of shear stress to identify grain sizes present in the coming precipiton and in the surficial stratigraphic layer that are in motion.
					excessTau_corr=tau-tauc_corr

					# Calculate the fraction of each grain size, i.e., ratio of heights.
					F=stratTSLp/numpy.sum(stratTSLp)

					for i in range(0,Nd):
						# For grain sizes that are available in the topmost storage layer:
						if excessTau_corr[i]>0:

							# Calculate the erosion flux.
							dote[i]=calculateErosionFlux(ke[i],tauc_corr[i],tau,F[i])

			if dt_temp==0:
				if xi[1]!=0:
					dotdCoarse[p]=-hsp[1]/xi[1]*(A/dy/dt)
				if xi[0]!=0:
					dotdFine[p]=-hsp[0]/xi[0]*(A/dy/dt)

				if F[1]!=0:
					doteCoarse[p]=dote[1]/F[1]
				if F[0]!=0:
					doteFine[p]=dote[0]/F[0]
				FCoarse[p]=F[1]
				xiCoarse[p]=xi[1]
				xiFine[p]=xi[0]

			# Calculate the sediment load variation in the precipiton.
			#dhsp=(xi*dote*(dt-dt_temp)/dx-hsp)*(1-numpy.exp(-dx/xi))

			dhsp[xi!=0]=(xi[xi!=0]*dote[xi!=0]*(dt-dt_temp)/dx-hsp[xi!=0])*(1-numpy.exp(-dx/xi[xi!=0]))
			dhsp[xi==0]=-hsp[xi==0]

			# Identify if grain sizes are deposited or eroded.
			# In case of deposition:
			if len(dhsp[dhsp<0])>0:

				#if p==0:
				#	print(t,'DEPOT')
				#if p==0 and Zpmoinsun[p]!=z[p]+numpy.sum(numpy.sum(stratp))/(1-porosity):
				#	print(t,p,Zpmoinsun[p],z[p]+numpy.sum(numpy.sum(stratp))/(1-porosity),z[p]+numpy.sum(numpy.sum(basin[p]))/(1-porosity),'DEPOT')

				# Set the change of sediment load in precipiton for the sizes eroded to zero.
				dhsp[dhsp>0]=0.0

				# Calculate the sediment load variation in the precipiton due its interaction with the topmost storage layer at the pixel p.
				dhsp_eff=dhsp

				# Special case with the "virtual" lake
				if p>=Nxriver:

					if Zp[p-1]>Zpmoinsun[p]:

						if(numpy.sum(-dhsp_eff)>(Zp[p-1]-Zpmoinsun[p])*(1-porosity)):
							dhsp_eff=(Zp[p-1]-Zpmoinsun[p])*(1-porosity)*dhsp_eff/numpy.sum(-dhsp_eff)
				
					else:
						dhsp_eff=numpy.zeros(Nd)

				#if numpy.sum(numpy.sum(stratp))==0:
				#	stratp=numpy.delete(stratp,0,0)

				## Create a layer made of this deposit and keep the surficial stratigraphic layer unchanged.
				#stratp=numpy.insert(stratp,[0],-dhsp_eff,axis=0)


				strat_below=numpy.zeros(Nd)


				if numpy.sum(numpy.sum(stratp[0,:]))>=minThick:
					agep=numpy.insert(agep,[0],t_start+t*dt,axis=0)

				elif numpy.sum(numpy.sum(stratp[0,:]))==0:
					agep[0]=t_start+t*dt

				if numpy.sum(stratp[0,:])<minThick:
					strat_below=numpy.copy(stratp[0,:])
					stratp=numpy.delete(stratp,0,0)

				stratp=numpy.insert(stratp,[0],-dhsp_eff+strat_below,axis=0)


				# Set the time elapsed since the precipiton is on the pixel p
				dt_temp=dt



			# In case of erosion only:
			elif (len(dhsp[dhsp>0])>0) and (len(dhsp[dhsp<0])==0):

				# Determine whether erosion of some grain sizes exceeds the heights available in the topmost storage layer.
				ratio=numpy.divide(dhsp[stratTSLp>0],stratTSLp[stratTSLp>0])


				# Identify the minimum portion of the topmost storage layer to erode.
				ratio_min=numpy.amin(ratio)

				# If all the grain sizes of the topmost storage layer are eroded and the heights to remove exceed the ones available.
				if ratio_min>=1:

					# Calculate the sediment load variation in the precipiton since the precipiton is on the pixel p.
					dhsp_eff=numpy.copy(stratTSLp)

					# Update the stratigraphy.
					# If the pixel stratigraphy is only made of one layer:
					if len(stratp)==1:

						# Set all heights to zero.
						stratTSLp=numpy.zeros(Nd)
						stratp[0,:]=stratTSLp

						agep[0]=0.

						# Set the time elapsed since the precipiton is on the pixel p.
						dt_temp=dt

					# If the pixel stratigraphy is made of at least two layers:
					else:

						countDeletions=countDeletions+1

						# Delete the topmost storage layer.
						stratp=numpy.delete(stratp,0,0)

						agep=numpy.delete(agep,0,0)

						# Calculate the longest time that was needed to erode one grain size from this layer.
						dt_layer=numpy.amax(-(dt-dt_temp)/dx*numpy.multiply(xi[stratTSLp>0],numpy.log(1-numpy.divide(stratTSLp[stratTSLp>0],numpy.multiply(xi[stratTSLp>0],dote[stratTSLp>0])*(dt-dt_temp)/dx-hsp[stratTSLp>0]))))

						# Set the time elapsed since the precipiton is on the pixel p.
						dt_temp=dt_temp+dt_layer

				# If at least one grain sizes of the topmost storage layer is eroded and the heights to remove does NOT exceed the one available.
				else:

					# Identify the finest grain size among the ones thave a ratio equal to the minimum one.
					d_eroded=dBinsCenters[stratTSLp>0]
					dlim=numpy.amin(d_eroded[ratio==ratio_min])

					if len(dlim[dlim!=0])!=1:
						dlim=dlim[0]

					# Calculate the additional thickness of the topmost storage layer where erosion can take place.
					hlim=alpha*dlim

					# Calculate the sediment load variation in the precipiton based on the limitng grain size.
					dhsp_lim=ratio_min*stratTSLp

					# Add the amount available on the thickness where erosion can occur based on the limiting grain size.
					dhsp_lim[hlim>=0.5*dBinsCenters]=dhsp_lim[hlim>=0.5*dBinsCenters]+(1-porosity)*F[hlim>=0.5*dBinsCenters]*alpha*dlim

					# Calculate the sediment load variation in the precipiton since the precipiton is on the pixel p.
					dhsp_eff=numpy.minimum(stratp[0,:],dhsp,dhsp_lim)

					# Update the topmost storage layer.
					#stratTSLp=stratTSLp-dhsp_eff
					#stratp[0,:]=stratTSLp

					stratTSLp=stratp[0,:]
					
					if numpy.amax(dhsp_eff[stratTSLp!=0]/stratTSLp[stratTSLp!=0])<1 and numpy.sum(stratTSLp)>=minThick:

						stratp=numpy.insert(stratp,[0],stratp[0,:],axis=0)
						agep=numpy.insert(agep,[0],agep[0],axis=0)

						stratTSLp[stratTSLp!=0]=numpy.amax(dhsp_eff[stratTSLp!=0]/stratTSLp[stratTSLp!=0])*stratTSLp[stratTSLp!=0]#(ratio_min+(1-porosity)*alpha*dlim/numpy.sum(stratp[0,:]))*stratp[0,:]
						stratp[0,:]=stratTSLp
						stratp[1,:]=stratp[1,:]-stratp[0,:]



					stratp[0,:]=stratp[0,:]-dhsp_eff

					stratTSLp=stratp[0,:]
					if len(stratTSLp[stratTSLp<0])!=0:
						#print(stratTSLp[stratTSLp<0])
						stratTSLp[stratTSLp<0]=0.0

					stratp[0,:]=stratTSLp


					if numpy.sum(stratp[0,:])==0 and len(stratp)>1:
						stratp=numpy.delete(stratp,0,0)
						agep=numpy.delete(agep,0,0)

					# Set the time elapsed since the precipiton is on the pixel p.
					dt_temp=dt


			# If no grains are eroded nor deposited:
			else:

				# Set the time elapsed since the precipiton is on the pixel p.
				dt_temp=dt

			# Updating the sediment load in precipiton
			hsp=hsp+dhsp_eff

			dhsp_tot=dhsp_tot+dhsp_eff

		#if p==0:
		#	print(t,p,tauc_corr,dhsp_tot,dote,ke,F,tau-tauc_corr)

		# Updating bed elevation
		Zp[p]=Zpmoinsun[p]-numpy.sum(dhsp_tot)/(1-porosity)


		if numpy.sum(hsp)>0:

			# Calculate the characteristic diameters, i.e., d50.
			d50p[p]=calculateD50(dBinsCenters,hsp)

		else:
			d50p[p]=None


		#if Zp[p]>Zp[p-1] and p>0:
		#	print(t,p,'PROBLEM - before',Zp[p]-Zp[p-1])
		##if p==0 and Zp[p]!=z[p]+numpy.sum(numpy.sum(stratp))/(1-porosity):
		##	print(t,numpy.sum(numpy.sum(basin[p]))-numpy.sum(dhsp_tot)==numpy.sum(numpy.sum(stratp)),numpy.sum(numpy.sum(basin[p]))-numpy.sum(dhsp_tot)-numpy.sum(numpy.sum(stratp)))


		## Issue: error propagation
		#if Zp[p]!=z[p]+numpy.sum(numpy.sum(stratp))/(1-porosity):
		#	#print('CORRECTED')
		#	#err=(Zp[p]-(z[p]+numpy.sum(numpy.sum(stratp))/(1-porosity)))*(1-porosity)
		#	# Correction:
		#	#print('err',t,p,err)
		#	Zp[p]=z[p]+numpy.sum(numpy.sum(stratp))/(1-porosity)

		#if Zp[p]>Zp[p-1] and p>0:
		#	print(t,p,'PROBLEM - mid',Zp[p]-Zp[p-1])

		#if Zp[p]>Zp[p-1] and p>0:
		#	err=Zp[p]-Zp[p-1]
		#	stratp[0]=stratp[0]*(1-err/numpy.sum(stratp[0]))
		#	Zp[p]=z[p]+numpy.sum(numpy.sum(stratp))/(1-porosity)

		#if Zp[p]>Zp[p-1] and p>0:
		#	print(t,p,'PROBLEM - after',Zp[p]-Zp[p-1])

		## Combine really thin topmost storage layer with second storage layer.
		#if (len(stratp)>1) and (numpy.sum(numpy.sum(stratp[0]))<1e-3):
		#	stratp[1]=stratp[1]+stratp[0]
		#	stratp=numpy.delete(stratp,0,0)

		# Merge layers of similar composition.
		# If there are at least two stratigraphic layers at pixel p:
		#if (numpy.sum(dhsp_tot)!=0) and len(stratp)>=2:
		if (numpy.sum(dhsp_tot)<0) and len(stratp)>=2 and (numpy.sum(numpy.sum(stratp[0,:]))>=minThick):

			#print(t,p,'merge',stratp,numpy.sum(stratp[0,:]),numpy.sum(stratp[1,:]))

			# Compare its grain size distribution to the one of the stratigraphic layer right below.
			RMSE=math.sqrt(numpy.sum(numpy.power(stratp[0,:]/numpy.sum(stratp[0,:])-stratp[1,:]/numpy.sum(stratp[1,:]),2))/Nd)

			# If the topmost storage layer and the layer below are significantly alike:
			if RMSE<=RMSE_lim/100:

				# Merge the topmost storage layer with the layer below.
				stratp[1,:]=stratp[0,:]+stratp[1,:]

				stratp=numpy.delete(stratp,0,0)
				agep=numpy.delete(agep,0,0)


		# Updating stratigraphy at pixel p in the basin.
		basin[p]=stratp

		chrono[p]=agep

		#if z[p]+numpy.sum(numpy.sum(basin[p]))!=Zp[p]:
			#print(t,p,z[p]+numpy.sum(numpy.sum(basin[p]))-Zp[p])




	if (t*dt)%1e3==0:


		#for p in range(0,Nx):

		#	if numpy.sum(numpy.sum(basin[p]))!=0 and len(basin[p])>2:
		#		strat=basin[p]
		#		stratNew=numpy.array([numpy.copy(strat[-1,:])])


		#		for i in range (1,len(strat)-1):

		#			if numpy.sum(strat[len(strat)-1-i,:])>minThick:
		#				stratNew=numpy.insert(stratNew,[0],strat[len(strat)-1-i,:],axis=0)

		#			else:
		#				stratNew[0,:]=stratNew[0,:]+strat[len(strat)-1-i,:]

		#		stratNew=numpy.insert(stratNew,[0],strat[0,:],axis=0)

		#		basin[p]=stratNew




		#print(countDeletions)

		filename='CaseM100_Qs01_'

		# SAVE DATA AS A VTK FILE
		directory='D:/PostDoc_Marine_Le_Minor/MANUSCRIPTS/Manuscript_1DStratigraphy/Codes_Python/'
		file=open(directory+filename+str(int(t*dt))+'.vtk','w')

		# Header
		file.write('# vtk DataFile Version 3.0 \n')
		file.write('Stratigraphic record \n')
		file.write('ASCII \n')
		file.write('\n')


		Nlayers=0
		voxels=[]

		ALLd50=[]
		ALLage=[]
		ALLFfine=[]

		ALLx=[]
		ALLz=[]
		d50surface=numpy.zeros(Nx)


		for p in range(0,Nx):

			strat=basin[p]
			agep=chrono[p]

			thickAbove=0

			if numpy.sum(numpy.sum(strat))!=0:
				ALLx=ALLx+[x[p].tolist()]
				ALLz=ALLz+[Zp[p].tolist()]

			if numpy.sum(numpy.sum(strat))==0:
				d50surface[p]=None

			count=0

			for i in range(0,len(strat)):

				d50_temp=0
				Ffine_temp=0
				thick_temp=0

				verts=numpy.zeros((8,3))



				if numpy.sum(strat[i,:])!=0 and numpy.sum(strat[i,:])>=1e-16:

					d50=calculateD50(dBinsCenters,strat[i,:])
					age=agep[i]
					Ffine=strat[i,0]/numpy.sum(strat[i,:])

					d50=(d50_temp*thick_temp+d50*numpy.sum(strat[i,:]))/(thick_temp+numpy.sum(strat[i,:]))
					ALLd50=ALLd50+[d50.tolist()]
					ALLage=ALLage+[age.tolist()]
					Ffine=(Ffine_temp*thick_temp+Ffine*numpy.sum(strat[i,:]))/(thick_temp+numpy.sum(strat[i,:]))
					ALLFfine=ALLFfine+[Ffine.tolist()]

					#if p==0 and i==0:
					#	print(strat[i,:],Ffine)

					thick=numpy.sum(strat[i,:])/(1-porosity)
					thick=round(thick,16)

					if i==len(strat)-1:
						verts=numpy.array([[x[p]-dx/2,0,z[p]-thickinit[p]],
							   [x[p]+dx/2,0,z[p]-thickinit[p]],
							   [x[p]-dx/2,W,z[p]-thickinit[p]],
							   [x[p]+dx/2,W,z[p]-thickinit[p]],
							   [x[p]-dx/2,0,Zp[p]-thickAbove],
							   [x[p]+dx/2,0,Zp[p]-thickAbove],
							   [x[p]-dx/2,W,Zp[p]-thickAbove],
							   [x[p]+dx/2,W,Zp[p]-thickAbove]])


					else:
						verts=numpy.array([[x[p]-dx/2,0,Zp[p]-thickAbove-thick],
							   [x[p]+dx/2,0,Zp[p]-thickAbove-thick],
							   [x[p]-dx/2,W,Zp[p]-thickAbove-thick],
							   [x[p]+dx/2,W,Zp[p]-thickAbove-thick],
							   [x[p]-dx/2,0,Zp[p]-thickAbove],
							   [x[p]+dx/2,0,Zp[p]-thickAbove],
							   [x[p]-dx/2,W,Zp[p]-thickAbove],
							   [x[p]+dx/2,W,Zp[p]-thickAbove]])


					thickAbove=thickAbove+thick
					d50_temp=0
					Ffine_temp=0
					thick_temp=0

					if Nlayers==0:
						voxels=verts
					else:
						voxels=numpy.append(voxels,verts,axis=0)

					if count==0:
						d50surface[p]=d50

					Nlayers=Nlayers+1
					count=count+1


				elif numpy.sum(strat[i,:])!=0 and numpy.sum(strat[i,:])<1e-16:
					d50_temp=(d50_temp*thick_temp+calculateD50(dBinsCenters,strat[i,:])*numpy.sum(strat[i,:]))/(thick_temp+numpy.sum(strat[i,:]))
					Ffine_temp=(Ffine_temp*thick_temp+strat[i,0]/numpy.sum(strat[i,:])*numpy.sum(strat[i,:]))/(thick_temp+numpy.sum(strat[i,:]))
					thick=numpy.sum(strat[i,:])/(1-porosity)
					thick_temp=thick_temp+thick
					thickAbove=thickAbove+thick

		voxels=numpy.unique(voxels,axis=0)



		file.write('DATASET UNSTRUCTURED_GRID \n')
		file.write('POINTS '+str(len(voxels))+' float \n')

		for i in range(0,len(voxels)):
			file.write(str(voxels[i,0])+' '+str(voxels[i,1])+' '+str(voxels[i,2])+' \n')

		file.write('\n')
		file.write('CELLS '+str(Nlayers)+' '+str(Nlayers*9)+' \n')



		for p in range(0,Nx):

			strat=basin[p]

			thickAbove=0

			for i in range(0,len(strat)):
				
				thick_temp=0

				count=0

				if numpy.sum(strat[i,:])!=0 and numpy.sum(strat[i,:])>=1e-16:

					thick=numpy.sum(strat[i,:])/(1-porosity)
					thick=round(thick,16)

					if i==len(strat)-1:
						verts=numpy.array([[x[p]-dx/2,0,z[p]-thickinit[p]],
							   [x[p]+dx/2,0,z[p]-thickinit[p]],
							   [x[p]-dx/2,W,z[p]-thickinit[p]],
							   [x[p]+dx/2,W,z[p]-thickinit[p]],
							   [x[p]-dx/2,0,Zp[p]-thickAbove],
							   [x[p]+dx/2,0,Zp[p]-thickAbove],
							   [x[p]-dx/2,W,Zp[p]-thickAbove],
							   [x[p]+dx/2,W,Zp[p]-thickAbove]])

					else:
						verts=numpy.array([[x[p]-dx/2,0,Zp[p]-thickAbove-thick],
							   [x[p]+dx/2,0,Zp[p]-thickAbove-thick],
							   [x[p]-dx/2,W,Zp[p]-thickAbove-thick],
							   [x[p]+dx/2,W,Zp[p]-thickAbove-thick],
							   [x[p]-dx/2,0,Zp[p]-thickAbove],
							   [x[p]+dx/2,0,Zp[p]-thickAbove],
							   [x[p]-dx/2,W,Zp[p]-thickAbove],
							   [x[p]+dx/2,W,Zp[p]-thickAbove]])

					thickAbove=thickAbove+thick
					thick_temp=0

					count=count+1

					vox0=numpy.where((voxels==(verts[0,0],verts[0,1],verts[0,2])).all(axis=1))
					vox0=int(vox0[0])
					vox1=numpy.where((voxels==(verts[1,0],verts[1,1],verts[1,2])).all(axis=1))
					vox1=int(vox1[0])
					vox2=numpy.where((voxels==(verts[2,0],verts[2,1],verts[2,2])).all(axis=1))
					vox2=int(vox2[0])
					vox3=numpy.where((voxels==(verts[3,0],verts[3,1],verts[3,2])).all(axis=1))
					vox3=int(vox3[0])
					vox4=numpy.where((voxels==(verts[4,0],verts[4,1],verts[4,2])).all(axis=1))
					vox4=int(vox4[0])
					vox5=numpy.where((voxels==(verts[5,0],verts[5,1],verts[5,2])).all(axis=1))
					vox5=int(vox5[0])
					vox6=numpy.where((voxels==(verts[6,0],verts[6,1],verts[6,2])).all(axis=1))
					vox6=int(vox6[0])
					vox7=numpy.where((voxels==(verts[7,0],verts[7,1],verts[7,2])).all(axis=1))
					vox7=int(vox7[0])

					file.write('8 '+str(vox0)+' '+str(vox1)+' '+str(vox2)+' '+str(vox3)+' '+str(vox4)+' '+str(vox5)+' '+str(vox6)+' '+str(vox7)+' \n')


				elif numpy.sum(strat[i,:])!=0 and numpy.sum(strat[i,:])<1e-16:
					thick=numpy.sum(strat[i,:])/(1-porosity)
					thick_temp=thick_temp+thick
					thickAbove=thickAbove+thick

				#elif numpy.sum(strat[i,:])!=0 and numpy.sum(strat[i,:])<1e-12:
				#	thick=numpy.sum(strat[i,:])/(1-porosity)
				#	thickAbove=thickAbove+thick

		file.write('\n')
		file.write('CELL_TYPES '+str(Nlayers)+' \n')




		for i in range(0,Nlayers):
			file.write('11 \n')

		file.write('\n')
		file.write('CELL_DATA '+str(Nlayers)+' \n')
		file.write('SCALARS age int \n')
		file.write('LOOKUP_TABLE default \n')


		allage=[]

		for i in range(0,Nlayers):
			allage=allage+[ALLage[i]]
			file.write(str(int(allage[i]))+' \n')

		file.write('\n')
		#file.write('CELL_DATA '+str(Nlayers)+' \n')
		file.write('SCALARS d50 float \n')
		file.write('LOOKUP_TABLE default \n')

		alld50=[]

		for i in range(0,Nlayers):
			alld50=alld50+[ALLd50[i]]
			file.write(str(round(alld50[i],16))+' \n')

		file.write('\n')
		#file.write('CELL_DATA '+str(Nlayers)+' \n')
		file.write('SCALARS Ffine float \n')
		file.write('LOOKUP_TABLE default \n')

		allFfine=[]

		for i in range(0,Nlayers):
			allFfine=allFfine+[ALLFfine[i]]
			file.write(str(round(allFfine[i],16))+' \n')

		file.close()

		# SAVE DATA AS A VTK FILE
		directory='D:/PostDoc_Marine_Le_Minor/MANUSCRIPTS/Manuscript_1DStratigraphy/Codes_Python/'
		file=open(directory+filename+str(int(t*dt))+'.csv','w')

		# Header
		writer=csv.writer(file)
		writer.writerow(['x','z','d50-surface','d50-precipiton','tau','tauc-coarse','tauc-fine','dotd-coarse','dotd-fine','dote/F-coarse','dote/F-fine','F-coarse','xi-coarse','xi-fine'])
		allx=[]
		allz=[]
		for p in range(0,len(ALLx)):
			allx=allx+[ALLx[p]]
			allz=allz+[ALLz[p]]
			writer.writerow([str(round(allx[p],16)),str(round(allz[p],16)),str(round(d50surface[p],16)),str(round(d50p[p],16)),str(round(taup[p],16)),str(round(taucCoarse[p],16)),str(round(taucFine[p],16)),str(round(dotdCoarse[p],16)),str(round(dotdFine[p],16)),str(round(doteCoarse[p],16)),str(round(doteFine[p],16)),str(round(FCoarse[p],16)),str(round(xiCoarse[p],16)),str(round(xiFine[p],16))])

		file.close
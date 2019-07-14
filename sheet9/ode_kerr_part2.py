#!/usr/bin/env python


#read in packages
#for basic math
import math

#numerical routines
from numpy import *

#for integration etc ..
from scipy.integrate import odeint
from scipy.integrate import solve_ivp


#for plotting
from pylab import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def initial(alpha,beta,r0,theta0,phi0,a):

    ###initial conditions for each ray
    ##alpha  =  x-coordinate in the image plane 
    ##beta   =  y-coordinate in the image plane
    ##r0     =  distance to the observer
    ##theta0 =  inclination of the observer
    ##phi0   =  azimuthal angle of the observer
    ##a      =  spin of black hole


    #####
    #calculate initial position
    #####
    ##correct for small angles -->aviod zeros
    if abs(theta0)<1e-8:
        theta0=sign(theta0)*1e-8
    if abs(phi0)<1e-8:
        phi0=sign(phi0)*1e-8

    ##define a**2
    a2 = a * a

    ##define r0**2
    r2 = r0*r0


    ##transformation from observers frame to black hole frame (cartesian coordinates)
    ##Equation 18 in Pu et al. 2016 / Eq. 11 from Ex4
    D=sqrt(r2+a2)*sin(theta0)-beta*cos(theta0)

    xprime     = D*cos(phi0) - alpha*sin(phi0)
    yprime     = D*sin(phi0)+alpha*cos(phi0)
    zprime     = r0*cos(theta0)+beta*sin(theta0)
    
    ##define w 
    w     = xprime*xprime+yprime*yprime+zprime*zprime-a2

    ##intial vector 6dim (r,theta,phi,t,p_r,p_theta)
    #will be returned by the sub-routine
    #y[0]=r
    #y[1]=theta
    #
    y0=np.zeros(6)

    ##convert cartesian to Boyer-Lindquist coordinates
    ##Equations 20 -22 in Pu et al. 2016
    ##radius
    y0[0] = sqrt((w+sqrt(w*w+(2.*a*zprime)*(2.*a*zprime)))/2.)   

    ##theta
    y0[1] = np.arccos(zprime/y0[0])

    ##phi
    y0[2] = np.arctan2(yprime,xprime)
                                   

    ##t --> set time to zero at the image plane (inverse ray-tracing!)
    y0[3]=0

    

    ########
    #calculate initial velocities
    ########

    ##to make programming easier
    r=y0[0]
    theta=y0[1]
    phi=y0[2]
    
    ##define auxiliary variables
    sigma = r*r+(a*cos(theta))*(a*cos(theta))
    R     = sqrt(a2+r*r)
    v     = -sin(theta0)*cos(phi)
    ###reflect velocities (inverse ray-tracing!)
    zdot=-1.

    ###initial velocities in Boyer-Lindquist
    ##Equations 24 -26 in Pu et al. 2016 / Eq 16 -18 in Ex4
    rdot0 = zdot*(-R*R*cos(theta0)*cos(theta)+r*R*v*sin(theta))/sigma         
    thetadot0 = zdot*(cos(theta0)*r*sin(theta)+R*v*cos(theta))/sigma         
    phidot0 = zdot*sin(theta0)*sin(phi)/(R*sin(theta))

    ##define auxiliary variables --> make programming easier
    ##angles
    sintheta=sin(theta)
    costheta=cos(theta)
    cos2 = costheta*costheta
    sin2 = sintheta*sintheta
    
    ##additional variables
    r2 = r * r
    delta = r2 - 2.0 * r + a2
    s1 = sigma - 2.0 * r

    ###get conserved four momentum
    ##Equations 4 and 5 in Pu et al. 2016 / Eq 2 and 3 in Ex4
    y0[4]= rdot0*sigma/delta
    y0[5]= thetadot0*sigma
    
    ##compute the square of the energy
    ##Equation 7 in Pu et al. 2016 / Eq 8 in Ex4
    energy2 = s1*(rdot0*rdot0/delta+thetadot0*thetadot0)+ delta*sin2*phidot0*phidot0
    energy =sqrt(energy2)
   
    
    #rescaled momentum by energy
    y0[4] = y0[4]/energy
    y0[5] = y0[5]/energy
    
    #set E = 1 and compute angular momentum
    ##Equation 8 in Pu et al. 2016 / Eq 9 in Ex4
    L = ((sigma*delta*phidot0-2.0*a*r*energy)*sin2/s1)/energy
    

    ##compute kappa
    ##Equation 15+ in Pu et al. 2016 /Eq 10 in Ex4
    kappa = y0[5]*y0[5]+a2*sin2+L*L/sin2

    ##return initial vector, energy, angular momentum and kappa
    ##as dictonary
    return {'y0': y0, 'E': energy, 'L':L, 'kappa':kappa}




def plot3d(x,y,z,dim,name,M,a):


    ##plot curve/points in 3d
    ##x,y,z    = coordinates to plot
    ##dim      = plot range [x0,x1,y0,y1,z0,z1]
    ##M        = mass of black hole
    ##a        = spin of black hole
    ##name     = file name

    ##use matplotlib and create a figure
    fig = plt.figure()
    ax = fig.gca(projection='3d',aspect='equal',facecolor='w')

    ##plot position
    plot([x],[y],[z],marker='o',color='k',ls='None')

    ##set labels and plot range
    ax.set_xlabel('X')
    ax.set_xlim3d(dim[0], dim[1])
    ax.set_ylabel('Y')
    ax.set_ylim3d(dim[2], dim[3])
    ax.set_zlabel('Z')
    ax.set_zlim3d(dim[4], dim[5])

    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)

    xx=np.outer(np.cos(u), np.sin(v))
    yy=np.outer(np.sin(u), np.sin(v))
    zz=np.outer(np.ones(np.size(u)), np.cos(v))

    #outer horizon radius
    rs=M+sqrt(M**2-a**2)

    ##static limit radius
    rsl=M+sqrt(M**2-a**2*cos(v)**2)

    xrs = rs * xx
    yrs = rs * yy
    zrs = rs * zz

    # Plot the surface
    xrsl = rsl * xx
    yrsl = rsl * yy
    zrsl = rsl * zz

    ##plot surfaces
    ax.plot_surface(xrs, yrs, zrs, color='k',alpha=1)
    ax.plot_surface(xrsl, yrsl, zrsl, color='b',alpha=0.2)
    
    fig.savefig('%s.pdf' %name, bbox_inches='tight', pad_inches = 0.06)
    show()



def oderay(t,y,a,kappa,E,L):

    ###initial conditions for each ray
    ##t      =  t or affine parameter 
    ##y      =  state vector (6-dim,r,theta,phi,t,p_r,p_theta)
    ##a      =  black hole spin
    ##kappa  =  kappa
    ##E      =  Energy conserved for each ray
    ##L      =  angular momentum conserved for each ray


    ##get vector variables -->easier to programm
    r = y[0]
    theta = y[1]
    phi = y[2]
    t = y[3]
    pr = y[4]
    pth = y[5]


    ##auxillary variables
    r2=r*r
    twor=2.0*r
    a2=a*a

    ##sin and cos theta
    sintheta=sin(theta)
    costheta=cos(theta)
    cos2=costheta*costheta
    sin2=sintheta*sintheta

    ##more auxillary variables from metric see below Eq.1 in Ex 4 
    sigma=r2+a2*cos2
    delta=r2-twor+a2
    sd=sigma*delta
    siginv=1.0/sigma
    bot=1.0/sd


    ##avoid small number
    if (abs(sintheta)<1e-8):
        sintheta=sign(sintheta)*1e-8
        sin2=1e-16

        
    ##set equation array for derivatives
    eqs=np.zeros((6))


    ##set RHS of equations eqs 4,5,9,10,14,15 in Pu et al and eqs.2 - 7 to ww in Ex 4
    #dr
    eqs[0]=-pr*delta*siginv

    #dtheta
    eqs[1]=-pth*siginv

    #dphi
    eqs[2]=-(twor*a+(sigma-twor)*L/sin2)*bot

    #dt
    eqs[3]=-(1.0+(twor*(r2+a2)-twor*a*L)*bot)

    #dpr
    eqs[4]=-(((r-1.0)*(-kappa)+twor*(r2 + a2)-2.0*a* L)*bot-2.0*pr*pr*(r - 1.0)*siginv)

    #dpth
    eqs[5]=-sintheta*costheta*(L*L/(sin2*sin2)- a2)*siginv

    #return the LHS of the equations
    return eqs
    


def hitrs(t,y):

    ##stop rays at the horizon
    rs=M+sqrt(M**2-a**2)
    
    return y[0]-(rs+1e-2)

def leavimg(t,y):

    ##stop rays after some distance
    return 12-abs(y[0])



def ray(M,a,r0,alpha,beta,theta0,phi0,tmin,tmax):

    ##combine all routines in one and compute a single ray
    ##M      = mass of black hole
    ##a      = spin of black hole
    ##r0     = distance from black hole 
    ##alpha  = x position on observers screen
    ##beta   = y position on observers screen
    ##theta0 = inclination of observer
    ##phi0   = azimuthal position of observer
    ##tmin   = start time of integration
    ##tmax   = final time of integration

    ##compute initial conditions
    initcond= initial(alpha,beta,r0,theta0,phi0,a)

    ##get conserved E,L and kappa
    kappa=initcond['kappa']
    E=initcond['E']
    L=initcond['L']

    ##set terminate conditions for integrator see scipy solve_ivp
    ##stop rays which hit the horizon
    hitrs.terminal=True

    ##stop rays which travelled to r>90M
    leavimg.terminal=True

    ##call ode solver
    ##for details see https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html#scipy.integrate.solve_ivp
    solver=solve_ivp(lambda t,y: oderay(t,y,a,kappa,E,L),t_span=[tmin,tmax],y0=initcond['y0'],method='Radau',events=(hitrs,leavimg),dense_output=True,atol=1e-8,rtol=1e-6)


    ##compute x,y,z in cartesian coordinates
    ##from r,theta and phi provided by the solver
    x=solver.y[0,:]*cos(solver.y[2,:])*sin(solver.y[1,:])
    y=solver.y[0,:]*sin(solver.y[2,:])*sin(solver.y[1,:])
    z=solver.y[0,:]*cos(solver.y[1,:])


    ##return solver output and coordinates
    return {'solver':solver, 'x':x, 'y':y, 'z':z}

    

def plot3dray(x,y,z,dim,name,M,a):


    ##plot curve/points in 3d
    ##x,y,z    = coordinates of trajectories
    ##dim      = plot range [x0,x1,y0,y1,z0,z1]
    ##M        = mass of black hole
    ##a        = spin of black hole
    ##name     = file name

    ##use matplotlib and create a figure
    fig = plt.figure()
    ax = fig.gca(projection='3d',facecolor='w') # aspect='equal'

    ##plot position
    for i in range(len(x)):
        plot(x[i],y[i],z[i],marker='None',color='r',ls='-')

    ##set labels and plot range
    ax.set_xlabel('X')
    ax.set_xlim3d(dim[0], dim[1])
    ax.set_ylabel('Y')
    ax.set_ylim3d(dim[2], dim[3])
    ax.set_zlabel('Z')
    ax.set_zlim3d(dim[4], dim[5])

    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)

    xx=np.outer(np.cos(u), np.sin(v))
    yy=np.outer(np.sin(u), np.sin(v))
    zz=np.outer(np.ones(np.size(u)), np.cos(v))

    #outer horizon radius
    rs=M+sqrt(M**2-a**2)

    ##static limit radius
    rsl=M+sqrt(M**2-a**2*cos(v)**2)

    xrs = rs * xx
    yrs = rs * yy
    zrs = rs * zz

    # Plot the surface
    xrsl = rsl * xx
    yrsl = rsl * yy
    zrsl = rsl * zz

    ##plot surfaces
    ax.plot_surface(xrs, yrs, zrs, color='k',alpha=1)
    ax.plot_surface(xrsl, yrsl, zrsl, color='b',alpha=0.2)
    
    fig.savefig('%s.pdf' %name, bbox_inches='tight', pad_inches = 0.06)
    show()



#####################  main  ##############################

##set some initial properties of the black hole
##mass of black hole
M=1

##spin of black hole
a=0.9


##set observers orientation
##inclination (in radian)
theta0=90*np.pi/180

##azimuthal position (in radian)
phi0=0*np.pi/180

##distance from black hole (in M)
r0=10

##inital position of ray in the image plane (observers frame)
alpha=6.475479
beta=0

##start and stop time of array integration
tstart=0
tstop=90

##call the integrator for current settings
ray1=ray(M,a,r0,alpha,beta,theta0,phi0,tstart,tstop)

##set dimensions of plot
dim=[-10,10,-10,10,-10,10]

##set name of output file
name='test3d'

##plot black hole and ray
niceplotray=plot3dray([ray1['x']],[ray1['y']],[ray1['z']],dim,name,M,a)

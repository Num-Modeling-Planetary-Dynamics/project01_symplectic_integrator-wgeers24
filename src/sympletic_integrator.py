
import numpy as np
import os
import matplotlib.pyplot as plt



#Initial Conditions

#Global constants in SI units
G= 6.6743E-11
m_Sun=1.9885E+30
mu= -G*(m_Sun)
m_Neptune= 102.409E+24
m_Pluto= 1.307E+22
Gmass=(m_Neptune + m_Pluto)
GMcb= -G*(m_Sun)
JD2S=86400
YR2S=np.longdouble(365.25*JD2S)
dt=5.5*JD2S
#Neptune Initial Conditions

x_Neptune=2.266223889656152E+08
y_Neptune=4.461841849742382E+09
z_Neptune=-9.710614902790475E+07
vx_Neptune=-5.447617894349518E+00
vy_Neptune= 3.044797889906726E-01
vz_Neptune= 1.187834180290304E-01

#Pluto Initial Conditions

x_Pluto= 1.540431633806554E+09
y_Pluto= 6.753704837972438E+09
z_Pluto=-1.168941164341141E+09
vx_Pluto=-3.741653796726164E+00
vy_Pluto= 3.621999156857085E-01
vz_Pluto= 1.046185526637238E+00



def xv2el (mu,x,y,z,vx,vy,vz):
#input position and velocity vectors in cartesian coordinates
    r_vec=np.array([x,y,z])
    r_mag=np.linalg.norm(r_vec)
    v_vec=np.array([vx,vy,vz])
    v_mag=np.linalg.norm(v_vec)
    h_vec=np.cross(r_vec,v_vec)
    h_mag=np.linalg.norm(h_vec)
    r_r_dot=np.vdot(r_vec,v_vec)
    r_dot=np.sign(r_r_dot)*np.sqrt((v_mag**2)-(h_mag**2/r_mag**2))


    a=((2/r_mag)-(v_mag**2/mu))**-1
    ecc=np.sqrt(1-(h_mag**2/(mu*a)))
    I=np.arccos(h_vec[2]/h_mag)

    sO=(np.sign(h_vec[2])*h_vec[0]/(h_mag*np.sin(I)))
    cO=(np.sign(h_vec[2])*h_vec[1])/(h_mag*np.sin(I))
    omega=np.arctan2(sO,cO)

    sof=z/(r_mag*np.sin(I))
    cof=((x/r_mag)+np.sin(omega)*sof*np.cos(I))/np.cos(omega)
    pf=np.arctan2(sof,cof)

    sf=((a*(1-ecc**2))/(h_mag*ecc))*(r_dot)
    cf=(1/ecc)*((a*(1-ecc**2)/r_mag)-1)
    f=np.arctan2(sf,cf)
    peri=pf-f



#output orbital elements
#Use equations from sxn 2.8 pg 52 and 53
    return a,ecc,I,omega,peri,f

def vh2vb (vhvec, mu, Gmass):
    #Converts heliocentric velocity vector to barycentric
    Gmtot= mu + np.sum(Gmass)
    vbcb = -np.sum(Gmass* vhvec,axis=0)/Gmtot
    vbvec= vhvec + vbcb
    #vhvec: heliocentric velocity vectors of all n bodies
    #Gmass: masses of n bodies in system
    #vb: barycentric velocity vectors of n bodies
    #mu: central body gravitational parameter
    #vbcv: barycentric velocity vector of central body
    return vbvec, vbcb

def vb2vh (vbvec, mu, Gmass):
    vbcb = -np.sum(Gmass * vbvec, axis=0)/ mu
    vhvec= vbvec - vbcb
    # vhvec: heliocentric velocity vectors of all n bodies
    # Gmass: masses of n bodies in system
    # vb: barycentric velocity vectors of n bodies
    # mu: central body gravitational parameter
    # vbcv: barycentric velocity vector of central body
    return vhvec



def Kep_drift (M,r_vec0,v_vec0,mu,dt):
    #Uses Kepler's equation for eccentric anomaly E
    # Equation 2.49 solve for E
    def Danby(M,ecc,accuracy=1E-14):

        def f(E):
            return E - ecc * np.sin(E) - M

        def fP(E):
            return 1 - ecc * np.cos(E)

        def fP2(E):
            return ecc * np.sin(E)

        def fP3(E):
            return ecc * np.cos(E)

        def D1(E):
            return -(f(E)) / (fP(E))

        def D2(E):
            return -(f(E)) / (fP(E) + 0.5 * (D1(E) * fP2(E)))

        def D3(E):
            return -(f(E)) / (fP(E) + 0.5 * (D2(E) * fP2(E)) + (1 / 6) * ((D2(E) ** 2) * fP3(E)))

        while True:
            Enew = E[-1] + D3(E[-1])
            E.append(Enew)
            if np.abs((E[-1] - E[-2])/E[-2])<accuracy:
                break
    #utilizes f and g functions
    #Danby method
    r_mag0=np.linalg.norm(r_vec0)
    v_mag2=np.vdot(v_vec0,v_vec0)
    h_vec0= np.cross(r_vec0, v_vec0)
    h_mag2= np.vdot(h_vec0,h_vec0)
    a = 1.0/(2.0/ r_mag0-v_mag2/mu)
    ecc= np.sqrt(1-h_mag2 / (mu*a))

    n= np.swrt(mu / a**3)
    E0 = np.where(ecc > np.finfo(np.float64).tiny, np.arccos(-(r_mag0 - a) / (a * ecc)), 0)

    if ecc < np.finfo(np.float64).tiny: #Uses M as E in 0 ecc orbits
        E0=0.0
    else:
        E0=np.arccos(-(r_mag0-a)/(a*ecc))

    if np.sign(np.vdot(r_vec0,v_vec0))<0.0:
        E0= 2*np.pi-E0

    M0=E0-ecc*np.sin(E0)

    M-M0+n*dt
    E=Danby(M,ecc)
    dE= E - E0

    f=a/r_mag0*(np.cos(dE)-1.0)+1.0
    g=dt+1.0/n*(np.sin(dE)-dE)

    r_vec=f * r_vec0 + g*v_vec0
    r_mag= np.linalg.norm(r_vec)
    fdot = -a**2 / (r_mag*r_mag0) *n * np.sin(dE)
    gdot= a/r_mag* (np.cos(dE)-1.0)+1.0

    v_vec= fdot* r_vec0 +gdot* v_vec0
    return r_vec, v_vec
    #Define fdot and gdot for vmag using eq 2.71 on pg 37
def kick ():
    #E_int step
    #uses the Hamiltonian for the interaction step
    #half-time step (solar drift/linear drift)
    def getacc_one(Gm, rvec, i):
        #Gm: G*mass values of each body
        #rvec: cartesian position vectors
        #i: indec of body to return the accelerations

    #acc: accelerations acting on body i due to all other bodies
        drvec= rvec - rvec[i, :]
        irij3= np.linalg.norm(drvec, axis=1)**3
        irij3[i]= 1
        irij3= Gm / irij3

        return np.sum(drvec.T * irij3, axis=1)

    n= rhvec.shape[0]
    acc = np.array([getacc_one(Gmass.flatten(), rhvec, i)for i in range(n)])
    vbvec += acc * dt
    return
def Sun_drift(Gmass,vbvec,dt):
    #updates barycenter of Sun
    pt = np.sum(Gmass *vbvec, axis=0) / param['GMcb']
    rhvec += pt*dt

def step(dt):
    #Advances simulation by dt
    #Use dt as 5
    dth=0.5*dt
    Sun_drift(Gmass,dth)
    kick(dth)
    Kep_drift(rvec0,vvec0,mu,dt)
    kick(dth)
    Sun_drift(dth)
def drift_one(mu,x,y,z,vx,vy,vz,dt):
    #inputs cartesian coordinates

    #Outputs new position and velocity cartesian coordinates

    def calc_energy_one(vhvec, rhvec, Gmass, mu):
        # Computes system energy at particular state of simulation
        # rhvec: heliocentric position vectors of n bodies
        # vhvec: heliocentric velocity vectors of n bodies
        # Gmass: masses of n bodies in system
        # mu: central body gravitational parameter
        # vbvec,vbcb: from vh2vb
        vbvec,vbcb= vh2vb(vhvec, mu, Gmass)
        vbmag2 = np.einsum("ij,ij->i", vbvec, vbvec)
        irh = 1.0 / np.linalg.norm(rhvec, axis=1)
        ke = 0.5 * (mu * np.vdot(vbcb, vbcb) + np.sum(Gmass * vbmag2)) / G

        def pe_one(Gm, r_vec, i):
            drvec = r_vec[i + 1:, :] - r_vec[i, :]
            irij = np.linalg.norm(drvec, axis=1)
            irij = Gm[i + 1:] * Gm[i] / irij
            return np.sum(drvec.T * irij, axis=1)

        n = rhvec.shape[0]
        pe = (-mu * np.sum(Gmass * irh) - np.sum([pe_one(Gmass.flatten(), rhvec, i) for i in range(n - 1)])) / G

        return ke + pe


#Make plots
fig, ax = plt.subplots(figsize=(8,6))
    data['phi'].plot(x="time",ax=ax)
    plt.savefig(os.path.join(os.pardir,"plots","resonance_angle.png"),dpi=300)
    plt.close()

    fig, ax = plt.subplots(figsize=(8,6))
    data['dE/E0'].plot(x="time",ax=ax)
    plt.savefig(os.path.join(os.pardir,"plots","energy.png"),dpi=300)
    plt.close()

def main():

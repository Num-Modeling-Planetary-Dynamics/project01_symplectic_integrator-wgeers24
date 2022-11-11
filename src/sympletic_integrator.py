
import numpy as np
import os
import matplotlib.pyplot as plt



#Initial Conditions 01/01/1900-01/01/2099

#Global constants in SI units
G= 6.6743E-11
m_Sun=1.9885E+30
m_Neptune= 102.409E+24
m_Pluto= 1.307E+22
m_planet=np.array(m_Neptune, m_Pluto)
mu=G*(m_Sun+m_planet)
Gmass=m_Neptune + m_Pluto
Gm= Gmass*G
GMcb= G*m_Sun
JD2S=86400
YR2S=np.longdouble(365.25*JD2S)
dt=5.5*YR2S
#Neptune Initial Conditions SI

x_Neptune=2.266223889656152E+08
y_Neptune=4.461841849742382E+09
z_Neptune=-9.710614902790475E+07
vx_Neptune=-5.447617894349518E+00
vy_Neptune= 3.044797889906726E-01
vz_Neptune= 1.187834180290304E-01

#Pluto Initial Conditions SI

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

    a = ((2 / r_mag) - (v_mag ** 2 / mu)) ** -1
    ecc = np.sqrt(1 - (h_mag ** 2 / (mu * a)))
    inc = np.arccos(h_vec[2] / h_mag)

    sO = (np.sign(h_vec[2]) * h_vec[0] / (h_mag * np.sin(inc)))
    cO = (np.sign(h_vec[2]) * h_vec[1]) / (h_mag * np.sin(inc))
    Omega = np.arctan2(sO, cO)

    sof = z / (r_mag * np.sin(inc))
    cof = ((x / r_mag) + np.sin(Omega) * sof * np.cos(inc)) / np.cos(Omega)
    of = np.arctan2(sof, cof)

    sf = ((a * (1 - ecc ** 2)) / (h_mag * ecc)) * (r_dot)
    cf = (1 / ecc) * ((a * (1 - ecc ** 2) / r_mag) - 1)

    f = np.arctan2(sf, cf)
    omega = of - f

    omega=np.mod(omega, 2*np.pi)

    varpi= np.mod(Omega+omega, 2*np.pi)

    E= np.where(ecc>np.finfo(np.float64).tiny, np.arccos(-(r_mag-a)/(a*ecc)),0)
    if np.sign(np.vdot(r_vec, v_vec))>0.0:
        E=2*np.pi-E

    M=E-ecc*np.sin(E)
    lam= np.mod(M+varpi, 2*np.pi)

    inc=np.rad2deg(inc)
    Omega=np.rad2deg(Omega)
    omega=np.rad2deg(omega)
    varpi=np.rad2deg(varpi)
    f= np.rad2deg(f)
    M= np.rad2deg(M)
    lam= np.rad2deg(lam)

    return a, ecc, inc, Omega, omega, varpi, f, M, lam

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
    #Treats n-body problem as 2-body for each body in system
    def Danby(M,ecc,accuracy=1E-14):
        E=[M]
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

        Maxloops=50
        for i in range(Maxloops):
            Enew = E[-1] + D3(E[-1])
            E.append(Enew)
            if np.abs((E[-1] - E[-2])/E[-2])<accuracy:
                return Enew

    r_mag0 = np.linalg.norm(r_vec0)
    v_mag2 = np.vdot(v_vec0, v_vec0)
    h_vec0 = np.cross(r_vec0, v_vec0)
    h_mag2 = np.vdot(h_vec0, h_vec0)
    a = 1.0 / (2.0 / r_mag0 - v_mag2 / mu)
    ecc = np.sqrt(1 - h_mag2 / (mu * a))

    n = np.sqrt(mu / a ** 3)
    E0 = np.where(ecc > np.finfo(np.float64).tiny, np.arccos(-(r_mag0 - a) / (a * ecc)), 0)

    if ecc < np.finfo(np.float64).tiny:  # Uses M as E in 0 ecc orbits
        E0 = 0.0
    else:
        E0 = np.arccos(-(r_mag0 - a) / (a * ecc))

    if np.sign(np.vdot(r_vec0, v_vec0)) < 0.0:
        E0 = 2 * np.pi - E0

    M0 = E0 - ecc * np.sin(E0)

    M - M0 + n * dt
    E = Danby(M, ecc)
    dE = E - E0

    f = a / r_mag0 * (np.cos(dE) - 1.0) + 1.0
    g = dt + 1.0 / n * (np.sin(dE) - dE)

    r_vec = f * r_vec0 + g * v_vec0
    r_mag = np.linalg.norm(r_vec)
    fdot = -a ** 2 / (r_mag * r_mag0) * n * np.sin(dE)
    gdot = a / r_mag * (np.cos(dE) - 1.0) + 1.0

    v_vec = fdot * r_vec0 + gdot * v_vec0
    return r_vec, v_vec
    #utilizes f and g functions
    #Danby method
    #Updates semi-major axis and eccentricity

    #Define fdot and gdot for vmag using eq 2.71 on pg 37
def kick (rhvec, vbvec, dt):
    #E_int step
    #uses the Hamiltonian for the interaction step
    #half-time step (solar drift/linear drift)
    #force of gravity between 2 bodies that aren't sun
    def getacc_one(Gm, r_vec, i):
        #Gm: G*mass values of each body
        #rvec: cartesian position vectors
        #i: indec of body to return the accelerations

    #acc: accelerations acting on body i due to all other bodies
        dr_vec= r_vec - r_vec[i, :]
        irij3= np.linalg.norm(dr_vec, axis=1)**3
        irij3[i]= 1
        irij3= Gm / irij3

        return np.sum(dr_vec.T * irij3, axis=1)

    n= rhvec.shape[0]
    acc = np.array([getacc_one(Gmass, rhvec, i)for i in range(n)])
    vbvec += acc * dt
    return vbvec
def Sun_drift(Gmass,rhvec,vbvec,GMcb,dt):
    #updates barycenter of Sun
    pt = np.sum(Gmass *vbvec, axis=0) / GMcb
    rhvec += pt*dt

def step(rvec0,vvec0,mu,Gmass,dt):
    #Advances simulation by dt
    #Use dt as 5
    dth=0.5*dt
    rvec, vvec= Sun_drift(Gmass, rhvec, vbvec, GMcb,dth)
    vvec = kick(rhvec, vbvec, dth)
    rvec, vvec= Kep_drift(M,rvec0,vvec0,mu,dt)
    vvec= kick(rhvec, vbvec, dth)
    rvec, vvec= Sun_drift(Gmass, rhvec, vbvec, GMcb, dth)
    return rvec, vvec
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

    def pe_one(Gm, r_vec, i): #Heart of integrator
        drvec = r_vec[i + 1:, :] - r_vec[i, :]
        irij = np.linalg.norm(drvec, axis=1)
        irij = Gm[i + 1:] * Gm[i] / irij
        return np.sum(drvec.T * irij, axis=1)

    n = rhvec.shape[0]
    pe = (-mu * np.sum(Gmass * irh) - np.sum([pe_one(Gmass.flatten(), rhvec, i) for i in range(n - 1)])) / G

    return ke + pe


#Make plots from help, make own using matplotlib
#fig, ax = plt.subplots(figsize=(8,6))
#data['phi'].plot(x="time",ax=ax)
#plt.savefig(os.path.join(os.pardir,"plots","resonance_angle.png"),dpi=300)
#plt.close()
def simulation(rhvec, vhvec, mu, Gmass, dt, tfinal):
    noutput=tfinal/dt

    for i in range(noutput):
        step()


#fig, ax = plt.subplots(figsize=(8,6))
#data['dE/E0'].plot(x="time",ax=ax)
#plt.savefig(os.path.join(os.pardir,"plots","energy.png"),dpi=300)
#plt.close()

#Plot of total change in system energy


if __name__=="__main__":
    dt = 5.0 * YR2S
    tfinal = 1e5 * YR2S
    rhvec = np.empty((2,3))
    vhvec = np.empty((2,3))
    rhvec[0,:] = np.array([x_Neptune, y_Neptune, z_Neptune])
    rhvec[1,:] = np.array([x_Pluto, y_Pluto, z_Pluto])
    vhvec[0,:] = np.array([vx_Neptune, vy_Neptune, vz_Neptune])
    vhvec[1,:] = np.array([vx_Pluto, vy_Pluto, vz_Pluto])
    Gmass = np.array([G*m_Neptune, G*m_Pluto])
    mu= Gmass + G*m_Sun
    rhhist= [np.copy(rhvec)]
    vhhist= [np.copy(vhvec)]

    elem_Neptune = [np.array(xv2el(rhvec[0,0], rhvec[0,1], rhvec[0,2], vhvec[0,0], vhvec[0,1], vhvec[0,2]))]
    elem_Pluto = [np.array(xv2el(rhvec[1,0], rhvec[1,1], rhvec[1,2], vhvec[1,0], vhvec[1,1], vhvec[1,2]))]
    time = np.arange(start=0.0, stop =tfinal, step=dt)
    vbvec = vh2vb(vhvec, mu, Gmass)
    elem_Neptune = np.empty((len(time)+1, 9))
    for t in time:

        rhvec, vhvec=step(rhvec, vhvec, mu, Gmass, dt)
        vhvec= vb2vh(vbvec, mu, Gmass)
        elem_Neptune.append(xv2el(rhvec[0,0], rhvec[0,1], rhvec[0,2], vhvec[0,0], vhvec[0,1], vhvec[0,2]))
        elem_Pluto.append(xv2el(rhvec[1,0], rhvec[1,1], rhvec[1,2], vhvec[1,0], vhvec[1,1], vhvec[1,2]))
        rhhist.append([np.copy(rhvec)])
        vhhist.append([np.copy(vhvec)])

    #rhvec_new, vhvec_new = simulation(rhvec, vhvec, mu, Gmass, dt, tfinal)

    plt.subplots(figsize=(8, 6))
    dE= E[-1]-E
    M_Neptune= elem_Neptune[7]
    M_Pluto = elem_Pluto[7]

    y_values = dE / E0
    x_values = [time]
    plt.plot(x_values, y_values)
    plt.ylabel('Change in Energy')
    plt.xlabel('Time (Year)')
    plt.title('Change in Energy of System')
    plt.savefig('Delta_Energy.pdf')
    # Plot of resonance angle vs time
    plt.subplots(figsize=(8, 6))
    y_values = res_angle
    x_values = [time]
    plt.plot(x_values, y_values)
    plt.ylabel('Resonance angle')
    plt.xlabel('Time')
    plt.title('Change in Resonance angle over time')
    plt.savefig('Delta_Res.pdf')





import numpy as np
import math
from decimal import Decimal
import time
import os
import sys

# ------Constantes

Ms = 1.98911 * 10**33  # g # Mass unit
Rs = 0.00465047  # Radio solar en UA
Mo = 1.15 * 10**26 / Ms  # g
Me = 5.972 * 10**27  # g
# C=(Mo/(3))**(1/3.)
# hillrad=C*amax*F # Hill´s radius
Mpla = Mo  # 0.055*Me/Ms # Masa de mercurio, se consideran planetas mayores a Mpla

Standt = 0.01  # Standard dt
LU = 5.0939 * 10**13  # cm length unit
UA = 14959787070000.0  # en cm
F = UA / LU  # Convertion Factor

##############CONVERTION FACTORS#############
g_to_Ms = 1. / Ms
kg_to_Ms = 1000. / Ms
yr_to_sec = 31556926.
cm_to_LU = 1. / LU
m_to_LU = 100. / LU

############COLLISION PARAMETERS#############
S = (3e7) * kg_to_Ms * yr_to_sec**2 / (m_to_LU)  # Impact Strength
vc = (363) * yr_to_sec * m_to_LU  # Sound Velocity at Basalt
den = (3e3) * kg_to_Ms / (m_to_LU)**3  # Density of the body, assume constant
ci = 0.7  # Restitution Coefficient
cii = 0.5  # Modified Restitution Coefficient
Vc = 55 * m_to_LU * yr_to_sec  # Maximum rebound velocity
Vcc = 2. * S / (vc * den)  # Calculated Maximum rebound velocity
cej = 3e6 * (m_to_LU * yr_to_sec)**(9. / 4.)  # Ejecta velocity coefficient
Kms = 1. / (1e7) * (1. / (m_to_LU * yr_to_sec)
                    )**2  # Mass excavated coefficient
phi = math.pi / 4.  # 0. #Angule orientation of debris velocities

###########JUPITER & SATURN PARAMETERS###########
Mjup = 1.898 * 10**30 / Ms
Msat = 5.683 * 10**29 / Ms


############################################################################
#################################FUNCTIONS##################################
############################################################################

# Delete elements of an array
def RMArray(Ind, MDIV, CDIV, MI, VI, RI, Bol1):  # MDIV:MDiv,CDIV:CDiv,MI:Mass,VI:Vel,RI:Pos
    if(len(Ind) != 0):
        for hy in range(0, len(Ind)):
            if(Bol1 == True):
                MDIV += MI[Ind[hy] - hy]
                CDIV += 1
            MI = np.delete(MI, Ind[hy] - hy, 0)  # Delete diverge masses
            VI = np.delete(VI, Ind[hy] - hy, 0)
            RI = np.delete(RI, Ind[hy] - hy, 0)

    return(MDIV, CDIV, MI, VI, RI)


# ------Main integrador

def mainIntegrator(systemName: str, folder_save: str):

    ############INITIAL CONDITIONS###########
    record = 100
    CDiv = 0  # Counter of diverge masses
    MDiv = 0.  # Total mass of diverge masses
    Coun100 = 0  # Count the times that pass 100 yr
    Coun10 = 0  # Count the times that pass 10 yr
    POrt = 0  # Swicth to calculate orbital data
    CPOrt = 0  # Counter to get the five points of the position
    Reg = (Mo / (4. * math.pi / 3. * den))**(1. / 3.)  # Fixed Regularitation
    RegDt = 0.0000000001  # Regularizacion para calcular el Dt, para evitar singularidades
    Col = 0  # Collisions and divergences
    Div = 0
    wdir = '.'  # Working directory, if working locally set to '.'
    wdir += '/'  # Adds slash to user-defined working directory
    odir = sys.argv[1]
    odir += '/'

    ############# READ INITIAL SYSTEM #################
    archivo = np.genfromtxt(systemName, delimiter=',')
    # print(systemName)
    M = np.array(archivo[:, 0])
    Rix = np.array(archivo[:, 1])
    Riy = np.array(archivo[:, 2])
    Vix = np.array(archivo[:, 3])
    Viy = np.array(archivo[:, 4])
    Mstar = archivo[0, 5]  # Mass of the Star
    # tf = archivo[0, 6]  # End time
    tf = 1e6  # this is hardcoded, comment this line if using tf from prepared InitialSystem file
    Mmin = (8 / 100) * Mo  # Merge - frin

    N = len(M)  # 's body
    ri = []
    vi = []
    par = []
    Mi = []
    t = 0.0  # Initial time

    for i in range(0, N):
        ri += [[Rix[i], Riy[i]]]
        vi += [[Vix[i], Viy[i]]]  # math.sqrt((Mo+1)/Ri[i][0])
        Mi += [M[i]]

    Mtot = np.sum(Mi)
    ri = np.array(ri)
    vi = np.array(vi)
    Mi = np.array(Mi)

    ############# END READ INITIAL SYSTEM #################
    folder_save_states = folder_save + '/'
    try:
        os.mkdir(wdir + folder_save_states)
    except:
        pass
        # print('')

    try:
    	os.mkdir(wdir + folder_save_states + systemName[11:-4])
    except:
    	pass
    
    Nmfile = 'Nmass.txt'
    outfile2 = open(wdir + folder_save_states + systemName[11:-4] + '/' + Nmfile, 'w')
    outfile2.write("%.5f , %d , %d , %d, %d, %.20f , %.20f , %.20f, %.5f %s " % (0.6, len(Mi), 0, len(Mi), len(np.where(Mi > Mpla)[0]), Mtot, Mtot, Mtot, Mstar, '\n'))
    # outfile2.write(systemName + '\n')
    outfile2.close()

    start = time.time()
    # print(systemName, 'Start', time.strftime("%c"))
    myk = 0

    while(t <= tf):
        Nend = len(Mi)
        NewThings = []
        RadMass = (Mi / (4. * math.pi / 3. * den)
                   )**(1. / 3.)  # Radio de las masas
        # e=2*np.amin(RadMass) #Regularitation 2 time de min radius
        # e=2*np.mean(RadMass) #Regularitation 2 time de mean radius
        # print(e)
        ###################COLLISION FRAGMENTATION CODE####################
        for i in range(0, Nend):
            for j in range(0, i + 1):
                if (i != j and Mi[i] != 0 and Mi[j] != 0):
                    rij = np.array([ri[i][0] - ri[j][0], ri[i][1] - ri[j][1]])
                    if np.linalg.norm(rij) <= (RadMass[i] + RadMass[j]):  # 160*
                        # print(t, "Colision", len(Mi))
                        ###################COLLISION CODE####################
                        # Impact energy due to the motion of bodies relative
                        E = 0.5 * (Mi[i] * Mi[j]) / (Mi[i] + Mi[j]) * \
                            ((vi[i][0] - vi[j][0])**2 +
                             (vi[i][1] - vi[j][1])**2)
                        E1 = S * Mi[i] / den  # Impact energy m1
                        E2 = S * Mi[j] / den  # Impact energy m2
                        Vrx = (vi[i][0] - vi[j][0])
                        Vry = (vi[i][1] - vi[j][1])
                        Vr = math.sqrt(Vrx**2 + Vry**2)
                        # Velocity of center of mass
                        Tmass = Mi[i] + Mi[j]  # Total mass
                        Vgx = (Mi[i] * vi[i][0] + Mi[j] * vi[j][0]) / \
                            Tmass  # x-velocity of center of mass
                        Vgy = (Mi[i] * vi[i][1] + Mi[j] * vi[j][1]) / \
                            Tmass  # y-velocity of center of mass
                        rcmx = (Mi[i] * ri[i][0] + Mi[j] * ri[j][0]) / \
                            Tmass  # CM Position
                        rcmy = (Mi[i] * ri[i][1] + Mi[j] * ri[j][1]) / \
                            Tmass  # CM Position
                        #####################REBOUND######################
                        # Very low relative velocity-No change of mass
                        # Relative velocity
                        # Para masa minima puede haber Rebound and Acrecion
                        if (Mi[i] < Mmin or Mi[j] < Mmin):
                            if Vr < Vc:
                                #print("Rebound min mass")
                                # Rebound Velocity
                                Vrebx = -ci * Vrx
                                Vreby = -ci * Vry
                                # New Velocities
                                vi[i][0] = Vgx - Mi[j] / Tmass * Vrebx
                                vi[i][1] = Vgy - Mi[j] / Tmass * Vreby
                                vi[j][0] = Vgx + Mi[i] / Tmass * Vrebx
                                vi[j][1] = Vgy + Mi[i] / Tmass * Vreby
                                rsc = (Mi[i] / (4. * math.pi / 3. * den)
                                       )**(1. / 3.)
                                # Velocidad de escape de ambos objetos
                                Vsc = math.sqrt(2 * Mi[i] / rsc)
                                rscpj = (
                                    Mi[j] / (4. * math.pi / 3. * den))**(1. / 3.)
                                Vscpj = math.sqrt(2 * Mi[j] / rscpj)
                                Vresx = (vi[i][0] - vi[j][0])
                                Vresy = (vi[i][1] - vi[j][1])
                                Vrscp = math.sqrt(Vresx**2 + Vresy**2)
                                if (Vrscp < Vscpj or Vrscp < Vsc):
                                    Mi[i] = 0
                                    Mi[j] = Tmass
                                    vi[j][0] = Vgx
                                    vi[j][1] = Vgy
                                    ri[j][0] = rcmx
                                    ri[j][1] = rcmy
                                # What happened when the velocity is less than the
                                # escape velocity
                            else:
                                #print("Accreted 1")
                                Mi[i] = 0
                                Mi[j] = Tmass
                                vi[j][0] = Vgx
                                vi[j][1] = Vgy
                                ri[j][0] = rcmx
                                ri[j][1] = rcmy
                        else:
                            if Vr < Vc:
                                # print("Rebound")
                                # Rebound Velocity
                                Vrebx = -ci * Vrx
                                Vreby = -ci * Vry
                                # New Velocities
                                vi[i][0] = Vgx - Mi[j] / Tmass * Vrebx
                                vi[i][1] = Vgy - Mi[j] / Tmass * Vreby
                                vi[j][0] = Vgx + Mi[i] / Tmass * Vrebx
                                vi[j][1] = Vgy + Mi[i] / Tmass * Vreby
                                rsc = (Mi[i] / (4. * math.pi / 3. * den)
                                       )**(1. / 3.)
                                Vsc = math.sqrt(2 * Mi[i] / rsc)
                                rscpj = (
                                    Mi[j] / (4. * math.pi / 3. * den))**(1. / 3.)
                                Vscpj = math.sqrt(2 * Mi[j] / rscpj)
                                Vresx = (vi[i][0] - vi[j][0])
                                Vresy = (vi[i][1] - vi[j][1])
                                Vrscp = math.sqrt(Vresx**2 + Vresy**2)
                                if (Vrscp < Vscpj or Vrscp < Vsc):
                                    Mi[i] = 0
                                    Mi[j] = Tmass
                                    vi[j][0] = Vgx
                                    vi[j][1] = Vgy
                                    ri[j][0] = rcmx
                                    ri[j][1] = rcmy
                                # What happened when the velocity is less than the
                                # escape velocity

                            ###########REBOUND WITH CRATER FORMATION###########
                            if (Vr > Vc and E < E1 and E < E2):
                                #print("Rebound with crater Formation")
                                # Rebound Velocity
                                Vrebx = -cii * Vrx
                                Vreby = -cii * Vry
                                # New Velocities
                                vi[i][0] = Vgx - Mi[j] / Tmass * Vrebx
                                vi[i][1] = Vgy - Mi[j] / Tmass * Vreby
                                vi[j][0] = Vgx + Mi[i] / Tmass * Vrebx
                                vi[j][1] = Vgy + Mi[i] / Tmass * Vreby
                                # escape velocity from m1
                                Ve1 = (6 * Mi[i]**2 / (math.pi * den))
                                # Amount of mass that escapes from m1
                                mej1 = cej * Kms * Mi[j] * \
                                    Vr**2 * (Ve1)**(-3. / 8.)
                                # escape velocity from m2
                                Ve2 = (6 * Mi[j]**2 / (math.pi * den))
                                # Amount of mass that escapes from m2
                                mej2 = cej * Kms * Mi[i] * \
                                    Vr**2 * (Ve2)**(-3. / 8.)
                                # What happened when the velocity is less than the
                                # escape velocity
                                Mi[i] = Mi[i] - mej1 + mej2
                                Mi[j] = Mi[j] + mej1 - mej2
                                rsc = (Mi[i] / (4. * math.pi / 3. * den)
                                       )**(1. / 3.)
                                Vsc = math.sqrt(2 * Mi[i] / rsc)
                                rscpj = (
                                    Mi[j] / (4. * math.pi / 3. * den))**(1. / 3.)
                                Vscpj = math.sqrt(2 * Mi[j] / rscpj)
                                Vresx = (vi[i][0] - vi[j][0])
                                Vresy = (vi[i][1] - vi[j][1])
                                Vrscp = math.sqrt(Vresx**2 + Vresy**2)
                                if (Vrscp < Vscpj or Vrscp < Vsc):
                                    Mi[i] = 0
                                    Mi[j] = Tmass
                                    vi[j][0] = Vgx
                                    vi[j][1] = Vgy
                                    ri[j][0] = rcmx
                                    ri[j][1] = rcmy

                            ###################FRAGMENTATION###################
                            if (E > E1 and E > E2):
                                #print("Fragmentation 5")
                                radius = (
                                    Tmass / (4 * math.pi / 3 * den))**(1. / 3.)
                                # (6.*Tmass**2/(math.pi*den))
                                Ve = math.sqrt(2 * Tmass / radius)
                                # Amount of mass that leaves the two-body
                                # system
                                mej = cej * Tmass * (Ve)**(-9. / 4.)
                                # Mass of the biggest fragment, !!!!Marginal
                                # gramentation??
                                mmax = 0.5 * Tmass * \
                                    (E / (0.5 * (E1 + E2)))**(-1.24)
                                alfa = 4. * mmax / Tmass
                                MFrag = [0, 0, 0, 0]
                                if alfa < 1:
                                    MFrag[0] = mej / 4.
                                    MFrag[1] = mej / 4.
                                    MFrag[2] = mej / 4.
                                    MFrag[3] = mej / 4.
                                else:
                                    MFrag[0] = mej * mmax / Tmass
                                    MFrag[1] = mmax / alfa
                                    MFrag[2] = (3 - alfa) / (2 * alfa) * mmax
                                    MFrag[3] = (3 - alfa) / (2 * alfa) * mmax
                                # initial condition of remaining core
                                coremass = Tmass - mej
                                # print(coremass)
                                vcmx = (Mi[i] * vi[i][0] + Mi[j] *
                                        vi[j][0]) / Tmass  # CM Velocity
                                vcmy = (Mi[i] * vi[i][1] + Mi[j] *
                                        vi[j][1]) / Tmass  # CM Velocity
                                R = (coremass / (4. * math.pi / 3. * den))**(1. / 3.)
                                # Escape velocity coremass
                                Vecorm = (2 * coremass / R)**(1. / 2.)
                                # if(R==0 or math.isnan(R)):
                                # print(R,coremass,mej,Tmass,Mi[i],M[j])
                                # Initial condition of new ragment
                                RFrag = []  # Fragment position
                                VFrag = []  # Fragment Velocity
                                Mi[i] = 0
                                Mi[j] = coremass
                                vi[j][0] = vcmx
                                vi[j][1] = vcmy
                                ri[j][0] = rcmx
                                ri[j][1] = rcmy
                                for s in range(0, 4):
                                    RFrag += [[rcmx + 4 * R * math.cos(
                                        phi + s * 0.5 * math.pi), rcmy + 4 * R * math.sin(phi + s * 0.5 * math.pi)]]
                                    VFrag += [[Vecorm * math.cos(phi + s * 0.5 * math.pi), Vecorm * math.sin(
                                        phi + s * 0.5 * math.pi)]]
                                NewThings += [[[MFrag[0], MFrag[1], MFrag[2], MFrag[3]], [RFrag[0], RFrag[1], RFrag[2], RFrag[
                                    3]], [VFrag[0], VFrag[1], VFrag[2], VFrag[3]]]]  # Masses to delete,New Masses,################
                            ###################END COLLISION CODE##############
        ###################COLLISION FRAGMENTATION CODE####################

        #####################REMOVING MASSLESS BODIES#############
        # Eliminar masas iguales a cero
        indcero = (Mi == 0)
        indcero = np.arange(len(Mi))[indcero]
        MDiv, CDiv, Mi, vi, ri = RMArray(
            indcero, MDiv, CDiv, Mi, vi, ri, False)

        if(len(NewThings) != 0):
            # print("Frag",len(NewThings))
            for ll in range(0, len(NewThings)):
                for ee in range(0, 4):
                    Mi = np.concatenate(
                        (Mi, np.array([NewThings[ll][0][ee]])), axis=0)
                    ri = np.concatenate(
                        (ri, np.array([[NewThings[ll][1][ee][0], NewThings[ll][1][ee][1]]])), axis=0)
                    vi = np.concatenate(
                        (vi, np.array([[NewThings[ll][2][ee][0], NewThings[ll][2][ee][1]]])), axis=0)
            # print(len(Mi),len(ri),len(vi))
            # print(Mi,ri)
        #####################END REMOVING MASSLESS BODIES#############

        ################CALCULATING dt's##################
        DelT = []
        for jk in range(0, len(Mi)):
            for jg in range(0, jk + 1):
                if jk != jg:
                    RelR = np.array(
                        [ri[jk][0] - ri[jg][0], ri[jk][1] - ri[jg][1]])
                    RelV = np.array(
                        [vi[jk][0] - vi[jg][0], vi[jk][1] - vi[jg][1]])
                    DelT += [abs((np.dot(RelR, RelR)) /
                                 (4 * np.dot(RelR, RelV)+RegDt))]
        # print(DelT)

        if DelT == []:
            print('All planets collapsed in', systemName)
            POrt=1
            Coun100+=1
            #Graficar todas las masas
            GN=len(Mi)+CDiv
            Gt=math.log10(t)

            Rx=ri[0:len(Mi),0:1]/F
            Ry=ri[0:len(Mi),1:2]/F
            ###################CHECKING ECCENTRICITY####################
            tol = 1e-7
            mu_param = Mstar
            h = ri[:, 0] * vi[:, 1] - ri[:, 1] * vi[:, 0]
            h2 = h * h

            r_norm = np.linalg.norm(ri, axis=1)
            v2 = np.sum(vi*vi, axis=1)

            # Calcular eje semimayor
            a_semi_major = mu_param*r_norm/(2*mu_param - r_norm*v2)
            
            # Calcular excentricidad
            e2 = 1 - h2 / mu_param * (1. / a_semi_major)
            e2[np.abs(e2) < tol] = 0
            
            e = np.sqrt(e2)
            e[e < tol] = 0
            
            a_semi_major[e >= 1] = r_norm[e >= 1]
            b_semi_minor = a_semi_major*np.sqrt(1 - e**2)

            theta = np.arctan((ri[:,1]*a_semi_major)/(ri[:,0]*b_semi_minor))
            ###################END CHECKING ECCENTRICITY################
            
            ###################PLANET RADIUS################
            earth_radius = 6.3781e8 # centimeters 
            radius_planet = (Mi/(4*np.pi/3*den))**(1./3.)
            radius_planet = radius_planet / (cm_to_LU * earth_radius)
            ###################END PLANET RADIUS#############
            
            ###################ENERGY################
            # especific_orbital_energy = v2/2 - mu_param/r_norm
            # print(especific_orbital_energy)
            ###################END ENERGY############

            outfile2 = open(wdir + folder_save_states + systemName[11:-4] + '/' + 'PltAllGraph%05d.txt'%Coun100,'w')
            orbital_elements = open(wdir + folder_save_states + systemName[11:-4] + '/' + 'OrbitalElements%05d.txt'%Coun100,'w')
            
            for d in range(0,len(Mi)):
                outfile2.write('%.5f , %.5f , %.5f %s'%(Mi[d]*Ms/Me,Rx[d],Ry[d],'\n'))
                orbital_elements.write('%.7f, %.7f, %.7f, %.7f, %.7f, %.7f %s '%(Mi[d]*Ms/Me, e[d], a_semi_major[d]/F, b_semi_minor[d]/F, theta[d], radius_planet[d], '\n'))
            outfile2.close()
            orbital_elements.close()
            #Graficar masa mayores a un %

            indcGra=np.where(Mi>Mpla)[0]
            GMi=Mi[indcGra]
            GNP=len(GMi)
            Gri=ri[indcGra]
            Mtot=(np.sum(Mi)+MDiv) #Masa total del sistema
            Mtot2=(np.sum(Mi)) #Masa total de los elmentos math.sin diverger
            Mtot3=(np.sum(GMi)) #Masa total de 80%
            outfile2 = open(wdir + folder_save_states + systemName[11:-4] + '/' + 'Nmass.txt','a')
            outfile2.write("%.5f , %d , %d , %d, %d, %.20f , %.20f , %.20f , %.5f %s " % (t,len(Mi)+CDiv,CDiv,len(Mi),len(GMi),Mtot,Mtot2,Mtot3,Mstar,'\n'))
            outfile2.close()
            break

        minDelT = np.nanmin(DelT)

        # print(minDelT,Standt)
        dt = np.nanmin(np.array([minDelT, Standt]))

        
        # if(dt!=0.001):
        #    print(dt)
        ################END CALCULATING dt's##################

        ################EVOLUTION#################
        # Leap Froog Integrator
        Fij = []
        # Positions at n+1/2
        ri += 0.5 * dt * vi
        # Calculate the forces at n+1/2
        for rr in range(0, len(Mi)):
            b = np.array([0., 0.])
            for ff in range(0, len(Mi)):
                if (rr != ff):
                    rij = np.array(
                        [ri[rr][0] - ri[ff][0], ri[rr][1] - ri[ff][1]])
                    # Remove Mo to get the aceleration
                    b += Mi[ff] / (np.linalg.norm(rij)**2 +
                                   Reg**2)**(3. / 2) * rij
            Fij += [b]  # Total aceleration over i-th particle
        Fij = np.array(Fij)
        # Velocities and Posictions at n+1
        #print(len(ri), rr, ff, len(Mi))

        vi += (-Fij - Mstar /
               (np.linalg.norm(ri, axis=1).reshape(len(ri), 1))**3 * ri) * dt

        ri += 0.5 * dt * vi
        t += dt
    	###############END EVOLUTION##################

        ###################CHECKING ECCENTRICITY & ENERGY####################
        tol = 1e-7 # Tolerance factor
        mu_param = Mstar # Gravitational parameter

        h = ri[:, 0] * vi[:, 1] - ri[:, 1] * vi[:, 0] # Angular momentum per unit mass
        h2 = h * h

        r_norm = np.linalg.norm(ri, axis=1)
        v2 = np.sum(vi*vi, axis=1)

        a_semi_major = mu_param*r_norm/(2*mu_param - r_norm*v2) # Semimajor axis: to save

        e2 = 1 - h2 / mu_param * (1. / a_semi_major)
        # e2 = 1 + h2/mu_param*(v2/mu_param - 2/r_norm)

        e2[np.abs(e2) < tol] = 0

        e = np.sqrt(e2)

        e[e < tol] = 0 # Eccentricity: To save

        a_semi_major[e >= 1] = r_norm[e >= 1]

        # especific_orbital_energy = v2/2 - mu_param/r_norm
        ###################END CHECKING ECCENTRICITY & ENERGY################

    	##############REMOVING BODIES##############
        # Remove SuperEccentric orbits
        super_eccentric_index = e - 1 >= -tol
        super_eccentric_index = np.arange(len(Mi))[super_eccentric_index]
        MDiv,CDiv,Mi,vi,ri= RMArray(super_eccentric_index,MDiv,CDiv,Mi,vi,ri,True)
        
        # Remove near bodies to star
        Ndist0 = np.linalg.norm(ri, axis=1)
        Nind0 = Ndist0 < 0.1 * F
        Nind0 = np.arange(len(Mi))[Nind0]
        MDiv, CDiv, Mi, vi, ri = RMArray(Nind0, MDiv, CDiv, Mi, vi, ri, True)

        # Remove diverge bodies Velocity
        Ndist1 = np.linalg.norm(ri, axis=1)
        Nvelp1 = np.linalg.norm(vi, axis=1)
        NTmass1 = np.sum(Mi)
        NVEES1 = np.power(2. * (Mstar) / Ndist1, 0.5)  # Velocidad de scape del sistema
        Nind1 = Ndist1 > 30. * F
        Nindrm1 = Nvelp1 >= NVEES1  # Objetos a remover
        Nindrm1_vel = np.arange(len(Mi))[Nindrm1]
        Nindrm1 = np.arange(len(Mi))[Nind1 & Nindrm1]
        MDiv, CDiv, Mi, vi, ri = RMArray(Nindrm1, MDiv, CDiv, Mi, vi, ri, True)
        MDiv, CDiv, Mi, vi, ri = RMArray(Nindrm1_vel, MDiv, CDiv, Mi, vi, ri, True)
        ##############END REMOVING BODIES##############

        ###############SAVING DATES#######################
        minimum_record = t % record
        if minimum_record < dt and t != dt: # 100.025 % 100, dt = 1e-3
            POrt=1
            Coun100+=1
            #Graficar todas las masas
            GN=len(Mi)+CDiv
            Gt=math.log10(t)

            Rx=ri[0:len(Mi),0:1]/F
            Ry=ri[0:len(Mi),1:2]/F

            b_semi_minor = a_semi_major*np.sqrt(1 - e**2)
        
            theta = np.arctan((ri[:,1]*a_semi_major)/(ri[:,0]*b_semi_minor))

            ###################PLANET RADIUS################
            earth_radius = 6.3781e8 # centimeters 
            radius_planet = (Mi/(4*np.pi/3*den))**(1./3.)
            radius_planet = radius_planet / (cm_to_LU * earth_radius)
            ###################END PLANET RADIUS#############
            
            outfile2 = open(wdir + folder_save_states + systemName[11:-4] + '/' + 'PltAllGraph%05d.txt'%Coun100,'w')
            orbital_elements = open(wdir + folder_save_states + systemName[11:-4] + '/' + 'OrbitalElements%05d.txt'%Coun100,'w')
            
            for d in range(0,len(Mi)):
                outfile2.write('%.5f , %.5f , %.5f %s'%(Mi[d]*Ms/Me,Rx[d],Ry[d],'\n'))
                orbital_elements.write('%.7f, %.7f, %.7f, %.7f, %.7f, %.7f %s '%(Mi[d]*Ms/Me, e[d], a_semi_major[d]/F, b_semi_minor[d]/F, theta[d], radius_planet[d], '\n'))
            outfile2.close()
            orbital_elements.close()
            #Graficar masa mayores a un %

            indcGra=np.where(Mi>Mpla)[0]
            GMi=Mi[indcGra]
            GNP=len(GMi)
            Gri=ri[indcGra]
            Mtot=(np.sum(Mi)+MDiv) #Masa total del sistema
            Mtot2=(np.sum(Mi)) #Masa total de los elmentos math.sin diverger
            Mtot3=(np.sum(GMi)) #Masa total de 80%
            outfile2 = open(wdir + folder_save_states + systemName[11:-4] + '/' + 'Nmass.txt','a')
            outfile2.write("%.5f , %d , %d , %d, %d, %.20f , %.20f , %.20f , %.5f %s " % (t,len(Mi)+CDiv,CDiv,len(Mi),len(GMi),Mtot,Mtot2,Mtot3,Mstar,'\n'))
            outfile2.close()
        ################END SAVING DATES##############
    # print('End', time.strftime("%c"))
    done = time.time()
    # print('Elapsed time: ' + str(int(done - start)) + ' s')
    Pfile = 'Parameters.txt'
    outP = open(wdir + folder_save_states + systemName[11:-4] + '/' + Pfile, 'w')
    outP.write("tf [yr] = %e %s " % (tf, '\n'))
    outP.write("dt [yr] = %e %s " % (dt, '\n'))
    outP.write('Elapsed time: ' + str(int(done - start)) + ' s')
    outP.close()

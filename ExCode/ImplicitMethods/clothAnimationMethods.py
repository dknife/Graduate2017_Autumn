from OpenGL.GLUT import *
from OpenGL.GL import *
from OpenGL.GLU import *

import numpy as np
import Cloth
import MatrixVector as MV
import math

nW = 5
nH = 5
clothWidth = 10
clothHeight = 10


Verts, Springs, L0 = Cloth.create(nW, nH, clothWidth, clothHeight)

Vel = np.zeros(shape=(len(Verts),3), dtype = np.float64)
Force = np.zeros(shape=(len(Verts),3), dtype = np.float64)
Mass = [10.0 for x in range(len(Verts))]
gravity = np.array([0, -9.8, 0])
K = 5000000.0
damp = 0.1
dt = 0.01

def flatten(x) :
    flattenArray = []
    for lists in x :
        for elements in lists :
            flattenArray.append(elements)
    return flattenArray

def computeForces(springs, verts) :
    global Vel, K
    force = np.zeros((len(verts),3), dtype=np.float64)
    for i in range(0, len(springs)):
        idx0, idx1 = springs[i][0], springs[i][1]
        p0  , p1   = verts[idx0], verts[idx1]
        L = np.linalg.norm(p1 - p0)
        dir = (p1 - p0) / L
        Vrel = Vel[idx0] - Vel[idx1]
        F = K * (L-L0[i])* dir - damp*Vrel  # spring force + damping force
        force[idx0] = force[idx0] + F
        force[idx1] = force[idx1] - F

    for i in range(0, len(verts)) :
        force[i] = force[i] + gravity*Mass[i]
    return force


def IntegrateJacobiIteration(locs, springs, force, vel, h) :
    global L0, Mass, gravity, K, nW

    Mii  = np.zeros(shape=(len(locs), 3, 3))
    Jxii = np.zeros(shape=(len(locs), 3, 3))
    Jvii = np.zeros(shape=(len(locs), 3, 3))
    Jxij = np.zeros(shape=(len(springs), 3, 3))
    Jvij = np.zeros(shape=(len(springs), 3, 3))

    for i in range(0, len(locs)) :
        Mii[i] = Mass[i] * np.identity(3)

    for s in range(0, len(springs)) :
        i,j = springs[s][0], springs[s][1]
        xji = locs[j]-locs[i]
        lt  = np.linalg.norm(xji)
        Jxij[s] = K*(((lt-L0[s])/lt)*np.identity(3) + (L0[s]/(lt**3.0))*np.outer(xji, xji))
        Jvij[s] = damp*np.identity(3)
        Jxii[i], Jxii[j] = Jxii[i] - Jxij[s], Jxii[j] - Jxij[s]
        Jvii[i], Jvii[j] = Jvii[i] - Jvij[s], Jvii[j] - Jvij[s]

    b = np.zeros(shape=(len(locs), 3), dtype = np.float64)
    deltaV = np.zeros(shape=(len(locs), 3), dtype=np.float64)
    for i in range(len(locs)) :
        b[i] = force[i]*h + h*h*Jxii[i].dot(vel[i])
    for s in range(len(springs)) :
        i,j = springs[s][0], springs[s][1]
        b[i] = b[i] + h * h * Jxij[s].dot(vel[j])
        b[j] = b[j] + h * h * Jxij[s].dot(vel[i])


    Minv = np.zeros(shape=(len(locs), 3, 3))
    for i in range(len(locs)) :
        Minv[i] = np.linalg.inv(Mii[i] - h * h * Jxii[i] - h * Jvii[i])
        #deltaV[i] = Minv[i].dot(b[i])
        if i is not 0 and i is not nW - 1:
            deltaV[i] = Minv[i].dot(b[i])
        else:
            deltaV[i] = np.zeros(shape=(3,))

    for iter in range(1) :
        for s in range(len(springs)):
            i, j = springs[s][0], springs[s][1]
            OffDiag = -h*h * Jxij[s] - h*Jvij[s]
            b[i] -= OffDiag.dot(deltaV[j])
            b[j] -= OffDiag.dot(deltaV[i])
        for i in range(len(locs)):
            if i is not 0 and i is not nW-1 :
                deltaV[i] = Minv[i].dot(b[i])
            else :
                deltaV[i] = np.zeros(shape=(3,))


    return flatten(deltaV)


def IntegrateImplicit(locs, springs, force, vel, h) :
    global L0, Mass, gravity, K, nW
    longF = MV.convertLongVector(force)
    bigM = np.zeros(shape=(len(locs) * 3, len(locs) * 3))
    Jx   = np.zeros(shape=(len(locs) * 3, len(locs) * 3))
    Jv   = np.zeros(shape=(len(locs) * 3, len(locs) * 3))
    for i in range(0, len(locs)) :
        MV.set33SubMatrix(bigM, i, i, Mass[i] * np.identity(3))

    Jxii = np.zeros(shape=(len(locs), 3, 3))
    Jvii = np.zeros(shape=(len(locs), 3, 3))

    for spring in range(0, len(springs)) :
        i,j = springs[spring][0], springs[spring][1]
        xji = locs[j]-locs[i]
        lt  = np.linalg.norm(xji)
        Jxij = K*(((lt-L0[spring])/lt)*np.identity(3) + (L0[spring]/(lt**3.0))*np.outer(xji, xji))
        Jvij = damp*np.identity(3)
        Jxii[i], Jxii[j] = Jxii[i] - Jxij, Jxii[j] - Jxij
        Jvii[i], Jvii[j] = Jvii[i] - Jvij, Jvii[j] - Jvij
        MV.set33SubMatrixSymetric(Jx, i, j, Jxij)
        MV.set33SubMatrixSymetric(Jv, i, j, Jvij)

    for i in range(0, len(locs)):
        MV.set33SubMatrix(Jx, i, i, Jxii[i])
        MV.set33SubMatrix(Jv, i, i, Jvii[i])

    bigM = bigM - h*h*Jx - h*Jv
    bigMinv = np.linalg.inv(bigM)
    longV = MV.convertLongVector(vel)
    deltaV = bigMinv.dot(longF*h + (h**2)*Jx.dot(longV))
    return deltaV

def compute_Int_f_dt(locs, springs, vel, h) :
    global L0, Mass, gravity, K, nW, damp

    phi  = np.zeros(shape=(len(locs), 3))

    for s in range(0, len(springs)) :
        i,j = springs[s][0], springs[s][1]
        xji, vji = locs[j] - locs[i], vel[j] - vel[i]
        lt  = np.linalg.norm(xji)
        mi, mj = Mass[i], Mass[j]
        M = mi*mj/(mi + mj)
        d = lt - L0[s]
        dddt = 0
        xjiHat = np.array([1,0,0])

        if lt > 0.0:
            xjiHat = xji / lt
            dldt = np.dot(xjiHat, vji)
        else:
            dldt = 0.0

        omega = math.sqrt(K/M) * math.sqrt(1.0 - damp**2.0/(4.0*K*M))
        alpha = damp/ (2.0*M)

        A = math.exp(alpha*3.141592/omega)*np.sqrt(d**2.0 + (dldt**2.0) * (mi * mj) / (K * (mi + mj)))

        dDA = 0.0
        arcsine = 0.0
        if A != 0.0:
            dDA = d / A;
            arcsine = np.arcsin(dDA)
        if d > 0 and dldt < 0:
            arcsine = np.pi - arcsine
        if d <= 0:
            if dldt > 0:
                arcsine = 2.0 * np.pi + arcsine
            else:
                arcsine = np.pi - arcsine

        T0 = arcsine / omega
        T1 = T0 + h

        Et = np.exp(-alpha * T0)
        Eh = np.exp(-alpha * h )
        coeff = -  A * K  /  (omega**2.0 + alpha**2.0)
        st = np.sin(omega * T0)
        ct = np.cos(omega * T0)
        sh = np.sin(omega * h)
        ch = np.cos(omega * h)

        phi_ij = coeff * Et * (
        st * (Eh * (alpha * ch - omega * sh) - alpha) + ct * (Eh * (alpha * sh + omega * ch) - omega))
        phi_ij = phi_ij / M
        phi[i] += (phi_ij*mj/((mi+mj))) * xjiHat
        phi[j] -= (phi_ij*mi/((mi+mj))) * xjiHat

    for i in range(len(locs)) :
        phi[i] += h*gravity

    return phi

def IntegrateHybrid(locs, springs, force, vel, h) :
    global L0, Mass, gravity, K, nW, damp

    Mii  = np.zeros(shape=(len(locs), 3, 3))
    Jxii = np.zeros(shape=(len(locs), 3, 3))
    Jxij = np.zeros(shape=(len(springs), 3, 3))


    for i in range(0, len(locs)) :
        Mii[i] = Mass[i] * np.identity(3)

    for s in range(0, len(springs)) :
        i,j = springs[s][0], springs[s][1]
        xji = locs[j]-locs[i]
        lt  = np.linalg.norm(xji)
        Jxij[s] = K*(((lt-L0[s])/lt)*np.identity(3) + (L0[s]/(lt**3.0))*np.outer(xji, xji))
        Jxii[i], Jxii[j] = Jxii[i] - Jxij[s], Jxii[j] - Jxij[s]


    b = np.zeros(shape=(len(locs), 3), dtype = np.float64)
    deltaV = np.zeros(shape=(len(locs), 3), dtype=np.float64)

    phi = compute_Int_f_dt(locs, springs, vel, h)

    for i in range(len(locs)) :
        b[i] = force[i]*h +  h*h*Jxii[i].dot(vel[i])
    for s in range(len(springs)) :
        i,j = springs[s][0], springs[s][1]
        b[i] = b[i] + h * h * Jxij[s].dot(vel[j])
        b[j] = b[j] + h * h * Jxij[s].dot(vel[i])


    for i in range(len(locs)) :
        b[i] = 0.5*b[i] + 0.5*phi[i]

    Minv = np.zeros(shape=(len(locs), 3, 3))
    for i in range(len(locs)) :
        Minv[i] = np.linalg.inv(Mii[i] - h * h * Jxii[i])
        deltaV[i] = Minv[i].dot(b[i])

    deltaV = flatten(deltaV)
    return deltaV

def IntegrateHarmonic(locs, springs, force, vel, h) :
    global L0, Mass, gravity, K, nW, damp

    deltaV = np.zeros(shape=(len(locs), 3), dtype=np.float64)
    phi = compute_Int_f_dt(locs, springs, vel, h)

    deltaV = flatten(phi)
    return deltaV

def IntegrateEuler(locs, springs, force, vel, h) :
    global L0, Mass, gravity, K, nW, damp

    deltaV = np.zeros(shape=(len(locs), 3), dtype=np.float64)
    for i in range(0, len(locs)):
        deltaV[i] = force[i]*h/Mass[i]

    deltaV = flatten(deltaV)
    return deltaV

def updateState(locs, vel, deltaV, h) :
    for i in range(0, len(locs)):
        vel[i] = vel[i] + deltaV[i*3:i*3+3]
        if i is not 0 and i is not nW-1:
            locs[i] = locs[i] + vel[i] * h
        else :
            vel[i] = np.zeros(shape=(3,))



def moveImplicit():
    global Verts, Springs, Vel, Force, Mass, K, damp, L0, dt, gravity
    # Force Initialization

    Force = computeForces(Springs, Verts)
    #deltaV = IntegrateEuler(Verts, Springs, Force, Vel, dt)
    #deltaV = IntegrateImplicit(Verts, Springs, Force, Vel, dt)
    #deltaV = IntegrateJacobiIteration(Verts, Springs, Force, Vel, dt)
    deltaV = IntegrateHybrid(Verts, Springs, Force, Vel, dt)
    #deltaV = IntegrateHarmonic(Verts, Springs, Force, Vel, dt)



    updateState(Verts, Vel, deltaV, dt)
    glutPostRedisplay()  ## Hey! Draw Function! Draw!



def draw():
    global clothWidth, clothHeight
    glClear(GL_COLOR_BUFFER_BIT)
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    gluPerspective(60, 1.0, 0.1, 1000)
    glMatrixMode(GL_MODELVIEW)
    glLoadIdentity()
    gluLookAt(clothWidth*0.5, 1.0 - clothHeight*0.5, clothHeight*1.5, clothWidth*0.5, 1.0 - clothHeight*0.5, 0.0, 0.0, 1.0, 0.0)

    for i in range(0, len(Springs)):
        Cloth.drawSpring(Verts[Springs[i][0]], Verts[Springs[i][1]])
    glFlush()

glutInit(sys.argv)
glutInitDisplayMode(GLUT_SINGLE | GLUT_RGBA)
glutInitWindowSize(250, 250)
glutInitWindowPosition(100, 100)
glutCreateWindow(b"OpenGL with Python")
glClearColor(0.5, 0.5, 0.5, 1)
glutDisplayFunc(draw)
glutIdleFunc(moveImplicit)
glutMainLoop()

# End of program
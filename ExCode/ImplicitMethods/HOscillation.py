from OpenGL.GLUT import *
from OpenGL.GL import *
from OpenGL.GLU import *
import Shader as sc

import numpy as np
import Cloth
import MatrixVector as MV
import math

nW = 35
nH = 35
clothWidth = 1
clothHeight = 1


Verts, Springs, L0 = Cloth.create(nW, nH, clothWidth, clothHeight)

Vel = np.zeros(shape=(len(Verts),3), dtype = np.float64)
Force = np.zeros(shape=(len(Verts),3), dtype = np.float64)
Mass = [2.0 for x in range(len(Verts))]
gravity = np.array([0, -9.8, 0])
K = 2200000.0
damp = 0.1
dt = 0.005



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
        phi[i] += (phi_ij*mj/((mi+mj))) * xjiHat
        phi[j] -= (phi_ij*mi/((mi+mj))) * xjiHat

    for i in range(len(locs)) :
        phi[i] += h*gravity/Mass[i]

    return phi

def Integrate(locs, springs, force, vel, h) :
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


    for i in range(0, len(locs)):
        vel[i] = vel[i] + deltaV[i]
        if i is not 0 and i is not nW-1:
            locs[i] = locs[i] + vel[i] * h
        else :
            vel[i] = np.zeros(shape=(3,))



def moveImplicit():
    global Verts, Springs, Vel, Force, Mass, K, damp, L0, dt, gravity
    # Force Initialization
    Force = computeForces(Springs, Verts)
    Integrate(Verts, Springs, Force, Vel, dt)
    glutPostRedisplay()  ## Hey! Draw Function! Draw!



def drawSpring(loc1, loc2):
    glColor3f(1, 1, 1)
    glBegin(GL_LINES)
    glVertex3f(loc1[0], loc1[1], loc1[2])
    glVertex3f(loc2[0], loc2[1], loc2[2])
    glEnd()


def flatten(x) :
    flattenArray = []
    for lists in x :
        for elements in lists :
            flattenArray.append(elements)
    return flattenArray

shader = None
buffer = None
def draw():
    global clothWidth, clothHeight, Verts, shader, buffer

    vertices = flatten(Verts)
    if shader == None :
        shader = sc.Shader()
    if buffer == None :
        buffer = glGenBuffers(1)
        glBindBuffer(GL_ARRAY_BUFFER, buffer)

    glBufferData(GL_ARRAY_BUFFER,
                 len(vertices) * 4,  # byte size
                 (ctypes.c_float * len(vertices))(*vertices),
                 GL_STATIC_DRAW)

    glClear(GL_COLOR_BUFFER_BIT)
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    gluPerspective(60, 1.0, 0.1, 1000)
    glMatrixMode(GL_MODELVIEW)
    glLoadIdentity()
    gluLookAt(0, 1.0 - clothHeight*0.5, clothHeight*1.5, clothWidth*0.5, 1.0 - clothHeight*0.5, 0.0, 0.0, 1.0, 0.0)

    shader.begin()
    glEnableClientState(GL_VERTEX_ARRAY);
    glBindBuffer(GL_ARRAY_BUFFER, buffer);
    glVertexPointer(3, GL_FLOAT, 0, None);
    glDrawArrays(GL_POINTS, 0, len(vertices));
    glDisableClientState(GL_VERTEX_ARRAY);
    for i in range(0, len(Springs)):
        drawSpring(Verts[Springs[i][0]], Verts[Springs[i][1]])
    shader.end()

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
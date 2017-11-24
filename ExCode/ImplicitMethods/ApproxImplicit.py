from OpenGL.GLUT import *
from OpenGL.GL import *
from OpenGL.GLU import *

import numpy as np
import Cloth
import MatrixVector as MV

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
    return force


def Integrate(locs, springs, force, vel, h) :
    global L0, Mass, gravity, K, nW

    for i in range(0, len(locs)) :
        force[i] = force[i] + gravity

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
        drawSpring(Verts[Springs[i][0]], Verts[Springs[i][1]])
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
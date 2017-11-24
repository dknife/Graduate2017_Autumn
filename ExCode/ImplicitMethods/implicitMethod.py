from OpenGL.GLUT import *
from OpenGL.GL import *
from OpenGL.GLU import *

import numpy as np
import Cloth
import MatrixVector as MV

nW = 5
nH = 5
clothWidth = 10
clothHeight = 10


Verts, Springs, L0 = Cloth.create(nW, nH, clothWidth, clothHeight)

Vel = np.zeros(shape=(len(Verts),3), dtype = np.float64)
Force = np.zeros(shape=(len(Verts),3), dtype = np.float64)
Mass = [1.0 for x in range(len(Verts))]
gravity = np.array([0, -9.8, 0])
K = 5000000.0
damp = 0.1
dt = 0.01



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


def Integrate(locs, force, vel, h) :
    global Springs, L0, Mass, gravity, K, nW

    for i in range(0, len(locs)) :
        force[i] = force[i] + gravity
    longF = MV.convertLongVector(force)

    bigM = np.zeros(shape=(len(locs) * 3, len(locs) * 3))
    Jx   = np.zeros(shape=(len(locs) * 3, len(locs) * 3))
    Jv   = np.zeros(shape=(len(locs) * 3, len(locs) * 3))

    for i in range(0, len(locs)) :
        MV.set33SubMatrix(bigM, i, i, Mass[i] * np.identity(3))

    Jxii = np.zeros(shape=(len(locs), 3, 3))
    Jvii = np.zeros(shape=(len(locs), 3, 3))

    for spring in range(0, len(Springs)) :
        i,j = Springs[spring][0], Springs[spring][1]
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
    Integrate(Verts, Force, Vel, dt)
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
from OpenGL.GLUT import *
from OpenGL.GL import *
from OpenGL.GLU import *
import numpy as np

Verts = np.array(
    [[0, 2, 0], [1, 2, 0], [2, 2, 0],
     [0, 2, 1], [1, 2, 1], [2, 2, 1],
     [0, 2, 2], [1, 2, 2], [2, 2, 2]], np.float64
)

Vel = []
Force = []
Mass = []

for i in range(0, len(Verts)):
    Vel.append([0.0, 0.0, 0.0])
    Force.append([0.0, 0.0, 0.0])
    Mass.append(1.0)

Vel = np.array(Vel)
Force = np.array(Force)

Springs = np.array(
    [[0, 1], [1, 2],
     [3, 4], [4, 5],
     [6, 7], [7, 8],
     [0, 3], [3, 6],
     [1, 4], [4, 7],
     [2, 5], [5, 8],
     [0,4], [1,3],
     [1,5], [2,4],
     [3,7], [4,6],
     [4,8], [5,7]], np.int32)
gravity = np.array([0, -9.8, 0])

K = 900.0
damp = 0.0
sqrt2 = np.sqrt(2.0)
L0 = [ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, sqrt2, sqrt2, sqrt2, sqrt2, sqrt2, sqrt2, sqrt2, sqrt2]

dt = 0.03

def drawAxes():
    glBegin(GL_LINES)
    glColor3f(1, 0, 0)
    glVertex3f(0, 0, 0)
    glVertex3f(1, 0, 0)
    glColor3f(0, 1, 0)
    glVertex3f(0, 0, 0)
    glVertex3f(0, 1, 0)
    glColor3f(0, 0, 1)
    glVertex3f(0, 0, 0)
    glVertex3f(0, 0, 1)
    glEnd()


def drawSpring(loc1, loc2):
    glColor3f(1, 1, 1)
    glPushMatrix()
    glTranslatef(loc1[0], loc1[1], loc1[2])
    glutWireSphere(0.1, 10, 10)
    glPopMatrix()
    glPushMatrix()
    glTranslatef(loc2[0], loc2[1], loc2[2])
    glutWireSphere(0.1, 10, 10)
    glPopMatrix()
    glBegin(GL_LINES)
    glVertex3f(loc1[0], loc1[1], loc1[2])
    glVertex3f(loc2[0], loc2[1], loc2[2])
    glEnd()

def computeForces(springs, verts) :
    force = np.zeros((len(verts),3), dtype=np.float64)
    for i in range(0, len(springs)):
        idx0 = springs[i][0]
        idx1 = springs[i][1]
        p0 = verts[idx0]
        p1 = verts[idx1]
        L = np.linalg.norm(p1 - p0)
        dir = (p1 - p0) / L
        d = L - L0[i]
        F = K * d * dir
        f_damp = -np.dot(Vel[idx0],dir)*damp*dir
        force[idx0] = force[idx0] + F + f_damp
        f_damp =  np.dot(Vel[idx1], -dir) * damp * dir
        force[idx1] = force[idx1] - F + f_damp
    return force

def copyVectors(verts, n):
    v = np.zeros((n, 3), dtype=np.float64)
    for i in range(n):
        v[i] = verts[i]

    return v

def updateLocation(locs, force, vel, h):
    for i in range(0, len(locs)):
        vel[i] = vel[i] + (force[i] / Mass[i] + gravity) * h
        if i is not 0 and i is not 2:
            locs[i] = locs[i] + vel[i] * h

def moveRK():
    global Verts, Springs, Vel, Force, Mass, K, damp, L0, dt, gravity

    v = copyVectors(Verts, len(Verts))
    vel1 = copyVectors(Vel, len(Vel))
    force1 = computeForces(Springs, v)

    v = copyVectors(Verts, len(Verts))
    vel2 = copyVectors(vel1, len(vel1))
    updateLocation(v, force1, vel2, 0.5*dt)
    force2 = computeForces(Springs, v)

    v = copyVectors(Verts, len(Verts))
    vel3 = copyVectors(vel1, len(vel2))
    updateLocation(v, force2, vel3, 0.5 * dt)
    force3 = computeForces(Springs, v)

    v = copyVectors(Verts, len(Verts))
    vel4 = copyVectors(vel3, len(vel3))
    updateLocation(v, force3, vel4, dt)
    force4 = computeForces(Springs, v)

    Force = (force1+2.0*force2+2.0*force3+force4)/6.0
    Vel = (vel1 + 2.0*vel2 + 2.0*vel3 + vel4)/6.0
    updateLocation(Verts, Force, Vel, dt)

    glutPostRedisplay()  ## Hey! Draw Function! Draw!

def moveMidPoint() :
    global Verts, Springs, Vel, Force, Mass, K, damp, L0, dt, gravity
    v = copyVectors(Verts, len(Verts))
    velTemp = copyVectors(Vel, len(Vel))
    Force = computeForces(Springs, v)
    updateLocation(v, Force, velTemp, 0.5*dt)
    Force = computeForces(Springs, v)
    updateLocation(Verts, Force, Vel, dt)
    glutPostRedisplay()  ## Hey! Draw Function! Draw!

def moveEuler():
    global Verts, Springs, Vel, Force, Mass, K, damp, L0, dt, gravity

    # Force Initialization
    for i in range(0, len(Force)):
        Force[i] = np.array([0, 0, 0])

    # Compute spring force for each spring and accumulate it for each vert
    for i in range(0, len(Springs)):
        idx0 = Springs[i][0]
        idx1 = Springs[i][1]
        p0 = Verts[idx0]
        p1 = Verts[idx1]
        L = np.linalg.norm(p1 - p0)
        dir = (p1 - p0) / L
        d = L - L0[i]
        F = K * d * dir
        f_damp = -np.dot(Vel[idx0],dir)*damp*dir
        Force[idx0] = Force[idx0] + F + f_damp
        f_damp =  np.dot(Vel[idx1], -dir) * damp * dir
        Force[idx1] = Force[idx1] - F + f_damp

    for i in range(0, len(Verts)):
        Vel[i] = Vel[i] + (Force[i]/ Mass[i] + gravity) * dt
        if i is not 0 and i is not 2:
            Verts[i] = Verts[i] + Vel[i] * dt

    glutPostRedisplay()  ## Hey! Draw Function! Draw!


def draw():
    global Verts, Springs
    glClear(GL_COLOR_BUFFER_BIT)
    glMatrixMode(GL_MODELVIEW)
    glLoadIdentity()
    gluLookAt(0.5, 2.0, 5.0, 0.0, 2.0, 0.0, 0.0, 1.0, 0.0)
    drawAxes()
    for i in range(0, len(Springs)):
        drawSpring(Verts[Springs[i][0]], Verts[Springs[i][1]])
    glFlush()


def GLinit():
    glClearColor(0.5, 0.5, 0.5, 1)

def reshape(w, h) :
    aspectRatio = w/h
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    gluPerspective(60, aspectRatio, 0.1, 1000)
    glViewport(0, 0, w, h)


glutInit(sys.argv)
glutInitDisplayMode(GLUT_SINGLE | GLUT_RGBA)
glutInitWindowSize(250, 250)
glutInitWindowPosition(100, 100)
glutCreateWindow(b"OpenGL with Python")
GLinit()
glutDisplayFunc(draw)
glutReshapeFunc(reshape)
#glutIdleFunc(moveEuler)
#glutIdleFunc(moveMidPoint)
glutIdleFunc(moveRK)
glutMainLoop()

# End of program
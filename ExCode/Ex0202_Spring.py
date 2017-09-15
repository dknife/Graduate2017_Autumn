from OpenGL.GLUT import *
from OpenGL.GL import *
from OpenGL.GLU import *
import numpy as np

Loc1 = np.array([0.0,1.0,0.0])
Loc2 = np.array([1.0,1.0,0.0])
Vel1 = np.array([0.01,0,0])
Vel2 = np.array([0,0,0])
gravity = np.array([0,-9.8,0])

mass1 = 1
mass2 = 1
#K = 29999.0
K = 1000.0
damp = 0.1
L0 = 1.0
dt = 0.005


def drawAxes():
    glBegin(GL_LINES)
    glColor3f(1,0,0)
    glVertex3f(0,0,0)
    glVertex3f(1,0,0)
    glColor3f(0, 1, 0)
    glVertex3f(0, 0, 0)
    glVertex3f(0, 1, 0)
    glColor3f(0, 0, 1)
    glVertex3f(0, 0, 0)
    glVertex3f(0, 0, 1)
    glEnd()

def drawSpring(loc1, loc2) :
    glColor3f(1,1,1)
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

def move() :
    global Loc1, Loc2, Vel1, Vel2, L0, K, mass1, mass2, dt, damp, gravity
    L = np.linalg.norm(Loc2-Loc1)
    dir = (Loc2 - Loc1) / L
    d = L - L0
    f1 = K * d * dir
    f_damp = - damp * Vel2
    #print(np.linalg.norm(f_damp))
    f2 = -f1 + f_damp
    a1 = gravity + f1 / mass1
    a2 = gravity + f2 / mass2


    Vel1 = Vel1 + a1 * dt
    #Loc1 = Loc1 + Vel1 * dt
    Vel2 = Vel2 + a2 * dt
    Loc2 = Loc2 + Vel2 * dt
    glutPostRedisplay()  ## Hey! Draw Function! Draw!

def draw():
    glClear(GL_COLOR_BUFFER_BIT)
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    gluPerspective(60, 1.0, 0.1, 1000);
    glMatrixMode(GL_MODELVIEW)
    glLoadIdentity()
    gluLookAt(0.5,2.0,5.0, 0.0,0.0,0.0, 0.0,1.0,0.0)
    drawAxes()
    drawSpring(Loc1, Loc2)
    glFlush()

def GLinit() :
    glClearColor(0.5, 0.5, 0.5, 1)

glutInit(sys.argv)
glutInitDisplayMode(GLUT_SINGLE | GLUT_RGBA)
glutInitWindowSize(250, 250)
glutInitWindowPosition(100, 100)
glutCreateWindow(b"OpenGL with Python")
GLinit()
glutDisplayFunc(draw)
glutIdleFunc(move)
glutMainLoop()

# End of program
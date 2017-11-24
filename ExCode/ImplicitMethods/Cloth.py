from OpenGL.GLUT import *
from OpenGL.GL import *
from OpenGL.GLU import *

import numpy as np
import math

def create(nW, nH, width, height) :
    dx = width/(nW-1)
    dz = height / (nH-1)
    dDiag = math.sqrt(dx**2.0 + dz**2.0)

    Verts = np.zeros(shape=(nW * nH, 3))
    Springs = np.zeros(shape=((nW - 1) * nH + (nH - 1) * nW + (nW - 1) * (nH - 1) * 2, 2), dtype=np.int32)

    for i in range(0, nW * nH):
        Verts[i] = [(i % nW) * dx, 1, (int(i / nW)) * dz]


    L0 = np.zeros(shape=(len(Springs),))

    for i in range(0, (nW - 1) * nH):
        row = int(i / (nW - 1))
        col = i % (nW - 1)
        Springs[i] = [row * nW + col, row * nW + col + 1]
        L0[i] = dx

    for i in range(0, (nH - 1) * nW):
        row = int(i / (nW))
        col = i % (nW)
        Springs[(nW - 1) * nH + i] = [row * nW + col, (row + 1) * nW + col]
        L0[(nW - 1) * nH + i] = dz

    for i in range(0, (nW - 1) * (nH - 1)):
        row = int(i / (nW - 1))
        col = i % (nW - 1)
        Springs[(nW - 1) * nH + (nH - 1) * nW + 2 * i] = [row * nW + col, (row + 1) * nW + col + 1]
        Springs[(nW - 1) * nH + (nH - 1) * nW + 2 * i + 1] = [row * nW + col + 1, (row + 1) * nW + col]
        L0[(nW - 1) * nH + (nH - 1) * nW + 2 * i] = dDiag
        L0[(nW - 1) * nH + (nH - 1) * nW + 2 * i + 1] = dDiag

    return Verts, Springs, L0

def drawSpring(loc1, loc2):
    glColor3f(1, 1, 1)
    glBegin(GL_LINES)
    glVertex3f(loc1[0], loc1[1], loc1[2])
    glVertex3f(loc2[0], loc2[1], loc2[2])
    glEnd()
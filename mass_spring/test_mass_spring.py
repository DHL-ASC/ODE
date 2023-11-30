import sys
sys.path.append('/Users/joachim/texjs/lva/ws2324/ScientificComputing/ASC-ODE/build/mass_spring')

from mass_spring import *
import ngsolve.bla


mss = MassSpringSystem3d()
mss.gravity = (0,0,-9.81)

mA = mss.Add (Mass(1, (1,0,0)))
mB = mss.Add (Mass(2, (2,0,0)))
f1 = mss.Add (Fix( (0,0,0)) )
mss.Add (Spring(1, 10, (f1, mA)))
mss.Add (Spring(1, 20, (mA, mB)))


print ("state = ", mss.GetState())

Simulate (mss, 0.1, 10)

print ("state = ", mss.GetState())


Simulate (mss, 0.1, 10)

print ("state = ", mss.GetState())

for m in mss.masses:
    print (m.mass, m.pos)

mss.masses[0].mass = 5

for m in mss.masses:
    print (m.mass, m.pos)

from cmath import pi
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
import vpython as vp
from vpython import sphere, vector, color, cylinder, rate

def wykres_potencjal(t,l,m,g,f,A,fi, fidot):
    w = 2*pi*f
    fi_check = []
    u_check = []
    dt = 0.01
    while t < 100:
        fiddot = -(A * w * 2 * sp.cos(t * w) + g) * sp.sin(fi) / l
        #fiddot = -g / l * sp.sin(fi) - 1 / 2 * (A * w / l) ** 2 * sp.sin(fi) * sp.cos(fi)
        fidot = fidot + fiddot * dt
        fi = fi + fidot * dt
        U = (-1) * g * m * l * sp.cos(fi) + (m * (A * 2 * pi * f * sp.sin(fi)) ** 2) / 4
        fi_check.append(fi*180/pi)
        u_check.append(U*1000)
        t = t + dt

    plt.plot(fi_check, u_check)
    plt.title("Efektywny potencjał wahadła w funkcji 'powolnej' składowej kąta wychylenia")
    plt.xlabel("kąt [stopnie]")
    plt.ylabel("potencjał efektywny U [mJ]")
    plt.grid()
    plt.ylim(-4, 4)
    plt.xlim(-300, 300)

    plt.show()


def wykres_fi(t, l, g, f, A, fi, fidot):
    w = 2*pi*f
    check = []
    time = []
    dt = 0.01
    while t < 50:
        #fiddot = -(A * w * 2 * sp.cos(t * w) + g) * sp.sin(fi) / l
        fiddot = -g/l*sp.sin(fi)-1/2*(A * w/l)**2 * sp.sin(fi)*sp.cos(fi)
        fidot = fidot + fiddot * dt
        fi = fi + fidot * dt
        check.append(fi*180/pi)
        time.append(t)
        t = t + dt

    plt.plot(time, check)
    plt.grid()
    plt.title("'Powolna' składowa kąta wychylenia w funkcji czasu")
    plt.xlabel("czas t[s]")
    plt.ylabel("kąt [stopnie]")
    plt.ylim(0,200)
    plt.xlim(0,50)

    plt.show()

def symulacja(t, l, g, f, A, fi, fidot):
    dt = 0.001
    w = 2 * pi * f
    x = l * sp.sin(fi)
    y = -A * sp.cos(w * t) - l * sp.cos(fi)
    node = sphere(pos=A * vector(0, sp.sin(w * t), 0), radius=0.0005, make_trial=True)
    mass = sphere(pos=vector(x, y, 0), radius=0.0005, color=color.yellow, make_trail=True)
    rod = cylinder(pos=node.pos, axis=mass.pos - node.pos, radius=0.00025)

    while t < 5:
        rate(10000)
        fiddot = -(A * w ** 2 * sp.cos(t * w) + g) * sp.sin(fi) / l
        fidot = fidot + fiddot * dt
        fi = fi + fidot * dt
        x = l * sp.sin(fi)
        y = -A * sp.cos(w * t) - l * sp.cos(fi)
        node.pos = A * vector(0, sp.sin(w * t), 0)
        mass.pos = vector(x, y, 0)
        rod.pos = node.pos
        rod.axis = mass.pos - node.pos
        t = t + dt

if __name__=="__main__":
    g = 9.8
    m = 0.005
    l = 0.045
    A = 0.016
    #f = 10.5
    t = 0
    dt = 0.01
    fi = pi - 0.1
    fidot = 0


    print("Wahadło Kapicy\nPodaj następujące wielkości:")
    # l = float(input("- dlugość ramienia = "))
    # A = float(input("- amplituda = "))
    fg = 1 / (2 * pi * A) * sp.sqrt(2 * g * l)
    print("Częstotliwość graniczna dla podanych parametrów =", fg)
    f = float(input("Podaj czestotliwosc = "))

    #fg = 1/(2*pi*A)*sqrt(2*g*l)



    #print(rownanie_ruchu(t, l,m,g,f,A,fi))
    #print(potencjal_efektywny(t,l,m,g,f,A,fi))
    wykres_fi(t, l, g, f, A, fi, fidot)
    wykres_potencjal(t, l, m, g, f, A, fi, fidot)
    symulacja(t,l,g,f,A,fi, fidot)

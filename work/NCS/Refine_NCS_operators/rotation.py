from __future__ import division
import numpy as nm

def Rx(theta):
  """
  theta is the angle in degree
  """
  th = nm.pi*theta/180.
  rx = nm.round(nm.matrix([[1,0,0],
                           [0,nm.cos(th),-nm.sin(th)],
                           [0,nm.sin(th),nm.cos(th)]]),5)
  return rx

def Ry(psi):
  """
  psi is the angle in degree
  """
  p = nm.pi*psi/180.
  ry = nm.round(nm.matrix([[nm.cos(p),0,nm.sin(p)],
                           [0,1,0],
                           [-nm.sin(p),0,nm.cos(p)]]),5)
  return ry

def Rz(phi):
  """
  phi is the angle in degree
  """
  p = nm.pi*phi/180.
  rz = nm.round(nm.matrix([[nm.cos(p),-nm.sin(p),0],
                           [nm.sin(p),nm.cos(p),0],
                           [0,0,1]]),5)
  return rz

def R(theta,psi,phi):
  return Rx(theta).dot(Ry(psi).dot(Rz(phi)))


theta=90
psi=90
phi=90
print nm.linalg.det(Rx(theta))
print
print nm.linalg.det(Ry(psi))
print
print nm.linalg.det(Rz(phi))
print
print Rx(90).dot(Rx(90))
print
print R(theta=theta,psi=psi,phi=phi)
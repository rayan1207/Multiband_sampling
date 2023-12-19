import numpy
import math
import os

headstr=" beta Remu Immu H kdim kx ky reW imW"
#5.0 0.0 2 3.14159 1.0472 0 0.628318

Remu=0.0

kdim=2
kx1=0*numpy.pi/8
ky1=0*numpy.pi/4
reW=0
b=(5,)
kgrid=4
tol=1e-8
Nk=numpy.pi/2.0
ANk=numpy.pi 
z=0.0

data=[]

list1=numpy.arange(0.01,1.021,0.02)
list2=numpy.arange(2.0,100,2)
list3=numpy.concatenate((list1,list2))

for Remu in (0.0,):#(-1.3,4.0,10.0):#(0,5.0001,10.2):
  for H in numpy.arange(0,3.01,3.2):
    for kx in numpy.arange(0,ANk+0.001,ANk*10):
      for ky in numpy.arange(0,ANk,ANk*10):
        for beta in b: #list3: #({5}):#list3: #numpy.arange(0.01, 1.0, 0.01):#({2}): #(.05,0.1,0.2, 0.3, 0.4, 0.5, 0.6, 0.75, 1, 2 ,3, 4, 5,8): #(2.7,2.8,3,3.5, 4,4.5, 5,5.5,6,8, 10,12, 13, 15, 20,30, 40, 80, 120): #numpy.arange(2, , 1): #(1,2,5,10,11,12,13,14):#numpy.arange(1, 6, 10):
          for n_mat in numpy.arange(0,2,1):
              for Immu in (0.5,1):#-10,11,1): #(0,.110, .01): #numpy.arange(0,5.5,.125+tol):
            #imW=numpy.pi/beta*(2*n_mat) #+1) #+1)
                imW=numpy.pi/beta*(2*n_mat+1) #0.02
            # imW=n_mat
                reW=0 #n_mat #n_mat #n_mat
                data.append((beta, Remu, Immu, H, kdim,numpy.pi/2,numpy.pi/2,0,imW))               
# for Remu in (0.0,):#(-1.3,4.0,10.0):#(0,5.0001,10.2):
#   for H in numpy.arange(0,2,10):
#     for kx in numpy.arange(ANk,ANk+tol,ANk):
#       for ky in numpy.arange(ANk,ANk+tol,ANk):
#         for beta in b: #list3: #({5}):#list3: #numpy.arange(0.01, 1.0, 0.01):#({2}): #(.05,0.1,0.2, 0.3, 0.4, 0.5, 0.6, 0.75, 1, 2 ,3, 4, 5,8): #(2.7,2.8,3,3.5, 4,4.5, 5,5.5,6,8, 10,12, 13, 15, 20,30, 40, 80, 120): #numpy.arange(2, , 1): #(1,2,5,10,11,12,13,14):#numpy.arange(1, 6, 10):
#           for n_mat in numpy.arange(0,0.01,10):#-10,11,1): #(0,.110, .01): #numpy.arange(0,5.5,.125+tol):
#             #imW=numpy.pi/beta*(2*n_mat) #+1) #+1)
#             imW=numpy.pi/beta*(2*n_mat+1) #0.02
#             # imW=n_mat
#             reW=0 #n_mat #n_mat #n_mat
#             data.append((beta, 0.01, Immu, H, kdim, kx,ky,0,0))
# for Remu in (0.0,):#(-1.3,4.0,10.0):#(0,5.0001,10.2):
#   for H in numpy.arange(0,2,10):
#     for kx in numpy.flip(numpy.arange(0,ANk+tol,ANk*100)):
#       for ky in numpy.flip(numpy.arange(ANk/kgrid,ANk+tol,ANk/kgrid)):
#         for beta in b: #list3: #({5}):#list3: #numpy.arange(0.01, 1.0, 0.01):#({2}): #(.05,0.1,0.2, 0.3, 0.4, 0.5, 0.6, 0.75, 1, 2 ,3, 4, 5,8): #(2.7,2.8,3,3.5, 4,4.5, 5,5.5,6,8, 10,12, 13, 15, 20,30, 40, 80, 120): #numpy.arange(2, , 1): #(1,2,5,10,11,12,13,14):#numpy.arange(1, 6, 10):
#           for n_mat in numpy.arange(0,0.01,10):#-10,11,1): #(0,.110, .01): #numpy.arange(0,5.5,.125+tol):
#             #imW=numpy.pi/beta*(2*n_mat) #+1) #+1)
#             imW=numpy.pi/beta*(2*n_mat+1) #0.02
#             # imW=n_mat
#             reW=0 #n_mat #n_mat #n_mat
#             data.append((beta, 0.01, Immu, H, kdim, ky,ky,0,0))
            
print(data)
numpy.savetxt("ext_vars.dat", data, fmt='%f', header=headstr)

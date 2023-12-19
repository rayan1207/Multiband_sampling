# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 18:05:50 2023

@author: Rayan
"""

import numpy as np
import matplotlib.pyplot as plt
axis_font = {'fontname':'Arial', 'size':'14'}
pi = np.pi
tol= 0.01
tp =0.5
tpp=0.0
ksize = 6
om = np.arange(-8,8+0.001,0.25)

qx,qy = np.meshgrid(np.arange(0,pi+tol,pi/ksize),np.arange(0,pi+tol,pi/ksize))
qx1= np.concatenate((qx[0,:],qx[1:-1,-1],np.flip(np.diag(qx))))
qy1= np.concatenate((qy[0,:],qy[1:-1,-1],np.flip(np.diag(qy))))



 
omega,qx = np.meshgrid(om,qx1)
omega,qy = np.meshgrid(om,qy1)
omega,l = np.meshgrid(om,range(len(qy1)))



def spectral(real,imag,omega,qx,qy,tp,tpp, i=0):
    e = -2*(1+tpp)*(np.cos(qx)+np.cos(qy))-tp
    den = (omega - e - real)**2 + (i+imag)**2
    num = imag+i
    return -num/(den*np.pi)


tp33_v0 = np.loadtxt('tp03_33_U2_v05.dat')
tp33_v15= np.loadtxt('tp03_33_U3_v075.dat')
tp11_v0 = np.loadtxt('tp03_11_U2_v05.dat')
tp11_v15= np.loadtxt('tp03_11_U3_v075.dat')


# print('omega i produce is \n')
# print(omega)
# print('omega produced is \n')
# print(np.reshape(tp33_v0[:,0] , (19, 65)))

size = ksize*3+1

def rotate(matrix):
    a=ksize*2+1
    b = ksize*3+2
    # Split the matrix into 17x65 and 9x65 parts
    part1 = matrix[:ksize*2+1, :]
    part2 = matrix[a:b, :]

    # Rotate the 9x65 part by 180 degrees
    rotated_part2 = np.flipud(part2)

    # Join both parts together vertically
    result_matrix = np.vstack((part1, rotated_part2))

    return result_matrix

def makegrid(data,tp,tpp):
    re=np.reshape(data[:,3] , (size, 65))
    imag=np.reshape(data[:,5] , (size, 65))
    omega=np.reshape(data[:,0] , (size, 65))
    # print(omega)
    qx =np.reshape(data[:,1] , (size, 65))
    qy =np.reshape(data[:,2] , (size, 65))
    return [re,rotate(spectral(re,imag,omega,qx,qy,tp,tpp))]




print(omega)


tp33_v15=makegrid(tp33_v15,-0.5,-0.3)
tp11_v15=makegrid(tp11_v15,0.5,0.3)


tp33_v0=makegrid(tp33_v0,-0.5,-0.3)
tp11_v0=makegrid(tp11_v0,0.5,0.3)




A=tp33_v15[1]+tp11_v15[1]
B= tp33_v0[1]+tp11_v0[1]

plt.figure(figsize=[12,8])
plt.subplot(1,2,2)
plt.pcolormesh(l,omega,A,vmin=-0.5,vmax=0.5)
plt.xticks([0, ksize, ksize*2, ksize*3], ['$[0,0]$',
                            '$[\pi,0]$', '$[\pi,\pi]$', '$[0,0]$'])


plt.text(4,-7.8,r'$t_{\perp}$='+str(tp) +' '+ r'$t^\prime_{\perp}$='+str(tpp)+' \n U=3.00, V=0.75',**axis_font,bbox=dict(facecolor='yellow', alpha=0.5))
plt.colorbar()
plt.subplot(1,2,1)
plt.ylabel(r'$A(k,\omega)$',**axis_font)
plt.pcolormesh(l,omega,B,vmin=-0.5,vmax=0.5)
plt.xticks([0, ksize, ksize*2, ksize*3], ['$[0,0]$',
                            '$[\pi,0]$', '$[\pi,\pi]$', '$[0,0]$'])
plt.text(4,-7.8,r'$t_{\perp}$='+str(tp) +' '+ r'$t^\prime_{\perp}$='+str(tpp)+' \n U=2.00, V=0.50',**axis_font,bbox=dict(facecolor='yellow', alpha=0.5))
plt.colorbar()
plt.show()


plt.plot(A[:,32],'x-',label='a')
plt.plot(B[:,32],'x-',label='b')
plt.xticks([0, ksize, ksize*2, ksize*3], ['$[0,0]$',
                            '$[\pi,0]$', '$[\pi,\pi]$', '$[0,0]$'])
plt.show()
           



# def dos(matrix):
#     num_rows = len(matrix)
#     num_cols = len(matrix[0])
#     col_sums = [0] * num_cols

#     for row in matrix:
#         for j in range(num_cols):
#             col_sums[j] += row[j]

#     return col_sums

# dos = dos(tp33_v0[1]+tp11_v0[1])


# plt.plot(om,dos)
# plt.show()


    


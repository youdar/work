import numpy as nm
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

 
def Rx(theta):
    return nm.matrix([[1,0,0],[0,nm.cos(theta),-nm.sin(theta)],[0,nm.sin(theta),nm.cos(theta)]])

def Ry(theta):
    return nm.matrix([[nm.cos(theta),0,nm.sin(theta)],[0,1,0],[-nm.sin(theta),0,nm.cos(theta)]])

def Rz(theta):
    return nm.matrix([[nm.cos(theta),-nm.sin(theta),0],[nm.sin(theta),nm.cos(theta),0],[0,0,1]])

def vec_to_list(x):
    #return map(round((x.transpose()).tolist()[0],5)
    return [round(y,6) for y in x]

def mat_to_list(m,t):
    #return map(round((x.transpose()).tolist()[0],5) [(m.flatten()).tolist()[0],[0,0,0]]
    return [[round(y,6) for y in (m.flatten()).tolist()[0]],t]

def collect_pos_data(x,datax,datay,dataz):
    if type(x) == nm.matrix:
        x = (x.flatten()).tolist()[0]
    elif type(x) == nm.array:  
        x = (x.flatten()).tolist()[0]
    datax.append(round(x[0],6))
    datay.append(round(x[1],6))
    dataz.append(round(x[2],6))
    
    
    return datax,datay,dataz

def run():
    '''
    Produce data for test in multimer class
    '''
    datax = []
    datay = []
    dataz = []     
    theta = nm.pi 
    # output data point
    transformation_data = [[[1,0,0,0,1,0,0,0,1],[0,0,0]]]  
    # Base point
    #x = [1,1,1]
    #datax,datay,dataz = collect_pos_data(x,datax,datay,dataz)
    # convert to matrix type
    #xv = nm.matrix(x).transpose()
    # start collecting data
    m = Rx(theta/2)
    #datax,datay,dataz = collect_pos_data(m*xv,datax,datay,dataz)
    transformation_data.append(mat_to_list(m,[0,0,0]))
    #
    m = Ry(theta/2)
    #datax,datay,dataz = collect_pos_data(m*xv,datax,datay,dataz)
    transformation_data.append(mat_to_list(m,[0,0,0]))  
    #
    m = Rz(theta/2)
    #datax,datay,dataz = collect_pos_data(m*xv,datax,datay,dataz)
    transformation_data.append(mat_to_list(m,[0,0,0]))    
    #
    m = Ry(theta/2)
    #datax,datay,dataz = collect_pos_data(m*xv,datax,datay,dataz)
    transformation_data.append(mat_to_list(m,[0,0,0]))  
    #
    m = Rz(theta/2)*Ry(theta/2)
    #datax,datay,dataz = collect_pos_data(m*xv,datax,datay,dataz)
    transformation_data.append(mat_to_list(m,[0,0,0]))     
    #
    m = Rz(theta/3)
    #datax,datay,dataz = collect_pos_data(m*xv,datax,datay,dataz)
    transformation_data.append(mat_to_list(m,[0,0,0]))     
    #
    m = Rz(theta/3)*Rz(theta/3)
    #datax,datay,dataz = collect_pos_data(m*xv,datax,datay,dataz)
    transformation_data.append(mat_to_list(m,[0,0,0]))    
    # Translation 
    t = [0,0.5,0]
    tv = nm.matrix(t).transpose()
    #datax,datay,dataz = collect_pos_data(xv+tv,datax,datay,dataz)
    transformation_data.append([[1,0,0,0,1,0,0,0,1],t])
    # Rotation + Translation
    t = [0,0,-1]
    tv = nm.matrix(t).transpose()
    m = Rz(theta/3)*Rz(theta/3) 
    #datax,datay,dataz = collect_pos_data(m*xv+tv,datax,datay,dataz)
    transformation_data.append(mat_to_list(m,t))    
    
    
    # Positions to test:
    atoms_xyz = [[1,1,1],
                 [94.618, -5.253, 91.582],
                 [62.395, 51.344, 80.786],
                 [39.954, 51.526 ,72.372]]
    #fig = plt.figure()
    #ax = Axes3D(fig)   
    #ax.plot(datax,datay,dataz,'.')
    #plt.xlabel('x')
    #plt.ylabel('y')      
    #plt.show() 
    
    # check that the we are not missing any data   
    #assert len(transformation_data) == len(datax),'Missing data: len(transformation_data)={}, len(datax)={} '\
           #.format(len(transformation_data),len(datax))
    
    
    print 'test pdb BIOMT data'
    print '='*50
    for i in range(len(transformation_data)):
        [d11,d12,d13] = transformation_data[i][0][:3]
        [d21,d22,d23] = transformation_data[i][0][3:6]
        [d31,d32,d33] = transformation_data[i][0][6:]
        [x,y,z] = transformation_data[i][1]
        #print 'REMARK 350   BIOMT1  {0:2} {1:9} {2:9}  {3:9}        {4:8} '.format(i+1,d11,d12,d13,datax[i])
        print 'REMARK 350   BIOMT1  %2i % 1.6f % 1.6f % 1.6f       % 1.6f ' % (i+1,d11,d12,d13,x)
        print 'REMARK 350   BIOMT2  %2i % 1.6f % 1.6f % 1.6f       % 1.6f ' % (i+1,d21,d22,d23,y)
        print 'REMARK 350   BIOMT3  %2i % 1.6f % 1.6f % 1.6f       % 1.6f ' % (i+1,d31,d32,d33,z)
    
    
    print '='*50
    #print 'ATOM    426  CA  LEU A   1      94.618  -5.253  91.582  1.00 87.10           C '
    xyz_out_list = []
    for xyz in atoms_xyz:
        for i in range(len(transformation_data)):
            xyz_in = nm.matrix(xyz).transpose()
            m = nm.matrix(transformation_data[i][0]).reshape(3,3)
            t = nm.matrix(transformation_data[i][1] ).transpose()
            xyz_out = vec_to_list(m*xyz_in + t)
            xyz_out_list.append(xyz_out)
            print 'ATOM    426  CA  LEU A   1      {: .3f}  {: .3f}  {: .3f}  1.00 87.10           C '.format(xyz_out[0],xyz_out[1],xyz_out[2])
            
    print xyz_out_list                                                                          
    #print '='*50
    #print 'expected coordinates after reconstruction'
    #print '='*50
    #for i in range(len(datax)):
        #print '% 1.6f % 1.6f % 1.6f'%(datax[i],datay[i],dataz[i])
        
        
    #f1 = 'REMARK 350   BIOMT1  12 -0.815628 -0.473258  0.332832        0.00000 '
    #f2 = 'REMARK 350   BIOMT2  12 -0.815628 -0.473258  0.332832        0.00000 '
    #f3 = 'REMARK 350   BIOMT3  12 -0.815628 -0.473258  0.332832        0.00000 '  
    # ATOM    426  CA  LEU A   1      94.618  -5.253  91.582  1.00 87.10           C  
    # ATOM   3733  CA  ARG B   2      62.395  51.344  80.786  1.00107.25           C  
    # HETATM 3754  C1  EDO A   3      39.954  51.526  72.372  0.33 60.93           C  
    
    

if __name__=="__main__":
    run()
    
    

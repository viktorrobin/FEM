from __future__ import print_function
import os
import platform
import numpy as np
import random
import shutil

#http://dunne.uni-hd.de/VisuSimple/documents/vtkfileformat.html
#http://people.sc.fsu.edu/~jburkardt/data/vtk/vtk.html

syst = platform.system()
base=os.path.dirname(os.path.abspath(__file__))


## Import parameters
File=base+'/Output/parameters.txt'
data=np.loadtxt(File,skiprows=0)

nelem=int(data[0])
nnode=int(data[1])
npoint=int(data[2])
nincs=int(data[3])+1
ntype=int(data[4])
ncrit=int(data[5])

#num_files=print(len([name for name in os.listdir('.') if os.path.isfile("/vtkmovie")]))
#path, dirs, files = os.walk("vtkmovie/").next()

#Name of the model
if ncrit==1:
    name_crit='Tresca'
elif ncrit==2:
    name_crit='Von Mises'
elif ncrit==3:
    name_crit='Mohr-Coulomb'
elif ncrit==4:
    name_crit='Drucker-Prager'
elif ncrit==5:
    name_crit='Modified Cam Clay'
elif ncrit==6:
    name_crit='Lime Treated Soil Model'


import sys
for iincs in range(0,nincs):

    sys.stdout.write("\r{0} %".format((float(iincs)/float(nincs+1.0))*100))
    sys.stdout.flush()

    f = open('vtkmovie/vtk_'+str(iincs)+'.vtk','w')

    print("# vtk DataFile Version 1.0", file=f)
    print("2D Unstructured Grid of Linear Triangles", file=f)
    print("ASCII", file=f)
    print("", file=f)
    print("DATASET UNSTRUCTURED_GRID", file=f)
    print("POINTS "+str(npoint)+" float", file=f)

    #Write coordinates


    # Import coordinates
    #File=base+'/Output/coordinates.txt'
    File=base+'/Output/output_increm/displacements/'+str(iincs)+'.dat'

    #data=np.loadtxt(File,skiprows=1)
    data=np.loadtxt(File,skiprows=0)
    #dat2=np.loadtxt(File2,skiprows=1)

    coord=[]
    xmin=9999.0
    xmax=-9999.0
    ymin=9999.0
    ymax=-9999.0

    for i in range(0,len(data)):
        coord.append([])
        #coord[i].append((float(data[i,0])))
        #coord[i].append(float(data[i,1]))

        coord[i].append( float(data[i,0])+float(data[i,2]) )
        coord[i].append( float(data[i,1])+float(data[i,3]) )

        #if (float(data[i,0]) < xmin):
        #    xmin = float(data[i,0])
        #
        #if (float(data[i,0]) > xmax):
        #    xmax = float(data[i,0])
        #
        #if (float(data[i,1]) < ymin):
        #    ymin = float(data[i,1])
        #
        #if (float(data[i,1]) > ymax):
        #    ymax = float(data[i,1])

    for i in range(0,len(coord)):
        print(str(coord[i][0])+"  "+str(coord[i][1])+"  "+"0.0", file=f)


    print("", file=f)

    #Connectivity table
    # Import connectivity table
    File=base+'/Output/conntable.txt'
    data=np.loadtxt(File,skiprows=0)
    lnods=[]

    if (nelem>1):
        for i in range(0,len(data)):
            lnods.append([])
            for j in range(0,nnode):
                lnods[i].append((int(data[i,j+1])-1))


    if (nelem==1):
        for i in range(0,nelem+1):
            lnods.append([])
            for j in range(0,nnode):
                lnods[i].append((int(data[j+1])-1))

    print("CELLS "+str(nelem)+" "+str(nelem*(nnode+1)), file=f)

    #Re-order nodes
    for i in range(0,len(lnods)):
        myorder=[0,2,4,6,1,3,5,7]
        lnods_temp = [ lnods[i][j] for j in myorder]
        lnods[i]=lnods_temp



    for i in lnods:
        var=str(i)
        var=var[0:-1]
        var = var.replace("[",str(nnode)+"  ")
        newstr = var.replace(",", " ")
        print(newstr, file=f)

    print("", file=f)

    print("CELL_TYPES "+str(nelem), file=f)
    if (nnode==8):
        type=23
    elif (nnode==4):
        type=9
    elif (nnode==4):
        type=28


    for i in range(0,nelem):
        print(type, file=f)


    #START PLOTING THE RESULTS: POINTS
    print("", file=f)
    print("POINT_DATA "+str(npoint), file=f)

    #-------------------------
    #        STRAINS
    #-------------------------

    File=base+'/Output/output_increm/strains/'+str(iincs)

    x=[]
    y=[]
    eps_xx=[]
    eps_yy=[]
    eps_xy=[]
    eps_zz=[]
    ps_max=[]
    epstn=[]

    data=np.loadtxt(File,skiprows=1)
    for i in range(0,len(data)):
        if data[i,0]==0:
            x.append((float(data[i,3])))
            y.append(float(data[i,4]))
            eps_xx.append(float(data[i,5]))
            eps_yy.append(float(data[i,6]))
            eps_xy.append(float(data[i,7]))
            eps_zz.append(float(data[i,8]))
            ps_max.append(float(data[i,9]))
            epstn.append(float(data[i,12]))

    strains=[eps_xx,eps_yy,eps_xy,eps_zz,ps_max,epstn]

    if ntype==3:
        names=['eps_rr','eps_zz','eps_rz','eps_tt','eps_max_principal','eps_p']
    else:
        names=['eps_xx','eps_yy','eps_xy','eps_zz','eps_max_principal','eps_p']
    compt=-1

    for res in strains:
        compt += 1
        print("SCALARS " +names[compt]+ " float", file=f)
        print("LOOKUP_TABLE default", file=f)
        for i in res:
            print(i, file=f)
        print("", file=f)

    #-------------------------
    #        STRESSES
    #-------------------------

    File=base+'/Output/output_increm/stresses/'+str(iincs)

    x=[]
    y=[]
    s_xx=[]
    s_yy=[]
    s_xy=[]
    s_zz=[]
    s_max=[]

    data=np.loadtxt(File,skiprows=1)
    for i in range(0,len(data)):
        if data[i,0]==0:
            x.append((float(data[i,3])))
            y.append(float(data[i,4]))
            s_xx.append(float(data[i,5]))
            s_yy.append(float(data[i,6]))
            s_xy.append(float(data[i,7]))
            s_zz.append(float(data[i,8]))
            s_max.append(float(data[i,9]))


    stresses=[s_xx,s_yy,s_xy,s_zz,s_max]
    if ntype==3:
        names=['s_rr','s_zz','s_rz','s_tt','s_max_principal']
    else:
        names=['s_xx','s_yy','s_xy','s_zz','s_max_principal']

    compt=-1

    for res in stresses:
        compt += 1
        print("SCALARS " +names[compt]+ " float", file=f)
        print("LOOKUP_TABLE default", file=f)
        for i in res:
            print(i, file=f)
        print("", file=f)

    #------------------------------
    #   DEVIATORIC STRESS SXX-SYY
    #------------------------------


    print("SCALARS " +names[compt]+ " float", file=f)
    print("LOOKUP_TABLE default", file=f)
    for i,j in zip(s_xx,s_yy):
        print(i, file=f)
    print("", file=f)


    #NODAL DISPLACEMENTS
    File=base+'/Output/output_increm/displacements/'+str(iincs)+'.dat'
    data=np.loadtxt(File,skiprows=0)

    ux=[]
    uy=[]
    u_norme=[]

    for i in range(0,len(data)):
        ux.append((float(data[i,2])))
        uy.append(float(data[i,3]))
        u_norme.append( np.sqrt( (float(data[i,2]))**2 + (float(data[i,3]))**2 ) )

    strains=[ux,uy,u_norme]
    names=["ux","uy","u_norme"]

    compt=-1

    for res in strains:
        compt += 1
        print("SCALARS " +names[compt]+ " float", file=f)
        print("LOOKUP_TABLE default", file=f)
        for i in res:
            print(i, file=f)

        print("", file=f)

    #NODAL FORCES
    File=base+'/Output/output_increm/reactions/'+str(iincs)+'.dat'
    data=np.loadtxt(File,skiprows=0)

    fx=[]
    fy=[]
    f_norme=[]

    for i in range(0,len(data)):
        fx.append((float(data[i,1])))
        fy.append(float(data[i,2]))
        f_norme.append( np.sqrt( (float(data[i,1]))**2 + (float(data[i,2]))**2 ) )

    strains=[fx,fy,f_norme]
    names=["Fx","Fy","F_norme"]

    compt=-1

    for res in strains:
        compt += 1
        print("SCALARS " +names[compt]+ " float", file=f)
        print("LOOKUP_TABLE default", file=f)
        for i in res:
            print(i, file=f)

        print("", file=f)

    if ncrit==5 or ncrit==6:
        #YIELD STRESS AT NODAL POINTS
        File=base+'/Output/output_increm/plastic_strains/'+str(iincs)
        data=np.loadtxt(File,skiprows=0)

        yield_node=[]
        epstnp=[]
        deviatoric_stress=[]
        effective_mean_stress=[]

        for i in range(0,len(data)):
            yield_node.append((float(data[i,5])))
            epstnp.append((float(data[i,6])))
            deviatoric_stress.append((float(data[i,7])))
            effective_mean_stress.append((float(data[i,8])))

        print("SCALARS yield_stress float", file=f)
        print("LOOKUP_TABLE default", file=f)
        for i in yield_node:
            print(i, file=f)

        print("", file=f)

        #Plastic volumetric deformations AT NODAL POINTS
        print("SCALARS epstn_p float", file=f)
        print("LOOKUP_TABLE default", file=f)
        for i in epstnp:
            print(i, file=f)

        print("", file=f)

        #Deviatoric stress AT NODAL POINTS
        print("SCALARS deviatoric_stress float", file=f)
        print("LOOKUP_TABLE default", file=f)
        for i in deviatoric_stress:
            print(i, file=f)

        print("", file=f)

        #effective mean stress AT NODAL POINTS
        print("SCALARS eff_mean_stress float", file=f)
        print("LOOKUP_TABLE default", file=f)
        for i in effective_mean_stress:
            print(i, file=f)

        print("", file=f)

        #SPECIFIC VOLUME AT NODAL POINTS
        File=base+'/Output/output_increm/specific_volume/'+str(iincs)+'.dat'
        data=np.loadtxt(File,skiprows=0)

        spe_vol=[]


        for i in range(0,len(data)):
            spe_vol.append((float(data[i,5])))

        print("SCALARS specific_volume float", file=f)
        print("LOOKUP_TABLE default", file=f)
        for i in spe_vol:
            print(i, file=f)

        print("", file=f)


    #START PLOTING THE RESULTS: VECTORS
    #DISPLACEMENTS VECTORS
    print("VECTORS displacements_vec float", file=f)
    #print("LOOKUP_TABLE default", file=f)

    for i,j in zip(ux,uy):
        alpha=np.sqrt( (i*i + j*j)*1.5)

        if alpha==0.0:
            alpha=1.0

        alpha=1.0

        print(str(i/alpha)+'  '+str(j/alpha)+'  '+'0.0',file=f)
    #    print(str(ux[i])+"  "+str(uy[i])+"  "+"0.0", file=f)


    f.close()

print("")
print("VTK FILES GENERATED")
    #shutil.copy2('vtkmovie/vtk_'+str(iincs)+'.txt', 'vtkmovie/vtk_'+str(iincs)+'.vtk')

#f = open('myfile.txt','w')
#
#
#f.write('hi there\n') # python will convert \n to os.linesep
#f.close() # you

#f = open('testvtk2.txt','w')
#
#print("# vtk DataFile Version 1.0", file=f)
#print("2D Unstructured Grid of Linear Triangles", file=f)
#print("ASCII", file=f)
#print("", file=f)
#print("DATASET UNSTRUCTURED_GRID", file=f)
#print("POINTS "+str(npoint+1)+" float", file=f)
#
##Write coordinates
#
#
## Import coordinates
#File=base+'/Output/coordinates.txt'
#data=np.loadtxt(File,skiprows=1)
#
#coord=[]
#xmin=9999.0
#xmax=-9999.0
#ymin=9999.0
#ymax=-9999.0
#
#for i in range(0,len(data)):
#    coord.append([])
#    coord[i].append((float(data[i,0])))
#    coord[i].append(float(data[i,1]))
#
#    if (float(data[i,0]) < xmin):
#        xmin = float(data[i,0])
#
#    if (float(data[i,0]) > xmax):
#        xmax = float(data[i,0])
#
#    if (float(data[i,1]) < ymin):
#        ymin = float(data[i,1])
#
#    if (float(data[i,1]) > ymax):
#        ymax = float(data[i,1])
#
#for i in range(0,len(coord)):
#    print(str(coord[i][0])+"  "+str(coord[i][1])+"  "+"0.0", file=f)
#
#print(str(1.0)+"  "+str(1.0)+"  "+"0.0", file=f)
#
#
#print("", file=f)
#
##Connectivity table
## Import connectivity table
#File=base+'/Output/conntable.txt'
#data=np.loadtxt(File,skiprows=0)
#lnods=[]
#
#if (nelem>1):
#    for i in range(0,len(data)):
#        lnods.append([])
#        for j in range(0,nnode):
#            lnods[i].append((int(data[i,j+1])-1))
#
#
#if (nelem==1):
#    for i in range(0,nelem+1):
#        lnods.append([])
#        for j in range(0,nnode):
#            lnods[i].append((int(data[j+1])-1))
#
#print("CELLS "+str(nelem+1)+" "+str(nelem*(nnode+1)+2), file=f)
#
##Re-order nodes
#for i in range(0,len(lnods)):
#    myorder=[0,2,4,6,1,3,5,7]
#    lnods_temp = [ lnods[i][j] for j in myorder]
#    lnods[i]=lnods_temp
#
#
#
#for i in lnods:
#    var=str(i)
#    var=var[0:-1]
#    var = var.replace("[",str(nnode)+"  ")
#    newstr = var.replace(",", " ")
#    print(newstr, file=f)
#
#print("1  21", file=f)
#
#print("", file=f)
#
#print("CELL_TYPES "+str(nelem+1), file=f)
#if (nnode==8):
#    type=23
#elif (nnode==4):
#    type=9
#for i in range(0,nelem):
#    print(type, file=f)
#
#print("1", file=f)
#print("", file=f)
#
##Write strains
#print("POINT_DATA "+str(npoint+1), file=f)
#print("SCALARS pressure float", file=f)
#print("LOOKUP_TABLE default", file=f)
#
##Plot strains
#File=base+'/Output/strains.txt'
#
#x=[]
#y=[]
#eps_xx=[]
#eps_yy=[]
#eps_xy=[]
#eps_zz=[]
#ps_max=[]
#
#uy=[]
#u_norme=[]
#
#data=np.loadtxt(File,skiprows=1)
#for i in range(0,len(data)):
#    x.append((float(data[i,3])))
#    y.append(float(data[i,4]))
#    eps_xx.append(float(data[i,5]))
#    eps_yy.append(float(data[i,6]))
#    eps_xy.append(float(data[i,7]))
#    eps_zz.append(float(data[i,8]))
#    ps_max.append(float(data[i,9]))
#
#
#for i in range(0,npoint):
#    print(random.random()*10, file=f)
#print(1000.0, file=f)
#
#f.close()
#
#
#import shutil
#shutil.copy2('testvtk2.txt', 'testvtk2.vtk')
#
##f = open('myfile.txt','w')
##
##
##f.write('hi there\n') # python will convert \n to os.linesep
##f.close() # you
#

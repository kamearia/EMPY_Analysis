from netgen.occ import * 
def FieldPlot(mesh, field, line, xcomp, ycomp, N): 
    import matplotlib.pylab as plt
    N=100
    xs=line.start
    xe=line.end
    dx=(1./N)*(xe-xs)

    xp=[]
    ypreal=[]
    ypimag=[]
    for n in range(N+1):
        xx=xs+(n*dx)
        pnt=mesh(xx[0], xx[1], xx[2])
        #print( "x= ", x, "  Bz= ", Bfield(pnt)[2].real, "  ", Bfield(pnt)[2].imag)
        xp.append(xx[xcomp])
        ypreal.append(field(pnt)[ycomp].real)
        ypimag.append(field(pnt)[ycomp].imag)

    plt.plot(xp, ypreal )  
    plt.plot(xp, ypimag ) 
    plt.xlabel("x["+str(xcomp)+"]")  # Add x-axis label
    plt.ylabel("Field["+str(ycomp)+"]")  # Add y-axis label
    plt.show()  
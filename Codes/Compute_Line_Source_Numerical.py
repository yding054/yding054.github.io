from Similarity_Wind import Similarity_Wind
from Reformulation_Sigmayz import Reformulation_Sigmayz
from Mean_Plume_Height import Mean_Plume_Height
from Solving_Equations_Iteratively import Solving_Equations_Iteratively
from scipy.special import erf
import numpy as np

def Compute_Line_Source_Numerical(Line_Source,Receptor,Met):
    xb,yb,xe,ye=Line_Source.Rotate_Line_Source(Met)
    xr,yr=Receptor.Rotate_Receptor(Met)
    
    theta=(270-Met.wdir-Met.Ndegree)*np.pi/180
    x0=xb+(yr-yb)*np.tan(theta)
    y0=yb+(xr-xb)/np.tan(theta)
    
    Conc=0
    
    if xb<xe and yb<ye:
        k_line=(ye-yb)/(xe-xb)
        k=(yr-yb)/(xr-xb)
        if k<k_line and xr>xb and xr<xe:
            xd=xr-x0
            R=xd
            
            z_old,z_new=Solving_Equations_Iteratively(Met.ustar,Met.L,Met.sigmav,\
                                                      xd,Line_Source.zs,0.053,1,100)
            
            U,Ue=Similarity_Wind(Met.ustar,z_new,0.053,Met.L,Met.sigmav)
            
            sigmay,sigmaz=Reformulation_Sigmayz(xd,Met.sigmav,Met.ustar,Met.L,Ue)
            
            f=2*Met.sigmav**2/Ue**2
            
            tb=(yb-yr)/(np.sqrt(2)*sigmay)
            

#            VERT=1/(np.sqrt(2*np.pi)*sigmaz)*(np.exp(-0.5*((Line_Source.zs-\
#                   Receptor.zr)/sigmaz)**2)+np.exp(-0.5*((Line_Source.zs+Receptor.zr)\
#            /sigmaz)**2))
            
            VERT=np.sqrt(2/np.pi)/sigmaz
            
#            Finite
            
            HORZ_pl=(1-erf(tb))/2
            HORZ_m=(y0-yb)/(2*np.pi*R)
            
            Conc=(1-f)*Line_Source.Emis/Ue/np.cos(theta)*VERT*HORZ_pl+\
            f*Line_Source.Emis/Ue/np.cos(theta)*VERT*HORZ_m
            
##            Infinite
#        
#            HORZ_pl=1
#            HORZ_m=(y0-yb)/(2*np.pi*R)
#            
#            Conc=(1-f)*Line_Source.Emis/Ue/np.cos(theta)*VERT*HORZ_pl+\
#            f*Line_Source.Emis/Ue/np.cos(theta)*VERT*HORZ_m
            
        elif k<k_line and xr>xe:
            xd=xr-x0
            R=xd
            
            z_old,z_new=Solving_Equations_Iteratively(Met.ustar,Met.L,Met.sigmav,\
                                                      xd,Line_Source.zs,0.053,1,100)
            
            U,Ue=Similarity_Wind(Met.ustar,z_new,0.053,Met.L,Met.sigmav)
            
            sigmay,sigmaz=Reformulation_Sigmayz(xd,Met.sigmav,Met.ustar,Met.L,Ue)
            
            f=2*Met.sigmav**2/Ue**2
            
            te=(ye-yr)/(np.sqrt(2)*sigmay)
            tb=(yb-yr)/(np.sqrt(2)*sigmay)
            

            VERT=1/(np.sqrt(2*np.pi)*sigmaz)*(np.exp(-0.5*((Line_Source.zs-\
                   Receptor.zr)/sigmaz)**2)+np.exp(-0.5*((Line_Source.zs+Receptor.zr)\
            /sigmaz)**2))
            
#            VERT=np.sqrt(2/np.pi)/sigmaz
            
#           Finite
            
            HORZ_pl=(erf(te)-erf(tb))/2
            HORZ_m=(y0-yb)/(2*np.pi*R)
            
            Conc=(1-f)*Line_Source.Emis/Ue/np.cos(theta)*VERT*HORZ_pl+\
            f*Line_Source.Emis/Ue/np.cos(theta)*VERT*HORZ_m
            
#            Infinite
        
#            HORZ_pl=1
#            HORZ_m=(y0-yb)/(2*np.pi*R)
#            
#            Conc=(1-f)*Line_Source.Emis/Ue/np.cos(theta)*VERT*HORZ_pl+\
#            f*Line_Source.Emis/Ue/np.cos(theta)*VERT*HORZ_m
        
#        elif k>k_line and xr>xe:
#            xd=xr-xe+(xe-xb)/2
#            
#            sigmay,sigmaz=Reformulation_Sigmayz(xd,Met)
#            
#            te=(ye-yr)/(np.sqrt(2)*sigmay*(xr-xe))
#            tb=(yb-yr)/(np.sqrt(2)*sigmay*(xr-xb))
#            
#            Conc=np.sqrt(2/np.pi)*Line_Source.Emis/Met.Ue/np.cos(theta)\
#            /sigmaz*(erf(te)-erf(tb))/2
#        
#        elif k>k_line and xr>xb and xr<xe:
#            xd=(xr-xb)/2
#            
#            sigmay,sigmaz=Reformulation_Sigmayz(xd,Met)
#            
#            tb=(yb-yr)/(np.sqrt(2)*sigmay*(xr-xb))
#            
#            if yr>ye:
#                Conc=np.sqrt(2/np.pi)*Line_Source.Emis/Met.Ue/np.cos(theta)\
#                /sigmaz*(-1-erf(tb))/2
#            else:
#                Conc=np.sqrt(2/np.pi)*Line_Source.Emis/Met.Ue/np.cos(theta)\
#                /sigmaz*(1-erf(tb))/2
        
#        
#    elif xe<xb and ye<yb:
#        
#        
#    elif xb<xe and yb>ye:
#        
#        
#    elif xe<xb and ye>yb:
        
        
    return Conc
    
    
from math import tau,pi,cos,sin,atan,acos,asin,sqrt,atan2,tan,tau,hypot
import numpy as np
from Gears.InvoluteGear2D import Gear2D,InvoluteGear2D,InvoluteGear2D3
#from InvoluteGear2D import InvoluteGear2D3

#class InvoluteGear(BaseGear):
class InvoluteGear():
    def __init__(self, m, z1, z2,**kwargs):
        self.P=None
        self.C=None
        self.E=None
        self.A=None
        self.alpha=kwargs.pop("alpha",20*pi/180)
        print(self.alpha*180/pi)
        self.x1=kwargs.pop("x1",0)
        self.x2=kwargs.pop("x2",0)
        self.hap1=kwargs.pop("hap1",1)
        self.hap2=kwargs.pop("hap2",1)
        self.hfp1=kwargs.pop("hfp1",1.25)
        self.hfp2=kwargs.pop("hfp2",1.25)
        self.w1=kwargs.pop("w1",0.3)  #Winkelgeschwindigkeit omega
        self.M1=kwargs.pop("M1",1)  #Moment
        self.mu=kwargs.pop("mu",0.2)  #Moment omega
        self.z1=z1
        self.m=m
        self.z2=z2
       
        self.alpha_w=self.alpha
        #super().__init__(a, z1, z2, hap, hfp, ha, hf, d, df, dw, db, da1, da2, df1, df2,**kwargs)
        self.calculate_parameter()



    def calculate_alpha_w(self):
        from gekko import Gekko
        m=Gekko(remote=True)
        a_w=m.Var(value=20*pi/180,lb=0,ub=60*pi/180)
        alpha=self.alpha
        x1,x2=self.x1,self.x2
        z1,z2=self.z1,self.z2
        rightSide=tan(alpha)-alpha+2*(tan(alpha)*(x1+x2))/(z1+z2)
        print(f"check for Vorzeichen alpha_w")
        m.Equation(m.tan(a_w)-a_w==rightSide)
        m.options.SOLVER=1
        m.solve(disp=False)
        self.alpha_w=a_w.value[0]

    def calculate_parameter(self):
        self.calculate_alpha_w()
        
        x1,x2=self.x1,self.x2
        hap1,hap2,hfp1,hfp2=self.hap1,self.hap2,self.hfp1,self.hfp2
        a_w,alpha=self.alpha_w,self.alpha
        z1,z2=self.z1,self.z2
        
        m=self.m
        self.a=cos(alpha)/cos(a_w)*m/2*(z1+z2)
        self.i=-z2/z1
        self.r1=(z1*m/2)
        self.r2=(z2*m/2)
        self.rw1=(m*z1/2*cos(alpha)/cos(a_w))
        self.rw2=(m*z2/2*cos(alpha)/cos(a_w))
        self.rb1=(m*z1/2*cos(alpha))
        self.rb2=(m*z2/2*cos(alpha))
        self.ra1=(m*z1/2+m*(hap1+x1))
        self.rf1=(m*(z1/2-hfp1+x1))

        # if z2<0:
        #     self.ra2=m*(abs(z2)/2+x2+hfp2)
        #     self.rf2=m*(abs(z2)/2+x2-hap2)
        # else:
        #     self.ra2=m*(z2+x2+hap2)
        #     self.rf2=m*(z2+x2-hfp2)
        self.ra2=m*(z2/2+x2+hap2)
        self.rf2=m*(z2/2+x2-hfp2)

        #check for undercut
        if min(abs(self.ra1),abs(self.rf1))>self.rb1:
            print(f"undercut in gear 1")
        if min(abs(self.ra2),abs(self.rf2))>self.rb2:
            print(f"undercut in gear 2")

        alpha_kopf=acos(self.rb1/self.ra1) #eingriffswinkel am Kopf
        alpha_fuss=acos(self.rb1/self.r1) #eingriffswinkel am Kopf
        inv=lambda a:tan(a)-a
        self.phi_fuss=inv(alpha_fuss)
        self.phi_kopf=inv(alpha_kopf)-self.phi_fuss
        

        print(f"phi_fuss={self.phi_fuss*180/pi:.1f}° vs phi_kopf={self.phi_kopf*180/pi:.1f}")
        self.calculate_points()

    def delta_i(self,r:float):
        alpha_i=acos(self.r1*cos(self.alpha)/r)
        inv=lambda x: tan(x)-x
        return inv(alpha_i)

    def calculate_points(self):
        self.C=np.array((0,self.rw1))

        A=self.calculate_A()
        self.A=A

        E=self.calculate_E()    
        self.E=E

        self.phi_A=atan2(A[1],A[0])-pi/2-self.phi_fuss*0
        self.phi_E=atan2(E[1],E[0])-pi/2+self.phi_kopf*0
        self.calculate_state(self.phi_A)

    def calculate_E(self):

        #Variante von A nach E
        rb1=self.rb1
        a_w=self.alpha_w
        ra1=self.ra1
        A0=np.array((sin(a_w),cos(a_w)))*rb1
        self.A0=A0
        A0E=sqrt(ra1**2-rb1**2)

        E=A0+np.array((-cos(a_w),sin(a_w)))*A0E
        return E

    def calculate_A(self):

        #Variante von E0 nach A -> Pythagoras
        rb2=self.rb2
        a_w=self.alpha_w
        ra2=self.ra2
        s=np.sign(self.z2)
        E0=-np.array((sin(a_w),cos(a_w)))*rb2
        # if s<1:
        #     print(f"s<1 {s}")
        #     E0=-np.array((sin(a_w),cos(a_w)))*rb2*s
        # else:
        #     E0=-np.array((sin(a_w),cos(a_w)))*rb2*s
        
        self.E0=E0+self.center2
        if abs(ra2)>=abs(rb2):
            E0A=sqrt(ra2**2-rb2**2)
        else:
            print(f"Undercut in AI Gear")
            E0A=0
        A=self.E0+np.array((s*cos(a_w),-sin(a_w)*s))*E0A
       
        return A


    @property
    def center2(self): return np.array((0,self.a))

    @property
    def e_alpha(self):
        '''Überdeckungsfaktor epsilon_alpha'''
        db1,db2=self.rb1*2,self.rb2*2
        da1=self.ra1*2
        da2=self.ra2*2
        alpha=self.alpha
        a_w=self.alpha_w
        denom=2*self.m*pi*cos(alpha)
        zaehler=sqrt(da1**2-db1**2)+sqrt(da2**2-db2**2)-(db1+db2)*tan(a_w)
        
        rb1,rb2=self.rb1,self.rb2
        ra1,ra2=self.ra1,self.ra2
        m=self.m
        T2A=sqrt(ra2**2-rb2**2)
        T1E=sqrt(ra1**2-rb1**2)
        T1T2=tan(a_w)*(rb1+rb2)
        denom=self.m*pi*cos(alpha)
        print(f"T1E={T1E:.2f} T2A={T2A:.2f} T1T2={T1T2:.2f}")
        ea=(T2A+T1E+T1T2)
        A,E=self.A,self.E
        T1,T2=self.A0,self.E0
        T2A_calculated=np.linalg.norm(T2-A)
        T1E_calculated=np.linalg.norm(T1-E)
        print(f"T2A_calc={T2A_calculated} T1E_calc={T1E_calculated}")
        AE=hypot(E[1]-A[1],E[0]-A[0])
        print(f"AE={AE}")
        print(f"ea={ea} ohne denom")
        return AE/denom
        denom=2*self.m*pi*cos(alpha)
        zaehler=sqrt(da1**2-db1**2)+sqrt(da2**2-db2**2)-(db1+db2)*tan(a_w)



        return zaehler/denom

    def calculate_P(self,phi:float=0):
        '''Berechnung des Eingriffspunkt'''
        a_w=self.alpha_w
        r_p=self.rw1*sin(pi/2-a_w)/sin(pi-(pi/2-a_w)+phi)
        P=np.array((-sin(phi),cos(-phi)))*r_p
        return P

    def norm_phi(self,phi:float):
        '''normalisiert den Winkel auf den Bereich des Eingriffs'''
        phi_A=self.phi_A
        phi_E=self.phi_E
        teilungswinkel=-phi_A+phi_E
        #print(f"teilungswinkel={teilungswinkel*180/pi:.1f}°\t phi={phi*180/pi:.1f} \t test={(phi%teilungswinkel)*180/pi:.1f}")
        phi=(phi-phi_A)%teilungswinkel+phi_A
        return phi
    

    def calculate_state(self,phi:float=0,**kwargs):
        w_1=kwargs.pop("w1",self.w1)
        
        mu=kwargs.pop("mu",self.mu)
        phi=self.norm_phi(phi)
        a_w=self.alpha_w
 

        r_p=self.rw1*sin(pi/2-a_w)/sin(pi-(pi/2-a_w)+phi)
        #print(f"rp={r_p:.2f}")
        P=self.calculate_P(phi)
        self.P=P
        
        Fu_length=self.M1/r_p
        Fn=np.array((-cos(a_w),sin(a_w)))*Fu_length/(cos(a_w))
        rho=np.sign(phi)*atan(mu)*np.sign(self.z2)
        F_res=np.linalg.norm(Fn)/cos(rho)*np.array((-cos(a_w+rho),sin(a_w+rho)))
        Fr2=F_res-Fn
        #print(f"Fu={Fu_length:.3f}\t Fn={np.linalg.norm(Fn):.2f} vs F_res={np.linalg.norm(F_res):.2f}")
        #print(f"F_length{F_length:.1f} vs F_res{F_res}")
        #Mr2
        s=np.sign(self.z2)
        P2=P*np.array([1,-s])-self.center2
        #print(f"P2={P2}\t P={P}")
        phi_P2=atan2(P2[1],P2[0])
        angle=s*tau/4+phi_P2
        #print(f"angle={angle*180/pi:.1f}° und phi_p2={phi_P2*180/pi:.0f}°")
        hebel=np.linalg.norm(P2)*sin(angle+(pi/2-a_w-rho))
        #print(f"rp={r_p:.1f} vs hebel={hebel:.1f} rp2={np.linalg.norm(P2)}")
        #print(f"phi={phi*180/pi:.1f} and psi= {1/self.i*phi} \thebel={hebel:.2f} rho={rho*180/pi:.1f}°")
        M_res=np.linalg.norm(F_res)*hebel*np.sign(self.i)
        #print(f"M2={M_res} vs {self.M1*self.i}")
        #sliding velocity
        vt1=np.array((-cos(phi),-sin(phi)))*w_1*r_p
        angle=phi/self.i
        vt2=-np.array((cos(angle),-sin(angle)))*w_1/self.i*np.linalg.norm(P2)
        #print(f"vt1={vt1}\t vt2={vt2} \t delta={vt1-vt2}")
        vg=vt1-vt2
        self.vt1=vt1
        self.vt2=vt2
        self.vg=vg
        self.phi=phi
        result=dict(P=P,F_res=F_res,Fr2=Fr2,Fn=Fn,M_res=M_res,vg=vg,vt1=vt1,vt2=vt2)
        return result
    
    


    def gearGeometry1(self,num:int=1000,**kwargs):
        gamma=kwargs.get("rotation_angle",pi/(2*self.z1)+pi/2)
        x,y=InvoluteGear2D3(z=self.z1,m=self.m,alpha=self.alpha,x=self.x1,rotation_angle=gamma,M=np.array((0,0)),hap=self.hap1,hfp=self.hfp1)
        return x,y
    
    def gearGeometry2(self,num:int=1000,**kwargs):
        gamma=kwargs.get("rotation_angle",pi/(2*self.z2)+pi/2)
        x,y=InvoluteGear2D3(z=self.z2,m=self.m,alpha=self.alpha,x=self.x2,rotation_angle=gamma,M=self.center2,hap=self.hap2,hfp=self.hfp2)
        return x,y
    
    def calculate_path_of_contact(self):
        A,E=self.A,self.E
        return np.array((A[0],E[0])),np.array((A[1],E[1]))
    
    def print_Parameter(self):
        print(f"z1={self.z1}")
        print(f"z2={self.z2}")
        print(f"a={self.a:2f}")
        print(f"r1={self.r1:.1f}")
        print(f"r2={self.r2:.1f}")
        print(f"rb1={self.rb1:.2f}")
        print(f"rb2={self.rb2:.2f}")
        print(f"ra1={self.ra1:.2f}")
        print(f"ra2={self.ra2:.2f}")
        print(f"rf1={self.rf1:.2f}")
        print(f"rf2={self.rf2:.2f}")
        print(f"rw1={self.rw1:.1f}")
        print(f"rw2={self.rw2:.1f}")
        print(f"alpha={self.alpha*180/pi:.1f}")
        print(f"alpha_w={self.alpha_w*180/pi:.1f}")
        print(f"x1={self.x1:.1f}")
        print(f"x2={self.x2:.1f}")
    

if __name__ =="__main__":
    import sys
    import os

    sys.path.append(os.path.dirname(os.path.abspath(__file__)))
    gear_pair = InvoluteGear(
        z1=19,
        z2=-30,
        a=62,
        alpha=20*pi/180,
        hap1=1,hap2=1,
        hfp1=1.25,hfp2=1.25,
        x1=0,x2=0
    )
    print(f"center2={gear_pair.center2}")
    result=gear_pair.calculate_state(-12*pi/180)
    for i, (k, v) in enumerate(result.items()):
        print(i, k, v)
    

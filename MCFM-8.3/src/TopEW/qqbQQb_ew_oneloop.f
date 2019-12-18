      subroutine qqbQQb_ew_oneloop(corr,s,beta,z)
      implicit none
      include 'types.f'
      include 'nf.f'
      include 'constants.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'masses.f'
      include 'epinv.f'
      include 'scale.f'
      include 'qcdcouple.f'
      include 'anomcoup.f'
      real(dp):: corr(nf),s,t,u,beta,z,alpha,sigma0,gvt,gat,
     .     xI1,xI2,xI3,xI4,D6,sw2,cw2,mw,mz,mh,born,f1,f2,
     .     rz,rw,rb,rh,yphi,ys,ini_corr(nf),dFD(nf),T3(nf),gvq(nf),
     .     gaq(nf),Qt,qa(5),bxew(nf,-1:0),bxqcd(nf,-1:0),BB(nf,-1:0),BQCDEW(5),fac
      real(dp) :: db0
      integer ep,i,fl
      real(dp):: vev, Lambdainv, gvt_smeft, gat_smeft, gw_smeft, gvt_smeft2, gat_smeft2, gw_smeft2, Cpq3, Cpu, voL2, voL4, c(1:7),gw2gvt,gw2gat,gatgvt,gw4      
      real(dp)::qa4,qa3,PreFac
      

      
      ep = 0
c      musq = (2._dp*mt)**2
      alpha = esq/fourpi
c      alpha = 1._dp/126.3_dp
!      as = gsq/fourpi

      mw = wmass
      mz = zmass
      mh = hmass
      
      sw2 = xw
! on-shell condition cross-check w/ Doreen
!      sw2 = 1._dp - mw**2/mz**2
      cw2 = 1._dp - sw2

      T3 = (/-0.5_dp,0.5_dp,-0.5_dp,0.5_dp,-0.5_dp/)

      gvq(1:5) = (T3(1:5)-2._dp*sw2*Q(1:5))/2._dp/sqrt(sw2*cw2)
      gaq(1:5) = T3(1:5)/2._dp/sqrt(sw2*cw2)
      gw = 0.5_dp/sqrt(2._dp*sw2)

      Qt=Q(2)
      gvt = gvq(2)
      gat = gaq(2)
      

      
      call smeft_coupl(gvt,gat,gw,gvt_smeft,gat_smeft,gvt_smeft2,gat_smeft2,gw_smeft2,voL2,voL4,Cpq3,Cpu,c)
      

      gw2gvt=(gw**2*((2*Cpq3 - Cpu)*Sqrt(cw2*sw2)*(voL2 + 2*Cpq3*voL4) + 4*cw2*gvt*sw2*(1 + 2*Cpq3*voL2 + Cpq3**2*voL4)))/(4._dp*cw2*sw2)
      gw2gat=(gw**2*((2*Cpq3 + Cpu)*Sqrt(cw2*sw2)*(voL2 + 2*Cpq3*voL4) + 4*cw2*gat*sw2*(1 + 2*Cpq3*voL2 + Cpq3**2*voL4)))/(4._dp*cw2*sw2)
      gatgvt=(16*cw2*gat*gvt*sw2 + 8*Cpq3*(gat + gvt)*Sqrt(cw2*sw2)*voL2 + 4*Cpq3**2*voL4 - Cpu*(4*(gat - gvt)*Sqrt(cw2*sw2)*voL2 + Cpu*voL4))/(16._dp*cw2*sw2)
      gw4=(1 + 4*Cpq3*voL2 + 6*Cpq3**2*voL4)*gw**4


!      as = 1._dp
      sigma0 = pi*as**2/8._dp*(Nc**2-1._dp)/Nc**2*beta/s
      


      t = -s*(1._dp-beta*z)/2._dp+mt**2
      u = -s*(1._dp+beta*z)/2._dp+mt**2

      rz = mz**2/s
      rb = mb**2/s
      rw = mw**2/s
      rh = mh**2/s

      yphi = mb**2/mw**2
      ys = 1._dp - beta**2 + 4._dp*rb
      
!-----BEGIN TILL (gvt->gvt_smeft, gat->gat_smeft)   
C --- set BB  = 0 for now
C      BB = 0._dp
      fac = pi*alpha*as*(Nc**2-1._dp)*8._dp*s/(s - mz**2)
      BB(:,0) = - fac*(
     .     + (2._dp - beta**2*(1._dp - z**2))*gvq(:)*gvt_smeft*c(1)
     .     + 2._dp*beta*z*gaq(:)*gat_smeft*c(2) )

!      write(6,*) 'BB: ', BB(1:2,0)
!      write(6,*) 'born: ', 4._dp*s/(s - mz**2)*(
!     .     + (2._dp - beta**2*(1._dp - z**2))*gvq(1:2)*gvt
!     .     - 2._dp*beta*z*gaq(1:2)*gat )
!      stop
!      BB(:,-1) = - 2._dp*fac*(-gvq(:)*gvt + 3._dp*beta*z*gaq(:)*gat)

      dFD = - alpha/8._dp/pi*((gvq**2+gaq**2)*f1(rz)+2._dp*gw**2*f1(rw))

      born = sigma0*(2._dp - beta**2 + beta**2*z**2)

      ini_corr = 2._dp*born*dFD
      
      PreFac=sigma0/8._dp

      if (vol2>0._dp) then  
      do fl=1,5
         BQCDEW(fl)=1/32._dp/Pi*beta/s*((32*alpha**2*Pi**2*(MZ**4*Qt**2*(s**2 + 2*s*t + 2*t**2)*Q(fl)**2 - 2*MZ**2*Qt*s*Q(fl)*(gat_smeft*s*(s + 2*t)*gaq(fl) + (s**2 + 2*s*t + 2*t**2)*(gvt_smeft*gvq(fl) + Qt*Q(fl))) + 
     -      s**2*((gat_smeft2 + gvt_smeft2)*(s**2 + 2*s*t + 2*t**2)*gaq(fl)**2 + 2*s*(s + 2*t)*gaq(fl)*(2*gatgvt*gvq(fl) + gat_smeft*Qt*Q(fl)) + 
     -         (s**2 + 2*s*t + 2*t**2)*((gat_smeft2 + gvt_smeft2)*gvq(fl)**2 + 2*gvt_smeft*Qt*gvq(fl)*Q(fl) + Qt**2*Q(fl)**2)) + 
     -      2*mt**4*(MZ**4*Qt**2*Q(fl)**2 - 2*MZ**2*Qt*s*Q(fl)*(gvt_smeft*gvq(fl) + Qt*Q(fl)) + s**2*((gat_smeft2 + gvt_smeft2)*(gaq(fl)**2 + gvq(fl)**2) + 2*gvt_smeft*Qt*gvq(fl)*Q(fl) + Qt**2*Q(fl)**2)) + 
     -      4*mt**2*(s**2*(gat_smeft*MZ**2*Qt*gaq(fl)*Q(fl) - s*(gat_smeft2*gaq(fl)**2 + 2*gatgvt*gaq(fl)*gvq(fl) + gat_smeft2*gvq(fl)**2 + gat_smeft*Qt*gaq(fl)*Q(fl))) - 
     -         t*(MZ**4*Qt**2*Q(fl)**2 - 2*MZ**2*Qt*s*Q(fl)*(gvt_smeft*gvq(fl) + Qt*Q(fl)) + s**2*((gat_smeft2 + gvt_smeft2)*(gaq(fl)**2 + gvq(fl)**2) + 2*gvt_smeft*Qt*gvq(fl)*Q(fl) + Qt**2*Q(fl)**2)))))/
     -  ((MZ**2 - s)**2*s**2))
      enddo
      BQCDEW(5)=BQCDEW(5)+1/32._dp/Pi*beta/s*(-32*alpha*as*Pi**2*(mt**6 + mt**4*(2*MW**2 + s - 2*t) + 2*MW**2*(s + t)**2 + mt**2*(t**2 - 2*MW**2*(s + 2*t)))*(1 + 2*Cpq3*vol2 + Cpq3**2*vol4))/(9._dp*MW**2*s*sw2*(MW**2 - t))

      BQCDEW(5)=BQCDEW(5)+1/32._dp/Pi*beta/s*((64*alpha**2*GW4*Pi**2*(mt**8 - 2*mt**6*t - 8*mt**2*MW**4*(s + t) + 4*MW**4*(s + t)**2 + mt**4*(4*MW**2*(MW**2 + s) + t**2)))/(MW**4*(MW**2 - t)**2))
      BQCDEW(5)=BQCDEW(5)+1/32._dp/Pi*beta/s*((-64*alpha**2*GW_smeft2*Pi**2*Qt*(mt**6 + mt**4*(2*MW**2 + s - 2*t) + 2*MW**2*(s + t)**2 + mt**2*(t**2 - 2*MW**2*(s + 2*t)))*Q(5))/(3._dp*MW**2*s*(MW**2 - t)))
      BQCDEW(5)=BQCDEW(5)+1/32._dp/Pi*beta/s*((64*alpha**2*Pi**2*(GW2gvt*(mt**6 + mt**4*(2*MW**2 + s - 2*t) + 2*MW**2*(s + t)**2 + mt**2*(t**2 - 2*MW**2*(s + 2*t))) + 
     -      GW2gat*(-mt**6 + 2*MW**2*(s + t)**2 + mt**4*(2*MW**2 + s + 2*t) - mt**2*(t**2 + MW**2*(6*s + 4*t))))*(gaq(5) + gvq(5)))/(3._dp*MW**2*(MZ**2 - s)*(MW**2 - t)))
      else 
      do fl=1,5
         BQCDEW(fl)=0._dp
      enddo
      end if

      if (1==1) then
      qa(1) = 
     .     (0.25_dp*alpha*sigma0*(-2._dp*(gat_smeft2*c(5) + gvt_smeft2*c(4))*(1._dp + z**2) 
     .     + 2._dp*((1._dp - beta**2)*(-3._dp*gat_smeft2*c(5) + gvt_smeft2*c(4)) 
     .     + 2._dp*(gat_smeft2*c(5) + gvt_smeft2*c(4))*rz)*s*(2._dp - beta**2*(1._dp 
     .     - z**2))*db0(mt**2,mz**2,mt**2) + (8._dp*(gat_smeft2*c(5) 
     .     + gvt_smeft2*c(4))*(1._dp + z**2)*(-xI1(mt**2,musq,ep) 
     .     + xI1(mz**2,musq,ep)))/((1._dp - beta**2)*s) 
     .     - 4._dp*(gat_smeft2*c(5)*(1._dp - 5._dp*z**2 - 3._dp*beta**2*(1._dp 
     .     - z**2)) - gvt_smeft2*c(4)*(3._dp + z**2 - beta**2*(1._dp - z**2)) 
     .     + ((gat_smeft2*c(5) + gvt_smeft2*c(4))*rz*(1._dp - beta**2 - 3._dp*z**2 
     .     + 7._dp*beta**2*z**2 + 2._dp*beta**4*(1._dp - z**2)))/(beta**2*
     .     (1._dp - beta**2)))*xI2(mt**2,mz**2,mt**2,musq,ep) 
     .     + 2._dp*(-4._dp*gvt_smeft2*c(4)*(1._dp + z**2) - (-3._dp*gat_smeft2*c(5) 
     .     + gvt_smeft2*c(4))*f2(z,beta) + (2._dp*(gat_smeft2*c(5) + gvt_smeft2*c(4))*rz*f2(z,beta))
     .     /beta**2)*xI2(s,mt**2,mt**2,musq,ep) + 2._dp*s*(-(1._dp 
     .     + beta**2)*gvt_smeft2*c(4)*(2._dp - beta**2*(1._dp - z**2)) 
     .     - gat_smeft2*c(5)*(2._dp - 3._dp*beta**2 + 5._dp*beta**2*z**2 
     .     + 3._dp*beta**4*(1._dp - z**2)) + (2._dp*(gat_smeft2*c(5) 
     .     + gvt_smeft2*c(4))*rz**2*f2(z,beta))/beta**2 + 4._dp*rz*(-gvt_smeft2*c(4)*
     .     (1._dp + z**2) + gat_smeft2*c(5)*f2(z,beta)))*xI3(mt**2,mt**2,s,mt**2,
     .     mz**2,mt**2,musq,ep)))/pi
     

       
!-----END TILL     

!-----BEGIN TILL (gw->gw_smeft)
      qa(2) = 
     .     (0.5_dp*alpha*gw_smeft2*c(6)*sigma0*(-2._dp*(1._dp + z**2) 
     .     + (-1._dp + beta**2 + 4._dp*(-rb + rw))*s*(2._dp - beta**2*(1._dp 
     .     - z**2))*db0(mt**2,mb**2,mw**2) + (8._dp*(1._dp 
     .     + z**2)*(-xI1(mb**2,musq,ep) + xI1(mw**2,musq,ep)))/((1._dp 
     .     - beta**2)*s) - (1._dp*(1._dp - 3._dp*z**2 - 2._dp*beta**4*(1._dp 
     .     - z**2) - 5._dp*beta**2*(1._dp + z**2) + (4._dp*(-rb + rw)*(1._dp 
     .     - beta**2 - 3._dp*z**2 + 7._dp*beta**2*z**2 + 2._dp*beta**4*
     .     (1._dp - z**2)))/(1._dp - beta**2))*xI2(mt**2,mb**2,mw**2,musq,
     .     ep))/beta**2 + ((1._dp - 3._dp*z**2 - 2._dp*beta**4*(1._dp 
     .     - z**2) - 5._dp*beta**2*(1._dp + z**2) + 4._dp*(-rb + rw)*
     .     f2(z,beta))*xI2(s,mb**2,mb**2,musq,ep))/beta**2 + (0.25_dp*s*
     .     (1._dp - 16._dp*beta**2 + beta**4 - 3._dp*z**2 - 4._dp*beta**2*
     .     z**2 - 11._dp*beta**4*z**2 - 2._dp*beta**6*(1._dp - z**2) 
     .     + 8._dp*rw*(1._dp - 3._dp*z**2 - 2._dp*beta**4*(1._dp - z**2) 
     .     - 3._dp*beta**2*(1._dp + z**2)) + 8._dp*(-(1._dp + beta**2)*rb 
     .     + 2._dp*(-rb + rw)**2)*f2(z,beta))*xI3(mt**2,mt**2,s,mb**2,
     .     mw**2,mb**2,musq,ep))/beta**2))/pi


!-----END TILL     

      qa(3) = (0.5_dp*alpha*gat_smeft2*mt**2*sigma0*(-2._dp*(1._dp + z**2) 
     .     + 4._dp*rz*s*(2._dp - beta**2*(1._dp - z**2))*db0(mt**2,mz**2,
     .     mt**2) - (8._dp*(1._dp + z**2)*(xI1(mt**2,musq,ep) 
     .     - xI1(mz**2,musq,ep)))/((1._dp - beta**2)*s) - (4._dp*rz*(1._dp 
     .     - beta**2 - 3._dp*z**2 + 7._dp*beta**2*z**2 + 2._dp*beta**4*
     .     (1._dp - z**2))*xI2(mt**2,mz**2,mt**2,musq,ep))/(beta**2*(1._dp 
     .     - beta**2)) + (2._dp*(beta**2*(1._dp + z**2) + 2._dp*rz*
     .     f2(z,beta))*xI2(s,mt**2,mt**2,musq,ep))/beta**2 + (4._dp*s*
     .     (beta**2*(1._dp - beta**2)*rz*(1._dp - z**2) 
     .     + rz**2*f2(z,beta))*xI3(mt**2,mt**2,s,mt**2,mz**2,mt**2,musq,
     .     ep))/beta**2))/(mz**2*pi)
      

      qa(4) = (0.5_dp*alpha*gw_smeft2*sigma0*((-0.25_dp*ys*(1._dp + z**2))/rw 
     .     + 0.125_dp*s*((-(1._dp - beta**2)**2)/rw + 8._dp*(1._dp 
     .     - beta**2)*yphi - 16._dp*rb*yphi + 4._dp*ys)*(2._dp - beta**2*
     .     (1._dp - z**2))*db0(mt**2,mb**2,mw**2) + (ys*(1._dp 
     .     + z**2)*(-xI1(mb**2,musq,ep) + xI1(mw**2,musq,ep)))/((1._dp 
     .     - beta**2)*mw**2) - (0.125_dp*(-16._dp*beta**2*
     .     (1._dp - beta**2)**2*yphi*(1._dp - z**2) - 16._dp*rb*yphi*(1._dp 
     .     - beta**2 - 3._dp*z**2 + 7._dp*beta**2*z**2 + 2._dp*beta**4*
     .     (1._dp - z**2)) + 4._dp*ys*(1._dp - beta**2 - 3._dp*z**2 
     .     + 7._dp*beta**2*z**2 + 2._dp*beta**4*(1._dp - z**2)) 
     .     + ((1._dp - beta**2)**2*(1._dp - 3._dp*z**2 - 2._dp*beta**4*(1._dp 
     .     - z**2) + 3._dp*beta**2*(1._dp + z**2)))/rw)*xI2(mt**2,mb**2,
     .     mw**2,musq,ep))/(beta**2*(1._dp - beta**2)) + (0.125_dp*(((1._dp 
     .     - beta**2)*(1._dp - 3._dp*z**2 - 2._dp*beta**4*(1._dp - z**2) 
     .     + 3._dp*beta**2*(1._dp + z**2)))/rw + 4._dp*(-2._dp*beta**2*yphi 
     .     - 4._dp*rb*yphi + ys)*f2(z,beta))*xI2(s,mb**2,mb**2,musq,ep))/
     .     beta**2 + (0.03125_dp*s*(8._dp*(1._dp - beta**2)**2*(1._dp 
     .     - 3._dp*z**2 + 2._dp*beta**2*(1._dp - z**2)) - 4._dp*(1._dp 
     .     - beta**2)*yphi*(1._dp - 11._dp*beta**2 - 3._dp*z**2 
     .     - 3._dp*beta**2*z**2 + 2._dp*beta**4*(1._dp - z**2)) 
     .     + ((1._dp - beta**2)**2*(1._dp - 7._dp*beta**2 - 3._dp*z**2 
     .     + beta**2*z**2 + 2._dp*beta**4*(1._dp - z**2)))/rw 
     .     - 16._dp*rb*yphi*(1._dp + beta**2 - 3._dp*z**2 
     .     + 9._dp*beta**2*z**2 + 2._dp*beta**4*(1._dp - z**2)) 
     .     + 16._dp*(-8._dp*rb**2 + 4._dp*rb**2*yphi + rw*ys)*f2(z,beta))*
     .     xI3(mt**2,mt**2,s,mb**2,mw**2,mb**2,musq,ep))/beta**2))/pi
     
      qa(5) = 
     .     (0.5_dp*alpha*gw**2*mt**2*sigma0*(-1._dp - z**2 + 2._dp*(-1._dp 
     .     + beta**2 + rh)*s*(2._dp - beta**2*(1._dp - z**2))*db0(mt**2,
     .     mt**2,mh**2) + (4._dp*(1._dp + z**2)*
     .     (xI1(mh**2,musq,ep) - xI1(mt**2,musq,ep)))/((1._dp - beta**2)*
     .     s) - 2._dp*(2._dp*(1._dp - beta**2)*(1._dp - z**2) + (rh*(1._dp 
     .     - beta**2 - 3._dp*z**2 + 7._dp*beta**2*z**2 + 2._dp*beta**4*
     .     (1._dp - z**2)))/(beta**2*(1._dp - beta**2)))*xI2(mt**2,mt**2,
     .     mh**2,musq,ep) + ((beta**2*(5._dp - 3._dp*z**2 - 4._dp*beta**2*
     .     (1._dp - z**2)) + 2._dp*rh*f2(z,beta))*xI2(s,mt**2,mt**2,musq,
     .     ep))/beta**2 + (2._dp*s*(3._dp*beta**2*(1._dp - beta**2)*rh*
     .     (1._dp - z**2) - beta**2*(1._dp - beta**2)*(2._dp - beta**2*
     .     (1._dp - z**2)) + rh**2*f2(z,beta))*xI3(mt**2,mt**2,s,mt**2,
     .     mh**2,mt**2,musq,ep))/beta**2))/(mw**2*pi)
       

!-----BEGIN TILL (gvt->gvt_smeft, gat->gat_smeft)
      do ep = -1,0
         bxew(:,ep) = 
     .     (-0.015625_dp*as*BB(:,0)*beta*(-(1._dp + beta*z)*
     .     xI3(0._dp,u,mt**2,0._dp,0._dp,mt**2,musq,ep) 
     .     + (1._dp - beta*z)*
     .     xI3(0._dp,t,mt**2,0._dp,0._dp,mt**2,musq,ep)))
     .     /(Nc**2*pi) + (alpha*sigma0*((-2._dp*beta**2*s*(1._dp - z**2)*
     .     (2._dp*gaq*gat_smeft*c(2)*(-rz**2 + beta*z + rz*(1._dp + beta*z)) 
     .     + gvq*gvt_smeft*c(1)*(3._dp - beta**2 + 2._dp*rz**2 + 2._dp*beta*z*(1._dp 
     .     + beta*z) - rz*(1._dp + beta**2 + 2._dp*beta*z)))*D6(0._dp,0._dp
     .     ,mt**2,mt**2,s,u,0._dp,0._dp,mz**2,mt**2,musq,ep))/
     .     ((1._dp - rz)**2*(1._dp + beta*z)) + (2._dp*beta**2*s*(1._dp 
     .     - z**2)*(2._dp*gaq*gat_smeft*c(2)*(rz**2 + beta*z - rz*(1._dp - beta*z)) 
     .     + gvq*gvt_smeft*c(1)*(3._dp - beta**2 + 2._dp*rz**2 - rz*(1._dp + beta**2 
     .     - 2._dp*beta*z) - 2._dp*beta*z*(1._dp - beta*z)))*D6(0._dp,0._dp,
     .     mt**2,mt**2,s,t,0._dp,0._dp,mz**2,mt**2,musq,ep))/
     .     ((1._dp - rz)**2*(1._dp - beta*z)) - (2._dp*(1._dp - beta**2)*
     .     (2._dp*beta*gaq*gat_smeft*c(2) + gvq*gvt_smeft*c(1)*z*(1._dp + 2._dp*beta**2 
     .     - beta**2*z**2))*xI2(mt**2,0._dp,mt**2,musq,ep))/(beta*(1._dp 
     .     - beta**2*z**2)) - (2._dp*(1._dp - beta**2)*(2._dp*beta*gaq*gat_smeft*c(2) 
     .     + gvq*gvt_smeft*c(1)*z*(1._dp + 2._dp*beta**2 - beta**2*z**2))*xI2(mt**2,
     .     mz**2,mt**2,musq,ep))/(beta*(1._dp - beta**2*z**2)) 
     .     + (4._dp*(beta*gaq*gat_smeft*c(2) + gvq*gvt_smeft*c(1)*z)*xI2(s,0._dp,mz**2,musq,ep))
     .     /beta - (2._dp*(-gaq*gat_smeft*c(2) + gvq*gvt_smeft)*c(1)*(1._dp - 2._dp*beta**2 
     .     + beta**2*z**2)*xI2(u,0._dp,mt**2,musq,ep))/(1._dp + beta*z) 
     .     + (2._dp*(gaq*gat_smeft*c(2) + gvq*gvt_smeft*c(1))*(1._dp - 2._dp*beta**2 
     .     + beta**2*z**2)*xI2(t,0._dp,mt**2,musq,ep))/(1._dp - beta*z) 
     .     - (1._dp*s*(2._dp*gaq*gat_smeft*c(2)*(1._dp - beta**2 + 2._dp*beta*z 
     .     - beta**3*z + 2._dp*beta**2*z**2 + beta**3*z**3 + rz*(-3._dp 
     .     + 3._dp*beta**2 - 3._dp*beta*z + beta**3*z - 2._dp*beta**2*z**2) 
     .     - rz**3*(1._dp - 2._dp*beta**2 + beta**2*z**2) + rz**2*(3._dp 
     .     - 4._dp*beta**2 + beta*z - 2._dp*beta**3*z + beta**2*z**2 
     .     + beta**3*z**3)) + gvq*gvt_smeft*c(1)*(2._dp - beta**2 + 4._dp*beta*z 
     .     - 2._dp*beta**3*z + 3._dp*beta**2*z**2 - beta**4*z**2 
     .     + 2._dp*beta**3*z**3 + beta**4*z**4 + 2._dp*rz**3*(1._dp 
     .     - 2._dp*beta**2 + beta**2*z**2) - rz**2*(4._dp - 5._dp*beta**2 
     .     - beta**4 - 2._dp*beta**3*z + beta**2*z**2 + beta**4*z**2 
     .     + 2._dp*beta**3*z**3) + beta*rz*(-4._dp*beta + beta**3 - 4._dp*z 
     .     - 2._dp*beta**3*z**2 + beta**3*z**4)))*xI3(0._dp,mt**2,u,0._dp,
     .     mz**2,mt**2,musq,ep))/((1._dp - rz)**2*(1._dp + beta*z)) 
     .     + (s*(2._dp*gaq*gat_smeft*c(2)*(-1._dp + beta**2 + 2._dp*beta*z - beta**3*z 
     .     - 2._dp*beta**2*z**2 + beta**3*z**3 + rz**3*(1._dp 
     .     - 2._dp*beta**2 + beta**2*z**2) + rz*(3._dp - 3._dp*beta**2 
     .     - 3._dp*beta*z + beta**3*z + 2._dp*beta**2*z**2) + rz**2*(-3._dp 
     .     + 4._dp*beta**2 + beta*z - 2._dp*beta**3*z - beta**2*z**2 
     .     + beta**3*z**3)) + gvq*gvt_smeft*c(1)*(2._dp - beta**2 - 4._dp*beta*z 
     .     + 2._dp*beta**3*z + 3._dp*beta**2*z**2 - beta**4*z**2 
     .     - 2._dp*beta**3*z**3 + beta**4*z**4 + 2._dp*rz**3*(1._dp 
     .     - 2._dp*beta**2 + beta**2*z**2) - rz**2*(4._dp - 5._dp*beta**2 
     .     - beta**4 + 2._dp*beta**3*z + beta**2*z**2 + beta**4*z**2 
     .     - 2._dp*beta**3*z**3) + beta*rz*(-4._dp*beta + beta**3 + 4._dp*z 
     .     - 2._dp*beta**3*z**2 + beta**3*z**4)))*xI3(0._dp,mt**2,t,0._dp,
     .     mz**2,mt**2,musq,ep))/((1._dp - rz)**2*(1._dp - beta*z)) 
     .     + (2._dp*(1._dp - beta**2)*s*(2._dp*beta*gaq*gat_smeft*c(2)*(rz**2 
     .     + beta**2*z**2 - rz*(1._dp - beta**2*z**2)) + gvq*gvt_smeft*c(1)*z*(1._dp 
     .     + 2._dp*beta**2 - beta**4 - beta**2*z**2 + beta**4*z**2 
     .     - rz**2*(-1._dp - 2._dp*beta**2 + beta**2*z**2) + rz*(-2._dp 
     .     - beta**4 + 2._dp*beta**2*z**2 + beta**4*z**2)))*xI3(s,mt**2,
     .     mt**2,0._dp,mz**2,mt**2,musq,ep))/(beta*(1._dp - rz)*(1._dp 
     .     - beta**2*z**2))))/pi

         bxqcd(:,ep) = 
     .     (-0.015625_dp*as*BB(:,0)*beta*(-(1._dp + beta*z)*
     .     xI3(0._dp,u,mt**2,0._dp,0._dp,mt**2,musq,ep)
     .     + (1._dp - beta*z)*
     .     xI3(0._dp,t,mt**2,0._dp,0._dp,mt**2,musq,ep)))
     .     /(Nc**2*pi) + (alpha*sigma0*((-beta**2*s*(1._dp - z**2)*
     .     (2._dp*beta*gaq*gat_smeft*c(2)*z + gvq*gvt_smeft*c(1)*(3._dp - beta**2 
     .     + 2._dp*beta*z*(1._dp + beta*z)))*D6(0._dp,0._dp,mt**2,mt**2,s,
     .     u,0._dp,0._dp,0._dp,mt**2,musq,ep))/(1._dp + beta*z) 
     .     + (beta**2*s*(1._dp - z**2)*(2._dp*beta*gaq*gat_smeft*c(2)*z 
     .     + gvq*gvt_smeft*c(1)*(3._dp - beta**2 - 2._dp*beta*z*(1._dp 
     .     - beta*z)))*D6(0._dp,0._dp,mt**2,mt**2,s,t,0._dp,0._dp,0._dp,
     .     mt**2,musq,ep))/(1._dp - beta*z) + (2._dp*(1._dp - beta**2)*
     .     (-beta*gaq*gat_smeft*c(2)*(1._dp + beta**2*z**2) + gvq*gvt_smeft*c(1)*z*(-1._dp 
     .     - 2._dp*beta**2 + beta**2*z**2))*
     .     xI2(mt**2,0._dp,mt**2,musq,ep))/(beta*(1._dp - beta**2*z**2)) 
     .     + (2._dp*(beta*gaq*gat_smeft*c(2) + gvq*gvt_smeft*c(1)*z)*xI2(s,0._dp,0._dp,musq,ep))
     .     /beta - (1._dp*(beta*gaq*gat_smeft*c(2)*(beta + z)*(1._dp - beta*z) 
     .     + gvq*gvt_smeft*c(1)*(1._dp - 2._dp*beta**2 + beta**2*z**2))*xI2(u,0._dp,
     .     mt**2,musq,ep))/(1._dp + beta*z) + ((-beta*gaq*gat_smeft*c(2)*(beta - z)*
     .     (1._dp + beta*z) + gvq*gvt_smeft*c(1)*(1._dp - 2._dp*beta**2 
     .     + beta**2*z**2))*xI2(t,0._dp,mt**2,musq,ep))/(1._dp - beta*z) 
     .     + ((1._dp - beta**2)*s*(beta*gaq*gat_smeft*c(2)*(1._dp + beta**2*z**2) 
     .     + gvq*gvt_smeft*c(1)*z*(1._dp + beta**2 + beta**2*(1._dp - beta**2)*(1._dp 
     .     - z**2)))*xI3(s,mt**2,mt**2,0._dp,0._dp,mt**2,musq,ep))/
     .     (beta*(1._dp - beta**2*z**2))))/(pi*(1._dp - rz))

           
      end do 

!-----END TILL

! check the coefficient prop to epinv
!      write(6,*) bxqcd(1:2,-1)*BB(1:2,0)/BB(1:2,-1)
!      write(6,*) as/32._dp/pi/Nc**2*beta/s*BB(1:2,0)*log((t-mt**2)
!     .     /(u-mt**2))
! check the finite piece from epinv/epinv
!      write(6,*) bxqcd(1:2,-1)
!      write(6,*) as/32._dp/pi/Nc**2*beta/s*BB(1:2,-1)*log((t-mt**2)
!     .     /(u-mt**2))
!      bxew(:,-1) = bxew(:,-1)*(BB(:,0)/BB(:,-1)*epinv + 1._dp)
!      bxqcd(:,-1) = bxqcd(:,-1)*(BB(:,0)/BB(:,-1)*epinv + 1._dp)

      bxew(:,-1) = bxew(:,-1)*epinv
      bxqcd(:,-1) = bxqcd(:,-1)*epinv

c--- use the value of "tevscale" (passed via anomcoup.f)
c--- as an anomalous top Yukawa coupling:
c--- g(top Y) = (tevscale) x g(SM, top Y)
c--- it therefore affects Higgs diagrams as the square
      qa(5)=qa(5)*tevscale**2
 

      corr = BQCDEW(:) + qa(1) + qa(2) + qa(3) + qa(4) + qa(5) 
     .     + bxew(:,0) + bxew(:,-1) + bxqcd(:,0) + bxqcd(:,-1)
      corr = corr + ini_corr
      endif

c      corr = BQCDEW(:)

      corr = corr/born

      end subroutine qqbQQb_ew_oneloop


      function f1(x)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'cplx.h'
      real(dp):: x,cli2,f1

      f1 = real(1._dp + 2._dp*(
     .     + (1._dp+log(cplx1(x)))*(2._dp*x+3._dp)
     .     - 2._dp*(1._dp+x)**2*(cli2(cplx1(1._dp+1._dp/x))
     .     -pi**2/6._dp)),dp)

      end function f1


      function f2(z,beta)
      implicit none
      include 'types.f'
      real(dp):: z,beta,f2

      f2 = 1._dp - 3._dp*z**2 - 2._dp*beta**2*(1._dp - z**2)

      end function f2


      function D6(p1sq,p2sq,p3sq,p4sq,s12,s23,m1sq,
     .     m2sq,m3sq,m4sq,musq,ep)
      implicit none
      include 'types.f'
      include 'constants.f'
      real(dp):: p1sq,p2sq,p3sq,p4sq,s12,s23,m1sq,m2sq,m3sq,
     .     m4sq,musq,xI3,xI4,p1dp2,p1dp3,p2dp3,m1,m2,m3,m4,D6
      integer ep
      external xI3,xI4

      p1dp2 = (s12 - p1sq - p2sq)/2._dp
      p2dp3 = (s23 - p2sq - p3sq)/2._dp
      p1dp3 = (p2sq + p4sq - s12 -s23)/2._dp
      m1 = sqrt(m1sq)
      m2 = sqrt(m2sq)
      m3 = sqrt(m3sq)
      m4 = sqrt(m4sq)

      D6 = 
     .     (0.25_dp*(-((-m1**2 + m2**2)*p1dp3*p2sq + p1dp2*((m1**2 
     .     - m2**2 - p1sq)*p2dp3 + p1dp3*(m2**2 - m3**2 + p2sq)) 
     .     - p1dp2**2*(m3**2 - m4**2 + 2._dp*p2dp3 + p3sq) 
     .     + p1sq*((-m2**2 + m3**2)*p2dp3 + p2sq*(m3**2 - m4**2 + p1dp3 
     .     + p2dp3 + p3sq)))*xI3(p1sq,p2sq,s12,m1**2,m2**2,m3**2,musq,
     .     ep) + (p1dp3**2*(m2**2 - m3**2 + p2sq) - p1dp3*(m1**2 - m2**2 
     .     - p1sq)*(p2dp3 + p2sq) - p1dp2**2*(m3**2 - m4**2 + 2._dp*p2dp3 
     .     + p3sq) + p1sq*(2._dp*p2dp3**2 + (m3**2 - m4**2 + p2dp3)*p2sq 
     .     - p2dp3*(m2**2 - 2._dp*m3**2 + m4**2 - p3sq) 
     .     + (-m2**2 + m3**2)*p3sq) + p1dp2*(p1dp3*(m2**2 - 2._dp*m3**2 
     .     + m4**2 - 2._dp*p2dp3 + p2sq - p3sq) + (m1**2 - m2**2 - p1sq)*
     .     (p2dp3 + p3sq)))*xI3(p1sq,s23,p4sq,m1**2,m2**2,m4**2,musq,ep) 
     .     - ((m1**2 - m2**2 - p1sq)*(p2dp3**2 - p2sq*p3sq) 
     .     + p1dp2*(-2._dp*p2dp3**2 + (m2**2 - m3**2 + p2sq)*p3sq 
     .     - p2dp3*(m3**2 - m4**2 + p3sq)) + p1dp3*((-m2**2 + m3**2)*
     .     p2dp3 + p2sq*(m3**2 - m4**2 + p2dp3 + p3sq)))*xI3(p2sq,p3sq,
     .     s23,m2**2,m3**2,m4**2,musq,ep) + (-m3**2*p1sq*p2dp3 
     .     + m4**2*p1sq*p2dp3 + m1**2*p2dp3**2 - m2**2*p2dp3**2 
     .     - p1sq*p2dp3**2 + p1dp3**2*(-m2**2 + m3**2 + p2sq) 
     .     + 2._dp*p1dp2**2*p3sq + m2**2*p1sq*p3sq - m3**2*p1sq*p3sq 
     .     - p1sq*p2dp3*p3sq - m1**2*p2sq*p3sq + m2**2*p2sq*p3sq 
     .     + p1dp2*(-2._dp*p2dp3**2 + (-m1**2 + 2._dp*m2**2 - m3**2 + p1sq 
     .     + p2sq)*p3sq - p2dp3*(m3**2 - m4**2 + p3sq) + p1dp3*(m3**2 
     .     - m4**2 - 2._dp*p2dp3 + p3sq)) + p1dp3*((m1**2 - 2._dp*m2**2 
     .     + m3**2 - p1sq)*p2dp3 + p2sq*(m3**2 - m4**2 + p2dp3 + p3sq))
     .     )*xI3(s12,p3sq,p4sq,m1**2,m3**2,m4**2,musq,ep) 
     .     + (-2._dp*m2**2*m3**2*p1sq*p2dp3 + 2._dp*m3**4*p1sq*p2dp3 
     .     + 2._dp*m2**2*m4**2*p1sq*p2dp3 - 2._dp*m3**2*m4**2*p1sq*p2dp3 
     .     - m1**4*p2dp3**2 + 2._dp*m1**2*m2**2*p2dp3**2 - m2**4*p2dp3**2 
     .     + 2._dp*m1**2*p1sq*p2dp3**2 - 2._dp*m2**2*p1sq*p2dp3**2 
     .     + 4._dp*m3**2*p1sq*p2dp3**2 - p1sq**2*p2dp3**2 
     .     + m3**4*p1sq*p2sq - 2._dp*m3**2*m4**2*p1sq*p2sq 
     .     + m4**4*p1sq*p2sq + 2._dp*m3**2*p1sq*p2dp3*p2sq 
     .     - 2._dp*m4**2*p1sq*p2dp3*p2sq - p1dp3**2*((m2**2 - m3**2)**2 
     .     - 2._dp*(m2**2 + m3**2)*p2sq + p2sq**2) + m2**4*p1sq*p3sq 
     .     - 2._dp*m2**2*m3**2*p1sq*p3sq + m3**4*p1sq*p3sq 
     .     - 2._dp*m2**2*p1sq*p2dp3*p3sq + 2._dp*m3**2*p1sq*p2dp3*p3sq 
     .     + m1**4*p2sq*p3sq - 2._dp*m1**2*m2**2*p2sq*p3sq 
     .     + m2**4*p2sq*p3sq - 2._dp*m1**2*p1sq*p2sq*p3sq 
     .     - 2._dp*m4**2*p1sq*p2sq*p3sq + p1sq**2*p2sq*p3sq 
     .     + 2._dp*p1sq*p2dp3*p2sq*p3sq + p1sq*p2sq**2*p3sq 
     .     + p1sq*p2sq*p3sq**2 - p1dp2**2*((m3**2 - m4**2)**2 
     .     + 4._dp*p2dp3**2 - 2._dp*(2._dp*m2**2 - m3**2 + m4**2)*p3sq 
     .     + p3sq**2 + 4._dp*p2dp3*(m3**2 - m4**2 + p3sq)) 
     .     - 2._dp*p1dp3*(m1**2 - m2**2 - p1sq)*((-m2**2 + m3**2)*p2dp3 
     .     + p2sq*(m3**2 - m4**2 + p2dp3 + p3sq)) + 2._dp*p1dp2*((m1**2 
     .     - m2**2 - p1sq)*(2._dp*p2dp3**2 - (m2**2 - m3**2 + p2sq)*p3sq 
     .     + p2dp3*(m3**2 - m4**2 + p3sq)) + p1dp3*(-2._dp*(m2**2 
     .     + m3**2)*p2dp3 + (m2**2 - m3**2)*(m3**2 - m4**2 + p3sq) 
     .     + p2sq*(m3**2 - m4**2 + 2._dp*p2dp3 + p3sq))))*xI4(p1sq,p2sq,
     .     p3sq,p4sq,s12,s23,m1**2,m2**2,m3**2,m4**2,musq,ep)))/
     .     (-2._dp*p1dp2*p1dp3*p2dp3 + p1dp3**2*p2sq + p1dp2**2*p3sq 
     .     + p1sq*(p2dp3**2 - p2sq*p3sq))

      D6 = -2._dp*pi*D6

      end function D6

      subroutine smeft_coupl(gvt_sm,gat_sm,gw_sm,gvt_smeft,gat_smeft,gvt_smeft2,gat_smeft2,gw_smeft2,voL2,voL4,Cpq3,Cpu,c)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'ewcouple.f'
      include 'masses.f'
      real(dp):: Cpq3i, Cpui, voL2i, voL4i
      common /smeftcpl/ Cpq3i, Cpui, vol2i, vol4i  
      real(dp):: mw,mz,gvt_sm,gat_sm,gw_sm,sw2,cw2
      integer i
      real(dp):: vev, Lambdainv, gvt_smeft, gat_smeft, gw_smeft, gvt_smeft2, gat_smeft2, gw_smeft2,voL2,voL4,Cpq3,Cpq1,Cpu,c(1:7)

      sw2 = xw
      cw2 = 1._dp-sw2      
      mw = wmass
      mz = zmass
            
      do i=1,7
       c(i)=1._dp
      enddo
      
      
      vev=1._dp/sqrt(2._dp*esq)*mw/gw
      Lambdainv=0.001_dp
      voL2=vol2i*(vev*Lambdainv)**2
      voL4=vol4i*voL2**2

      if (vol2i>0._dp) then
       Cpq3=Cpq3i/vol2
       Cpq1=-Cpq3
       Cpu=Cpui/vol2
      else
       Cpq3=0._dp
       Cpq1=0._dp
       Cpu=0._dp
      end if

      gvt_smeft=gvt_sm-voL2*0.5_dp*(Cpq1-Cpq3+Cpu)*0.5_dp/sqrt(sw2*cw2)
      gat_smeft=gat_sm-voL2*0.5_dp*(Cpq1-Cpq3-Cpu)*0.5_dp/sqrt(sw2*cw2)
      gw_smeft=gw_sm*(1._dp+voL2*Cpq3)  

      gvt_smeft2=((Cpq1 - Cpq3 + Cpu)**2*voL4 - 
     -    8._dp*(Cpq1 - Cpq3 + Cpu)*Sqrt(cw2*sw2)*voL2*gvt_sm + 
     -    16._dp*cw2*sw2*gvt_sm**2)/(16._dp*cw2*sw2)

      gat_smeft2=((-Cpq1 + Cpq3 + Cpu)**2*voL4 + 
     -    8._dp*gat_sm*((-Cpq1 + Cpq3 + Cpu)*Sqrt(cw2*sw2)*voL2 + 
     -       2._dp*cw2*sw2*gat_sm))/(16._dp*cw2*sw2)

      gw_smeft2=(1._dp + 2._dp*Cpq3*voL2 + Cpq3**2*voL4)*gw_sm**2
      
      
      end subroutine smeft_coupl

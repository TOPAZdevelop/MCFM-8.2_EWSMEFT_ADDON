      subroutine ggQQb_ew_oneloop(vew,s,beta,z)
C --- the analytical result is taken from [arXiv:hep-ph/0610335] 
C --- by J. H. Kuehn, A. Scharf and P. Uwer
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'ewcouple.f'
      include 'masses.f'
      include 'scale.f'
      include 'anomcoup.f'
      real(dp):: vew,s,t,beta,z,vrts(5),slf(5),vrt(5),bx(5),alpha,
     .     genfacs,genfacself,genfacvert,genfacbox,genfacschi,sigma0,
     .     chifacs,chifac,gvt,gat,xI1,xI2,xI3,xI4,sw2,cw2,
     .     mw,mz,mh,born,aemmz,trih,trizx
      real(dp) :: db0
      integer ep,i
      real(dp):: gvt_smeft, gat_smeft, gw_smeft, 
     .     gvt_smeft2, gat_smeft2, gw_smeft2, vol2, vol4, Cpq3, 
     .     Cpu, c(1:7)            
      common/em/aemmz
      real(dp)::trizx2,vrts3,vrts4,bx3,bx4,vrt3,vrt4,slf3,slf4,PreFac

      ep = 0
c      musq = (2._dp*mt)**2
c      alpha = 1._dp/126.3_dp
c--- alpha -> standard MCFM value
      alpha=aemmz
      sw2 = xw


      mw = wmass
      mz = zmass
      mh = hmass

c      sw2 = 1._dp - mw**2/mz**2
      cw2 = 1._dp-sw2

      gvt = (0.5_dp-4._dp*sw2/3._dp)/2._dp/sqrt(sw2*cw2)
      gat = 0.5_dp/2._dp/sqrt(sw2*cw2)
      gw = 0.5_dp/sqrt(2._dp*sw2)
      
!-----BEGIN TILL (SMEFT couplings)
      
      call smeft_coupl(gvt,gat,gw,gvt_smeft,gat_smeft,gvt_smeft2,gat_smeft2,gw_smeft2,voL2,voL4,Cpq3,Cpu,c)
     
!-----END TILL      

      genfacs = Nc*mt**2/s**2*z**2/(1._dp-beta**2*z**2)
      genfacself = (2._dp - Nc**2*(1._dp-beta*z))/Nc/(1._dp-beta**2*z**2)
      genfacvert = (2._dp - Nc**2*(1._dp-beta*z))/Nc/(1._dp+beta*z)**2
      genfacbox = (Nc**2*(1._dp-beta*z)-2._dp)/Nc/(1._dp-beta**2*z**2)
      genfacschi = Nc*mt**4/s**2*z**2/(1._dp-beta**2*z**2)
      chifacs = gat**2
      chifac = gat**2/mz**2
c --- sigma0 is jacobian additional to matrix element square
c --- sigma0 = gsq**2*beta/64._dp/pi/s/(Nc**2-1._dp)
      sigma0 = 1._dp
c      sigma0 = 1._dp/16._dp

      t = -s*(1._dp+beta*z)/2._dp+mt**2
      
C --- leading order differential cross section, use to normalize the result
C --- in the ref. (II.1)
      born = sigma0*(Nc**2*(1._dp+beta**2*z**2)-2._dp)/
     .     Nc/(1._dp-beta**2*z**2)**2*(1._dp-beta**4*z**4
     .     +2._dp*beta**2*(1._dp-beta**2)*(1._dp-z**2))
C --- the differential cross section is normalized when 
C --- sigma0 = (pi^2*alpha_s^2)/2 = gsq^2/32

C --- add triangle t,b contribution (II.16-17)
      trih = 
     .     + alpha/pi*sigma0*mt**2/mw**2/sw2*tevscale*
     .     beta**2/(1._dp - beta**2*z**2)/(s - mh**2)*
     .     ( + tevscale*
     .     ( + mt**2*(s - 4._dp*mt**2)*
     .     xI3(0._dp,0._dp,s,mt**2,mt**2,mt**2,musq,ep) - 2._dp*mt**2)
     .     + mb**2*(s - 4._dp*mb**2)*
     .     xI3(0._dp,0._dp,s,mb**2,mb**2,mb**2,musq,ep) - 2._dp*mb**2)

      
      trizx = 
     .     + 16._dp*alpha/pi*sigma0*gat_smeft2*mt**2/mz**2
     .     /(1._dp-beta**2*z**2)*(
     .     + mt**2*xI3(0._dp,0._dp,s,mt**2,mt**2,mt**2,musq,ep)
     .     )


C --- In the ref. (II.18-22), and Appendix B (B.1-15).
      vrts(1) = 
     .     (-2._dp*alpha*genfacs*sigma0*
     .     (beta**2*s*(2._dp*(gat_smeft2*c(5) + gvt_smeft2*c(4))*mz**2 
     .     + (1._dp - beta**2)*(-3._dp*gat_smeft2*c(5) + gvt_smeft2*c(4))*s)*
     .     db0(mt**2,mt**2,mz**2) + 2._dp*(2._dp*(gat_smeft2*c(5) + gvt_smeft2*c(4))*mz**2 
     .     - beta**2*(-3._dp*gat_smeft2*c(5) + gvt_smeft2*c(4))*s)*( - xI2(mt**2,
     .     mt**2,mz**2,musq,ep) + xI2(s,mt**2,mt**2,musq,ep)) + (4._dp*
     .     (gat_smeft2*c(5) + gvt_smeft2*c(4))*mz**4 + 8._dp*beta**2*gat_smeft2*c(5)*mz**2*s
     .     - beta**2*(gat_smeft2*c(5) + gvt_smeft2*c(4) + beta**2*(-3._dp*gat_smeft2*c(5) + gvt_smeft2*c(4))
     .     )*s**2)*xI3(s,mt**2,mt**2,mt**2,mt**2,mz**2,musq,ep)))/pi  


      vrts(2) = 
     .     (-8._dp*alpha*genfacs*gw_smeft2*c(6)*sigma0*
     .     (-0.25_dp*beta**2*s*(4._dp*(mb**2 - mw**2) + (1._dp 
     .     - beta**2)*s)*db0(mt**2,mb**2,mw**2) + 0.5_dp*(4._dp*
     .     (-mb**2 + mw**2) + (1._dp + beta**2)*s)*
     .     (-xI2(mt**2,mb**2,mw**2,musq,ep) 
     .     + xI2(s,mb**2,mb**2,musq,ep)) 
     .     + 0.125_dp*(-4._dp*beta**2*s**2 + (4._dp*(-mb**2 + mw**2) 
     .     + (1._dp + beta**2)*s)**2)*
     .     xI3(s,mt**2,mt**2,mb**2,mb**2,mw**2,musq,ep)))/pi


      vrts(3) = 
     .     (-8._dp*alpha*chifacs*genfacschi*sigma0*
     .     (beta**2*s*(1 + 4*Cpq3*vol2 + 2*Cpu*vol2 + (2*Cpq3 + Cpu)**2*vol4)*DB0(MT**2,MT**2,MZ**2) - 
     -  2*(1 + 4*Cpq3*vol2 + 2*Cpu*vol2 + (2*Cpq3 + Cpu)**2*vol4)*xI2(MT**2,MT**2,MZ**2,musq,ep) + 
     -  2*(1 + 4*Cpq3*vol2 + 2*Cpu*vol2 + (2*Cpq3 + Cpu)**2*vol4)*xI2(s,MT**2,MT**2,musq,ep) + 
     -  (2*MZ**2 + beta**2*s)*(1 + 4*Cpq3*vol2 + 2*Cpu*vol2 + (2*Cpq3 + Cpu)**2*vol4)*
     -   xI3(MT**2,MT**2,s,MT**2,MZ**2,MT**2,musq,ep)))/pi
     
      vrts(4) = 
     .     -(-8._dp*alpha*genfacs*gw**2*sigma0*
     .     ((beta**2*(-1 + beta**2)*s**2*(4*MW**2 + (-1 + beta**2)*s)*(1 + 2*Cpq3*vol2 + Cpq3**2*vol4)*
     -     DB0(MT**2,MB**2,MW**2))/32._dp - ((-1 + beta**2)*s*(4*MW**2 + s + beta**2*s)*
     -     (1 + 2*Cpq3*vol2 + Cpq3**2*vol4)*xI2(MT**2,0._dp,MW**2,musq,ep))/16._dp + 
     -  ((-1 + beta**2)*s*(4*MW**2 + s + beta**2*s)*(1 + 2*Cpq3*vol2 + Cpq3**2*vol4)*xI2(s,0._dp,0._dp,musq,ep))/
     -   16._dp + ((-1 + beta**2)*s*(4*MW**2 + (-1 + beta)**2*s)*(4*MW**2 + (1 + beta)**2*s)*
     -     (1 + 2*Cpq3*vol2 + Cpq3**2*vol4)*xI3(MT**2,MT**2,s,0._dp,MW**2,0._dp,musq,ep))/64._dp))
     .     /(mw**2*pi)



      vrts(5) = 
     .     (4._dp*alpha*genfacs*gw**2*mt**2*sigma0*
     .     (beta**2*s*(-mh**2 + (1._dp - beta**2)*s)*db0(mt**2
     .     ,mt**2,mh**2) + 2._dp*(mh**2 + beta**2*s)*(xI2(mt**2,mt**2,
     .     mh**2,musq,ep) - xI2(s,mt**2,mt**2,musq,ep)) + 
     .     (-2._dp*mh**4 - 3._dp*beta**2*mh**2*s + beta**2*(1._dp 
     .     - beta**2)*s**2)*xI3(s,mt**2,mt**2,mt**2,mt**2,mh**2
     .     ,musq,ep)))/(mw**2*pi)



      slf(1) = 
     .     (0.125_dp*alpha*genfacself*sigma0*
     .     ((2._dp*(2._dp*(gat_smeft2*c(5) + gvt_smeft2*c(4))*mz**2 + 
     .     (1._dp - beta**2)*(-3._dp*gat_smeft2*c(5) + gvt_smeft2*c(4))*s)*
     .     (1._dp - beta**4*z**4 + (1._dp - beta**2)*
     .     (2._dp*beta**2 - beta*z - 3._dp*beta**2*z**2))*
     .     db0(mt**2,mt**2,mz**2))/(1._dp + beta*z) + 
     .     (16._dp*(gat_smeft2*c(5) + gvt_smeft2*c(4))*(1._dp - beta**4*z**4 + 
     .     beta**2*(1._dp - beta**2)*(1._dp - 3._dp*z**2))*
     .     (-xI1(mt**2,musq,ep) + xI1(mz**2,musq,ep)))/ 
     .     ((1._dp - beta**2)*s*(1._dp + beta**2 + 2._dp*beta*z))  
     .     - (4._dp*((1._dp - beta**2)*(-3._dp*gat_smeft2*c(5) + gvt_smeft2*c(4))*
     .     (1._dp + beta**2 - 2._dp*beta**4 -
     .     beta**2*(2._dp - 3._dp*beta**2)*z**2 - beta**4*z**4) 
     .     + (2._dp*(gat_smeft2*c(5) + gvt_smeft2*c(4))*mz**2*(2._dp + 2._dp*beta**2 - 
     .     5._dp*beta**4 + 2._dp*beta**6 + beta**3*(3._dp - 2._dp*beta**2)*z 
     .     - 3._dp*beta**2*(2._dp - 3._dp*beta**2 + beta**4)*z**2 - 
     .     3._dp*beta**3*(1._dp - beta**2)*z**3 - beta**4*(2._dp 
     .     - beta**2)*z**4 - beta**5*z**5))/((1._dp 
     .     - beta**2)*s))*xI2(mt**2,mt**2,mz**2,musq,ep))/(1._dp 
     .     + beta*z)**2 + (4._dp*((2._dp*beta**2*(gat_smeft2*c(5) + gvt_smeft2*c(4))*mz**2*
     .     (1._dp - z**2)*(2._dp + beta**2 - 2._dp*beta**4 + 
     .     (3._dp*beta - 2._dp*beta**3)*z + beta**4*z**2 + beta**3*z**3))
     .     /s + gvt_smeft2*c(4)*(2._dp + 2._dp*beta**2 - 4._dp*beta**4 - 
     .     beta**6 + 2._dp*beta**8 + beta*(4._dp + 2._dp*beta**2 - 
     .     8._dp*beta**4 + 4._dp*beta**6)*z + beta**2*(-4._dp + 
     .     7._dp*beta**2 + beta**4 - 3._dp*beta**6)*z**2 - beta**3*
     .     (10._dp - 16._dp*beta**2 + 6._dp*beta**4)*z**3 + beta**4*(-5._dp 
     .     + 3._dp*beta**2 + beta**4)*z**4 - beta**5*(4._dp 
     .     - 2._dp*beta**2)*z**5 - beta**6*z**6) - gat_smeft2*c(5)*
     .     (2._dp + 2._dp*beta**2 - 8._dp*beta**4 - 3._dp*beta**6 
     .     + 6._dp*beta**8 + 2._dp*beta*(2._dp - beta**2 
     .     - 8._dp*beta**4 + 6._dp*beta**6)*z - beta**2*(4._dp 
     .     - 5._dp*beta**2 - 7._dp*beta**4 + 9._dp*beta**6)*z**2 
     .     - 6._dp*beta**3*(1._dp - 4._dp*beta**2 + 3._dp*beta**4)*z**3 
     .     + beta**4*(1._dp - 3._dp*beta**2 + 3._dp*beta**4)*z**4 
     .     - beta**5*(4._dp - 6._dp*beta**2)*z**5 + beta**6*z**6))*
     .     xI2(t,mt**2,mz**2,musq,ep))/((1._dp + beta*z)**2*(1._dp 
     .     + beta**2 + 2._dp*beta*z))))/pi


      slf(2) = 
     .     (alpha*genfacself*gw_smeft2*c(6)*sigma0*
     .     ((0.25_dp*(-4._dp*mb**2 + 4._dp*mw**2 - s + beta**2*s)*
     .     (1._dp - beta**4*z**4 + beta*(1._dp - beta**2)*
     .     (2._dp*beta - z - 3._dp*beta*z**2))*
     .     db0(mt**2,mb**2,mw**2))/(1._dp + beta*z) + (4._dp*(1._dp - 
     .     beta**4*z**4 + beta**2*(1._dp - beta**2)*(1._dp - 
     .     3._dp*z**2))*(-xI1(mb**2,musq,ep) + xI1(mw**2,musq,ep)))
     .     /((1._dp - beta**2)*s*(1._dp + beta**2 + 2._dp*beta*z)) + 
     .     (0.5_dp*((-4._dp*(mb**2 - mw**2)*(-2._dp - 2._dp*beta**2 + 
     .     5._dp*beta**4 - 2._dp*beta**6 - beta**3*(3._dp - 
     .     2._dp*beta**2)*z + 3._dp*beta**2*(2._dp - 3._dp*beta**2 + 
     .     beta**4)*z**2 + 3._dp*beta**3*(1._dp - beta**2)*z**3 + 
     .     beta**4*(2._dp - beta**2)*z**4 + beta**5*z**5))/(1._dp - 
     .     beta**2) - beta**2*s*(1._dp - z**2)*(2._dp + 
     .     beta**2 - 2._dp*beta**4 + beta*(3._dp - 2._dp*beta**2)*z + 
     .     beta**3*z**2*(beta + z)))*xI2(mt**2,mb**2,mw**2,musq,ep))/
     .     (s*(1._dp + beta*z)**2) + (0.5_dp*beta**2*(1._dp - z**2)*
     .     (2._dp + beta**2 - 2._dp*beta**4 + beta*(3._dp - 2._dp*beta**2)*z 
     .     + beta**3*z**2*(beta + z))*(4._dp*(-mb**2 + mw**2) + 
     .     s*(1._dp + beta**2 + 2._dp*beta*z))*xI2(t,mb**2,mw**2,musq,ep))
     .     /(s*(1._dp + beta*z)**2*(1._dp + beta**2 + 2._dp*beta*z))))/pi


      slf(3) = 
     .     (2._dp*alpha*chifac*genfacself*sigma0*
     .     (-(s*(-3 + beta*z)*(4*beta*(-1 + beta**2)*Cpq3*vol2*z + 
     -        4*Cpq3**2*vol4*(1 - beta*z + beta**3*z - beta**2*z**2) + 
     -        Cpu*(2*beta*(-1 + beta**2)*vol2*z + Cpu*vol4*(1 - beta*z + beta**3*z - beta**2*z**2))))/16._dp
     -    + ((-1 + beta**2)*MZ**2*s*(1 + 2*Cpu*vol2 + 4*Cpq3**2*vol4 + Cpu**2*vol4 + 
     -       4*Cpq3*(vol2 + Cpu*vol4))*(-1 + beta*z - beta**3*z + beta**2*(-2 + 3*z**2) + 
     -       beta**4*(2 - 3*z**2 + z**4))*DB0(MT**2,MT**2,MZ**2))/(8 + 8*beta*z) + 
     -  ((2 + 4*Cpq3*vol2 + 2*Cpu*vol2 + 4*Cpq3**2*vol4 + Cpu**2*vol4)*
     -     (-1 + beta**2*(-1 + 3*z**2) + beta**4*(1 - 3*z**2 + z**4))*xI1(MT**2,musq,ep))/
     -   (4._dp*(1 + beta**2 + 2*beta*z)) - ((2 + 4*Cpq3*vol2 + 2*Cpu*vol2 + 4*Cpq3**2*vol4 + Cpu**2*vol4)*
     -     (-1 + beta**2*(-1 + 3*z**2) + beta**4*(1 - 3*z**2 + z**4))*xI1(MZ**2,musq,ep))/
     -   (4._dp*(1 + beta**2 + 2*beta*z)) + (MZ**2*
     -     (-2 - 6*Cpq3*vol2 - 3*Cpu*vol2 - 4*Cpq3**2*vol4 - 4*Cpq3*Cpu*vol4 - Cpu**2*vol4 + 
     -       3*beta**3*(1 + 2*Cpq3*vol2 + Cpu*vol2)*z*(-1 + z**2) - 
     -       beta**6*(1 + 2*Cpu*vol2 + 4*Cpq3**2*vol4 + Cpu**2*vol4 + 4*Cpq3*(vol2 + Cpu*vol4))*
     -        (2 - 3*z**2 + z**4) + beta**5*(1 + 2*Cpq3*vol2 + Cpu*vol2)*z*(2 - 3*z**2 + z**4) + 
     -       2*beta**2*(-1 + 3*z**2 + 4*Cpq3**2*vol4*z**2 + Cpu**2*vol4*z**2 + Cpu*vol2*(-1 + 4*z**2) + 
     -          Cpq3*(4*Cpu*vol4*z**2 + vol2*(-2 + 8*z**2))) + 
     -       beta**4*(5 - 9*z**2 + 2*z**4 + 4*Cpq3**2*vol4*(3 - 5*z**2 + z**4) + 
     -          Cpu**2*vol4*(3 - 5*z**2 + z**4) + Cpu*vol2*(8 - 14*z**2 + 3*z**4) + 
     -          2*Cpq3*(2*Cpu*vol4*(3 - 5*z**2 + z**4) + vol2*(8 - 14*z**2 + 3*z**4))))*
     -     xI2(MT**2,MT**2,MZ**2,musq,ep))/(4._dp*(1 + beta*z)**2) - 
     -  ((1 - beta*z)*(s*(1 + 2*Cpq3*vol2 + Cpu*vol2) + 
     -       2*MZ**2*(Cpu*vol2 + 2*Cpq3*(vol2 + 2*Cpu*vol4)) + 
     -       beta*(s*(2 + 6*Cpq3*vol2 + 3*Cpu*vol2) + 4*MZ**2*(Cpu*vol2 + 2*Cpq3*(vol2 + 2*Cpu*vol4)))*
     -        z + beta**8*(2*MZ**2*(1 + 2*Cpu*vol2 + 4*Cpq3**2*vol4 + Cpu**2*vol4 + 
     -             4*Cpq3*(vol2 + Cpu*vol4))*(2 - 3*z**2 + z**4) + 
     -          s*z**2*(1 - 3*z**2 + z**4 + 2*Cpq3*vol2*(2 - 3*z**2 + z**4) + 
     -             Cpu*vol2*(2 - 3*z**2 + z**4) + 4*Cpq3**2*vol4*(2 - 3*z**2 + z**4) + 
     -             Cpu**2*vol4*(2 - 3*z**2 + z**4))) + 
     -       beta**7*z*(2*MZ**2*(1 + 6*Cpq3*vol2 + 3*Cpu*vol2 + 8*Cpq3**2*vol4 + 8*Cpq3*Cpu*vol4 + 
     -             2*Cpu**2*vol4)*(2 - 3*z**2 + z**4) + 
     -          s*(2*(1 - 3*z**2 + z**4) + 2*Cpq3*vol2*(4 - 5*z**2 + 2*z**4) + 
     -             Cpu*vol2*(4 - 5*z**2 + 2*z**4) + 4*Cpq3**2*vol4*(4 - 4*z**2 - z**4 + z**6) + 
     -             Cpu**2*vol4*(4 - 4*z**2 - z**4 + z**6))) + 
     -       beta**3*z*(s*(-6*z**2 + 28*Cpq3**2*vol4*(-1 + z**2) + 7*Cpu**2*vol4*(-1 + z**2) - 
     -             2*Cpq3*vol2*(4 + z**2) - Cpu*vol2*(4 + z**2)) + 
     -          2*MZ**2*(3 - 3*z**2 + Cpu*vol2*(3 - 7*z**2) + 8*Cpq3**2*vol4*(-1 + z**2) + 
     -             2*Cpu**2*vol4*(-1 + z**2) - 2*Cpq3*(8*Cpu*vol4*z**2 + vol2*(-3 + 7*z**2)))) - 
     -       beta**2*(s*(2*z**2 + Cpq3*vol2*(2 - 4*z**2) - 8*Cpq3**2*vol4*(-1 + z**2) - 
     -             2*Cpu**2*vol4*(-1 + z**2) + Cpu*(vol2 - 2*vol2*z**2)) + 
     -          2*MZ**2*(2*(-1 + z**2) + Cpu*vol2*(-3 + 4*z**2) + 
     -             Cpq3*(4*Cpu*vol4*(-1 + 2*z**2) + vol2*(-6 + 8*z**2)))) + 
     -       beta**6*(s*(1 - 5*z**2 + 7*z**4 - z**6 - 2*Cpq3*vol2*(-2 + 3*z**2 - 5*z**4 + z**6) - 
     -             Cpu*vol2*(-2 + 3*z**2 - 5*z**4 + z**6) + 8*Cpq3**2*vol4*(1 - 2*z**4 + z**6) + 
     -             2*Cpu**2*vol4*(1 - 2*z**4 + z**6)) + 
     -          2*MZ**2*(-3 + 4*z**2 - z**4 + 4*Cpq3**2*vol4*(-1 + z**2)**3 + 
     -             Cpu**2*vol4*(-1 + z**2)**3 - Cpu*vol2*(4 - 6*z**2 + z**4) - 
     -             2*Cpq3*(2*Cpu*vol4*(1 - 2*z**2) + vol2*(4 - 6*z**2 + z**4)))) + 
     -       beta**4*(s*(-2 + 6*z**2 - 4*z**4 - 2*Cpq3*vol2*(2 + z**2 + 2*z**4) - 
     -             Cpu*vol2*(2 + z**2 + 2*z**4) + 4*Cpq3**2*vol4*(-1 - 7*z**2 + 8*z**4) + 
     -             Cpu**2*vol4*(-1 - 7*z**2 + 8*z**4)) + 
     -          2*MZ**2*(-1 + z**2 - Cpu*vol2*(-2 + z**2)**2 + 4*Cpq3**2*vol4*(-2 - z**2 + 3*z**4) + 
     -             Cpu**2*vol4*(-2 - z**2 + 3*z**4) - 
     -             2*Cpq3*(vol2*(-2 + z**2)**2 + 2*Cpu*vol4*(3 - 3*z**2 + z**4)))) - 
     -       beta**5*z*(s*(-16*Cpq3**2*vol4*z**2*(-1 + z**2) - 4*Cpu**2*vol4*z**2*(-1 + z**2) + 
     -             2*(2 - 6*z**2 + z**4) + 2*Cpq3*vol2*(3 - 6*z**2 + 2*z**4) + 
     -             Cpu*vol2*(3 - 6*z**2 + 2*z**4)) + 
     -          2*MZ**2*(5 - 6*z**2 + z**4 - 16*Cpq3**2*vol4*(-1 + z**2) - 4*Cpu**2*vol4*(-1 + z**2) + 
     -             Cpu*vol2*(11 - 16*z**2 + 3*z**4) + 
     -             Cpq3*(8*Cpu*vol4*(3 - 5*z**2 + z**4) + vol2*(22 - 32*z**2 + 6*z**4)))))*
     -     xI2(MT**2 - (s*(1 + beta*z))/2._dp,MT**2,MZ**2,musq,ep))/
     -   (8._dp*(-1 + beta*z)*(1 + beta*z)**2*(1 + beta**2 + 2*beta*z))))/pi

      slf(4) = 
     .     (alpha*genfacself*gw**2*sigma0*
     .     (-(Cpq3*s*(-3 + beta*z)*(2*beta*(-1 + beta**2)*vol2*z + 
     -        Cpq3*vol4*(1 - beta*z + beta**3*z - beta**2*z**2)))/8._dp + 
     -  ((-1 + beta**2)*s*(4*MW**2 + (-1 + beta**2)*s)*(1 + 2*Cpq3*vol2 + Cpq3**2*vol4)*
     -     (-1 + beta*z - beta**3*z + beta**2*(-2 + 3*z**2) + beta**4*(2 - 3*z**2 + z**4))*
     -     DB0(MT**2,MB**2,MW**2))/(32._dp*(1 + beta*z)) + 
     -  ((1 + beta**2*(1 - 3*z**2) - beta**4*(1 - 3*z**2 + z**4))*xI1(MW**2,musq,ep))/
     -   (2._dp*(1 + beta**2 + 2*beta*z)) + ((1 - beta*z)*
     -     (4*MW**2*(2 + 2*Cpq3*vol2 + 3*beta**3*(-1 + Cpq3**2*vol4)*z*(-1 + z**2) + 
     -          beta**6*(1 + 2*Cpq3*vol2 + Cpq3**2*vol4)*(2 - 3*z**2 + z**4) + 
     -          beta**5*(-1 + Cpq3**2*vol4)*z*(2 - 3*z**2 + z**4) + 
     -          2*beta**2*(1 - 3*z**2 - 2*Cpq3*vol2*z**2 + Cpq3**2*vol4*(-1 + z**2)) + 
     -          beta**4*(-5 + 9*z**2 - 2*z**4 + Cpq3**2*vol4*(-1 + z**2) - 
     -             2*Cpq3*vol2*(3 - 5*z**2 + z**4))) + 
     -       (-1 + beta**2)*s*(-2*Cpq3*(vol2 + Cpq3*vol4) + 
     -          3*beta**3*(1 + 4*Cpq3*vol2 + 3*Cpq3**2*vol4)*z*(-1 + z**2) + 
     -          beta**6*(1 + 2*Cpq3*vol2 + Cpq3**2*vol4)*(2 - 3*z**2 + z**4) + 
     -          beta**5*(1 + 4*Cpq3*vol2 + 3*Cpq3**2*vol4)*z*(2 - 3*z**2 + z**4) + 
     -          2*beta**2*(-1 + z**2 + 2*Cpq3*vol2*(-2 + 3*z**2) + Cpq3**2*vol4*(-3 + 5*z**2)) + 
     -          beta**4*(-1 + z**2 + 2*Cpq3*vol2*(1 - 3*z**2 + z**4) + 
     -             Cpq3**2*vol4*(3 - 7*z**2 + 2*z**4))))*xI2(MT**2,0._dp,MW**2,musq,ep))/
     -   (16._dp*(-1 + beta*z)*(1 + beta*z)**2) + 
     -  ((4*MW**2 + s*(1 + beta**2 + 2*beta*z))*
     -     (2*Cpq3*vol2 + 4*beta*Cpq3*vol2*z + 
     -       beta**8*(1 + 2*Cpq3*vol2 + Cpq3**2*vol4)*(2 - 3*z**2 + z**4) + 
     -       beta**7*(1 + 4*Cpq3*vol2 + 3*Cpq3**2*vol4)*z*(2 - 3*z**2 + z**4) + 
     -       beta**3*z*(3 - 3*z**2 - 8*Cpq3*vol2*z**2 + 7*Cpq3**2*vol4*(-1 + z**2)) + 
     -       2*beta**2*(1 - z**2 + Cpq3**2*vol4*(-1 + z**2) + Cpq3*(vol2 - 2*vol2*z**2)) + 
     -       beta**5*z*(-5 + 6*z**2 - z**4 - 4*Cpq3*vol2*(3 - 5*z**2 + z**4) + 
     -          Cpq3**2*vol4*(-3 + 2*z**2 + z**4)) + 
     -       beta**4*(-1 + z**2 - 2*Cpq3*vol2*(3 - 3*z**2 + z**4) + 
     -          3*Cpq3**2*vol4*(-1 - z**2 + 2*z**4)) + 
     -       beta**6*(-3 + 4*z**2 - z**4 + 2*Cpq3*vol2*(-1 + 2*z**2) + 
     -          Cpq3**2*vol4*(1 + 2*z**2 - 5*z**4 + 2*z**6)))*
     -     xI2(MT**2 - (s*(1 + beta*z))/2._dp,0._dp,MW**2,musq,ep))/
     -   (16._dp*(1 + beta*z)**2*(1 + beta**2 + 2*beta*z))))/(mw**2*pi)



      slf(5) = 
     .     (alpha*genfacself*gw**2*sigma0*
     .     ((-0.125_dp*(1._dp - beta**2)*s*(-mh**2 
     .     + (1._dp - beta**2)*s)*(1._dp - beta**4*z**4 
     .     + beta*(1._dp - beta**2)*(2._dp*beta - z 
     .     - 3._dp*beta*z**2))*db0(mt**2,mt**2,mh**2))/(1._dp + beta*z) 
     .     + (0.5_dp*(1._dp - beta**4*z**4 + beta**2*(1._dp 
     .     - beta**2)*(1._dp - 3._dp*z**2))*(xI1(mh**2,musq,ep) 
     .     - xI1(mt**2,musq,ep)))/(1._dp + beta**2 + 2._dp*beta*z) 
     .     + (0.25_dp*((1._dp - beta**2)**2*s*(1._dp + beta**2 
     .     - 2._dp*beta**4 - beta**2*(2._dp - 3._dp*beta**2)*z**2 
     .     - beta**4*z**4) + mh**2*(-2._dp - 2._dp*beta**2 
     .     + 5._dp*beta**4 - 2._dp*beta**6 - beta**3*(3._dp 
     .     - 2._dp*beta**2)*z + 3._dp*beta**2*(2._dp - 3._dp*beta**2 
     .     + beta**4)*z**2 + 3._dp*beta**3*(1._dp - beta**2)*z**3 
     .     + beta**4*(2._dp - beta**2)*z**4 + beta**5*z**5))*
     .     xI2(mt**2,mt**2,mh**2,musq,ep))/(1._dp + beta*z)**2 - 
     .     (0.125_dp*(1._dp - beta**2)*(s*(1._dp + beta**2 
     .     - 5._dp*beta**4 - 2._dp*beta**6 + 4._dp*beta**8 
     .     + 2._dp*beta*(1._dp - beta**2 - 5._dp*beta**4 
     .     + 4._dp*beta**6)*z - beta**2*(2._dp - 2._dp*beta**2 
     .     - 5._dp*beta**4 + 6._dp*beta**6)*z**2 - 2._dp*beta**3*(1._dp 
     .     - 7._dp*beta**2 + 6._dp*beta**4)*z**3 + beta**4*(2._dp 
     .     - 3._dp*beta**2 + 2._dp*beta**4)*z**4 - 2._dp*beta**5*(1._dp 
     .     - 2._dp*beta**2)*z**5 + beta**6*z**6) - 2._dp*beta**2*mh**2*
     .     (1._dp - z**2)*(2._dp + beta**2 - 2._dp*beta**4 + 
     .     beta*(3._dp - 2._dp*beta**2)*z + beta**3*z**2*(beta + z)))*
     .     xI2(t,mt**2,mh**2,musq,ep))/((1._dp + beta*z)**2*
     .     (1._dp + beta**2 + 2._dp*beta*z))))/(mw**2*pi)



      vrt(1) = 
     .     (0.125_dp*alpha*genfacvert*sigma0*
     .     (4._dp*(1._dp - beta**2)*(gat_smeft2*c(5) + gvt_smeft2*c(4)) 
     .     - (4._dp*(2._dp*(gat_smeft2*c(5) + gvt_smeft2*c(4))*mz**2 + 
     .     (1._dp - beta**2)*(-3._dp*gat_smeft2*c(5) + gvt_smeft2*c(4))*s)*
     .     (1._dp - beta**4*z**4 + beta*(1._dp - beta**2)*
     .     (2._dp*beta - z - 3._dp*beta*z**2))*db0(mt**2,mt**2,mz**2))/
     .     (1._dp - beta*z) + (16._dp*(gat_smeft2*c(5) + gvt_smeft2*c(4))*(1._dp 
     .     + 2._dp*beta**2 - beta**4 + beta*(1._dp + 4._dp*beta**2 
     .     - 3._dp*beta**4)*z - 2._dp*beta**4*z**4*(1._dp + beta*z) 
     .     - 2._dp*(1._dp - beta**2)*(2._dp*beta**2*z**2 
     .     + 3._dp*beta**3*z**3))*(xI1(mt**2,musq,ep) 
     .     - xI1(mz**2,musq,ep)))/((1._dp - beta**2)*s*
     .     (1._dp - beta*z)*(1._dp + beta**2 + 2._dp*beta*z)) 
     .     + (4._dp*((2._dp*(gat_smeft2*c(5) + gvt_smeft2*c(4))*mz**2*(1._dp + 8._dp*beta**2 
     .     - 11._dp*beta**4 + 4._dp*beta**6 + beta*(1._dp + 4._dp*beta**2 
     .     - 3._dp*beta**4)*z - 2._dp*beta**2*(5._dp - 8._dp*beta**2 
     .     + 3._dp*beta**4)*z**2 - 6._dp*beta**3*(1._dp - beta**2)*
     .     z**3 - 2._dp*beta**4*(2._dp - beta**2)*z**4 
     .     - 2._dp*beta**5*z**5))/s + (1._dp - beta**2)**2*
     .     (gvt_smeft2*c(4)*(1._dp + beta**2 - 4._dp*beta**4 
     .     - beta*(1._dp - beta**2)*z - 2._dp*beta**2*(1._dp 
     .     - 3._dp*beta**2)*z**2 - 2._dp*beta**4*z**4) + gat_smeft2*c(5)*(1._dp 
     .     - 7._dp*beta**2 + 12._dp*beta**4 - beta*(1._dp 
     .     - beta**2)*z + 6._dp*beta**2*(1._dp - 3._dp*beta**2)*z**2 
     .     + 6._dp*beta**4*z**4)))*xI2(mt**2,mt**2,mz**2,musq,ep))/
     .     ((1._dp - beta**2)*(1._dp - beta**2*z**2)) + 
     .     (4._dp*((2._dp*(gat_smeft2*c(5) + gvt_smeft2*c(4))*mz**2*(1._dp - 4._dp*beta**2 
     .     - 3._dp*beta**4 + 4._dp*beta**6 + beta*(1._dp - 8._dp*beta**2 
     .     + 5._dp*beta**4)*z + 2._dp*beta**2*(1._dp + 2._dp*beta**2 
     .     - 3._dp*beta**4)*z**2 + 6._dp*beta**3*(1._dp - beta**2)*
     .     z**3 + 2._dp*beta**6*z**4 + 2._dp*beta**5*z**5))/s 
     .     - gvt_smeft2*c(4)*(2._dp + 3._dp*beta**2 - 6._dp*beta**4 
     .     - beta**6 + 4._dp*beta**8 + beta*(3._dp + 7._dp*beta**2 
     .     - 13._dp*beta**4 + 7._dp*beta**6)*z - beta**2*(7._dp 
     .     - 18._dp*beta**2 + 3._dp*beta**4 + 6._dp*beta**6)*z**2 
     .     - beta**3*(14._dp - 26._dp*beta**2 + 12._dp*beta**4)*z**3 
     .     - beta**4*(10._dp - 6._dp*beta**2 - 2._dp*beta**4)*z**4 
     .     - beta**5*(8._dp - 4._dp*beta**2)*z**5 
     .     - 2._dp*beta**6*z**6) + gat_smeft2*c(5)*(-2._dp + 5._dp*beta**2 
     .     - 10._dp*beta**4 - 7._dp*beta**6 + 12._dp*beta**8 
     .     - beta*(3._dp - 9._dp*beta**2 + 35._dp*beta**4 
     .     - 25._dp*beta**6)*z - beta**2*(1._dp - 6._dp*beta**2 
     .     - 11._dp*beta**4 + 18._dp*beta**6)*z**2 - 2._dp*beta**3*(1._dp 
     .     - 19._dp*beta**2 + 18._dp*beta**4)*z**3 + 2._dp*beta**4*(1._dp 
     .     - 3._dp*beta**2 + 3._dp*beta**4)*z**4 - 4._dp*beta**5*(2._dp 
     .     - 3._dp*beta**2)*z**5 + 2._dp*beta**6*z**6))*
     .     xI2(t,mt**2,mz**2,musq,ep))/((1._dp + beta**2 + 2._dp*beta*z)*
     .     (1._dp - beta**2*z**2)) + 2._dp*(1._dp - beta**2)**2*
     .     (gat_smeft2*c(5) + gvt_smeft2*c(4))*s*
     .     xI3(0._dp,mt**2,t,mt**2,mt**2,mz**2,musq,ep)))/pi


      vrt(2) = 
     .     (alpha*genfacvert*gw_smeft2*c(6)*sigma0*
     .     (1._dp - beta**2 + (0.5_dp*(4._dp*(mb**2 - mw**2) + 
     .     (1._dp - beta**2)*s)*(1._dp - beta**4*z**4 + 
     .     beta*(1._dp - beta**2)*(2._dp*beta - z 
     .     - 3._dp*beta*z**2))*db0(mt**2,mb**2,mw**2))/
     .     (1._dp - beta*z) + (4._dp*(1._dp + 2._dp*beta**2 
     .     - beta**4 + beta*(1._dp + 4._dp*beta**2 - 3._dp*beta**4)*z 
     .     - 2._dp*beta**4*z**4*(1._dp + beta*z) 
     .     - 2._dp*(1._dp - beta**2)*(2._dp*beta**2*z**2 
     .     + 3._dp*beta**3*z**3))*(xI1(mb**2,musq,ep) 
     .     - xI1(mw**2,musq,ep)))/((1._dp - beta**2)*s*
     .     (1._dp - beta*z)*(1._dp + beta**2 + 2._dp*beta*z)) + 
     .     (0.5_dp*((1._dp - beta**2)*s*(3._dp + 3._dp*beta**4 
     .     - 4._dp*beta**6 - beta*(1._dp - 8._dp*beta**2 
     .     + 5._dp*beta**4)*z - 6._dp*beta**2*(1._dp - beta**4)*z**2 
     .     - 6._dp*beta**3*(1._dp - beta**2)*z**3 - 2._dp*beta**6*z**4 
     .     - 2._dp*beta**5*z**5) + 4._dp*(-mb**2 + mw**2)*(1._dp 
     .     + 8._dp*beta**2 - 11._dp*beta**4 + 4._dp*beta**6 + beta*(1._dp 
     .     + 4._dp*beta**2 - 3._dp*beta**4)*z - 2._dp*beta**2*(5._dp 
     .     - 8._dp*beta**2 + 3._dp*beta**4)*z**2 - 6._dp*beta**3*(1._dp 
     .     - beta**2)*z**3 - 2._dp*beta**4*(2._dp 
     .     - beta**2)*z**4 - 2._dp*beta**5*z**5))*
     .     xI2(mt**2,mb**2,mw**2,musq,ep))/((1._dp - beta**2)*s*
     .     (1._dp - beta**2*z**2)) + (0.5_dp*(-s*(1._dp 
     .     + beta**2 + 2._dp*beta*z)*(3._dp + 3._dp*beta**4 - 4._dp*beta**6 
     .     - beta*(1._dp - 8._dp*beta**2 + 5._dp*beta**4)*z 
     .     - 6._dp*beta**2*(1._dp - beta**4)*z**2 - 6._dp*beta**3*
     .     (1._dp - beta**2)*z**3 - 2._dp*beta**6*z**4 
     .     - 2._dp*beta**5*z**5) + 4._dp*(-mb**2 + mw**2)*
     .     (1._dp - 4._dp*beta**2 - 3._dp*beta**4 + 4._dp*beta**6 
     .     + beta*(1._dp - 8._dp*beta**2 + 5._dp*beta**4)*z 
     .     + 2._dp*beta**2*(1._dp + 2._dp*beta**2 - 3._dp*beta**4)*z**2 
     .     + 6._dp*beta**3*(1._dp - beta**2)*z**3 
     .     + 2._dp*beta**6*z**4 + 2._dp*beta**5*z**5))*
     .     xI2(t,mb**2,mw**2,musq,ep))/(s*(1._dp + beta**2 
     .     + 2._dp*beta*z)*(1._dp - beta**2*z**2)) 
     .     + 2._dp*(1._dp - beta**2)*mb**2*
     .     xI3(0._dp,mt**2,t,mb**2,mb**2,mw**2,musq,ep)))/pi



      vrt(3) = 
     .     (2._dp*alpha*chifac*genfacvert*sigma0*
     .     ((s*(-3 - 6*Cpu*vol2 + 20*Cpq3**2*vol4 + 5*Cpu**2*vol4 - 12*Cpq3*(vol2 + Cpu*vol4) - 
     -       beta*(-3 + 15*(2*Cpq3 + Cpu)*vol2 + 2*(4*Cpq3**2 + Cpu**2)*vol4)*z + 
     -       beta**2*(6 + 18*Cpq3*vol2 + 9*Cpu*vol2 + 8*Cpq3**2*vol4 + 12*Cpq3*Cpu*vol4 + 
     -          2*Cpu**2*vol4 - 3*(3*(2*Cpq3 + Cpu)*vol2 + 4*(4*Cpq3**2 - Cpq3*Cpu + Cpu**2)*vol4)*z**2)
     -        + beta**5*z*(3 + 6*Cpq3*vol2 + 3*Cpu*vol2 + 4*Cpq3**2*vol4 + Cpu**2*vol4 - 
     -          3*(4*Cpq3*vol2 + 4*Cpq3**2*vol4 + Cpu*(2*vol2 + Cpu*vol4))*z**2) - 
     -       beta**6*(4*Cpq3**2 + Cpu**2)*vol4*(2 - 3*z**2 + z**4) + 
     -       beta**4*(-3 - 6*Cpq3*vol2 - 3*Cpu*vol2 + 16*Cpq3**2*vol4 + 4*Cpu**2*vol4 - 
     -          3*(2*Cpq3 + Cpu)*(-3*vol2 + (2*Cpq3 + Cpu)*vol4)*z**2 + 4*(4*Cpq3**2 + Cpu**2)*vol4*z**4)
     -         + beta**3*z*(-6 + 4*Cpq3**2*vol4*(7 - 3*z**2) + 12*Cpq3*vol2*(2 + z**2) + 
     -          Cpu*(Cpu*vol4*(7 - 3*z**2) + 6*vol2*(2 + z**2)))))/(24._dp*(-1 + beta*z)) + 
     -  ((-1 + beta**2)*MZ**2*s*(1 + 4*Cpq3*vol2 + 2*Cpu*vol2 + (2*Cpq3 + Cpu)**2*vol4)*
     -     (-1 + beta*(z + beta*(-2 - beta*z + 3*z**2 + beta**2*(2 - 3*z**2 + z**4))))*
     -     DB0(MT**2,MT**2,MZ**2))/(-4 + 4*beta*z) + 
     -  ((-3 - 6*beta**2 + 3*beta**4 - 6*Cpq3*vol2 - 12*beta**2*Cpq3*vol2 + 6*beta**4*Cpq3*vol2 - 
     -       3*Cpu*vol2 - 6*beta**2*Cpu*vol2 + 3*beta**4*Cpu*vol2 - 8*Cpq3**2*vol4 + 
     -       12*beta**4*Cpq3**2*vol4 - 8*beta**6*Cpq3**2*vol4 - 2*Cpu**2*vol4 + 3*beta**4*Cpu**2*vol4 - 
     -       2*beta**6*Cpu**2*vol4 + beta*(3*(-1 - 4*beta**2 + 3*beta**4)*(1 + 2*Cpq3*vol2 + Cpu*vol2) + 
     -          (-2 + beta**2)*(4*Cpq3**2 + Cpu**2)*vol4)*z + 
     -       beta**2*(-1 + beta**2)*(-12*(1 + 2*Cpq3*vol2 + Cpu*vol2) + 
     -          (-4 + 3*beta**2)*(4*Cpq3**2 + Cpu**2)*vol4)*z**2 - 
     -       3*beta**3*(-1 + beta**2)*(6 + 6*(2*Cpq3 + Cpu)*vol2 + (4*Cpq3**2 + Cpu**2)*vol4)*z**3 + 
     -       beta**4*(6 + 6*(2*Cpq3 + Cpu)*vol2 - (-2 + beta**2)*(4*Cpq3**2 + Cpu**2)*vol4)*z**4 + 
     -       beta**5*(6 + 6*(2*Cpq3 + Cpu)*vol2 + (4*Cpq3**2 + Cpu**2)*vol4)*z**5)*xI1(MT**2,musq,ep))/
     -   (6._dp*(-1 + beta*z)*(1 + beta**2 + 2*beta*z)) + 
     -  (((beta**2*(-1 + z**2)*(1 + 2*Cpq3*vol2 + 8*Cpq3**2*vol4 + Cpu*(vol2 + 2*Cpu*vol4) + 
     -            3*beta*(4*Cpq3**2 + Cpu**2)*vol4*z + 
     -            beta**4*(1 + 2*Cpq3*vol2 + 4*Cpq3**2*vol4 + Cpu*(vol2 + Cpu*vol4))*(-2 + z**2) + 
     -            beta**3*(4*Cpq3**2 + Cpu**2)*vol4*z*(-2 + z**2) + 
     -            beta**2*(1 + 2*Cpq3*vol2 + Cpu*vol2 + 4*Cpq3**2*vol4 + Cpu**2*vol4 - 
     -               (1 + 2*Cpq3*vol2 + Cpu*vol2)*z**2)))/(1 + beta**2 + 2*beta*z) - 
     -       (1 + 2*Cpq3*vol2 + 4*Cpq3**2*vol4 + Cpu*(vol2 + Cpu*vol4))*
     -        (-1 + beta*(z + beta*(-2 - beta*z + 3*z**2 + beta**2*(2 - 3*z**2 + z**4)))))*
     -     xI1(MZ**2,musq,ep))/(2._dp*(-1 + beta*z)) - 
     -  ((-1 + beta**2)*(4*Cpq3**2 + Cpu**2)*s*vol4*
     -     (-1 + beta*(z + beta*(-2 - beta*z + 3*z**2 + beta**2*(2 - 3*z**2 + z**4))))*
     -     xI2(0._dp,MT**2,MT**2,musq,ep))/(24._dp*(-1 + beta*z)) + 
     -  ((-4*MZ**2*(1 + beta*z) + 4*MZ**2*(1 + beta*z)*
     -        (-((2*Cpq3 + Cpu)*vol2) + beta*(1 + 2*Cpq3*vol2 + Cpu*vol2)*z - 
     -          beta**3*(1 + 2*Cpq3*vol2 + Cpu*vol2)*z + 
     -          beta**2*(1 + 2*Cpq3*vol2 + Cpu*vol2)*(-2 + 3*z**2) + 
     -          beta**4*(1 + 2*Cpq3*vol2 + Cpu*vol2)*(2 - 3*z**2 + z**4)) + 
     -       (1 - beta**2)*((-1 + beta**2)*s*(1 + 2*Cpq3*vol2 + Cpu*vol2)*(-1 + beta*z)*
     -           (1 + beta**2 + 2*beta*z) + 
     -          2*MZ**2*(1 + 2*Cpq3*vol2 + Cpu*vol2 - beta*(1 + 2*Cpq3*vol2 + Cpu*vol2)*z + 
     -             beta**3*(1 + 2*Cpq3*vol2 + Cpu*vol2)*z + 
     -             2*beta**4*(1 + 4*Cpq3*vol2 + 2*Cpu*vol2 + (2*Cpq3 + Cpu)**2*vol4)*
     -              (2 - 3*z**2 + z**4) + beta**2*
     -              (-3 + 2*z**2 + (2*Cpq3 + Cpu)*
     -                 (2*(2*Cpq3 + Cpu)*vol4*(-1 + z**2) + vol2*(-5 + 4*z**2))))))*
     -     xI2(MT**2,MT**2,MZ**2,musq,ep))/(8._dp*(-1 + beta**2*z**2)) - 
     -  ((beta*s*((-1 + beta)*beta*(1 + beta)*(1 + beta**4) + 
     -          (beta + beta**3 - 3*beta**5 + beta**7)*(2*Cpq3 + Cpu)*vol2 - 
     -          2*beta*(-2 - beta**2 + 2*beta**4)*(4*Cpq3**2 + Cpu**2)*vol4 - 
     -          (-1 + 8*beta**2 - 8*beta**4 + beta**8 + 
     -             (1 - 6*beta**4 + 4*beta**6 + beta**8)*(2*Cpq3 + Cpu)*vol2 + 
     -             2*beta**2*(-7 + 4*beta**4)*(4*Cpq3**2 + Cpu**2)*vol4)*z - 
     -          beta*(-3 + 15*beta**2 - 19*beta**4 + 7*beta**6 + 
     -             (5 + beta**2 - 15*beta**4 + 9*beta**6)*(2*Cpq3 + Cpu)*vol2 + 
     -             2*(2 - 7*beta**2 + 2*beta**6)*(4*Cpq3**2 + Cpu**2)*vol4)*z**2 + 
     -          2*beta**2*(3*(-1 + beta**2)**2 + 2*(-1 + beta**4)*(2*Cpq3 + Cpu)*vol2 + 
     -             (-7 + 4*(beta**2 + beta**4))*(4*Cpq3**2 + Cpu**2)*vol4)*z**3 + 
     -          2*beta**3*(4 - 7*beta**2 + 3*beta**4 + (2 - 5*beta**2 + 3*beta**4)*(2*Cpq3 + Cpu)*vol2 + 
     -             (-8 + 4*beta**2 + 3*beta**4)*(4*Cpq3**2 + Cpu**2)*vol4)*z**4 + 
     -          2*beta**4*(-2*(-1 + beta**2)*(1 + 2*Cpq3*vol2 + Cpu*vol2) + 
     -             (-4 + beta**2)*(4*Cpq3**2 + Cpu**2)*vol4)*z**5 - 
     -          2*beta**5*((-1 + beta**2)*(1 + 2*Cpq3*vol2 + Cpu*vol2) + 
     -             (2 + beta**2)*(4*Cpq3**2 + Cpu**2)*vol4)*z**6 - 
     -          2*beta**6*(4*Cpq3**2 + Cpu**2)*vol4*z**7) - 
     -       2*MZ**2*(-1 - (2*Cpq3 + Cpu)*vol2 - 2*(4*Cpq3**2 + Cpu**2)*vol4 - 
     -          beta*(1 + 2*Cpq3*vol2 + 16*Cpq3**2*vol4 + Cpu*(vol2 + 4*Cpu*vol4))*z + 
     -          beta**2*(5 + 14*Cpq3*vol2 + 7*Cpu*vol2 + 8*Cpq3*Cpu*vol4 + 
     -             2*(-1 - 2*(2*Cpq3 + Cpu)*vol2 + (-2*Cpq3 + Cpu)**2*vol4)*z**2) + 
     -          beta**3*z*(9 + 26*Cpq3*vol2 + 13*Cpu*vol2 + 16*Cpq3*Cpu*vol4 + 
     -             2*(-3 - 5*(2*Cpq3 + Cpu)*vol2 + 4*(4*Cpq3**2 - 2*Cpq3*Cpu + Cpu**2)*vol4)*z**2) + 
     -          2*beta**8*(1 + 4*Cpq3*vol2 + 2*Cpu*vol2 + (2*Cpq3 + Cpu)**2*vol4)*(2 - 3*z**2 + z**4) + 
     -          beta**5*z*(-13 - 50*Cpq3*vol2 - 25*Cpu*vol2 - 32*Cpq3**2*vol4 - 48*Cpq3*Cpu*vol4 - 
     -             8*Cpu**2*vol4 + 4*(3 + 7*(2*Cpq3 + Cpu)*vol2 + 
     -                (4*Cpq3**2 + 16*Cpq3*Cpu + Cpu**2)*vol4)*z**2 - 
     -             2*(1 + 6*Cpq3*vol2 + 3*Cpu*vol2 + 8*Cpq3*Cpu*vol4)*z**4) + 
     -          beta**7*z*(5 + (2*Cpq3 + Cpu)*(13*vol2 + 8*(2*Cpq3 + Cpu)*vol4) - 6*z**2 - 
     -             6*(2*Cpq3 + Cpu)*(3*vol2 + 2*(2*Cpq3 + Cpu)*vol4)*z**2 + 
     -             2*(1 + 6*Cpq3*vol2 + 3*Cpu*vol2 + 2*(2*Cpq3 + Cpu)**2*vol4)*z**4) + 
     -          beta**6*(-7 - (2*Cpq3 + Cpu)*(9*vol2 + 2*(2*Cpq3 + Cpu)*vol4) + 10*z**2 + 
     -             4*(3*(2*Cpq3 + Cpu)*vol2 + (4*Cpq3**2 + 2*Cpq3*Cpu + Cpu**2)*vol4)*z**2 - 
     -             2*(1 + 2*Cpq3*vol2 + 12*Cpq3**2*vol4 + Cpu*(vol2 + 3*Cpu*vol4))*z**4 + 
     -             2*(4*Cpq3**2 + Cpu**2)*vol4*z**6) + 
     -          beta**4*(-1 - 2*z**2 + Cpu*vol2*(-5 + 4*z**2 - 2*z**4) + 
     -             8*Cpq3**2*vol4*(-1 - z**2 + 3*z**4) + 2*Cpu**2*vol4*(-1 - z**2 + 3*z**4) - 
     -             2*Cpq3*(4*Cpu*vol4*(2 - 3*z**2 + z**4) + vol2*(5 - 4*z**2 + 2*z**4)))))*
     -     xI2(MT**2 - (s*(1 + beta*z))/2._dp,MT**2,MZ**2,musq,ep))/
     -   (8._dp*(-1 + beta*z)*(1 + beta*z)*(1 + beta**2 + 2*beta*z)) - 
     -  ((-1 + beta**2)**2*s**2*(1 + 2*Cpq3*vol2 + Cpu*vol2)*(1 + beta**2 + 2*beta*z)*
     -     xI3(0._dp,MT**2,MT**2 - (s*(1 + beta*z))/2._dp,MT**2,MT**2,MZ**2,musq,ep))/16._dp))/pi


     

      vrt(4) = 
     .     (alpha*genfacvert*gw**2*sigma0*
     .     ((s*(-1 - 2*Cpq3*vol2 + 5*Cpq3**2*vol4 + beta*(1 - 3*Cpq3*(4*vol2 + Cpq3*vol4))*z - 
     -       beta**5*z*(-1 + Cpq3**2*vol4 + 2*Cpq3*(2*vol2 + Cpq3*vol4)*z**2) + 
     -       2*beta**2*(1 + Cpq3*vol2 - Cpq3*(3*vol2 + 5*Cpq3*vol4)*z**2) + 
     -       2*beta**3*z*(-1 - Cpq3**2*vol4*(-4 + z**2) + 2*Cpq3*vol2*(3 + z**2)) + 
     -       beta**4*(-1 + 6*Cpq3*vol2*z**2 + Cpq3**2*vol4*(1 + 2*(z**2 + z**4)))))/(-8 + 8*beta*z) + 
     -  ((-1 + beta**2)*s*(4*MW**2 + (-1 + beta**2)*s)*(1 + 2*Cpq3*vol2 + Cpq3**2*vol4)*
     -     (-1 + beta*(z + beta*(-2 - beta*z + 3*z**2 + beta**2*(2 - 3*z**2 + z**4))))*DB0(MT**2,MB**2,MW**2)
     -     )/(16._dp*(-1 + beta*z)) + ((1 + Cpq3**2*vol4 + 
     -       beta*(-(beta*(-2 + beta**2*(1 + Cpq3**2*vol4))) + 
     -          (1 + 4*beta**2 - 3*beta**4 + (-1 + beta**2)**2*Cpq3**2*vol4)*z + 
     -          2*beta*(-1 + beta**2)*(2 + Cpq3**2*vol4)*z**2 + 6*beta**2*(-1 + beta**2)*z**3 - 
     -          2*beta**3*z**4 - 2*beta**4*z**5))*xI1(MW**2,musq,ep))/
     -   (2._dp*(-1 + beta*z)*(1 + beta**2 + 2*beta*z)) + 
     -  (((1 - beta**2)*((-1 + beta**2)*s*(-3 - Cpq3*(4*vol2 + Cpq3*vol4) + beta*(-1 + Cpq3**2*vol4)*z + 
     -             beta**3*(z - Cpq3**2*vol4*z) + 
     -             beta**2*(3 + 4*Cpq3*vol2 + Cpq3**2*vol4)*(-1 + 2*z**2) + 
     -             2*beta**4*(1 + 2*Cpq3*vol2 + Cpq3**2*vol4)*(2 - 3*z**2 + z**4)) + 
     -          4*MW**2*(1 - Cpq3**2*vol4 + beta*(-1 + Cpq3**2*vol4)*z + beta**3*(z - Cpq3**2*vol4*z) + 
     -             beta**2*(-3 - 4*Cpq3*vol2 - Cpq3**2*vol4 + 2*(1 + 2*Cpq3*vol2 + Cpq3**2*vol4)*z**2) + 
     -             2*beta**4*(1 + 2*Cpq3*vol2 + Cpq3**2*vol4)*(2 - 3*z**2 + z**4))) - 
     -       2*(4*MW**2*(-1 + Cpq3**2*vol4) + (-1 + beta**2)*s*(1 + 4*Cpq3*vol2 + 3*Cpq3**2*vol4))*
     -        (1 + beta*z)*(-1 + beta*(z + beta*(-2 - beta*z + 3*z**2 + beta**2*(2 - 3*z**2 + z**4)))))*
     -     xI2(MT**2,0._dp,MW**2,musq,ep))/(16._dp*(-1 + beta**2*z**2)) + 
     -  ((4*MW**2*(-1 - 3*Cpq3**2*vol4 - beta*(z + 7*Cpq3**2*vol4*z) + 
     -          2*beta**8*(1 + 2*Cpq3*vol2 + Cpq3**2*vol4)*(2 - 3*z**2 + z**4) + 
     -          beta**2*(5 - 2*z**2 - 4*Cpq3*vol2*(-1 + z**2) + Cpq3**2*vol4*(-5 + 6*z**2)) + 
     -          beta**3*z*(9 - 6*z**2 - 8*Cpq3*vol2*(-1 + z**2) + Cpq3**2*vol4*(-9 + 22*z**2)) + 
     -          beta**5*z*(-13 + 12*z**2 - 2*z**4 - 8*Cpq3*vol2*(3 - 4*z**2 + z**4) + 
     -             Cpq3**2*vol4*(-3 - 4*z**2 + 2*z**4)) + 
     -          beta**4*(-1 - 2*z**2 - 4*Cpq3*vol2*(2 - 3*z**2 + z**4) + 
     -             Cpq3**2*vol4*(-3 - 2*z**2 + 12*z**4)) + 
     -          beta**6*(-7 + 10*z**2 - 2*z**4 + 4*Cpq3*vol2*(-1 + z**2) + 
     -             Cpq3**2*vol4*(3 - 2*z**2 - 10*z**4 + 4*z**6)) + 
     -          beta**7*z*(5 - 6*z**2 + 2*z**4 + 8*Cpq3*vol2*(2 - 3*z**2 + z**4) + 
     -             Cpq3**2*vol4*(11 + 6*z**2*(-3 + z**2)))) + 
     -       s*(1 + beta**2 + 2*beta*z)*(-1 + 4*Cpq3*vol2 + Cpq3**2*vol4 + 
     -          beta*(-1 + 8*Cpq3*vol2 + Cpq3**2*vol4)*z + 
     -          beta**2*(5 + 4*Cpq3*vol2 - 5*Cpq3**2*vol4 + 2*(-1 - 4*Cpq3*vol2 + Cpq3**2*vol4)*z**2) + 
     -          beta**3*z*(9 - 17*Cpq3**2*vol4 + 2*(-3 - 8*Cpq3*vol2 + 7*Cpq3**2*vol4)*z**2) + 
     -          2*beta**8*(1 + 2*Cpq3*vol2 + Cpq3**2*vol4)*(2 - 3*z**2 + z**4) + 
     -          beta**5*z*(-13 + 12*z**2 - 2*z**4 - 8*Cpq3*vol2*(3 - 5*z**2 + z**4) + 
     -             Cpq3**2*vol4*(-3 + 4*z**2 + 2*z**4)) + 
     -          beta**4*(-1 - 2*z**2 - 4*Cpq3*vol2*(3 - 3*z**2 + z**4) + 
     -             Cpq3**2*vol4*(-7 - 2*z**2 + 12*z**4)) + 
     -          beta**6*(-7 + 10*z**2 - 2*z**4 + 4*Cpq3*vol2*(-1 + 2*z**2) + 
     -             Cpq3**2*vol4*(3 + 2*z**2 - 10*z**4 + 4*z**6)) + 
     -          beta**7*z*(5 - 6*z**2 + 2*z**4 + 8*Cpq3*vol2*(2 - 3*z**2 + z**4) + 
     -             Cpq3**2*vol4*(11 + 6*z**2*(-3 + z**2)))))*
     -     xI2(MT**2 - (s*(1 + beta*z))/2._dp,0._dp,MW**2,musq,ep))/
     -   (16._dp*(-1 + beta*z)*(1 + beta*z)*(1 + beta**2 + 2*beta*z))))/(mw**2*pi)


      vrt(5) = 
     .     (alpha*genfacvert*gw**2*sigma0*
     .     (0.125_dp*(1._dp - beta**2)**2*s + (0.25_dp*(1._dp 
     .     - beta**2)*s*(-mh**2 + (1._dp - beta**2)*s)*
     .     (1._dp - beta**4*z**4 + beta*(1._dp - beta**2)*
     .     (2._dp*beta - z - 3._dp*beta*z**2))*
     .     db0(mt**2,mt**2,mh**2))/(1._dp - beta*z) 
     .     + (0.5_dp*(1._dp + 2._dp*beta**2 - beta**4 
     .     + beta*(1._dp + 4._dp*beta**2 - 3._dp*beta**4)*z 
     .     - 2._dp*beta**4*z**4*(1._dp + beta*z) 
     .     - 2._dp*(1._dp - beta**2)*(2._dp*beta**2*z**2 
     .     + 3._dp*beta**3*z**3))*(-xI1(mh**2,musq,ep) 
     .     + xI1(mt**2,musq,ep)))/((1._dp - beta*z)*
     .     (1._dp + beta**2 + 2._dp*beta*z)) + (0.125_dp*(-(1._dp 
     .     - beta**2)**2*s*(1._dp + 5._dp*beta**2 - 8._dp*beta**4 
     .     + beta*(1._dp - beta**2)*z - 6._dp*beta**2*(1._dp 
     .     - 2._dp*beta**2)*z**2 - 4._dp*beta**4*z**4) 
     .     + 2._dp*mh**2*(1._dp + 8._dp*beta**2 - 11._dp*beta**4 
     .     + 4._dp*beta**6 + beta*(1._dp + 4._dp*beta**2 - 3._dp*beta**4)*z 
     .     - 2._dp*beta**2*(5._dp - 8._dp*beta**2 + 3._dp*beta**4)*z**2 
     .     - 6._dp*beta**3*(1._dp - beta**2)*z**3 - 2._dp*beta**4*
     .     (2._dp - beta**2)*z**4 - 2._dp*beta**5*z**5))*
     .     xI2(mt**2,mt**2,mh**2,musq,ep))/(1._dp - beta**2*z**2) 
     .     + (0.125_dp*(1._dp - beta**2)*(2._dp*mh**2*(1._dp 
     .     - 4._dp*beta**2 - 3._dp*beta**4 + 4._dp*beta**6 + beta*(1._dp 
     .     - 8._dp*beta**2 + 5._dp*beta**4)*z + 2._dp*beta**2*(1._dp 
     .     + 2._dp*beta**2 - 3._dp*beta**4)*z**2 + 6._dp*beta**3*(1._dp 
     .     - beta**2)*z**3 + 2._dp*beta**6*z**4 + 2._dp*beta**5*z**5) 
     .     + beta*s*(beta*(3._dp - 8._dp*beta**2 - 5._dp*beta**4 
     .     + 8._dp*beta**6) + (1._dp + beta**2 - 23._dp*beta**4 
     .     + 17._dp*beta**6)*z - beta*(1._dp - 11._dp*beta**4 
     .     + 12._dp*beta**6)*z**2 - 2._dp*beta**2*(1._dp - 13._dp*beta**2 
     .     + 12._dp*beta**4)*z**3 + 2._dp*beta**3*(2._dp - 3._dp*beta**2 
     .     + 2._dp*beta**4)*z**4 - 4._dp*beta**4*(1._dp 
     .     - 2._dp*beta**2)*z**5 + 2._dp*beta**5*z**6))*
     .     xI2(t,mt**2,mh**2,musq,ep))/((1._dp + beta**2 + 2._dp*beta*z)*
     .     (1._dp - beta**2*z**2)) 
     .     + 0.0625_dp*(1._dp - beta**2)**2*s**2*(3._dp 
     .     - beta**2 + 2._dp*beta*z)*
     .     xI3(0._dp,mt**2,t,mt**2,mt**2,mh**2,musq,ep)))/(mw**2*pi)

      bx(1) = 
     .     (0.125_dp*alpha*genfacbox*sigma0*
     .     (-2._dp*beta*(gat_smeft2*c(5) + gvt_smeft2*c(4))*z*(1._dp - z**2) + 
     .     (8._dp*beta*(gat_smeft2*c(5) + gvt_smeft2*c(4))*(1._dp - z**2)*
     .     (4._dp*beta - z - 2._dp*beta*z**2)*(xI1(mt**2,musq,ep) 
     .     - xI1(mz**2,musq,ep)))/(s*(1._dp + beta**2 
     .     + 2._dp*beta*z)) + (4._dp*((-(gat_smeft2*c(5) + gvt_smeft2*c(4))*mz**2*
     .     (-6._dp*beta**3 - (5._dp - 8._dp*beta**2 + 4._dp*beta**4)*z 
     .     - beta*(5._dp - 12._dp*beta**2)*z**2 + (3._dp 
     .     - 4._dp*beta**2 + 2._dp*beta**4)*z**3 + beta*(3._dp 
     .     - 4._dp*beta**2)*z**4))/(beta*s) - gat_smeft2*c(5)*(3._dp*beta**2 
     .     - 5._dp*beta**4 - 3._dp*beta*(3._dp - 5._dp*beta**2 
     .     + 2._dp*beta**4)*z - (1._dp + 8._dp*beta**2 
     .     - 9._dp*beta**4)*z**2 + beta*(4._dp - 7._dp*beta**2 
     .     + 3._dp*beta**4)*z**3 + beta**2*(5._dp - 3._dp*beta**2)*z**4) 
     .     + gvt_smeft2*c(4)*(5._dp*beta**2 - 3._dp*beta**4 + beta*(1._dp + beta**2 
     .     - 2._dp*beta**4)*z + (1._dp - 4._dp*beta**2 + 3._dp*beta**4)*z**2 
     .     - beta**3*(1._dp - beta**2)*z**3 - beta**2*
     .     (1._dp + beta**2)*z**4))*xI2(mt**2,mt**2,mz**2,musq,ep))/
     .     (1._dp + beta*z) - 2._dp*((2._dp*(gat_smeft2*c(5) + gvt_smeft2*c(4))*mz**2*z*(5._dp 
     .     - 4._dp*beta**2 - (3._dp - 2._dp*beta**2)*z**2))/(beta*s) 
     .     - gat_smeft2*c(5)*(2._dp - 3._dp*beta*(5._dp - 4._dp*beta**2)*z 
     .     - 2._dp*z**2 + 3._dp*beta*(3._dp - 2._dp*beta**2)*z**3) 
     .     - gvt_smeft2*c(4)*(2._dp + beta*(1._dp - 4._dp*beta**2)*z 
     .     - 2._dp*z**2 + beta*(1._dp + 2._dp*beta**2)*z**3))*
     .     xI2(s,mt**2,mt**2,musq,ep) + (4._dp*((2._dp*beta*(gat_smeft2*c(5) 
     .     + gvt_smeft2*c(4))*mz**2*(beta + z)*(1._dp - 3._dp*beta**2 + beta*(1._dp 
     .     - 2._dp*beta**2)*z + 2._dp*beta**2*z**2 + beta**3*z**3))/s 
     .     - gvt_smeft2*c(4)*(1._dp + 8._dp*beta**2 - 3._dp*beta**6 
     .     + beta*(4._dp + 14._dp*beta**2 - 11._dp*beta**4 
     .     - 2._dp*beta**6)*z - beta**2*(2._dp + beta**2 
     .     + 3._dp*beta**4)*z**2 - beta**3*(11._dp - 6._dp*beta**2 
     .     - beta**4)*z**3 - 2._dp*beta**4*(1._dp 
     .     - beta**2)*z**4 - beta**5*z**5) - gat_smeft2*c(5)*
     .     (1._dp + 5._dp*beta**6 + beta*(4._dp - 10._dp*beta**2 
     .     + 5._dp*beta**4 + 6._dp*beta**6)*z + beta**2*(2._dp 
     .     - 17._dp*beta**2 + 9._dp*beta**4)*z**2 + beta**3*(1._dp 
     .     - 2._dp*beta**2 - 3._dp*beta**4)*z**3 + 6._dp*beta**4*(1._dp 
     .     - beta**2)*z**4 - beta**5*z**5))*
     .     xI2(t,mt**2,mz**2,musq,ep))/((1._dp + beta*z)*(1._dp 
     .     + beta**2 + 2._dp*beta*z)) + 2._dp*((8._dp*(gat_smeft2*c(5) 
     .     + gvt_smeft2*c(4))*mz**4)/s + 4._dp*mz**2*(-gat_smeft2*c(5)*(1._dp 
     .     - 4._dp*beta**2 - beta*z) + gvt_smeft2*c(4)*(3._dp + beta*z)) 
     .     + s*(gat_smeft2*c(5)*(4._dp - 6._dp*beta**2 + 7._dp*beta**4 
     .     - 2._dp*beta*(2._dp - 3._dp*beta**2)*z + beta**2*z**2) 
     .     + gvt_smeft2*c(4)*(4._dp + 2._dp*beta**2 - beta**4 
     .     + 2._dp*beta*(2._dp - beta**2)*z + beta**2*z**2)))*
     .     xI3(0._dp,0._dp,s,mt**2,mt**2,mt**2,musq,ep) - 2._dp*((8._dp* 
     .     (gat_smeft2*c(5) + gvt_smeft2*c(4))*mz**4*(1._dp + beta*z))/s + 4._dp*mz**2*(1._dp 
     .     + beta*z)*(-gat_smeft2*c(5) + 4._dp*beta**2*gat_smeft2*c(5) + 3._dp*gvt_smeft2*c(4) 
     .     + beta*(gat_smeft2*c(5) + gvt_smeft2*c(4))*z) + s*(gat_smeft2*c(5)*(4._dp - 6._dp*beta**2 
     .     + 7._dp*beta**4 - beta*(1._dp - 2._dp*beta**2 
     .     - 6._dp*beta**4)*z - 3._dp*beta**2*(1._dp - 2._dp*beta**2)*z**2 
     .     + beta**3*z**3) + gvt_smeft2*c(4)*(4._dp + 2._dp*beta**2 - beta**4 
     .     + beta*(7._dp + 2._dp*beta**2 - 2._dp*beta**4)*z + beta**2*(5._dp 
     .     - 2._dp*beta**2)*z**2 + beta**3*z**3)))*xI3(0._dp,mt**2,t,mt**2
     .     ,mt**2,mz**2,musq,ep) + 2._dp*((-2._dp*(gat_smeft2*c(5) + gvt_smeft2*c(4))*
     .     mz**4*z*(5._dp - 8._dp*beta**2 - 3._dp*z**2 
     .     + 2._dp*beta**2*z**2))/(beta*s) + beta*s*(gvt_smeft2*c(4)*(beta*(4._dp 
     .     - beta**2) + 4._dp*(1._dp + beta**2 - beta**4)*z 
     .     - beta**3*z**2 + beta**2*z**3 + beta**4*z**3) 
     .     - gat_smeft2*c(5)*(-3._dp*beta**3 - 4._dp*(1._dp - 3._dp*beta**2 
     .     + 3._dp*beta**4)*z + beta*(4._dp - 3._dp*beta**2)*z**2 
     .     - beta**2*(5._dp - 3._dp*beta**2)*z**3)) + 2._dp*mz**2*
     .     (gvt_smeft2*c(4)*(1._dp + beta**2 + 4._dp*beta*z - (1._dp 
     .     - beta**2)*z**2 + 2._dp*beta*z**3) + gat_smeft2*c(5)*(1._dp 
     .     + beta**2 - 4._dp*beta*(3._dp - 4._dp*beta**2)*z - (1._dp 
     .     - beta**2)*z**2 + 2._dp*beta*(3._dp - 2._dp*beta**2)*z**3))
     .     )*xI3(s,mt**2,mt**2,mt**2,mt**2,mz**2,musq,ep) + ((16._dp*
     .     (gat_smeft2*c(5) + gvt_smeft2*c(4))*mz**6)/s + 8._dp*mz**4*(-gat_smeft2*c(5)*(1._dp 
     .     - 5._dp*beta**2 - 2._dp*beta*z) + gvt_smeft2*c(4)*(3._dp + beta**2 
     .     + 2._dp*beta*z)) + 2._dp*mz**2*s*(gvt_smeft2*c(4)*(4._dp + 9._dp*beta**2 
     .     - 2._dp*beta**4 + 10._dp*beta*z + beta**2*(2._dp + beta**2)*z**2
     .     ) + gat_smeft2*c(5)*(4._dp - 7._dp*beta**2 + 14._dp*beta**4 - 2._dp*beta*
     .     (3._dp - 8._dp*beta**2)*z + beta**2*(2._dp + beta**2)*z**2)) 
     .     + s**2*(-gat_smeft2*c(5)*(3._dp - 7._dp*beta**2 + 5._dp*beta**4 
     .     - 6._dp*beta**6 - beta*(3._dp - 8._dp*beta**2 
     .     + 12._dp*beta**4)*z + 3._dp*beta**2*(1._dp - beta**2 
     .     - beta**4)*z**2 - beta**3*z**3) + gvt_smeft2*c(4)*(1._dp 
     .     + 7._dp*beta**2 - beta**4 - 2._dp*beta**6 + beta*(3._dp 
     .     + 8._dp*beta**2 - 4._dp*beta**4)*z + beta**2*(1._dp 
     .     + 3._dp*beta**2 - beta**4)*z**2 + beta**3*z**3)))*
     .     xI4(0._dp,0._dp,mt**2,mt**2,s,t,mt**2,mt**2,mt**2,mz**2,
     .     musq,ep)))/pi


      bx(2) = 
     .     (0.5_dp*alpha*genfacbox*gw_smeft2*c(6)*sigma0*
     .     (-beta*z*(1._dp - z**2) + (4._dp*beta*(1._dp 
     .     - z**2)*(4._dp*beta - z - 2._dp*beta*z**2)*
     .     (xI1(mb**2,musq,ep) - xI1(mw**2,musq,ep)))/(s*(1._dp 
     .     + beta**2 + 2._dp*beta*z)) + (0.5_dp*(4._dp*(-mb**2 
     .     + mw**2)*(6._dp*beta**3 + (5._dp - 8._dp*beta**2 + 4._dp*beta**4
     .     )*z + beta*(5._dp - 12._dp*beta**2)*z**2 - (3._dp 
     .     - 4._dp*beta**2 + 2._dp*beta**4)*z**3 - beta*(3._dp 
     .     - 4._dp*beta**2)*z**4) + s*(2._dp*beta**3*(5._dp - beta**2) 
     .     - beta*(3._dp + 5._dp*beta**2)*z**4 + (1._dp 
     .     - beta**2)*((5._dp + 12._dp*beta**2 - 4._dp*beta**4)*z 
     .     + 9._dp*beta*z**2 - (3._dp + 4._dp*beta**2 
     .     - 2._dp*beta**4)*z**3)))*xI2(mt**2,mb**2,mw**2,musq,ep))/
     .     (beta*s*(1._dp + beta*z)) - (0.5_dp*(4._dp*(-mb**2 
     .     + mw**2)*z*(5._dp - 4._dp*beta**2 - (3._dp - 2._dp*beta**2
     .     )*z**2) - s*(4._dp*beta - (5._dp + 5._dp*beta**2 
     .     - 4._dp*beta**4)*z - 4._dp*beta*z**2 + (3._dp + 5._dp*beta**2 
     .     - 2._dp*beta**4)*z**3))*xI2(s,mb**2,mb**2,musq,ep))/(beta*s) + 
     .     ((4._dp*beta*(-mb**2 + mw**2)*(beta + z)*(1._dp 
     .     - 3._dp*beta**2 + beta*(1._dp - 2._dp*beta**2)*z 
     .     + 2._dp*beta**2*z**2 + beta**3*z**3) - s*(1._dp + beta**2 
     .     + 2._dp*beta*z)*(2._dp + 5._dp*beta**2 - beta**4 + beta*
     .     (3._dp - 6._dp*beta**2 + 2._dp*beta**4)*z - beta**2*(7._dp 
     .     - 2._dp*beta**2)*z**2 + beta**3*(2._dp - beta**2)*z**3 
     .     - beta**4*z**4))*xI2(t,mb**2,mw**2,musq,ep))/(s*(1._dp 
     .     + beta*z)*(1._dp + beta**2 + 2._dp*beta*z)) + (0.5_dp*(16._dp*
     .     (-mb**2 + mw**2)**2 - 8._dp*beta*mb**2*s*(2._dp*beta + z) 
     .     + 8._dp*mw**2*s*(2._dp + beta**2 + beta*z) + s**2*(7._dp 
     .     + 2._dp*beta**2 + beta**4 + 2._dp*beta*(1._dp + beta**2)*z 
     .     + 2._dp*beta**2*z**2))*
     .     xI3(0._dp,0._dp,s,mb**2,mb**2,mb**2,musq,ep))
     .     /s - (0.5_dp*(16._dp*(-mb**2 + mw**2)**2*(1._dp + beta*z) 
     .     - 8._dp*beta*mb**2*s*(beta + z)*(2._dp + beta*z) + 8._dp*mw**2*
     .     s*(1._dp + beta*z)*(2._dp + beta**2 + beta*z) + s**2*(1._dp 
     .     + beta*z)*(7._dp + 2._dp*beta**2 + beta**4 + 2._dp*beta*(1._dp 
     .     + beta**2)*z + 2._dp*beta**2*z**2))*xI3(0._dp,mt**2,t,mb**2
     .     ,mb**2,mw**2,musq,ep))/s + (0.125_dp*(-16._dp*(-mb**2 
     .     + mw**2)**2*z*(5._dp - 8._dp*beta**2 - (3._dp 
     .     - 2._dp*beta**2)*z**2) - 8._dp*mb**2*s*(2._dp*beta*(1._dp 
     .     + beta**2) - (5._dp - 8._dp*beta**2)*(1._dp + beta**2)*z 
     .     - 2._dp*beta*(1._dp - beta**2)*z**2 + (3._dp 
     .     - 2._dp*beta**2)*(1._dp + beta**2)*z**3) + 8._dp*mw**2*s*
     .     (2._dp*beta*(1._dp + beta**2) - (5._dp - 5._dp*beta**2 
     .     - 8._dp*beta**4)*z - 2._dp*beta*(1._dp - beta**2)*z**2 + 
     .     (3._dp + 3._dp*beta**2 - 2._dp*beta**4)*z**3) - s**2*
     .     (-4._dp*beta*(1._dp + 4._dp*beta**2 + beta**4) + (5._dp 
     .     - 30._dp*beta**2 + beta**4 - 8._dp*beta**6)*z + 4._dp*beta*
     .     (1._dp + 2._dp*beta**2 - beta**4)*z**2 - (3._dp 
     .     + 4._dp*beta**2 + 11._dp*beta**4 - 2._dp*beta**6)*z**3))*
     .     xI3(s,mt**2,mt**2,mb**2,mb**2,mw**2,musq,ep))/(beta*s) 
     .     + (0.125_dp*(64._dp*(-mb**2 + mw**2)**3 + s**3*(1._dp 
     .     + beta**2 + 2._dp*beta*z)*(7._dp + 2._dp*beta**2 + beta**4 
     .     + 2._dp*beta*(1._dp + beta**2)*z + 2._dp*beta**2*z**2) 
     .     + 8._dp*s*(10._dp*mw**4 + 2._dp*(-mb**2 + mw**2)**2*
     .     (3._dp*beta**2 + 4._dp*beta*z) - 4._dp*mb**2*mw**2*(3._dp 
     .     + beta**2*z**2) + 2._dp*mb**4*(1._dp + 2._dp*beta**2*z**2)) 
     .     + 4._dp*s**2*(mw**2*(11._dp + 8._dp*beta**2 + 3._dp*beta**4 
     .     + 4._dp*beta*(3._dp + 2._dp*beta**2)*z + 6._dp*beta**2*z**2) 
     .     - mb**2*(11._dp - 4._dp*beta**2 + 3._dp*beta**4 
     .     + 8._dp*beta*(1._dp + beta**2)*z + 2._dp*beta**2*
     .     (6._dp + beta**2)*z**2)))*
     .     xI4(0._dp,0._dp,mt**2,mt**2,s,t,mb**2,mb**2,mb**2,mw**2,
     .     musq,ep))/s))/pi


      bx(3) = 
     .     (0.5_dp*alpha*chifac*genfacbox*mt**2*sigma0*
     .     (
     .     (3*(4*Cpq3**2 + Cpu**2)*vol4 + beta*z*
     -      (1 - 6*(2*Cpq3 + Cpu)*vol2 - 2*(4*Cpq3**2 + 2*Cpq3*Cpu + Cpu**2)*vol4 - 
     -        (1 + 4*Cpq3*vol2 + 2*Cpu*vol2 + (2*Cpq3 + Cpu)**2*vol4)*z**2) + 
     -     beta**5*z*(-2*(2*Cpq3 + Cpu)*(vol2 + (2*Cpq3 + Cpu)*vol4) - 
     -        (4*Cpq3*vol2 + 4*Cpq3**2*vol4 + Cpu*(2*vol2 + Cpu*vol4))*z**2) + 
     -     beta**2*(-2*(2*Cpq3 + Cpu)*(vol2 + (2*Cpq3 + Cpu)*vol4) + 
     -        (1 - 3*(-2*Cpq3 + Cpu)**2*vol4)*z**2 - 
     -        (1 + 4*Cpq3*vol2 + 2*Cpu*vol2 + (2*Cpq3 + Cpu)**2*vol4)*z**4) + 
     -     beta**4*(2*(2*Cpq3 + Cpu)*(vol2 + (2*Cpq3 + Cpu)*vol4) - 
     -        (1 + (4*Cpq3**2 + 12*Cpq3*Cpu + Cpu**2)*vol4)*z**2 + 
     -        (1 + 8*Cpq3**2*vol4 + 4*Cpq3*(vol2 + Cpu*vol4) + 2*Cpu*(vol2 + Cpu*vol4))*z**4) + 
     -     beta**3*z*(-1 + 24*Cpq3**2*vol4 + 6*Cpu**2*vol4 + z**2 + 4*Cpu*vol2*(2 + z**2) + 
     -        4*Cpq3*(2*vol2*(2 + z**2) + Cpu*vol4*(3 + z**2))))/((-1 + beta**2)*(1 + beta*z)) + 
     -  (4*(-4*beta**2 + 4*beta**4 - 12*beta**2*Cpq3*vol2 + 16*beta**4*Cpq3*vol2 - 6*beta**2*Cpu*vol2 + 
     -       8*beta**4*Cpu*vol2 - 4*Cpq3**2*vol4 - 12*beta**2*Cpq3**2*vol4 + 20*beta**4*Cpq3**2*vol4 - 
     -       8*beta**2*Cpq3*Cpu*vol4 + 16*beta**4*Cpq3*Cpu*vol4 - Cpu**2*vol4 - 3*beta**2*Cpu**2*vol4 + 
     -       5*beta**4*Cpu**2*vol4 + beta*(1 + 2*Cpu*vol2 + 4*Cpq3*(vol2 + Cpu*vol4) - 
     -          beta**2*(5 + 8*Cpq3*vol2 + 4*Cpu*vol2 - 4*Cpq3*Cpu*vol4) + 
     -          beta**4*(4 + 8*Cpq3*vol2 + 4*Cpu*vol2 + 4*Cpq3**2*vol4 + Cpu**2*vol4))*z - 
     -       beta**2*(-1 + beta**2)*(7 + 24*Cpq3*vol2 + 12*Cpu*vol2 + 
     -          4*(8*Cpq3**2 + 5*Cpq3*Cpu + 2*Cpu**2)*vol4)*z**2 - 
     -       beta*(-1 + beta**2)*(-1 + 3*beta**2*
     -           (2 + 4*Cpq3*vol2 + 2*Cpu*vol2 + 4*Cpq3**2*vol4 + Cpu**2*vol4) - 
     -          (2*Cpq3 + Cpu)*(2*vol2 + (2*Cpq3 + Cpu)*vol4))*z**3 + 
     -       beta**2*(-3 - 3*(2*Cpq3 + Cpu)*(2*vol2 + (2*Cpq3 + Cpu)*vol4) + 
     -          beta**2*(3 + 8*Cpq3*vol2 + 4*Cpu*vol2 + 2*(4*Cpq3**2 + 2*Cpq3*Cpu + Cpu**2)*vol4))*z**4
     -        + beta**3*(-2 + beta**2*(2 + 4*Cpq3*vol2 + 2*Cpu*vol2 + 4*Cpq3**2*vol4 + Cpu**2*vol4) - 
     -          2*(2*Cpq3 + Cpu)*(2*vol2 + (2*Cpq3 + Cpu)*vol4))*z**5)*xI1(MT**2,musq,ep))/
     -   ((-1 + beta**2)*s*(1 + beta*z)*(1 + beta**2 + 2*beta*z)) + 
     -  (4*(4*beta**2 - 4*beta**4 + 12*beta**2*Cpq3*vol2 - 16*beta**4*Cpq3*vol2 + 6*beta**2*Cpu*vol2 - 
     -       8*beta**4*Cpu*vol2 + 4*Cpq3**2*vol4 + 12*beta**2*Cpq3**2*vol4 - 20*beta**4*Cpq3**2*vol4 + 
     -       8*beta**2*Cpq3*Cpu*vol4 - 16*beta**4*Cpq3*Cpu*vol4 + Cpu**2*vol4 + 3*beta**2*Cpu**2*vol4 - 
     -       5*beta**4*Cpu**2*vol4 - beta*(1 + 2*Cpu*vol2 + 4*Cpq3*(vol2 + Cpu*vol4) - 
     -          beta**2*(5 + 8*Cpq3*vol2 + 4*Cpu*vol2 - 4*Cpq3*Cpu*vol4) + 
     -          beta**4*(4 + 8*Cpq3*vol2 + 4*Cpu*vol2 + 4*Cpq3**2*vol4 + Cpu**2*vol4))*z + 
     -       beta**2*(-1 + beta**2)*(7 + 24*Cpq3*vol2 + 12*Cpu*vol2 + 
     -          4*(8*Cpq3**2 + 5*Cpq3*Cpu + 2*Cpu**2)*vol4)*z**2 + 
     -       beta*(-1 + beta**2)*(-1 + 3*beta**2*
     -           (2 + 4*Cpq3*vol2 + 2*Cpu*vol2 + 4*Cpq3**2*vol4 + Cpu**2*vol4) - 
     -          (2*Cpq3 + Cpu)*(2*vol2 + (2*Cpq3 + Cpu)*vol4))*z**3 + 
     -       beta**2*(3 + 12*Cpq3*vol2 + 6*Cpu*vol2 + 3*(2*Cpq3 + Cpu)**2*vol4 - 
     -          beta**2*(3 + 8*Cpq3*vol2 + 4*Cpu*vol2 + 2*(4*Cpq3**2 + 2*Cpq3*Cpu + Cpu**2)*vol4))*z**4
     -        + beta**3*(2 + 8*Cpq3*vol2 + 4*Cpu*vol2 + 2*(2*Cpq3 + Cpu)**2*vol4 - 
     -          beta**2*(2 + 4*Cpq3*vol2 + 2*Cpu*vol2 + 4*Cpq3**2*vol4 + Cpu**2*vol4))*z**5)*
     -     xI1(MZ**2,musq,ep))/((-1 + beta**2)*s*(1 + beta*z)*(1 + beta**2 + 2*beta*z)) + 
     -  (2*((MZ**2*(-2*beta**7*(1 + 4*Cpq3*vol2 + 2*Cpu*vol2 + (2*Cpq3 + Cpu)**2*vol4)*z**2*
     -             (-2 + z**2) + (1 + 4*Cpq3*vol2 + 2*Cpu*vol2 + (2*Cpq3 + Cpu)**2*vol4)*z*
     -             (-5 + 3*z**2) + 2*beta*(1 + 4*Cpq3*vol2 + 2*Cpu*vol2 + (2*Cpq3 + Cpu)**2*vol4)*z**2*
     -             (-5 + 3*z**2) + 2*beta**3*
     -             (-3 + (2*Cpq3 + Cpu)*(-2*vol2 + (2*Cpq3 + Cpu)*vol4) + 15*z**2 + 
     -               (2*Cpq3 + Cpu)*(26*vol2 + 11*(2*Cpq3 + Cpu)*vol4)*z**2 - 
     -               7*(1 + 4*Cpq3*vol2 + 2*Cpu*vol2 + (2*Cpq3 + Cpu)**2*vol4)*z**4) + 
     -            beta**4*z*(-18 - 2*(2*Cpq3 + Cpu)*(17*vol2 + 8*(2*Cpq3 + Cpu)*vol4) + 
     -               (23 + (2*Cpq3 + Cpu)*(40*vol2 + 17*(2*Cpq3 + Cpu)*vol4))*z**2 - 
     -               7*(1 + 4*Cpq3*vol2 + 2*Cpu*vol2 + (2*Cpq3 + Cpu)**2*vol4)*z**4) + 
     -            beta**2*z*(13 + (2*Cpq3 + Cpu)*(28*vol2 + 15*(2*Cpq3 + Cpu)*vol4) - 12*z**2 - 
     -               12*(2*Cpq3 + Cpu)*(2*vol2 + (2*Cpq3 + Cpu)*vol4)*z**2 + 
     -               3*(1 + 4*Cpq3*vol2 + 2*Cpu*vol2 + (2*Cpq3 + Cpu)**2*vol4)*z**4) + 
     -            2*beta**6*z*(5 + (2*Cpq3 + Cpu)*(9*vol2 + 4*(2*Cpq3 + Cpu)*vol4) - 7*z**2 - 
     -               (2*Cpq3 + Cpu)*(11*vol2 + 4*(2*Cpq3 + Cpu)*vol4)*z**2 + 
     -               (2 + 6*Cpq3*vol2 + 3*Cpu*vol2 + (2*Cpq3 + Cpu)**2*vol4)*z**4) + 
     -            2*beta**5*(3 + 6*Cpq3*vol2 + 3*Cpu*vol2 - 
     -               4*(3 + (2*Cpq3 + Cpu)*(5*vol2 + 2*(2*Cpq3 + Cpu)*vol4))*z**2 + 
     -               (5 + (2*Cpq3 + Cpu)*(9*vol2 + 4*(2*Cpq3 + Cpu)*vol4))*z**4)))/beta - 
     -       (-1 + beta**2)**2*s*(-2 + beta**3*(1 + 2*Cpq3*vol2 + Cpu*vol2)*z + z**2 + 
     -          beta**2*(1 + (2*Cpq3 + Cpu)*(3*vol2 + 2*(2*Cpq3 + Cpu)*vol4) - z**2 - 
     -             (2*Cpq3 + Cpu)*(4*vol2 + 3*(2*Cpq3 + Cpu)*vol4)*z**2 + 
     -             (1 + 4*Cpq3*vol2 + 2*Cpu*vol2 + (2*Cpq3 + Cpu)**2*vol4)*z**4) + 
     -          (2*Cpq3 + Cpu)*((2*Cpq3 + Cpu)*vol4*(-1 + z**2) + vol2*(-3 + 2*z**2)) + 
     -          beta*z*(-3 + 2*z**2 + (2*Cpq3 + Cpu)*
     -              (2*(2*Cpq3 + Cpu)*vol4*(-1 + z**2) + vol2*(-5 + 4*z**2)))))*
     -     xI2(MT**2,MT**2,MZ**2,musq,ep))/((-1 + beta**2)*s*(1 + beta*z)**2) - 
     -  ((1 + 4*Cpq3*vol2 + 2*Cpu*vol2 + (2*Cpq3 + Cpu)**2*vol4)*
     -     (2*MZ**2*z*(5 - 3*z**2) + 2*beta*s*(-1 + z**2) + 
     -       beta**2*z*(4*MZ**2*(-2 + z**2) + s*(-1 + z**2)))*xI2(s,MT**2,MT**2,musq,ep))/(beta*s) + 
     -  (2*(-(s*(1 + 2*Cpq3*vol2 + Cpu*vol2)) - 2*(4*Cpq3**2 + Cpu**2)*MZ**2*vol4 + 
     -       beta**9*s*(1 + 2*Cpq3*vol2 + Cpu*vol2)*z - 
     -       beta*(s*(3 + 4*Cpq3*vol2 + 2*Cpu*vol2) + 
     -          2*MZ**2*(1 + 6*Cpq3*vol2 + 3*Cpu*vol2 + 4*(4*Cpq3**2 + 2*Cpq3*Cpu + Cpu**2)*vol4))*z + 
     -       beta**3*z*(2*MZ**2*(2 + 4*Cpq3*vol2 - 8*Cpq3**2*vol4 + 2*Cpu*(vol2 - Cpu*vol4) + 
     -             3*(-1 - (2*Cpq3 + Cpu)*vol2 + 2*(4*Cpq3**2 + Cpu**2)*vol4)*z**2) + 
     -          s*(17 + 44*Cpq3*vol2 + 22*Cpu*vol2 + 8*Cpq3**2*vol4 + 36*Cpq3*Cpu*vol4 + 2*Cpu**2*vol4 - 
     -             (9 + 13*(2*Cpq3 + Cpu)*vol2 + 2*(4*Cpq3**2 + 18*Cpq3*Cpu + Cpu**2)*vol4)*z**2)) + 
     -       beta**4*(-(s*(4 + 18*Cpq3*vol2 + 9*Cpu*vol2 + 24*Cpq3**2*vol4 + 20*Cpq3*Cpu*vol4 + 
     -               6*Cpu**2*vol4)) + 4*s*
     -           (5 + 6*(2*Cpq3 + Cpu)*vol2 + (4*Cpq3**2 + 11*Cpq3*Cpu + Cpu**2)*vol4)*z**2 + 
     -          2*s*(-5 - 7*(2*Cpq3 + Cpu)*vol2 + (4*Cpq3**2 - 12*Cpq3*Cpu + Cpu**2)*vol4)*z**4 + 
     -          2*MZ**2*(4 + 14*Cpq3*vol2 + 7*Cpu*vol2 + 16*Cpq3**2*vol4 + 12*Cpq3*Cpu*vol4 + 
     -             4*Cpu**2*vol4 + 2*(2 + 7*(2*Cpq3 + Cpu)*vol2 + (6*Cpq3 + Cpu)*(2*Cpq3 + 3*Cpu)*vol4)*
     -              z**2 - (3 + 7*(2*Cpq3 + Cpu)*vol2 + 16*Cpq3*Cpu*vol4)*z**4)) + 
     -       beta**5*z*(2*MZ**2*(2*(2 + 5*(2*Cpq3 + Cpu)*vol2 + 
     -                4*(4*Cpq3**2 + 3*Cpq3*Cpu + Cpu**2)*vol4) + 
     -             2*(1 - 8*Cpq3**2*vol4 + 2*Cpu*(vol2 - Cpu*vol4) + 4*Cpq3*(vol2 + Cpu*vol4))*z**2 - 
     -             (1 + 6*Cpq3*vol2 + 3*Cpu*vol2 + 8*Cpq3*Cpu*vol4)*z**4) + 
     -          s*(-17 - (2*Cpq3 + Cpu)*(29*vol2 + 13*(2*Cpq3 + Cpu)*vol4) + 
     -             2*(8 + 26*Cpq3*vol2 + 13*Cpu*vol2 + 24*Cpq3**2*vol4 + 32*Cpq3*Cpu*vol4 + 
     -                6*Cpu**2*vol4)*z**2 + 
     -             (-5 - 8*(2*Cpq3 + Cpu)*vol2 + (4*Cpq3**2 - 12*Cpq3*Cpu + Cpu**2)*vol4)*z**4)) + 
     -       beta**7*z*(2*MZ**2*(-5 - (2*Cpq3 + Cpu)*(9*vol2 + 4*(2*Cpq3 + Cpu)*vol4) + z**2 - 
     -             (2*Cpq3 + Cpu)*(vol2 + 2*(2*Cpq3 + Cpu)*vol4)*z**2 + 
     -             (1 + 6*Cpq3*vol2 + 3*Cpu*vol2 + 2*(2*Cpq3 + Cpu)**2*vol4)*z**4) + 
     -          s*(2*(1 + 16*Cpq3**2*vol4 + 8*Cpq3*(vol2 + Cpu*vol4) + 4*Cpu*(vol2 + Cpu*vol4)) - 
     -             (7 + 13*(2*Cpq3 + Cpu)*vol2 + (44*Cpq3**2 + 28*Cpq3*Cpu + 11*Cpu**2)*vol4)*z**2 + 
     -             (5 + 8*(2*Cpq3 + Cpu)*vol2 + 2*(4*Cpq3**2 + 6*Cpq3*Cpu + Cpu**2)*vol4)*z**4 + 
     -             (4*Cpq3**2 + Cpu**2)*vol4*z**6)) + 
     -       beta**6*(2*MZ**2*(-3 - 3*(2*Cpq3 + Cpu)*vol2 - 
     -             (20*Cpq3**2*vol4 + 12*Cpq3*(vol2 + 2*Cpu*vol4) + Cpu*(6*vol2 + 5*Cpu*vol4))*z**2 + 
     -             (2 + 5*(2*Cpq3 + Cpu)*vol2 + 12*Cpq3*Cpu*vol4)*z**4 + (4*Cpq3**2 + Cpu**2)*vol4*z**6)
     -           + s*(-1 - 2*Cpq3*vol2 - Cpu*vol2 + 4*Cpq3**2*vol4 - 4*Cpq3*Cpu*vol4 + Cpu**2*vol4 - 
     -             2*(9 + 2*(2*Cpq3 + Cpu)*(5*vol2 + (2*Cpq3 + Cpu)*vol4))*z**2 + 
     -             (12 + 15*(2*Cpq3 + Cpu)*vol2 + (4*Cpq3**2 + 20*Cpq3*Cpu + Cpu**2)*vol4)*z**4 + 
     -             (-1 - (2*Cpq3 + Cpu)*vol2 + 2*(4*Cpq3**2 + Cpu**2)*vol4)*z**6)) + 
     -       beta**8*(2*MZ**2*(1 + 4*Cpq3*vol2 + 2*Cpu*vol2 + (2*Cpq3 + Cpu)**2*vol4)*z**2*(-2 + z**2) + 
     -          s*(1 + (2*Cpq3 + Cpu)*(3*vol2 + 2*(2*Cpq3 + Cpu)*vol4) + 3*z**2 + 
     -             ((2*Cpq3 + Cpu)*vol2 - (4*Cpq3**2 + 12*Cpq3*Cpu + Cpu**2)*vol4)*z**2 - 
     -             (2 + 8*Cpq3**2*vol4 + 2*Cpq3*(vol2 - 2*Cpu*vol4) + Cpu*(vol2 + 2*Cpu*vol4))*z**4 + 
     -             (1 + 2*Cpq3*vol2 + 4*Cpq3**2*vol4 + Cpu*(vol2 + Cpu*vol4))*z**6)) - 
     -       beta**2*(s*(-5 - 8*(2*Cpq3 + Cpu)*vol2 - 2*(4*Cpq3**2 + 8*Cpq3*Cpu + Cpu**2)*vol4 + 
     -             (5 + 5*(2*Cpq3 + Cpu)*vol2 + 2*(4*Cpq3**2 + 8*Cpq3*Cpu + Cpu**2)*vol4)*z**2) + 
     -          2*MZ**2*(1 + 16*Cpq3**2*vol4 + 2*z**2 + 4*Cpu*(vol2 + Cpu*vol4 + vol2*z**2) + 
     -             4*Cpq3*(2*vol2 + 3*Cpu*vol4 + 2*(vol2 + Cpu*vol4)*z**2))))*
     -     xI2(MT**2 - (s*(1 + beta*z))/2._dp,MT**2,MZ**2,musq,ep))/
     -   ((-1 + beta**2)*s*(1 + beta*z)**2*(1 + beta**2 + 2*beta*z)) + 
     -  ((1 + 4*Cpq3*vol2 + 2*Cpu*vol2 + (2*Cpq3 + Cpu)**2*vol4)*
     -     (8*MZ**4 + 4*beta*MZ**2*s*(beta + z) + s**2*(2 + beta**4 + 2*beta*z + beta**2*(-2 + z**2)))*
     -     xI3(0._dp,0._dp,s,MT**2,MT**2,MT**2,musq,ep))/s - 
     -  ((-1 + beta**2*z**2)*(8*MZ**4*(1 + 4*Cpq3*vol2 + 2*Cpu*vol2 + (2*Cpq3 + Cpu)**2*vol4)*
     -        (1 + beta*z)**2 + 4*beta*s*(1 + 4*Cpq3*vol2 + 2*Cpu*vol2 + (2*Cpq3 + Cpu)**2*vol4)*
     -        (beta + z)*(MZ + beta*MZ*z)**2 + 
     -       s**2*(2 + 6*Cpq3*vol2 + 3*Cpu*vol2 + (2*Cpq3 + Cpu)**2*vol4 + 
     -          beta**5*(1 + 2*Cpq3*vol2 + Cpu*vol2)*z + 
     -          beta*(5 + (2*Cpq3 + Cpu)*(9*vol2 + 4*(2*Cpq3 + Cpu)*vol4))*z + 
     -          beta**4*(1 + (2*Cpq3 + Cpu)*(3*vol2 + 2*(2*Cpq3 + Cpu)*vol4) - 
     -             2*(2*Cpq3 + Cpu)*(vol2 + (2*Cpq3 + Cpu)*vol4)*z**2 + 
     -             (1 + 4*Cpq3*vol2 + 2*Cpu*vol2 + (2*Cpq3 + Cpu)**2*vol4)*z**4) + 
     -          2*beta**3*z*(-1 + 2*z**2 + 
     -             (2*Cpq3 + Cpu)*(-vol2 + 2*(2*vol2 + (2*Cpq3 + Cpu)*vol4)*z**2)) + 
     -          2*beta**2*(-1 + 3*z**2 + (2*Cpq3 + Cpu)*
     -              ((2*Cpq3 + Cpu)*vol4*(-1 + 4*z**2) + vol2*(-2 + 7*z**2)))))*
     -     xI3(0._dp,MT**2,MT**2 - (s*(1 + beta*z))/2._dp,MT**2,MT**2,MZ**2,musq,ep))/
     -   (s*(-1 + beta*z)*(1 + beta*z)**2) + 
     -  ((1 + 4*Cpq3*vol2 + 2*Cpu*vol2 + (2*Cpq3 + Cpu)**2*vol4)*
     -     (-2*beta*MZ**2*s*(-1 + z**2) + beta**3*s*(2*MZ**2 + s)*(1 + z**2) + 2*MZ**4*z*(-5 + 3*z**2) + 
     -       beta**4*s*z*(s - 2*MZ**2*(-4 + z**2)) + 
     -       beta**2*z*(s**2 - 4*MZ**4*(-4 + z**2) + 2*MZ**2*s*(-2 + z**2)))*
     -     xI3(MT**2,MT**2,s,MT**2,MZ**2,MT**2,musq,ep))/(beta*s) + 
     -  ((1 + 4*Cpq3*vol2 + 2*Cpu*vol2 + (2*Cpq3 + Cpu)**2*vol4)*
     -     (16*MZ**6 + 16*beta*MZ**4*s*(beta + z) + beta*s**3*(beta + z)*(1 + beta*z)**2 + 
     -       2*MZ**2*s**2*(2 + beta*(2*z + beta*(-1 + 4*beta*z + 2*z**2 + beta**2*(2 + z**2)))))*
     -     xI4(0._dp,0._dp,MT**2,MT**2,s,MT**2 - (s*(1 + beta*z))/2._dp,MT**2,MT**2,MT**2,MZ**2,musq,ep))/(2._dp*s)
     .     ))/pi


   

      bx(4) = 
     .     (0.25_dp*alpha*genfacbox*gw**2*sigma0*
     .     ((s*(-6*Cpq3**2*vol4 + 2*beta**5*Cpq3*z*(2*vol2*(1 + z**2) + Cpq3*vol4*(2 + z**2)) + 
     -       beta**3*z*(1 - z**2 + Cpq3**2*vol4*(-13 + z**2) - 6*Cpq3*vol2*(3 + z**2)) + 
     -       beta*z*(-1 + z**2 + Cpq3**2*vol4*(5 + z**2) + 2*Cpq3*vol2*(7 + z**2)) + 
     -       beta**4*(z**2 - z**4 + Cpq3**2*vol4*(-4 + z**2 - 3*z**4) - 2*Cpq3*vol2*(2 + z**2 + z**4)) + 
     -       beta**2*(-z**2 + z**4 + 2*Cpq3*vol2*(2 + z**2 + z**4) + Cpq3**2*vol4*(4 + 7*z**2 + z**4))))/
     -   (4 + 4*beta*z) + ((-2*Cpq3**2*vol4 - 2*beta**2*(2 + 2*Cpq3*vol2 + Cpq3**2*vol4) + 
     -       beta**4*(4 + 8*Cpq3*vol2 + 6*Cpq3**2*vol4) + 
     -       beta*(1 + 2*Cpq3*vol2 + beta**2*(-5 + 4*beta**2 + 2*Cpq3*vol2) + 
     -          (-1 + 5*beta**2 - 2*beta**4)*Cpq3**2*vol4)*z - 
     -       beta**2*(-1 + beta**2)*(7 + 10*Cpq3*vol2 + 9*Cpq3**2*vol4)*z**2 - 
     -       beta*(-1 + beta**2)*(-1 + 6*beta**2 - Cpq3*(2*vol2 + Cpq3*vol4))*z**3 + 
     -       beta**2*(-3 - 6*Cpq3*vol2 + beta**2*(3 + 2*Cpq3*vol2) + (-3 + beta**2)*Cpq3**2*vol4)*z**4 + 
     -       2*beta**3*(-1 + beta**2 - Cpq3*(2*vol2 + Cpq3*vol4))*z**5)*xI1(MW**2,musq,ep))/
     -   ((1 + beta*z)*(1 + beta**2 + 2*beta*z)) - 
     -  ((8*beta*s*(1 + Cpq3*vol2) - 2*beta*(20*MW**2 + 7*s)*(1 + 2*Cpq3*vol2 + Cpq3**2*vol4)*z**2 + 
     -       6*beta*(4*MW**2 + s)*(1 + 2*Cpq3*vol2 + Cpq3**2*vol4)*z**4 - 
     -       2*beta**9*s*(1 + 2*Cpq3*vol2 + Cpq3**2*vol4)*z**2*(-2 + z**2) + 
     -       (4*MW**2 + s)*(1 + 2*Cpq3*vol2 + Cpq3**2*vol4)*z*(-5 + 3*z**2) + 
     -       beta**2*z*(2*(s*(11 + 20*Cpq3*vol2 + 9*Cpq3**2*vol4) + 
     -             MW**2*(26 + 60*Cpq3*vol2 + 34*Cpq3**2*vol4)) - 
     -          (48*MW**2 + 19*s)*(1 + 2*Cpq3*vol2 + Cpq3**2*vol4)*z**2 + 
     -          3*(4*MW**2 + s)*(1 + 2*Cpq3*vol2 + Cpq3**2*vol4)*z**4) - 
     -       2*beta**8*s*z*(-1 + z**2 + Cpq3*
     -           (Cpq3*vol4 - (4*vol2 + 5*Cpq3*vol4)*z**2 + 2*(vol2 + Cpq3*vol4)*z**4)) - 
     -       2*beta**3*(4*MW**2*(3 - 2*Cpq3*vol2 - 5*Cpq3**2*vol4 - 
     -             (15 + 22*Cpq3*vol2 + 7*Cpq3**2*vol4)*z**2 + 7*(1 + 2*Cpq3*vol2 + Cpq3**2*vol4)*z**4)
     -           + s*(9*(1 + 2*Cpq3*vol2 + Cpq3**2*vol4) - 4*(4 + 11*Cpq3*vol2 + 7*Cpq3**2*vol4)*z**2 + 
     -             8*(1 + 2*Cpq3*vol2 + Cpq3**2*vol4)*z**4)) - 
     -       beta**4*z*(4*MW**2*(2*(9 + 16*Cpq3*vol2 + 7*Cpq3**2*vol4) - 
     -             (23 + 34*Cpq3*vol2 + 11*Cpq3**2*vol4)*z**2 + 7*(1 + 2*Cpq3*vol2 + Cpq3**2*vol4)*z**4)
     -           + 3*s*(9 - 9*z**2 + 2*z**4 + Cpq3**2*vol4*(9 - 13*z**2 + 2*z**4) + 
     -             2*Cpq3*vol2*(9 - 11*z**2 + 2*z**4))) + 
     -       beta**6*z*(8*MW**2*(5 - 7*z**2 + 2*z**4 - Cpq3**2*vol4*(-3 + z**2) + 
     -             2*Cpq3*vol2*(-2 + z**2)**2) + 
     -          s*(8 - 9*z**2 + 3*z**4 + 2*Cpq3*vol2*(12 - 21*z**2 + 5*z**4) + 
     -             Cpq3**2*vol4*(16 - 33*z**2 + 7*z**4))) - 
     -       2*beta**7*(4*MW**2*(1 + 2*Cpq3*vol2 + Cpq3**2*vol4)*z**2*(-2 + z**2) + 
     -          s*(1 + 2*z**2 + 2*Cpq3*vol2*(-2 + z**2)**2 + Cpq3**2*vol4*(7 + 2*z**2*(-5 + z**2)))) + 
     -       2*beta**5*(4*MW**2*(3 - 3*Cpq3**2*vol4 - 4*(3 + 4*Cpq3*vol2 + Cpq3**2*vol4)*z**2 + 
     -             (5 + 8*Cpq3*vol2 + 3*Cpq3**2*vol4)*z**4) + 
     -          s*(6 - 9*z**2 + 6*z**4 + Cpq3**2*vol4*(16 - 33*z**2 + 8*z**4) + 
     -             2*Cpq3*vol2*(11 + 7*z**2*(-3 + z**2)))))*xI2(MT**2,0._dp,MW**2,musq,ep))/
     -   (8._dp*beta*(1 + beta*z)**2) + ((-1 + beta**2)*(1 + 2*Cpq3*vol2 + Cpq3**2*vol4)*
     -     (2*beta**4*s*z*(-2 + z**2) + 4*beta*s*(-1 + z**2) - (4*MW**2 + s)*z*(-5 + 3*z**2) + 
     -       beta**2*z*(8*MW**2*(-2 + z**2) + 3*s*(-1 + z**2)))*xI2(s,0._dp,0._dp,musq,ep))/(8._dp*beta) + 
     -  ((2*s + 8*Cpq3**2*MW**2*vol4 + beta*
     -        (s*(7 + Cpq3**2*vol4) + 4*MW**2*(1 + 4*Cpq3*vol2 + 7*Cpq3**2*vol4))*z - 
     -       beta**10*s*(1 + 2*Cpq3*vol2 + Cpq3**2*vol4)*z**2*(-2 + z**2) + 
     -       beta**9*s*z*(1 - Cpq3**2*vol4 + (5 + 16*Cpq3*vol2 + 11*Cpq3**2*vol4)*z**2 - 
     -          (3 + 8*Cpq3*vol2 + 5*Cpq3**2*vol4)*z**4) + 
     -       beta**2*(-(s*(5 + 10*Cpq3*vol2 + Cpq3**2*vol4)) + 
     -          2*s*(5 + 6*Cpq3*vol2 + 3*Cpq3**2*vol4)*z**2 + 
     -          4*MW**2*(1 + 6*Cpq3*vol2 + 7*Cpq3**2*vol4 + 2*(1 + 2*Cpq3*vol2 - Cpq3**2*vol4)*z**2)) + 
     -       beta**3*z*(4*MW**2*(-2 + 3*z**2 + 3*Cpq3**2*vol4*(2 - 5*z**2)) + 
     -          s*(-21 + 13*z**2 + Cpq3**2*vol4*(5 + 3*z**2) + 4*Cpq3*vol2*(-7 + 9*z**2))) + 
     -       beta**4*(-4*MW**2*(4 + 6*Cpq3*vol2 + 4*Cpq3**2*vol4 + 
     -             4*(1 + 5*Cpq3*vol2 + 2*Cpq3**2*vol4)*z**2 + (-3 - 8*Cpq3*vol2 + 3*Cpq3**2*vol4)*z**4)
     -           + s*(-1 - 24*z**2 + 13*z**4 + 8*Cpq3*vol2*(1 - 2*z**2)**2 + 
     -             Cpq3**2*vol4*(9 + 4*z**2 - 13*z**4))) + 
     -       beta**5*z*(4*MW**2*(-4 - 2*z**2 + z**4 - Cpq3**2*vol4*(12 - 10*z**2 + z**4) + 
     -             4*Cpq3*vol2*(-3 - z**2 + z**4)) + 
     -          s*(8 - 17*z**2 + 7*z**4 + Cpq3**2*vol4*(24 - 31*z**2 - 7*z**4) + 
     -             4*Cpq3*vol2*(9 - 16*z**2 + 5*z**4))) + 
     -       beta**8*(-4*MW**2*(1 + 2*Cpq3*vol2 + Cpq3**2*vol4)*z**2*(-2 + z**2) + 
     -          s*(-1 + 6*z**2 - z**4 - 2*z**6 + Cpq3**2*vol4*(-7 + 8*z**2 + 15*z**4 - 8*z**6) + 
     -             8*Cpq3*vol2*(-1 + 2*z**2 + z**4 - z**6))) + 
     -       beta**7*z*(4*MW**2*(5 - z**2*(1 + z**2) + Cpq3**2*vol4*(3 + 5*z**2 - 3*z**4) + 
     -             4*Cpq3*vol2*(2 + z**2 - z**4)) - 
     -          s*(-5 + z**2 + 4*z**4 + 4*Cpq3*vol2*(2 - 3*z**2 + 3*z**4) + 
     -             Cpq3**2*vol4*(17 - 21*z**2 + 4*z**6))) + 
     -       beta**6*(s*(5 + 6*z**2 - 11*z**4 + 2*z**6 + 2*Cpq3*vol2*(5 - 19*z**4 + 4*z**6) - 
     -             Cpq3**2*vol4*(-3 + 8*z**2 + 5*z**4 + 4*z**6)) + 
     -          4*MW**2*(3 - 2*z**4 - 6*Cpq3*vol2*z**2*(-2 + z**2) + 
     -             Cpq3**2*vol4*(-3 + 2*z**2*(5 + z**2 - z**4)))))*
     -     xI2(MT**2 - (s*(1 + beta*z))/2._dp,0._dp,MW**2,musq,ep))/
     -   (4._dp*(1 + beta*z)**2*(1 + beta**2 + 2*beta*z)) - 
     -  ((-1 + beta**2)*(1 + 2*Cpq3*vol2 + Cpq3**2*vol4)*
     -     (16*MW**4 + 8*beta*MW**2*s*(beta + z) + 
     -       s**2*(3 + beta*(2*z + beta*(-2 + beta**2 + 2*beta*z + 2*z**2))))*xI3(0._dp,0._dp,s,0._dp,0._dp,0._dp,musq,ep))/
     -   8._dp + ((-1 + beta**2)*(1 + 2*Cpq3*vol2 + Cpq3**2*vol4)*(1 + beta*z)*
     -     (16*MW**4 + 8*beta*MW**2*s*(beta + z) + 
     -       s**2*(3 + beta*(2*z + beta*(-2 + beta**2 + 2*beta*z + 2*z**2))))*
     -     xI3(0._dp,MT**2,MT**2 - (s*(1 + beta*z))/2._dp,0._dp,0._dp,MW**2,musq,ep))/8._dp + 
     -  ((-1 + beta**2)*(1 + 2*Cpq3*vol2 + Cpq3**2*vol4)*
     -     ((4*MW**2 + s)**2*z*(5 - 3*z**2) + 2*beta**6*s**2*z*(-4 + z**2) + 
     -       4*beta*s*(4*MW**2 + s)*(-1 + z**2) - 4*beta**5*s**2*(1 + z**2) + 
     -       beta**4*s*z*(16*MW**2*(-4 + z**2) - 3*s*(-3 + z**2)) - 
     -       8*beta**3*s*(s*z**2 + 2*MW**2*(1 + z**2)) + 
     -       2*beta**2*z*(16*MW**4*(-4 + z**2) + 4*MW**2*s*(-1 + z**2) + s**2*(-11 + 2*z**2)))*
     -     xI3(MT**2,MT**2,s,0._dp,MW**2,0._dp,musq,ep))/(32._dp*beta) - 
     -  ((-1 + beta**2)*(1 + 2*Cpq3*vol2 + Cpq3**2*vol4)*(4*MW**2 + s*(1 + beta**2 + 2*beta*z))*
     -     (16*MW**4 + 8*beta*MW**2*s*(beta + z) + 
     -       s**2*(3 + beta*(2*z + beta*(-2 + beta**2 + 2*beta*z + 2*z**2))))*
     -     xI4(0._dp,0._dp,MT**2,MT**2,s,MT**2 - (s*(1 + beta*z))/2._dp,0._dp,0._dp,0._dp,MW**2,musq,ep))/32._dp))/(mw**2*pi)



      bx(5) = 
     .     (-0.25_dp*alpha*genfacbox*gw**2*mt**2*sigma0*
     .     (beta*z*(1._dp - z**2) - (4._dp*beta*(1._dp - z**2)*
     .     (4._dp*beta - z - 2._dp*beta*z**2)*(-xI1(mh**2,musq,
     .     ep) + xI1(mt**2,musq,ep)))/(s*(1._dp + beta**2 + 2._dp*beta*z)) 
     .     + ((2._dp*beta*(1._dp - beta**2)*s*(2._dp + 3._dp*beta**2 
     .     - beta*(3._dp - 4._dp*beta**2)*z - (1._dp 
     .     + 6._dp*beta**2)*z**2 + beta*(1._dp - 2._dp*beta**2)*z**3 
     .     + 2._dp*beta**2*z**4) + 2._dp*mh**2*(-6._dp*beta**3 - (5._dp 
     .     - 8._dp*beta**2 + 4._dp*beta**4)*z - beta*(5._dp 
     .     - 12._dp*beta**2)*z**2 + (3._dp - 4._dp*beta**2 + 2._dp*beta**4
     .     )*z**3 + beta*(3._dp - 4._dp*beta**2)*z**4))*xI2(mt**2,mt**2,
     .     mh**2,musq,ep))/(beta*s*(1._dp + beta*z)) + ((2._dp*mh**2*z*
     .     (5._dp - 4._dp*beta**2 - (3._dp - 2._dp*beta**2)*z**2) 
     .     - beta*s*(2._dp - beta*(7._dp - 8._dp*beta**2)*z 
     .     - 2._dp*z**2 + beta*(3._dp - 4._dp*beta**2)*z**3))*
     .     xI2(s,mt**2,mt**2,musq,ep))/(beta*s) + ((-4._dp*beta*mh**2*
     .     (beta + z)*(1._dp - 3._dp*beta**2 + beta*(1._dp - 2._dp*beta**2
     .     )*z + 2._dp*beta**2*z**2 + beta**3*z**3) - 2._dp*s*(1._dp 
     .     - 3._dp*beta**6 + beta*(2._dp + 2._dp*beta**2 - 5._dp*beta**4 
     .     - 4._dp*beta**6)*z + beta**2*(1._dp + 3._dp*beta**2 
     .     - 6._dp*beta**4)*z**2 + 2._dp*beta**5*(1._dp + beta**2)*z**3 
     .     + 4._dp*beta**6*z**4 + beta**5*z**5))*xI2(t,mt**2,mh**2,musq,
     .     ep))/(s*(1._dp + beta*z)*(1._dp + beta**2 + 2._dp*beta*z)) - 
     .     (1._dp*(8._dp*mh**4 - 4._dp*mh**2*s*(2._dp - 3._dp*beta**2 
     .     - beta*z) + s**2*(6._dp - 10._dp*beta**2 + 5._dp*beta**4 
     .     - 2._dp*beta*(1._dp - 2._dp*beta**2)*z + beta**2*z**2))*
     .     xI3(0._dp,0._dp,s,mt**2,mt**2,mt**2,musq,ep))/s + ((s**2*(6._dp 
     .     - 10._dp*beta**2 + 5._dp*beta**4 + beta*(3._dp - 4._dp*beta**2 
     .     + 4._dp*beta**4)*z - beta**2*(1._dp - 4._dp*beta**2)*z**2 
     .     + beta**3*z**3) + (1._dp + beta*z)*(8._dp*mh**4 - 4._dp*mh**2*s*
     .     (2._dp - 3._dp*beta**2 - beta*z)))*xI3(0._dp,mt**2,t,mt**2
     .     ,mt**2,mh**2,musq,ep))/s + ((2._dp*mh**4*z*(5._dp 
     .     - 8._dp*beta**2 - (3._dp - 2._dp*beta**2)*z**2) 
     .     - 2._dp*beta*mh**2*s*(1._dp + beta**2 - 2._dp*beta*(5._dp 
     .     - 6._dp*beta**2)*z - (1._dp - beta**2)*z**2 
     .     + 3._dp*beta*(1._dp - beta**2)*z**3) + beta**2*s**2*
     .     (beta - 2._dp*beta**3 - (7._dp - 13._dp*beta**2 
     .     + 8._dp*beta**4)*z + beta*(1._dp - 2._dp*beta**2)*z**2 
     .     - 2._dp*beta**2*(1._dp - beta**2)*z**3))*
     .     xI3(s,mt**2,mt**2,mt**2,mt**2,mh**2,musq,ep))/(beta*s) 
     .     + (0.5_dp*(-16._dp*mh**6 + 16._dp*mh**4*s*(1._dp - 2._dp*beta**2 
     .     - beta*z) - 2._dp*mh**2*s**2*(6._dp - 13._dp*beta**2 
     .     + 10._dp*beta**4 - 6._dp*beta*(1._dp - 2._dp*beta**2)*z 
     .     + beta**2*(2._dp + beta**2)*z**2) - s**3*(2._dp + beta**2 
     .     - 6._dp*beta**4 + 4._dp*beta**6 + beta*(5._dp - 10._dp*beta**2 
     .     + 8._dp*beta**4)*z + beta**4*(1._dp + 2._dp*beta**2)*z**2 
     .     + beta**3*z**3))*xI4(0._dp,0._dp,mt**2,mt**2,s,t,mt**2,mt**2,
     .     mt**2,mh**2,musq,ep))/s))/(mw**2*pi)

c--- use the value of "tevscale" (passed via anomcoup.f)
c--- as an anomalous top Yukawa coupling:
c--- g(top Y) = (tevscale) x g(SM, top Y)
c--- it therefore affects Higgs diagrams as the square
      vrts(5)=vrts(5)*tevscale**2
      slf(5) =slf(5) *tevscale**2
      vrt(5) =vrt(5) *tevscale**2
      bx(5)  =bx(5)  *tevscale**2

      vew = 
     .     + (vrts(1) + vrts(2) + vrts(3) + vrts(4) + vrts(5))/2._dp
     .     + slf(1) + slf(2) + slf(3) + slf(4) + slf(5)
     .     + vrt(1) + vrt(2) + vrt(3) + vrt(4) + vrt(5)
     .     + bx(1) + bx(2) + bx(3) + bx(4) + bx(5)

      vew = vew + (trih + trizx)/2._dp

      vew = vew/born

      end subroutine ggQQb_ew_oneloop


       

*   HASEL
      subroutine HASEL(along,alat,elav,bear,f,cha,nos,out,typo,tol,m,l)
      implicit real*8 (a-h,o-z)
      real*4 yp(301,15,2)
      character type*21,typo*21,ch1,cl1*21,cl2*21,cha,app
      dimension x(11),ytmp(11),xx(11),xtmp(11),out(11,14),typo(14)
     &,v(3),ryx(3,3)
      common /dat/foE,hmE,ymE,foF1,hmF1,ymF1,foF2,hmF2,ymF2,re,model
      common /grid/blhla,blhlo,rmg,dla,dlo,drg,dtz,yp,nr,nla,nlo
      common /mag/fhs,fr,xm,ym,zm,px,py,pz,ir
      re=6378.135
      scale=45./atan(1.)
      n=9
      noz=nos
      fr=f
      if(m.lt.0)then
      app='s'
      else
      app='L'
      endif
      model=abs(m)
      if(model.eq.7) dtz=1.
      if(tol.lt.0.0) then
      dtw=dtz
      dtz=0.0
      endif
      tolx=abs(tol)
***** dipole moment = (pxt,pyt,pzt) and dipole location = (xm,ym,zm) ***
      dlat=77.75
      dlong=295.75
      pxt=cos(dlat/scale)*cos(dlong/scale)
      pyt=cos(dlat/scale)*sin(dlong/scale)
      pzt=sin(dlat/scale)
      xm=-434.66
      ym=199.47
      zm=80.25
       fhs=2.8*.31*re**3
      px=fhs*pxt/fr
      py=fhs*pyt/fr
      pz=fhs*pzt/fr
      if(cha.ne.'y') fhs=0.0
***************************************************************************
      trhla=blhla+real(nla-1)*dla
      trhlo=blhlo+real(nlo-1)*dlo
      nhop=abs(nos)
       if(cha.eq.'n') then
       ira=0
       irb=0
       else
       ira=1
       irb=2
       if(abs(nos).eq.2) nhop=3
       if(abs(nos).ge.3) nhop=7
       endif
      nos=0
****************************************************************************
      do 909 irr=ira,irb
***** ray type loop ********************************************************
       if(irr.eq.0) then
       ch1='N'
       ir=0
       else if(irr.eq.1) then 
       ch1='0' 
       ir=1
       else if(irr.eq.2) then
       ch1='X'
       ir=-1
       endif
      dr=sin(elav/scale)
      dt=-cos(elav/scale)*cos(bear/scale)
      dp=cos(elav/scale)*sin(bear/scale)
      phi=atan(1.)*along/45.
      theta=2.*atan(1.)-atan(1.)*alat/45.
      x(1)=re*sin(theta)*cos(phi)
      x(2)=re*sin(theta)*sin(phi)
      x(3)=re*cos(theta)
      x1=x(1)
      x2=x(2)
      x3=x(3)
      x(4)=dr*sin(theta)*cos(phi)+dt*cos(theta)*cos(phi)-dp*sin(phi)
      x(5)=dr*sin(theta)*sin(phi)+dt*cos(theta)*sin(phi)+dp*cos(phi)
      x(6)=dr*cos(theta)-dt*sin(theta)
      t=0.0
      ita=2
      istep=1
      type='           '
      x(7)=0.0
      x(8)=0.0
      x(9)=0.0
***********************************************************************************
      do 979 ihop=1,nhop
***** hop loop ********************************************************************
      zmax=0.0
      sc=sqrt(x(4)**2+x(5)**2+x(6)**2)
      x(4)=x(4)/sc
      x(5)=x(5)/sc
      x(6)=x(6)/sc
       if(cha.eq.'y'.and.ihop.gt.1) then
        if(ihop.eq.2) then
          do k=1,n
           xtmp(k)=x(k)
          enddo
         ir=1
         ch1='0'
        else
         if(ihop.eq.3.or.ihop.eq.4.or.ihop.eq.6) ir=-1
         if(ihop.eq.5.or.ihop.eq.7) ir=1
          if(ihop.eq.3) then
            do k=1,n
             ytmp(k)=x(k)
             x(k)=xtmp(k)
            enddo
           cl1=type
           ita1=itb+1
          endif
          if(ihop.eq.4) then
            do k=1,n
             xtmp(k)=x(k)
            enddo
           cl2=type
           ita2=itb+1
          endif
           if(ir.eq.1) ch1='o'
           if(ir.eq.-1) ch1='X'
          if(ihop.eq.4.or.ihop.eq.5) then 
            do k=1,n
             x(k)=ytmp(k)
            enddo
           type=cl1
           ita=ita1
          endif
          if(ihop.eq.6.or.ihop.eq.7) then
            do k=1,n
             x(k)=xtmp(k)
            enddo
           type=cl2
           ita=ita2
          endif
        endif
       endif
      zad=0.0
      zed=0.0
      itb=ita+1
      type(ita:itb+1)='S'//ch1//' '
      itx=2
       if(ch1.eq.'N') then
       type(itb:itb)=' '
       itb=itb-1
       itx=itx-1
       endif
      hmax=50.
      hmin=.1
	h=hmax
****************************************	****
***** start of calculations for an individual hop *************************************
 999  istep=istep+1
       do i=1,6
       xx(i)=x(i)
       enddo 
      pl=x(7)
      dopx=x(8)
      pp=x(9)
      zbd=zad
      zad=zed
      xlong=rlong
      xlat=rlat
***** ODE solver**********************************************************************

*     h=hmax
       call rkf(t,x,n,h,tolx,hmin,hmax) 
*      x(1) 到x（3）是当前射线所在位置的直角坐标
*      x(7)到 x(9)是 同一位置的群路径，多普勒平移和相速度值
***************************************************************************************
      red=sqrt(x(1)**2+x(2)**2+x(3)**2)
      zedl=zed
      zed=red-re
      rlat=90.-acos(x(3)/red)*scale
       if(abs(x(1)).lt.1.0d-10) then
       rlong=atan2(x(2),1.0d-10)*scale
       else
       rlong=atan2(x(2),x(1))*scale
       endif
      if(rlong.lt.0.0) rlong=rlong+360.0
        if(l.gt.0) then
        call magnetic(x,v,ryx,ha,dip,dcl)
        vu=x(4)*v(1)+x(5)*v(2)+x(6)*v(3)
        xt=sqrt(abs(x(4)**2+x(5)**2+x(6)**2))
        angle=acos(.9999999*vu/(ha*xt))*scale
       endif
      if((rlat.gt.trhla.or.rlat.lt.blhla).and.model.ne.7) go to 909
      if((rlong.gt.trhlo.or.rlong.lt.blhlo)
     $.and.model.ne.7.and.noz.gt.0)  go to 909
****************************************************************************************
*     在这个解决方案中，zed表示射线的高度，rlong 表示射线的经度，rlat表示射线的纬度
***** label for highest reflection layer of current hop **************************
       if(zad.gt.max(zbd,zed,zmax)) then
        if(model.eq.10) then
          call layerparm (xlong,xlat,h1,h2,h3,f1,f2,f3)
          hmE=h1
          hmF1=h2
          hmF2=h3
        endif
        if(zad.le.hmE) then
         itb=ita+1
         type(ita:itb+1)='E'//ch1//' '
         itx=2
        else if(zad.le.hmF1) then
         type(ita:itb)='F1'//ch1
         itb=ita+2 
         itx=3
        else if(zad.le.hmF2) then
         itb=ita+2 
         type(ita:itb)='F2'//ch1
         itx=3
        else
         itb=ita+1
         type(ita:itb+1)='S'//ch1//' '
         itx=2
        endif
        if(ch1.eq.'N') then
         type(itb:itb)=' '
         itb=itb-1
         itx=itx-1
        endif
       endif
******************************************************************************
      zmax=max(zad,zmax)
      if(zed.ge.600.) go to 909
      if(zed.gt.0.) go to 999
      if(ihop.eq.1.or.cha.ne.'y') ita=ita+itx
***** interpolation for values at ground level********************************
      pop=zedl/(zedl-zed)
       do i=1,6
       x(i)=xx(i)+pop*(x(i)-xx(i))
       enddo
      xl=0.0
      xxl=0.0
      cosa=0.0
       do i=1,3
       xl=xl+x(i)**2
       xxl=xxl+x(i+3)**2
       cosa=cosa+x(i)*x(i+3)
        enddo
      elavf=abs(90.-45.*acos(cosa/sqrt(xl*xxl))/atan(1.))
      sd=sqrt((x(1)-x1)**2+(x(2)-x2)**2+(x(3)-x3)**2)
      sd=2.*re*asin(.5*sd/re)
      pl=pl+pop*(x(7)-pl)
      pp=pp+pop*(x(9)-pp)
      dopx=dopx+pop*(x(8)-dopx)
      x(7)=pl
      x(8)=dopx
      x(9)=pp
      dopx=1.0d6*dopx
      xlong=xlong+pop*(rlong-xlong)
      xlat=xlat+pop*(rlat-xlat)
      nos=nos+1
      out(1,nos)=dopx
      out(2,nos)=sd
      out(3,nos)=pl
      out(4,nos)=pp
      out(5,nos)=zmax
      out(6,nos)=elav
      out(7,nos)=elavf
      out(8,nos)=xlong
      out(9,nos)=xlat
      out(10,nos)=bear
      out(11,nos)=fr
      typo(nos)=type
      ss=x(1)**2+x(2)**2+x(3)**2
      xr=x(4)*x(1)+x(5)*x(2)+x(6)*x(3)
      x(4)=x(4)-2.*xr*x(1)/ss
      x(5)=x(5)-2.*xr*x(2)/ss
      x(6)=x(6)-2.*xr*x(3)/ss
      if(l.gt.0) then
      zod=0.0
      write(l,fmt='(4f15.7)') zod,xlong,xlat,angle
      endif
 979  continue
 909  continue
      if(tol.lt.0.0) dtz=dtw
      return
      end

* LAYERPARM 
* 在某一特定点的层参数
*     输入双精度，elong 表示采样点的经度，elat表示纬度
*     输出也是双精度，h1表示E层的高度，h2表示F1层的高度，
*     h3表示F2层的高度，f1表示E层的临界频率
*     f2表示F1层的临界频率，f3表示F2层的临界频率
***************************************************************************************
***** subroutines block*********************************************************
* ELDENX
***** common block ***************************************************************
* GRID (defined in HASEL)
**********************************************************************************
      subroutine layerparm (elong,elat,h1,h2,h3,f1,f2,f3)
      implicit real*8 (a-h,o-z)
       real*4 yp(301,15,2)
      common /grid/blhla,blhlo,rmg,dla,dlo,drg,dtz,yp,nr,nla,nlo 
      re=6378.135
      r0=rmg
      e0=eldenx(elong,elat,r0,d10,dthet,dphi,1)
      r1=r0+drg
      e1=eldenx(elong,elat,r1,d11,dthet,dphi,1)
      h0=rmg-re
      h1=0.0
      h2=0.0
      h3=0.0
       ea=0.0
       eb=0.0
      ec=0.0
        do i=3,nr
       r2=rmg+drg*(i-1)
       e2=eldenx(elong,elat,r2,d12,dthet,dphi,1)
       d21=(d12-d10)/(r2-r0)
       hh=0.0
       ee=0.0
       hz=0.0
        if(d11*d10.le.1.d-33.and.abs(d11-d10).gt.1.d-33) then
***** check for maximum***************************************************************
         hh=(d11*r0-d10*r1)/(d11-d10)-re
         ee=(d11*e0-d10*e1)/(d11-d10)
          if(hh.gt.100.0) then
            if(hh.lt.140.) then
              h1=hh
              ea=ee
            else
               if(hh.lt.230.) then
              h2=hh
              eb=ee
               else
              h3=hh
              ec=ee
               endif
            endif
          endif
        else
***** check for point of inflection ***********************
          if(i.gt.3) then
            if(r1-re.lt.230.) then
              if(d21*d20.le.1.d-33.and.abs(d21-d20).gt.1.d-33) then
              hz=(d21*r0-d20*r1)/(d21-d20)-re
              ez=(d21*e0-d20*e1)/(d21-d20)
               if(hz.gt.100.0.and.hz.lt.230.) then
                 if(hz.ge.140.) then
                   if(h2.lt.1.) then
                      h2=hz
                      eb=ez
                   endif
                 else
                   if(h1.lt.1.) then   
                       h1=hz
                       ea=ez
                   endif
                 endif
               endif
              endif
            endif
          endif
        endif
       d10=d11
       d11=d12
       d20=d21
       e0=e1
       e1=e2
       r0=r1
       r1=r2
       enddo
       if(h1.lt.h0) h1=h0
      if(h2.lt.h0) h2=h1
      if(h3.lt.h0) h3=h2
       if(ec.lt.1.) then
       ec=e2
       h3=r2-re
       endif
       f1=sqrt(abs(ea*80.6d-6))
       f2=sqrt(abs(ea*80.6d-6))
       f3=sqrt(abs(ea*80.6d-6))
      return
      end
***** ELECTRON 
*     电子密度和它的梯度
*     输入为双精度，其中x（*）表示地球上采样点的直角坐标，输出也是双精度输出，
*     en表示电子密度（每立方厘米的电子），dnx（*）表示笛卡尔衍生的电子密度

***** subroutines required*********************************************************
* ELDENX,ELDEN,LAYERX and ELECTRONR***************************************************
***** common blocks******************************************************************
* GRID and DAT (defined in HASEL)

      subroutine electron(x,en,dnx)
      implicit real*8 (a-h,o-z)
      dimension x(3),dnx(4)
       real*4 yp(301,15,2)
       common /grid/blhla,blhlo,rmg,dla,dlo,drg,dtz,yp,nr,nla,nlo
       common /dat/foE,hmE,ymE,foF1,hmF1,ymF1,foF2,hmF2,ymF2,re,model
      r=sqrt(x(1)**2+x(2)**2+x(3)**2)
      scale=45./atan(1.)
      theta=acos(x(3)/r)
      phi=atan2(x(2),x(1))
      elat=90.-scale*theta
      elong=scale*phi
      if(elong.lt.0.0) elong=elong+360.
***** 参数 da应该低于最小角度的20%
***** scale (in degrees).
      da=.002
       if(model.eq.10) then
       en=eldenx(elong,elat,r,diffr,dthet,dphi,1)
       dthet=scale*dthet
       dphi=scale*dphi
        if(dtz.gt.1.d-20) then 
         enp=eldenx(elong,elat,r,d1,d2,d3,-2)
         enm=en
        endif
       else
        if(model.eq.7) then 
***** 参数dr应该低于最小的高度的20%***************************************
***** scale (in kilometers).
        dr=.2
        t=0.0
        enpp=elden(t,elong+da,elat,r)
        enpm=elden(t,elong-da,elat,r)
        entp=elden(t,elong,elat+da,r)
        entm=elden(t,elong,elat-da,r)
         if(dtz.gt.1.d-20) then
          enp=elden(t+.5d0*dtz,elong,elat,r)
          enm=elden(t-.5d0*dtz,elong,elat,r)
         endif
       diffr=.5*(elden(t,elong,elat,r+dr)-elden(t,elong,elat,r-dr))/dr
       else
        call layerx(elong+da,elat,1)
        call electronr(r,enpp,dpp)
        call layerx(elong-da,elat,1)
        call electronr(r,enpm,dpm)
        call layerx(elong,elat+da,1)
        call electronr(r,entp,dtp)
        call layerx(elong,elat-da,1)
        call electronr(r,entm,dtm)
         if(dtz.gt.1.d-20) then
          call layerx(elong,elat,-2)
          call electronr(r,enp,dq)
          enm=.25*(enpp+enpm+entp+entm)
         endif
         diffr=.25*(dpp+dpm+dtp+dtm)
        endif
       en=.25*(enpp+enpm+entp+entm)
       dphi=.5*scale*(enpp-enpm)/da
       dthet=-.5*scale*(entp-entm)/da
       endif
       if(dtz.lt.1.d-20) then
       dtime=dphi/(240.*scale)
       else
       dtime=(enp-enm)/dtz
       endif
      dnx(1)=(diffr*x(1)+dthet*cos(theta)*cos(phi)
     &-dphi*sin(phi)/sin(theta))/r
      dnx(2)=(diffr*x(2)+dthet*cos(theta)*sin(phi)
     &+dphi*cos(phi)/sin(theta))/r
      dnx(3)=(diffr*x(3)-dthet*sin(theta))/r
      dnx(4)=dtime
      return
      end
***** LAYERX**********************************************************************
* 从网格值中得到层参数
*     输入双精度，elong表示采样点的经度，elat表示纬度，
*     iuu为整数，表示时间段的标签 （1或者2）
***** subroutine required*********************************************************
* TERP
***** common blocks***************************************************************
* GRID and DAT (defined in HASEL)************************************************
*   注意 输出被放置在DAT中*******************************************************

      subroutine layerx(elong,elat,iuu)
      implicit real*8 (a-h,o-z)
       real*4 yp(301,15,2)
      dimension yz(4),yy(9)
       common /grid/blhla,blhlo,rmg,dla,dlo,drg,dtz,yp,nr,nla,nlo
       common /dat/yy,re,model
      lip=abs(iuu)
      dlong=elong-blhlo
      if(dlong.gt.360.) dlong=dlong-360.
      la=min(max(2,int((elat-blhla)/dla)+1),nla-2) 
      lo=int(dlong/dlo)+1
      lo=lo-int((lo-1)/nlo)*nlo
      if(lo.lt.1) lo=lo+nlo
      dax=(elat-blhla)/dla-real(la-1)
      dox=dlong/dlo-real(lo-1)
      n0=(la-2)*nlo
      n1=n0+nlo
      n2=n1+nlo
      n3=n2+nlo
       do i=1,9
        do k=0,3
         kk=lo-1+k
         if(kk.lt.1) kk=kk+nlo
         if(kk.gt.nlo) kk=kk-nlo
         x0=yp(i,n0+kk,lip)
         x1=yp(i,n1+kk,lip)
         x2=yp(i,n2+kk,lip)
         x3=yp(i,n3+kk,lip)
         yz(k+1)=terp(dax,x0,x1,x2,x3)
        enddo
       yy(i)=terp(dox,yz(1),yz(2),yz(3),yz(4))
        enddo
       return
       end
*ELDENX
*从格点值中得到电子浓度和梯度
      real*8 function eldenx(elong,elat,r,diffr,dthet,dphi,iuu)
      implicit real*8 (a-h,o-z)
      real*4 yp(301,15,2)
      dimension yz(4),yzt(4),e(4),dt(4),dp(4)
       common /grid/blhla,blhlo,rmg,dla,dlo,drg,dtz,yp,nr,nla,nlo
      lip=abs(iuu)
      if(r.lt.rmg.or.r.gt.rmg+(nr-1)*drg)then
	    eldenx=0.0
	    diffr=0.0
          dthet=0.0
		dphi=0.0
	else
	    dlong=elong-blhlo
	    if(dlong.gt.360.) dlong=dlong-360.
	    la=min(max(2,int((elat-blhla)/dla+1)),nla-2)
	    lo=int(dlong/dlo)+1
	    lo=lo-int((lo-1)/nlo)*nlo
	    if(lo.lt.1) lo=lo+nlo
	    lr=min(max(1,int((r-rmg)/drg)+1),nr-2)
	    dax=(elat-blhla)/dla-real(la-1)
	    dox=dlong/dlo-real(lo-1)
	    drx=(r-rmg)/drg-real(lr-1)
	    n0=(la-2)*nlo
	    n1=n0+nlo
	    n2=n1+nlo
	    n3=n2+nlo
	    if(lr.lt.2)then
			nstart=0
	    else
			nstart=-1
	    endif
	    do i=nstart,2
			do k=0,3
				kk=lo-1+k
				if(kk.lt.1) kk=kk+nlo
				if(kk.gt.nlo) kk=kk-nlo
				x0=yp(lr+i,n0+kk,lip)
				x1=yp(lr+i,n1+kk,lip)
				x2=yp(lr+i,n2+kk,lip)
				x3=yp(lr+i,n3+kk,lip)
				if(iuu.gt.0) yzt(k+1)=dterp(dax,x0,x1,x2,x3)
	            yz(k+1)=terp(dax,x0,x1,x2,x3)
			enddo
			yy=terp(dox,yz(1),yz(2),yz(3),yz(4))
	        if(iuu.gt.0) then
	            yyt=terp(dox,yzt(1),yzt(2),yzt(3),yzt(4))
	            yyp=dterp(dox,yz(1),yz(2),yz(3),yz(4))
	        else
	            yyt=0
                  yyp=0
	        endif
	        e(i+2)=yy
              dt(i+2)=yyt
              dp(i+2)=yyp
	    enddo
	    if(lr.lt.2) then
	        e(1)=0.0
              dt(1)=0.0
              dp(1)=0.0
          endif
	    eldenx=terp(drx,e(1),e(2),e(3),e(4))
	    if(eldenx.lt.0.0) then
	        eldenx=0.0
              diffr=0.0
              dthet=0.0
              dphi=0.0
	    else
	        if(iuu.gt.0) then
	            diffr=dterp(drx,e(1),e(2),e(3),e(4))
	            dthet=-terp(drx,dt(1),dt(2),dt(3),dt(4))
	            dphi=terp(drx,dp(1),dp(2),dp(3),dp(4))
	        else
                  diffr=0.0
                  dthet=0.0
                  dphi=0.0
              endif
          endif
      endif
      dthet=dthet/dla
      dphi=dphi/dlo
      diffr=diffr/drg
      return
      end
     	        
              
              
*TERP
*三次插值函数
      real*8 function terp(d2,p0,p1,p2,p3)
      implicit real*8 (a-h,o-z)
      character app
      common /approx/app
      if(app.eq.'L') then
          d1=d2+1
          d3=d2-1
          d4=d2-2
          terp=(p3*d1*d2*d3-p0*d2*d3*d4)/6.+(p1*d1*d3*d4-p2*d1*d2*d4)/2.
      else
          terp=p1+.5*(p2-p0)*d2+(p0-2.5*p1+2.*p2-.5*p3)*d2**2
     $+(.5*p3-1.5*p2+1.5*p1-.5*p0)*d2**3
      endif
      return
      end
 
*DTERP
*三次插值函数的导数
      real*8 function dterp(d2,p0,p1,p2,p3)
      implicit real*8 (a-h,o-z)
      character app
      common /approx/app
      if(app.eq.'L') then
          d1=d2+1
          d3=d2-1
          d4=d2-2
          dterp=(p3*(d2*d3+d1*d3+d1*d2)-p0*(d3*d4+d2*d4+d2*d3))/6.
     $+(p1*(d3*d4+d1*d4+d1*d3)-p2*(d2*d4+d1*d4+d1*d2))/2.
      else
          dterp=.5*(p2-p0)+(2.*p0-5.*p1+4.*p2-p3)*d2
     $+(1.5*p3-4.5*p2+4.5*p1-1.5*p0)*d2
      endif
      return
      end               
 
*DET
*三维数组
      real*8 function det(a1,a2,a3)
	implicit real*8 (a-h,o-z)
	dimension a1(3),a2(3),a3(3)
	det=a1(1)*(a2(2)*a3(3)-a3(2)*a2(3))-a2(1)*(a1(2)*a3(3)
     $-a1(3)*a3(2))+a3(1)*(a1(2)*a2(3)-a1(3)*a2(2))
	return
	end

*ELECTRONR
*卡普曼层模型的电子浓度和梯度(foE等中使用了字母o而非数字0)
      subroutine electronr(r,en,diffr)
	implicit real*8 (a-h,o-z)
	dimension a1(3),a2(3),a3(3),b(3)
	common /dat/foE,hmE,ymE,foF1,hmF1,ymF1,foF2,hmF2,ymF2,re,model 
	if(r.lt.re) then
	    eN=0.0
	    diffr=0.0
	else
	    if(max(foE,foF1).lt.1.d-7) then
	        xN=foF2**2/80.6d-6
	        eN=xN*chap(1.d0,1.4142d0*(r-re-hmF2)/ymF2)
	        diffr=1.4142*xN*dchap(1.d0,1.4142d0*(r-re-hmF2)/ymF2)/ymF2
	    else
	        a1(1)=1.
	        a2(1)=chap(.5d0,2.d0*(hmE-hmF1)/ymF1)
	        a3(1)=chap(1.d0,1.4142d0*(hmE-hmF2)/ymF2)
	        a1(2)=chap(.5d0,2.d0*(hmF1-hmE)/ymE)
	        a2(2)=1.
	        a3(2)=chap(1.d0,1.4142d0*(hmF1-hmF2)/ymF2)
	        a1(3)=chap(.5d0,2.d0*(hmF2-hmE)/ymE)
	        a2(3)=chap(.5d0,2.d0*(hmF2-hmF1)/ymF1)
	        a3(3)=1
	        b(1)=foE**2/80.6d-6
	        b(2)=foF1**2/80.6d-6
	        b(3)=foF2**2/80.6d-6
	        de=det(a1,a2,a3)
	        c1=det(b,a2,a3)/de
	        c2=det(a1,b,a3)/de
	        c3=det(a1,a2,b)/de
	        eN=c1*chap(.5d0,2.d0*(r-re-hmE)/ymE)
     $+c2*chap(.5d0,2.d0*(r-re-hmF1)/ymF1) 
     $+c3*chap(1.d0,1.4142d0*(r-re-hmF2)/ymF2)
	        diffr=2.*c1*dchap(.5d0,2.d0*(r-re-hmE)/ymE)
     $+2.*c2*chap(.5d0,2.d0*(r-re-hmF1)/ymF1)
     $+1.4142*c3*chap(1.d0,1.4142d0*(r-re-hmF2)/ymF2)
	    endif
	endif
	return
	end


*CHAP
*卡曼层函数
      real*8 function chap(C,x)
	implicit real*8 (a-h,o-z)
	xx=exp(-x)
	if(abs(xx).gt.1.d30) then
	    chap=0
	else
	    chap=exp(C*(1.d0-x-xx))
	endif
	return
	end

*DCHAP
*卡普曼层函数的导数
      real*8 function dchap(C,x)
	implicit real*8 (a-h,o-z)
	xx=exp(-x)
	if(abs(xx).gt.1.d30) then
	    dchap=0
	else
	    dchap=C*(exp(-x)-1.d0)*exp(C*(1.d0-x-xx))
	endif
	return
	end

********Page 24,25 *************
********** MAGNETIC *************
	subroutine magnetic(x,v,ryx,ha,dip,dcl)
	implicit real*8 (a-h,o-z)
	dimension x(3),v(3),ryx(3,3)
	common /mag/fhs,fr,xm,ym,zm,px,py,pz,ir
	xmt=x(1)-xm
	ymt=x(2)-ym
	zmt=x(3)-zm
	pm=xmt*px+ymt*py+zmt*pz
	rr=xmt*xmt+ymt*ymt+zmt*zmt
      r=sqrt(rr)
	rp=1./(rr*r)
	v(1)=(px-3.*pm*xmt/rr)*rp
	v(2)=(py-3.*pm*ymt/rr)*rp
	v(3)=(pz-3.*pm*zmt/rr)*rp
	ra=sqrt(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
	ha=sqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))
	vr=v(1)*x(1)+v(2)*x(2)+v(3)*x(3)
	dip=-45.*asin(vr/(ha*ra))/atan(1.)
	dcl=45.*atan((v(2)*x(1)-v(1)*x(2))/(v(3)*ra-vr*x(3)))/atan(1.)
	ryx(1,1)=(-3.*pm/rr+6.*pm*xmt*xmt/(rr*rr)-3.*px*xmt/rr)*rp
     &-3.*xmt*(px-3.*pm*xmt/rr)*rp/rr
	ryx(1,2)=(-3.*py*xmt/rr+6.*pm*ymt*xmt/(rr*rr))*rp
     &-3.*ymt*(px-3.*pm*xmt/rr)*rp/rr
	ryx(1,3)=(-3.*pz*xmt/rr+6.*pm*zmt*xmt/(rr*rr))*rp
     &-3.*zmt*(px-3.*pm*xmt/rr)*rp/rr
	ryx(2,2)=(-3.*pm/rr+6.*pm*ymt*ymt/(rr*rr)-3.*py*ymt/rr)*rp
     &-3.*ymt*(py-3.*pm*ymt/rr)*rp/rr
	ryx(2,3)=(-3.*pz*ymt/rr+6.*pm*zmt*ymt/(rr*rr))*rp
     &-3.*zmt*(py-3.*pm*ymt/rr)*rp/rr
	ryx(3,3)=(-3.*pm/rr+6.*pm*zmt*zmt/(rr*rr)-3.*pz*zmt/rr)*rp
     &-3.*zmt*(pz-3.*pm*zmt/rr)*rp/rr
	ryx(2,1)=ryx(1,2)
	ryx(3,2)=ryx(2,3)
	ryx(3,1)=ryx(1,3)
	return
	end

*********FUNC************

	subroutine func(vec,f,n)
	implicit real*8 (a-h,o-z)
	dimension f(n),vec(n),v(3),ryx(3,3),dnx(4)
	common /mag/fhs,fr,xm,ym,zm,px,py,pz,ir
	call electron(vec,en,dnx)
	rxt=8.06d-5*dnx(4)/(fr*fr)
	rx=8.06d-5*en/(fr*fr)
	vv=vec(4)**2+vec(5)**2+vec(6)**2
	if(abs(fhs).gt.1.0e-20) then
		call magnetic(vec,v,ryx,ha,dip,dcl)
		vu=(vec(4)*v(1)+vec(5)*v(2)+vec(6)*v(3))/ha
		ry=ha
		ry2=ry*ry
		ga=-vu**2
		gg=-ga/vv
		al=1.-rx-ry2
		be=.5*ry*(1.-ga)
		fac=sqrt(abs(be**2-al*ga))
		if(ir.eq.1) then
			rp=-ga/(be+fac)
			rj=4.*fac
		else
			rp=-(be+fac)/al
			rj=-4.*fac
		endif
		rp2=rp**2
		rk=-2.*rx*vu*(rp*ry-1.)
		rl=ry*(1.-ry2)*rp2-2.*(1.-rx-ry2)*rp-ry
		rm=2.*rx*rp*(rp*ry-1.)

		do i=1,3
		    f(i)=rj*vec(i+3)-rk*v(i)
		    dxx=8.06e-5*dnx(i)/(fr*fr)
		    f(i+3)=rl*dxx+(rk*vec(4)+rm*v(1))*ryx(1,i)
     &+(rk*vec(5)+rm*v(2))*ryx(2,i)+(rk*vec(6)+rm*v(3))*ryx(3,i)
	    enddo
		rx1=1.-rx
		a=al+rx*ry2*gg
		b=rx1*al-.5*rx*ry2*(1.-gg)
		c=rx1*(rx1**2-ry2)
		emu=(b+ir*sqrt(abs(b**2-a*c)))/a
		eml=emu*a-b
		ems=sqrt(abs(emu))
		if(abs(eml).gt.1.0d-9) then
			aa=1.-1.5*rx-1.5*ry2+2.*rx*ry2*gg
			bb=rx1*(1.-3.*rx)-2.*ry2+1.5*rx*ry2*(1.+gg)
			cc=-1.5*rx*rx1**2-ry2*(.5-rx)
			emu=(aa*emu*ems-bb*ems+cc/ems)/eml
			ax=-1.+gg*ry2
			bx=-2.*rx1+.5*ry2*(1.+gg)
			cx=-3.*rx1**2+ry2
			emx=.25*(-ax*ems**4+2.*bx*ems**2-cx)/(a*ems**3-b*ems)
		else
			emu=1./sqrt(abs(1.-rx))
			emx=-.5*emu
		endif
	else
		ems=sqrt(abs(1.-rx))
		do i=1,3
			f(i)=vec(i+3)
			dxx=8.06e-5*dnx(i)/(fr*fr)
			f(i+3)=-.5*dxx
		enddo
		emu=1./ems
		emx=-.5*emu
	endif
	cosa=(vec(4)*f(1)+vec(5)*f(2)+vec(6)*f(3))/sqrt(vv)
	f(7)=emu*cosa
	f(8)=-cosa*fr*rxt*emx/3.e5
	f(9)=ems*cosa
	if(abs(fhs).gt.1.0e-20) then
		do i=1,n
			f(i)=real(ir)*f(i)
		enddo
	endif
	return
	end





*RKF 
*龙格库塔法解微分方程
      subroutine rkf(t,y,n,h,tol,hmin,hmax)
	implicit real*8 (a-h,o-z)
	dimension a1(11),a2(11),a3(11),a4(11),a5(11),a6(11),y(n),yt(11)
	iend=0
99	iend=iend+1
	call func(y,a1,n)
	alp=0.0
	do i=1,n
	    alp=alp+a1(i)**2
	    yt(i)=y(i)+.25*a1(i)*h
	enddo
	alp=max(sqrt(abs(alp)),1.0d-32)
	call func(yt,a2,n)
	do i=1,n
	    yt(i)=y(i)+(0.09375*a1(i)+.28125*a2(i))*h
	enddo
	call func(yt,a3,n)
	do i=1,n
	    yt(i)=y(i)+(0.87938097406d0*a1(i)-3.2771961766d0*a2(i)
     $+3.32089212563d0*a3(i))*h
	enddo
	call func(yt,a4,n)
	do i=1,n
	    yt(i)=y(i)+(2.0324074074d0*a1(i)-8.*a2(i)
     $+7.17348927875d0*a3(i)-.20589668616d0*a4(i))*h
	enddo
	call func(yt,a5,n)
	do i=1,n
	    yt(i)=y(i)+(-.29629629629d0*a1(i)+2.*a2(i)
     $-1.38167641326d0*a3(i)+.45279270955d0*a4(i)-.275*a5(i))*h
	enddo
	call func(yt,a6,n)
	erx=0.0
	do i=1,3
	    erz=.277777777777778d-2*a1(i)-.299415204678363d-1*a3(i)
     $-.291998936735779d-1*a4(i)+.2d-1*a5(i)+.363636363636364d-1*a6(i)
	    erx=erx+erz**2
	enddo
	hold=h
	if(erx.gt.0.001*tol) then
	    s=.84*sqrt(sqrt(abs(tol/erx)))
	    h=s*h
	else
	    h=2.*h
	endif
	if(h.lt.hmin/alp) h=hmin/alp
	if(h.gt.hmax/alp) h=hmax/alp
	if(erx.gt.tol.and.iend.lt.8) then
	    go to 99
	else
	    do i=1,n
	        y(i)=y(i)+(.11851851852*a1(i)+.51898635478d0*a3(i)
     $+.50613149034d0*a4(i)-.18*a5(i)+.03636363636*a6(i))*hold
	    enddo
	    t=t+hold
	    return
	endif
	end

********ELDEN**********
***DAT (defined in HASEL)*******
	real*8 function elden(t,elong,elat,r)
	implicit real*8 (a-h,o-z)
	common /dat/foE,hmE,ymE,foF1,hmF1,ymF1,foF2,hmF2,ymF2,re,model
	foF2=7.

	rb=re+180
	rm=rb+150
      hmF2=150
      ymF2=180
*	vel=.00
*	x=(r-re-hmF2-vel*t)/ymF2
*	elden=(foF2**2/80.6d-6)*chap(1.d0,1.4142d0*x)+0.*elong*elat
      EN=(foF2**2/80.6d-6)
      if(r.le.rb)then
          elden=0
	    return
      else 
          elden=(foF2**2/80.6d-6)*(1-((r-rm)/(rm-rb))**2*(rb/r)**2)
	    return
	endif
	end


	real*8 function output()
	implicit real*8 (a-h,o-z)
	character*21 typo(14),typo2(14)
      real*8 out(11,14),out2(11,14),along1,alat1,elav1,bear1,f1,cha1,
     $tol1
	real*4 yp(301,15,2),elde(301),elde1(151),elde2(150)
	integer*4 m1,l1
      common /dat/foE,hmE,ymE,foF1,hmF1,ymF1,foF2,hmF2,ymF2,re,model
      common /grid/blhla,blhlo,rmg,dla,dlo,drg,dtz,yp,nr,nla,nlo
      common /mag/fhs,fr,xm,ym,zm,px,py,pz,ir

      re=6378.135
      blhla=40.d0
	blhlo=40.d0
	rmg=re+180.d0
	dla=5.d0
	dlo=5.d0
	drg=.5d0
	nla=5
	nlo=3
	nr=300/1+1
      dtz=0.0
*      data elde/0,	40720.0394595130,	79969.8435564882,	
*     &117755.985417885 ,	
*     &154085.008388501,	188963.426181327,	222397.723027055, 	
*     &254394.353822725,	284959.744279542,	314100.291069843, 	
*     &341822.361973237,	368132.296021916,	393036.403645144 ,	
*     &416540.966812926,	438652.239178870,	459376.446222237 ,	
*     &478719.785389192,	496688.426233256,	513288.510554969 ,	
*     &528526.152540766,	542407.438901066,	554938.429007595 ,	
*     &566125.155029921,	575973.622071240,	584489.808303385 ,	
*     &591679.665101087,	597549.117175475,	602104.062706837 ,	
*     &605350.373476634,	607293.894998766,	607940.446650124/

      data elde1/0,	8261.94323553753,	16464.8659508962,	
     & 24608.8210184525,	32693.8612625425,	40720.0394595130,
     & 48687.4083377684,	56596.0205778202,	64445.9288123350,	
     & 72237.1856261835,	79969.8435564882,	87643.9550926725,	
     & 95259.5726765079,	102816.748702164,	110315.535516254,	
     & 117755.985417885,	125138.150658707,	132462.083442956,	
     & 139727.835927510,	146935.460221927,	154085.008388501,	
     & 161176.532442306,	168210.084351246,	175185.716036099,	
     & 182103.479370568,	188963.426181327,	195765.608248070,	
     & 202510.077303555,	209196.885033657,	215826.083077409,	
     & 222397.723027055,	228911.856428092,	235368.534779324,	
     & 241767.809532900,	248109.732094369,	254394.353822725,	
     & 260621.726030451,	266791.899983568,	272904.926901684,	
     & 278960.857958036,	284959.744279542,	290901.636946844,	
     & 296786.586994356,	302614.645410309,	308385.863136802,	
     & 314100.291069843,	319757.980059399,	325358.980909442,	
     & 330903.344377993,	336391.121177170,   341822.361973237,	
     & 347197.117386644,	352515.437992079,	357777.374318509,
     & 362982.976849232,	368132.296021916,	373225.382228652,	
     & 378262.285815993,	383243.057085007,	388167.746291315,	
     & 393036.403645144,	397849.079311367,	402605.823409552,	
     & 407306.686014007,	411951.717153824,	416540.966812926,	
     & 421074.484930112,	425552.321399103,	429974.526068584,	
     & 434341.148742255,	438652.239178870,	442907.847092288,	
     & 447108.022151514,	451252.813980744,	455342.272159413,	
     & 459376.446222237,	463355.385659262,	467279.139915902,
     & 471147.758392989,	474961.290446820,	478719.785389192,	
     & 482423.292487458,	486071.860964564,	489665.539999097,	
     & 493204.378725327,	496688.426233256,	500117.731568655,	
     & 503492.343733115,	506812.311684090,	510077.684334938,
     & 513288.510554969,	516444.839169486,	519546.718959831,	
     & 522594.198663430, 	525587.326973834,	528526.152540766,	
     & 531410.723970160,	534241.089824214,	537017.298621422,	
     & 539739.398836629,   542407.438901066,	545021.467202399,
     & 547581.532084770,	550087.681848840,	552539.964751837,	
     & 554938.429007595,	557283.122786596,	559574.094216019,	
     & 561811.391379779,	563995.062318574,	566125.155029921,	
     & 568201.717468208,	570224.797544731,	572194.443127741,
     & 574110.702042482,	575973.622071240,	577783.250953381,	
     & 579539.636385398,	581242.826020949,	582892.867470903,	
     & 584489.808303385,	586033.696043812,	587524.578174942,	
     & 588962.502136912,	590347.515327284,	591679.665101087,	
     & 592958.998770854,	594185.563606675,	595359.406836229,	
     & 596480.575644831,	597549.117175475,	598565.078528873,	
     & 599528.506763502,	600439.448895641,	601297.951899415,	
     & 602104.062706837,	602857.828207854,	603559.295250379,	
     & 604208.510640344,	604805.521141734,	605350.373476634,	
     & 605843.114325264,	606283.790326030,	606672.448075558,	
     & 607009.134128738,	607293.894998766,	607526.777157186,	
     & 607707.827033930,	607837.091017360,	607914.615454309,	
     & 607940.446650124/
      data elde2/
     & 4138.35248919233,	12370.7788525891,	20544.2111380034,	
     & 28658.7021937848,	36714.3048202693,	44711.0717698261,
     & 52649.0557469072,	60528.3094080956,	68348.8853621538,	
     & 76110.8361700725,	83814.2143451183,	91459.0723528831,	
     & 99045.4626113310,	106573.437490847,	114043.049314287,	
     & 121454.350357022,	128807.392846989,	136102.228964739,	
     & 143338.910843485,	150517.490569146,	157638.020180401,	
     & 164700.551668732,	171705.136978474,	178651.828006863,	
     & 185540.676604080,	192371.734573304,	199145.053670754,	
     & 205860.685605742,	212518.682040715,	219119.094591306,	
     & 225661.974826381,	232147.374268082,	238575.344391880,	
     & 244945.936626620,	251259.202354567,	257515.192911452,	
     & 263713.959586523,	269855.553622587,	275940.026216063,	
     & 281967.428517021,	287937.811629235,	293851.226610227,	
     & 299707.724471315,	305507.356177658,	311250.172648302,	
     & 316936.224756231,	322565.563328407,	328138.239145821,	
     & 333654.302943538,	339113.805410741,   344516.797190783,	
     & 349863.328881227,	355153.451033893,	360387.214154911,	
     & 365564.668704755,	370685.865098301,	375750.853704865,	
     & 380759.684848250,	385712.408806797,	390609.075813422,	
     & 395449.736055672,	400234.439675760,	404963.236770618,	
     & 409636.177391942,	414253.311546232,	418814.689194845,	
     & 423320.360254034,	427770.374594997,	432164.782043921,	
     & 436503.632382027,	440786.975345616,   445014.860626113,	
     & 449187.337870115,	453304.456679431,	457366.266611131,	
     & 461372.817177591,	465324.157846533,	469220.338041078,	
     & 473061.407139782,	476847.414476687,	480578.409341364,	
     & 484254.440978955,	487875.558590222,	491441.811331588,	
     & 494953.248315182,	498409.918608887,	501811.871236377,	
     & 505159.155177170,	508451.819366666,	511689.912696193,	
     & 514873.484013053,	518002.582120562,	521077.255778099,	
     & 524097.553701147,	527063.524561337,	529975.216986493,	
     & 532832.679560675,	535635.960824223,	538385.109273802,	
     & 541080.173362444,   543721.201499592,	546308.242051145,	
     & 548841.343339498,	551320.553643592,	553745.921198950,	
     & 556117.494197726,   558435.320788746,	560699.449077551,	
     & 562909.927126441,	565066.802954519,	567170.124537731,	
     & 569219.939808916,	571216.296657840,	573159.242931246,	
     & 575048.826432893,	576885.094923603,	578668.096121298,	
     & 580397.877701050,	582074.487295116,	583697.972492989,	
     & 585268.380841434,	586785.759844532,	588250.156963727,	
     & 589661.619617862,	591020.195183228,	592325.930993600,	
     & 593578.874340285,	594779.072472161,	595926.572595720,	
     & 597021.421875113,	598063.667432186,   599053.356346528,	
     & 599990.535655513,	600875.252354336,	601707.553396063,	
     & 602487.485691667,	603215.096110074,	603890.431478200,	
     & 604513.538580999,	605084.464161501,	605603.254920853,	
     & 606069.957518363,	606484.618571541,	606847.284656140,	
     & 607158.002306199,	607416.818014082,	607623.778230521,	
     & 607778.929364659,	607882.317784088,	607933.989814894/

      do i=1,150
	    elde(2*i-1)=elde1(i)
	    elde(2*i)=elde2(i)
	end do
	elde(301)=elde1(151)


      do j=1,nr
	    do i=1,nla*nlo     
              yp(j,i,1)=elde(j)
              yp(j,i,2)=elde(j)
	    end do
	end do

	along1=45.d0
	alat1=45.d0
	elav1=20.d0
	bear1=3.d0
      f1=14.d0
      nos=1
	tol=1d-16
	m1=10
      l1=7


*     
      call HASEL(along1,alat1,elav1,bear1,f1,'y',nos,out,typo,tol1,
     & m1,l1)
*	print*,out
*	call HASEL(45.d0,45.d0,30.d0,10.d0,15.d0,
*     $'y',nos,out2,typo2,1d-6,7,0)
*	j=0,output
	output=out(2,1)
	return
	end 
	
	subroutine mexFunction(nlhs,plhs,nrhs,prhs)
	integer plhs(*),prhs(*),nlhs,nrhs
	integer MXCREATEDOUBLEMATRIX,mxGetPr,mxGetM,mxGetN,mxlsNumeric
	integer ml,nl,m2,n2,size
	real*8 x(100),Y(100),Z(100)
	integer x_pr,y_pr,z_pr
	ml=mxGetM(prhs(1))
	nl=mxGetN(prhs(1))
	m2=mxGetM(prhs(2))
	n2=mxGetN(prhs(2))
	size = ml*n1
	plhs(1)=MXCREATEDOUBLEMATRIX(ml,n1,0)
	x_pr=mxGetPr(prhs(1))
	y_pr=mxGetPr(prhs(2))
	z_pr=mxGetPr(plhs(1))
	call mxCopyPtrToReal8(x_pr,X,size)
	call mxCopyPtrToReal8(y_pr,Y,size)
	z=output()
	call mxCopyReal8ToPtr(z,z_pr,size)
	return
	end
	


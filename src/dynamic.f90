module sh_dynamic
  use sh_constants
  use sh_io
  use sh_parameters
  use sh_hamiltonian
  use sh_random
  
  implicit none
  contains
  
  !=============================================!
  != init coordinate and Normal mode velocitie =!
  !=============================================!
  subroutine init_coordinate_velocity(xx,vv)
    implicit none

    real(kind=dp)::xx(1:nfreem),vv(1:nfreem)
    !real(kind=dp),external::gaussian_random_number

    do ifreem=1,nfreem
      !xx(ifreem)=gaussian_random_number(0.0d0,dsqrt(kb*temp/k))
      !vv(ifreem)=gaussian_random_number(0.0d0,dsqrt(kb*temp/mass))
      xx(ifreem)=gaussian_random_number(0.0d0,dsqrt(kb*temp/(womiga(ifreem))**2))
      vv(ifreem)=gaussian_random_number(0.0d0,dsqrt(kb*temp))
    enddo
    
  end subroutine init_coordinate_velocity

  !========================================================!
  != init dynamical variable                              =!
  !========================================================!
  !=incluse inintial state in adiabatic and diabatic .    =! 
  subroutine init_dynamical_variable(xx,ee,pp,cc_elec,ww_elec,cc_hole,ww_hole)

    implicit none

    !integer ibasis
    real(kind=dp):: xx(1:nfreem),ee(1:nbasis),pp(1:nbasis,1:nbasis),flagr,flagd
    complex(kind=dpc):: cc_elec(1:nbasis),ww_elec(1:nbasis)
    complex(kind=dpc):: cc_hole(1:nbasis),ww_hole(1:nbasis)
    !integer          :: icenter_index
    
    call Get_init_WFstat(init_elec_K,init_elec_band,init_elec_WF)
    call Get_init_WFstat(init_hole_K,init_hole_band,init_hole_WF)
    
    cc_elec=  cmplx_0
    cc_hole=  cmplx_0
    
    icenter_elec=(Rcenter(3)-1)*na2site*na1site*num_wann+&
                (Rcenter(2)-1)*na1site*num_wann+(Rcenter(1)-1)*num_wann+init_elec_WF
    icenter_hole=(Rcenter(3)-1)*na2site*na1site*num_wann+&
                (Rcenter(2)-1)*na1site*num_wann+(Rcenter(1)-1)*num_wann+init_hole_WF
                
    cc_elec(icenter_elec)=1.0d0
    cc_hole(icenter_hole)=1.0d0
    !call Coulomb(Relec,elec_WF,Rhole,hole_WF,E_ehExtion)
    
    call calculate_eigen_energy_state(xx,ee,pp)
    call Add_Coulomb(Rcenter,init_elec_WF,Rcenter,init_hole_WF,ee)
    call convert_diabatic_adiabatic(pp,cc_elec,ww_elec)
    call convert_diabatic_adiabatic(pp,cc_hole,ww_hole)
    call Get_surface(icenter_elec,isurface_elec)
    call Get_surface(icenter_hole,isurface_hole)
    
    contains
    subroutine Get_surface(icenter,isurface)
      implicit none
      integer ::  icenter,isurface
      
      call random_number(flagr)
      flagd=0.0d0
      do ibasis=1,nbasis
        flagd=flagd+pp(icenter,ibasis)**2
        if(flagr.le.flagd) then
          isurface=ibasis
          exit
        endif
      enddo
      
    end subroutine Get_surface
    
  endsubroutine init_dynamical_variable
  
  subroutine Get_init_WFstat(elec_K,elec_band,elec_WF)
    implicit none
    
    integer::elec_K,elec_band,elec_WF
    integer::i
    real(kind=dp),allocatable::elecKBproj(:)
    real::sumproj,flagr
    logical :: LnotfindWF
    
    allocate(elecKBproj(0:num_wann))
    elecKBproj(1:num_wann) = bands_projs(:,elec_band,elec_K)
    elecKBproj(0) = 0.0
    sumproj = 0.0
    do i=1,num_wann
      elecKBproj(i) = elecKBproj(i) +sumproj
      sumproj = sumproj + elecKBproj(i)
    enddo
    
    call random_number(flagr)
    LnotfindWF = .TRUE.
    elec_WF = 1
    do while(LnotfindWF)
      if(flagr > elecKBproj(elec_WF-1) .and. flagr<= elecKBproj(elec_WF) ) then
        LnotfindWF = .False.
      else
        elec_WF = elec_WF + 1
      endif
    enddo
    
    deallocate(elecKBproj)
    
  end subroutine Get_init_WFstat
  
  subroutine Add_Coulomb(Relec,elec_WF,Rhole,hole_WF,ee)
    use f95_precision
    use blas95
    implicit none
    integer::Relec(3),Rhole(3),elec_WF,hole_WF
    real(kind=dp)::ee(1,nbasis)
    real(kind=dp)::Relec_hole(3)
    real(kind=dp)::Reh
    real(kind=dp)::E_ehExtion
    real(kind=dp)::a(3,3)
    real(kind=dp)::x(3)
    real(kind=dp)::y(3)
    
    x = Relec-Rhole
    a = real_lattice
    call gemv(a,x,y,trans='T')
    
    Relec_hole(:) = y(:)+Rwann(:,elec_WF)-Rwann(:,hole_WF)
    Reh = sqrt(Sum(Relec_hole**2))
    E_ehExtion = -(elem_charge_SI**2)/(fopieps0*epsr*Reh)
    ee = ee +E_ehExtion
    
  end subroutine Add_Coulomb
    
  !subroutine Get_surface()
  !=======================================================================!
  != rk4 method to obtain coordinate and velocitie after a time interval =!
  !=======================================================================!
  != ref: http://en.wikipedia.org/wiki/runge_kutta_methods               =!
  !=======================================================================!

  subroutine rk4_nuclei(pp,xx,vv,tt)
    implicit none

    real(kind=dp):: pp(1:nbasis,1:nbasis),tt
    real(kind=dp):: xx(1:nfreem)
    real(kind=dp):: vv(1:nfreem)
    real(kind=dp):: tt2,tt6
    real(kind=dp):: xx0(1:nfreem),dx1(1:nfreem),dx2(1:nfreem),dx3(1:nfreem),dx4(1:nfreem)
    real(kind=dp):: vv0(1:nfreem),dv1(1:nfreem),dv2(1:nfreem),dv3(1:nfreem),dv4(1:nfreem)

    tt2=tt/2.0d0; tt6=tt/6.0d0
    call derivs_nuclei(pp,xx,vv,dx1,dv1)
    xx0=xx+tt2*dx1; vv0=vv+tt2*dv1
    call derivs_nuclei(pp,xx0,vv0,dx2,dv2)
    xx0=xx+tt2*dx2; vv0=vv+tt2*dv2
    call derivs_nuclei(pp,xx0,vv0,dx3,dv3)
    xx0=xx+tt*dx3; vv0=vv+tt*dv3
    call derivs_nuclei(pp,xx0,vv0,dx4,dv4)
    xx=xx+tt6*(dx1+2.0d0*dx2+2.0d0*dx3+dx4)
    vv=vv+tt6*(dv1+2.0d0*dv2+2.0d0*dv3+dv4)
  endsubroutine
  
  !====================================================!
  != calculate derivative of coordinate and velocitie =!
  !====================================================!
  != ref: notebook page 630                           =!
  !====================================================!

  subroutine derivs_nuclei(pp,xx,vv,dx,dv)
  implicit none

    !integer ibasis
    real(kind=dp)::pp(1:nbasis,1:nbasis),xx(1:nfreem),vv(1:nfreem),dx(1:nfreem),dv(1:nfreem)
    real(kind=dp)::dEa_dQf
    do ifreem=1,nfreem
      dEa_dQf=0.0
      do ibasis=1,nbasis
        dEa_dQf=dEa_dQf + hep(ibasis,ibasis,ifreem)*pp(ibasis,isurface_elec)**2-&
                          hep(ibasis,ibasis,ifreem)*pp(ibasis,isurface_hole)**2
      end do
      !dv(ifreem)=(-k(ifreem)*xx(ifreem)-dEa_dpf)/mass(ifreem)-gamma*vv(ifreem)
      dv(ifreem)= -(womiga(ifreem)**2)*xx(ifreem)- dEa_dQf -gamma*vv(ifreem)
      dx(ifreem)=vv(ifreem)
    enddo
  endsubroutine derivs_nuclei
  
  subroutine derivs_electron_diabatic(xx,cc,dc,llhole)
    implicit none

    integer ::jbasis
    logical :: llhole
    real(kind=dp) xx(1:nfreem)
    complex(kind=dpc) cc(1:nbasis),dc(1:nbasis)
    dc(:)=(0.0,0.0)
    
    !call set_HHt(xx)
    HHt=HH0
    do ifreem=1,nfreem
      HHt(:,:)=HHt(:,:)+xx(ifreem)*Hep(:,:,ifreem)
    enddo
	
    do ibasis=1,nbasis
      do jbasis=1,nbasis
        dc(ibasis)=dc(ibasis)+cc(jbasis)*HHt(ibasis,jbasis)
      end do
      if(llhole) then      
        dc(ibasis) = dc(ibasis)*(cmplx_i)/hbar_SI
      else
        dc(ibasis) = dc(ibasis)*(-cmplx_i)/hbar_SI
      endif
    enddo
  endsubroutine derivs_electron_diabatic

  !===========================================================!
  != rk4 method to obtain wavefunction after a time interval =!
  !===========================================================!
  != ref: http://en.wikipedia.org/wiki/runge_kutta_methods   =!
  !===========================================================!

  subroutine rk4_electron_diabatic(xx,cc,tt,llhole)
    implicit none

    real(kind=dp):: tt,tt2,tt6
    real(kind=dp):: xx(1:nfreem)
    complex(kind=dpc):: cc(1:nbasis),cc0(1:nbasis),dc1(1:nbasis),&
                        dc2(1:nbasis),dc3(1:nbasis),dc4(1:nbasis)
    logical ::  llhole
    tt2=tt/2.0d0; tt6=tt/6.0d0
    call derivs_electron_diabatic(xx,cc,dc1,llhole)
    cc0=cc+tt2*dc1
    call derivs_electron_diabatic(xx,cc0,dc2,llhole)
    cc0=cc+tt2*dc2
    call derivs_electron_diabatic(xx,cc0,dc3,llhole)
    cc0=cc+tt*dc3
    call derivs_electron_diabatic(xx,cc0,dc4,llhole)
    cc=cc+tt6*(dc1+2.0d0*dc2+2.0d0*dc3+dc4)
    
  endsubroutine rk4_electron_diabatic
  
  !=================================!
  != calculate hopping probability =!
  !=================================!
  != ref: notebook page 631        =!
  !=================================!

  subroutine calculate_hopping_probability(ww,vv,dd,tt,gg,gg1,isurface,llhole)
    implicit none

    integer ::jbasis,isurface
    logical ::llhole
    real(kind=dp):: vv(1:nfreem),dd(1:nbasis,1:nbasis,1:nfreem),&
                    gg(1:nbasis),gg1(1:nbasis),tt,sumvd
    complex(kind=dpc):: ww(1:nbasis)

    gg=0.0d0
    gg1=0.0d0
    do ibasis=1,nbasis
      if(ibasis /= isurface) then
        sumvd=0.0d0
        do ifreem=1,nfreem
          sumvd=sumvd+vv(ifreem)*dd(isurface,ibasis,ifreem)
        enddo
        !!when is the dynamic of hole ,the dijk while be dijk*(-1)
        if(llhole) sumvd = -sumvd
        
        gg(ibasis)=2.0d0*tt*real(conjg(ww(isurface))*ww(ibasis))*&
                    sumvd/real(conjg(ww(isurface))*ww(isurface))
        gg1(ibasis)=gg(ibasis)
        if(gg(ibasis) < 0.0d0) gg(ibasis)=0.0d0
      
      endif
    enddo
    
  endsubroutine calculate_hopping_probability
  
  !===========================!
  != nonadiabatic transition =!
  !===========================!
  != ref: notebook page 635  =!
  !===========================!
  
  subroutine calculate_sumg(sumg0,sumg1,w0,w,isurface,minde,llhole)
    implicit none
      real(kind=dp)::sumg0,sumg1,minde
      logical :: llhole
      integer :: isurface
      complex(kind=dpc) :: w0(nbasis),w(nbasis)
      
      sumg0=(abs(w0(isurface))**2-abs(w(isurface))**2)/abs(w0(isurface))**2
      if(llhole) then
        sumg1=sum(g1_hole)
      else 
        sumg1=sum(g1_elec)
      endif
      
      if(isurface == 1) then
        minde= e0(isurface+1)-e0(isurface)
      elseif(isurface == nbasis) then
        minde= e0(isurface)-e0(isurface-1)
      elseif((e0(isurface+1)-e0(isurface)) < (e0(isurface)-e0(isurface-1))) then
        minde=(e0(isurface+1)-e0(isurface))
      else
        minde= e0(isurface)-e0(isurface-1)
      endif    
      if(llhole) minde = -minde
      
  end subroutine calculate_sumg
  
  subroutine calculate_g(isurface,e0,sumg0,g1,g)
    implicit none
    integer :: isurface
    real(kind=dp) :: e0(nbasis),sumg0,g(nbasis),g1(nbasis)
    
      if(isurface == 1) then
        ibasis=isurface+1
      elseif(isurface == nbasis) then
        ibasis=isurface-1
      elseif((e0(isurface+1)-e0(isurface)) < (e0(isurface)-e0(isurface-1))) then
        ibasis=isurface+1
      else
        ibasis=isurface-1
      endif
      g(ibasis)=sumg0-(sum(g1)-g1(ibasis))
      if(g(ibasis) < 0.0d0) g(ibasis)=0.0d0
      if(sum(g) > 1.0d0) g=g/sum(g) 
      
  end subroutine calculate_g
  
  subroutine nonadiabatic_transition(ee,pp,dd,isurface,gg,ww,vv,cc,llhole)
  implicit none

    !integer ibasis,jbasis
    integer :: isurface
    real(kind=dp)::ee(1:nbasis),pp(1:nbasis,1:nbasis),gg(1:nbasis),&
                   dd(1:nbasis,1:nbasis,1:nfreem),vv(1:nfreem)
    real(kind=dp)::sumvd,sumdd,sumgg,flagr,flagd
    complex(kind=dpc):: ww(1:nbasis),cc(1:nbasis)
    logical :: llhole

    call more_random()
    call random_number(flagr)
    sumgg=0.0d0
    do ibasis=1,nbasis
      if(ibasis /= isurface) then
        sumgg=sumgg+gg(ibasis)
        if(flagr < sumgg) then
          sumvd=0.0d0
          sumdd=0.0d0
          do ifreem=1,nfreem
            !sumvd=sumvd+vv(ifreem)*dd(isurface,ibasis,ifreem)
            !sumdd=sumdd+dd(isurface,ibasis,ifreem)**2
            sumvd=sumvd+vv(ifreem)*dd(isurface,ibasis,ifreem)
            sumdd=sumdd+dd(isurface,ibasis,ifreem)**2
          enddo
          if(llhole) sumvd = - sumvd
          
          !flagd=1.0d0+2.0d0*(ee(isurface)-ee(ibasis))*sumdd/mass/sumvd**2
          if(llhole) then
            flagd=1.0d0-2.0d0*(ee(isurface)-ee(ibasis))*sumdd/sumvd**2
          else
            flagd=1.0d0+2.0d0*(ee(isurface)-ee(ibasis))*sumdd/sumvd**2
          end if
          
          if(flagd >= 0.0d0) then
            flagd=(sumvd/sumdd)*(-1.0d0+dsqrt(flagd))
            do ifreem=1,nfreem
              vv(ifreem)=vv(ifreem)+flagd*dd(isurface,ibasis,ifreem)
            enddo
            isurface=ibasis
          endif

          exit
        endif
      endif
    enddo
  endsubroutine nonadiabatic_transition
  
  !===============================================!
  != add bath effect to coordinate and velocitie =!
  !===============================================!
  != ref: notebook page 462 and 638              =!
  !===============================================!

  subroutine add_bath_effect(ee,dd,pp,tt,xx,vv)
  implicit none

    integer ik1site,ik2site
    real(kind=dp):: ee(1:nbasis),dd(1:nbasis,1:nbasis,1:nfreem),&
                    pp(1:nbasis,1:nbasis),tt,xx(1:nfreem),vv(1:nfreem)
    real(kind=dp):: sigmar
    real(kind=dp):: kk,r1,r2,r3,r4,z1,z2,z3,z4
    !real(kind=dp),external::gaussian_random_number_fast

    !sigmar=dsqrt(2.0d0*kb*temp*gamma*tt/mass) !gamma -> 1/t
    sigmar=dsqrt(2.0d0*kb*temp*gamma*tt)  !dangwei yu dp/dt   (mv) yiyang.
    do ifreem=1,nfreem
      kk=womiga(ifreem)**2
      do ibasis=1,nbasis
        !calculate d2Ei/dx(ifreem)2 
        if(ibasis /= isurface_elec) then
          do ik1site=1,nbasis
            do ik2site=1,nbasis
              kk=kk+2.0d0*dd(ibasis,isurface_elec,ifreem)*pp(ik1site,isurface_elec)*&
                  pp(ik2site,ibasis)*hep(ik1site,ik2site,ifreem)
            enddo
          enddo  
        endif
        if(ibasis /= isurface_hole) then
          do ik1site=1,nbasis
            do ik2site=1,nbasis
              kk=kk-2.0d0*dd(ibasis,isurface_hole,ifreem)*pp(ik1site,isurface_hole)*&
                  pp(ik2site,ibasis)*hep(ik1site,ik2site,ifreem)
            enddo
          enddo  
        endif
        
      end do

      r1=gaussian_random_number_fast(0.0d0,sigmar)
      r2=gaussian_random_number_fast(0.0d0,sigmar)
      r3=gaussian_random_number_fast(0.0d0,sigmar)
      r4=gaussian_random_number_fast(0.0d0,sigmar)
      z1=r1                                                     !dp/dt
      z2=tt*(r1/2.0d0+r2/sqrt3/2.0d0)                           !dp
      z3=tt**2*(r1/6.0d0+r2*sqrt3/12.0d0+r3/sqrt5/12.0d0)       !dpdt
      z4=tt**3*(r1/24.0d0+r2*sqrt3/40.0d0+r3/sqrt5/24.0d0+r4/sqrt7/120.0d0) !dp*(dt**2)
      xx(ifreem)=xx(ifreem)+(z2-gamma*z3+(-womiga(ifreem)**2+gamma**2)*z4)
      vv(ifreem)=vv(ifreem)+(z1-gamma*z2+(-womiga(ifreem)**2+gamma**2)*z3+&
                  (2.0d0*gamma*(womiga(ifreem)**2) - gamma**3 )*z4 )
    
    enddo
    
  endsubroutine add_bath_effect
  
  subroutine saveresult()
  implicit none
    character(len=maxlen) ::  pes_name,csit_name,wsit_name,xsit_name,psit_name,ksit_name
    integer               ::  pes_unit,csit_unit,wsit_unit,xsit_unit,psit_unit,ksit_unit
    character(len=maxlen) ::  inf_name,msd_name,msds_name,ipr_name
    integer               ::  inf_unit,msd_unit,msds_unit,ipr_unit
    
    pes_unit = io_file_unit()
    pes_name = 'pes.out'
    call open_file(pes_name,pes_unit)
    do iaver=1,1
      do isnap=1,nsnap
        write(pes_unit,'(999999e12.5)') dt*nstep*isnap,(pes(ibasis,isnap,iaver),ibasis=-1,nbasis)
      enddo
    enddo
    call close_file(pes_name,pes_unit)
   
    inf_unit = io_file_unit()
    inf_name = 'inf_elec.out'
    call open_file(inf_name,inf_unit)
    do iaver=1,1
      do isnap=1,nsnap
        write(inf_unit,'(999999e12.5)') dt*nstep*isnap,(inf_elec(ibasis,isnap,iaver),ibasis=1,3)
      enddo
    enddo
    call close_file(inf_name,inf_unit)
    
    inf_unit = io_file_unit()
    inf_name = 'inf_hole.out'
    call open_file(inf_name,inf_unit)
    do iaver=1,1
      do isnap=1,nsnap
        write(inf_unit,'(999999e12.5)') dt*nstep*isnap,(inf_hole(ibasis,isnap,iaver),ibasis=1,3)
      enddo
    enddo
    call close_file(inf_name,inf_unit)
    
    csit_unit = io_file_unit()
    csit_name = 'csit_elec.out'
    call open_file(csit_name,csit_unit)
    do isnap=1,nsnap
      write(csit_unit,'(999999e12.5)') dt*nstep*isnap,(csit_elec(ibasis,isnap),ibasis=1,nbasis)
    enddo
    call close_file(csit_name,csit_unit)
    
    csit_unit = io_file_unit()
    csit_name = 'csit_hole.out'
    call open_file(csit_name,csit_unit)
    do isnap=1,nsnap
      write(csit_unit,'(999999e12.5)') dt*nstep*isnap,(csit_hole(ibasis,isnap),ibasis=1,nbasis)
    enddo
    call close_file(csit_name,csit_unit)
    
    wsit_unit = io_file_unit()
    wsit_name = 'wsit_elec.out'
    call open_file(wsit_name,wsit_unit)
    do isnap=1,nsnap
      write(wsit_unit,'(999999e12.5)') dt*nstep*isnap,(wsit_elec(ibasis,isnap),ibasis=1,nbasis)
    enddo
    call close_file(wsit_name,wsit_unit)

    wsit_unit = io_file_unit()
    wsit_name = 'wsit_hole.out'
    call open_file(wsit_name,wsit_unit)
    do isnap=1,nsnap
      write(wsit_unit,'(999999e12.5)') dt*nstep*isnap,(wsit_hole(ibasis,isnap),ibasis=1,nbasis)
    enddo
    call close_file(wsit_name,wsit_unit)
    
    psit_unit = io_file_unit()
    psit_name = 'psit_elec.out'
    call open_file(psit_name,psit_unit)
    do isnap=1,nsnap
      write(psit_unit,'(999999e12.5)') dt*nstep*isnap,(psit_elec(ibasis,isnap),ibasis=1,nbasis)
    enddo
    call close_file(psit_name,psit_unit)
    
    psit_unit = io_file_unit()
    psit_name = 'psit_hole.out'
    call open_file(psit_name,psit_unit)
    do isnap=1,nsnap
      write(psit_unit,'(999999e12.5)') dt*nstep*isnap,(psit_hole(ibasis,isnap),ibasis=1,nbasis)
    enddo
    call close_file(psit_name,psit_unit)
    
    xsit_unit = io_file_unit()
    xsit_name = 'xsit.out'
    call open_file(xsit_name,xsit_unit)
    do isnap=1,nsnap
      write(xsit_unit,'(999999e12.5)') dt*nstep*isnap,(xsit(ibasis,isnap),ibasis=1,nbasis)
    enddo
    call close_file(xsit_name,xsit_unit)
    
    ksit_unit = io_file_unit()
    ksit_name = 'ksit.out'
    call open_file(ksit_name,ksit_unit)
    do isnap=1,nsnap
      write(ksit_unit,'(999999e12.5)') dt*nstep*isnap,(ksit(ibasis,isnap),ibasis=1,nbasis),kb*temp/2*au2ev
    enddo
    call  close_file(ksit_name,ksit_unit)
    
    msd_unit = io_file_unit()
    msd_name = 'msd.out'
    call open_file(msd_name,msd_unit)
    do isnap=1,nsnap
      write(msd_unit,'(999999e12.5)') dt*nstep*isnap*au2fs,msd(isnap)
    enddo
    call close_file(msd_name,msd_unit)
   
    ipr_unit = io_file_unit()
    ipr_name = 'ipr_elec.out'
    call open_file(ipr_name,ipr_unit)
    do isnap=1,nsnap
      write(ipr_unit,'(999999e12.5)') dt*nstep*isnap,ipr_elec(isnap)
    enddo
    call close_file(ipr_name,ipr_unit)
    
    ipr_unit = io_file_unit()
    ipr_name = 'ipr_hole.out'
    call open_file(ipr_name,ipr_unit)
    do isnap=1,nsnap
      write(ipr_unit,'(999999e12.5)') dt*nstep*isnap,ipr_hole(isnap)
    enddo
    call close_file(ipr_name,ipr_unit)
    
    msds_unit = io_file_unit()
    msds_name = 'msds.out'
    call open_file(msds_name,msds_unit)
    do isnap=1,nsnap
      write(msds_unit,'(999999e12.5)') dt*nstep*isnap*au2fs,(msds(isnap,iaver),iaver=1,naver)
    enddo
    call close_file(msds_name,msds_unit)
  
  end subroutine saveresult
  
end module
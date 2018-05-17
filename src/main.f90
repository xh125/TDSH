!=====================================================================================================!
!= this program is used to simulate charge transport with the self-consistent surface hopping method =!
!=====================================================================================================!
!= of 2(or 1 and 3)-dimensional materials stack with local and unlocal electron-phonon couplings and =!
!=====================================================================================================!
!=system-bath interactions by xiehua at department of physic, USTC;xh125@mail.ustc.edu.cn            =!     
!=====================================================================================================!
!= Last updata 2018-5.16 Version:0.2.5                                                                !   
!                                                                                                     !                                                                                                     !
!=====================================================================================================!

!===========!
!= program =!
!===========!

program SCSH
  use omp_lib
  use sh_constants
  use sh_parameters
  use sh_io
  use sh_utility
  use sh_hamiltonian
  use sh_random
  use sh_dynamic
  
  implicit none

  !===============!
  != preparation =!
  !===============!
	
  time0=io_time()
  stdout = io_file_unit()
  open(unit=stdout,file="SCSH.out")
  call io_date(cdate,ctime)
  write(stdout,*) 'SCSH :Execution started on ',cdate,' at ',ctime
  
  call read_parameters ()
  call treat_parameters()
  call set_HH0()
  call set_Hep()
  call init_random_seed()  
  
  !==========================!
  != loop over realizations =!
  !==========================!
  do iaver=1,naver
    write(stdout,'(a,i4.4,a)') '###### iaver=',iaver,' ######'
    !==================!
    != initialization =!
    !==================!

    call init_coordinate_velocity(Q,Vq)
    call init_dynamical_variable(Q,e,p,c_elec,w_elec,c_hole,w_hole)
    call calculate_nonadiabatic_coupling(e,p,d)
    Q0=Q; Vq0=Vq; e0=e; p0=p; d0=d; w0_elec=w_elec; w0_hole=w_hole

    !=======================!
    != loop over snapshots =!
    !=======================!

    do isnap=1,nsnap
      do istep=1,nstep

        !==========================!
        != update x,v,c,e,p,d,w,g =!
        !==========================!

        call rk4_nuclei(p0,Q,Vq,dt)
        
        call rk4_electron_diabatic(Q0,c_elec,dt,Llhole=.False.)
        call rk4_electron_diabatic(Q0,c_hole,dt,Llhole=.TRUE. )
        call calculate_eigen_energy_state(Q,e,p)
        call calculate_nonadiabatic_coupling(e,p,d)
        call convert_diabatic_adiabatic(p,c_elec,w_elec)
        call convert_diabatic_adiabatic(p,c_hole,w_hole)
        call calculate_hopping_probability(w0_elec,Vq0,d0,dt,g_elec,g1_elec,isurface_elec,Llhole=.False.)
        call calculate_hopping_probability(w0_hole,Vq0,d0,dt,g_hole,g1_hole,isurface_hole,Llhole=.True.)

        !===============================!
        != calculate sumg0,sumg1,minde =!
        !===============================!
        call calculate_sumg(sumg0_elec,sumg1_elec,w0_elec,w_elec,isurface_elec,minde_elec,llhole=.False.)
        call calculate_sumg(sumg0_hole,sumg1_hole,w0_hole,w_hole,isurface_hole,minde_hole,llhole=.TRUE.)

        !===================================!
        != change potential energy surface =!
        !===================================!
        call calculate_g(isurface_elec,e0,sumg0_elec,g1_elec,g_elec)
        call calculate_g(isurface_hole,e0,sumg0_hole,g1_hole,g_hole)

        call nonadiabatic_transition(e0,p0,d0,isurface_elec,g_elec,w_elec,Vq,c_elec,llhole=.False.)
        call nonadiabatic_transition(e0,p0,d0,isurface_hole,g_hole,w_hole,Vq,c_hole,llhole=.TRUE.)
        
        !===================!
        != add bath effect =!
        !===================!

        call add_bath_effect(e0,d0,p0,dt,Q,Vq)

        !============================!
        != reset dynamical variable =!
        !============================!

        Q0=Q; Vq0=Vq; E0=E; p0=p; d0=d; w0_elec=w_elec;w0_hole=w_hole
      enddo

      !=====================!
      != store information =!
      !=====================!

      pes(-1,isnap,iaver)=e(isurface_elec)
      pes( 0,isnap,iaver)=e(isurface_hole)
      inf_elec(1,isnap,iaver)=sumg0_elec
      inf_elec(2,isnap,iaver)=sumg1_elec
      inf_elec(3,isnap,iaver)=minde_elec
      inf_hole(1,isnap,iaver)=sumg0_hole
      inf_hole(2,isnap,iaver)=sumg1_hole
      inf_hole(3,isnap,iaver)=minde_hole      
      flagd_elec=0.0d0
      flagd_hole=0.0d0
      do ibasis=1,nbasis
        csit_elec(ibasis,isnap)=csit_elec(ibasis,isnap)+abs(c_elec(ibasis))**2
        csit_hole(ibasis,isnap)=csit_hole(ibasis,isnap)+abs(c_hole(ibasis))**2
        wsit_elec(ibasis,isnap)=wsit_elec(ibasis,isnap)+abs(w_elec(ibasis))**2
        wsit_hole(ibasis,isnap)=wsit_hole(ibasis,isnap)+abs(w_hole(ibasis))**2
        psit_elec(ibasis,isnap)=psit_elec(ibasis,isnap)+p(ibasis,isurface_elec)**2
        psit_hole(ibasis,isnap)=psit_hole(ibasis,isnap)+p(ibasis,isurface_hole)**2
        !ksit(ifreem,isnap)=ksit(ifreem,isnap)+0.5d0*mass(ifreem)*v(ifreem)**2
        pes(ibasis,isnap,iaver)=e(ibasis)
        !msds(isnap,iaver)=msds(isnap,iaver)+p(ibasis,isurface)**2*(ibasis-icenter_index)**2
        flagd_elec=flagd_elec+p(ibasis,isurface_elec)**4
        flagd_hole=flagd_hole+p(ibasis,isurface_hole)**4
      enddo
      
      do ifreem =1 , nfreem
        !!平均核构型
        xsit(ifreem,isnap)=xsit(ifreem,isnap)+q(ifreem)
        !!平均动能
        ksit(ifreem,isnap)=ksit(ifreem,isnap)+0.5d0*Vq(ifreem)**2
      enddo
      
      ipr_elec(isnap)=ipr_elec(isnap)+1/flagd_elec
      ipr_hole(isnap)=ipr_hole(isnap)+1/flagd_hole
      !msd(isnap)=msd(isnap)+msds(isnap,iaver)
    enddo
  
  enddo
  
  csit_elec=csit_elec/naver
  csit_hole=csit_hole/naver
  wsit_elec=wsit_elec/naver
  wsit_hole=wsit_hole/naver
  psit_elec=psit_elec/naver
  psit_hole=psit_hole/naver
  xsit=xsit/naver
  ksit=ksit/naver
  !msd=msd/naver
  ipr_elec=ipr_elec/naver
  ipr_hole=ipr_hole/naver

  !====================!
  != save information =!
  !====================!
  call saveresult()
  call cpu_time(t1)
  write(6,'(a,f10.2,a)') 'total time is',(t1-t0)/3600,'hours'
endprogram



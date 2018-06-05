!=====================================================================================================!
!= this program is --*used to simulate exciton transport with the self-consistent surface hopping method=!
!=====================================================================================================!
!= of 2dimensional materials stack with local and unlocal electron-phonon couplings and              =!
!=====================================================================================================!
!= elec-hole coulomb interactions and system interaction with environment by system-bath interactions=!     
!=====================================================================================================!
!= Last updata 2018-5.24 Version:0.2.6                                                                !   
!= Developed by xiehua at department of physic, USTC;xh125@mail.ustc.edu.cn                          =!
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
	
  time0   =io_time()
  stdout  = io_file_unit()
  open(unit=stdout,file="SCSH.out")
  call io_date(cdate,ctime)
  write(stdout,*) 'SCSH :Execution started on ',cdate,' at ',ctime
  
  call read_parameters ()
  call read_TB_parameters(num_wann,nfreem,HmnR_Tij_0,HmnR_Tij_ep)
  call init_random_seed()  
  
  !==========================!
  != loop over realizations =!
  !==========================!
  do iaver=1,naver
    write(stdout,'(a,i4.4,a)') '###### iaver=',iaver,' ######'
    !==================!
    != initialization =!
    !==================!
    
    !!得到简正坐标的初始位置和速度
    call init_coordinate_velocity(Q,Vq) 
    !!得到体系初始时刻的状态
    call init_dynamical_variable(Q,e_elec,p_elec,e_hole,p_hole,c_elec,w_elec,c_hole,w_hole)
    !!计算电子和空穴的非绝热耦合项大小，以及绝热表象下 dEij_dQq
    call calculate_nonadiabatic_coupling(e_elec,p_elec,d_elec)
    call calculate_nonadiabatic_coupling(e_hole,p_hole,d_hole)
    Q0=Q; Vq0=Vq
    e0_elec=e_elec; p0_elec=p_elec; d0_elec=d_elec; w0_elec=w_elec
    e0_hole=e_hole; p0_hole=p_hole; d0_hole=d_hole; w0_hole=w_hole
    !=======================!
    != loop over snapshots =!
    !=======================!

    do isnap=1,nsnap
      do istep=1,nstep

        !==========================!
        != update x,v,c,e,p,d,w,g =!
        !==========================!
        !rk4方法数组计算核的动力学，得到新的Q
        call rk4_nuclei(Q,Vq,p0_elec,p0_hole,d0_elec,d0_hole,dt)
        !rk4方法计算电子空穴在透热表象下的动力学，
        !得到dt时间后新的c_elec,n_ecle和c_hole,n_hole
        call rk4_electron_diabatic(Q0,c_elec,n_elec,dt,h_elec)
        call rk4_electron_diabatic(Q0,c_hole,n_hole,dt,h_hole)
        !在新的简正坐标Q和新的电荷密度n_elec和空穴密度n_hole下
        !对角化哈密顿量得到电子和空穴的能量本征态与本征值
        call calculate_eigen_energy_state(Q,n_elec,n_hole)
        !call calculate_eigen_energy_state(Q,e_hole,p_hole)
        !计算电子和空穴的非绝热耦合项
        call calculate_nonadiabatic_coupling(e_elec,p_elec,d_elec)
        call calculate_nonadiabatic_coupling(e_hole,p_hole,d_hole)
        !call calculate_nonadiabatic_coupling(e,p,d)
        !将透热表象下的电子和空穴的态转换到绝热表象下
        call convert_diabatic_adiabatic(p_elec,c_elec,w_elec)
        call convert_diabatic_adiabatic(p_hole,c_hole,w_hole)
        !计算电子和空穴的跃迁几率
        call calculate_hopping_probability(w0_elec,Vq0,d0_elec,dt,g_elec,g1_elec,isurface_elec)
        call calculate_hopping_probability(w0_hole,Vq0,d0_hole,dt,g_hole,g1_hole,isurface_hole)

        !===============================!
        != calculate sumg0,sumg1,minde =!
        !===============================!
        sumg1_elec=SUM(g1_elec)
        sumg1_hole=SUM(g1_hole)
        call calculate_sumg(sumg0_elec,E0_elec,w0_elec,w_elec,isurface_elec,minde_elec)
        call calculate_sumg(sumg0_hole,E0_hole,w0_hole,w_hole,isurface_hole,minde_hole)

        !===================================!
        != change potential energy surface =!
        !===================================!
        call calculate_g(isurface_elec,e0_elec,sumg0_elec,g1_elec,g_elec)
        call calculate_g(isurface_hole,e0_hole,sumg0_hole,g1_hole,g_hole)
        !!feedback
        call nonadiabatic_transition(e0_elec,p0_elec,d0_elec,isurface_elec,g_elec,w_elec,Vq,c_elec,llhole=.False.)
        call nonadiabatic_transition(e0_hole,p0_hole,d0_hole,isurface_hole,g_hole,w_hole,Vq,c_hole,llhole=.TRUE.)
        
        !===================!
        != add bath effect =!
        !===================!
        
        call add_bath_effect(d0_elec,p0_elec,d0_hole,p0_hole,dt,Q,Vq)

        !============================!
        != reset dynamical variable =!
        !============================!

        Q0=Q; Vq0=Vq; E0_elec=E_elec; p0_elec=p_elec; d0_elec=d_elec; w0_elec=w_elec;w0_hole=w_hole
        E0_hole=E_hole; p0_hole=p_hole; d0_hole=d_hole
      enddo

      !=====================!
      != store information =!
      !=====================!

      pes_exciton(-1,isnap,iaver)=e_elec(isurface_elec)
      pes_exciton( 0,isnap,iaver)=e_hole(isurface_hole)
      !pes_exciton(1:na1site*na2site*Num_occupied,isnap,iaver)=&
          !pes_hole(1:na1site*na2site*Num_occupied,isnap,iaver)
      !pes_exciton(na1site*na2site*Num_occupied+1:nbasis,isnap,iaver)=&
          !pes_elec(na1site*na2site*Num_occupied+1:nbasis,isnap,iaver)
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
        psit_elec(ibasis,isnap)=psit_elec(ibasis,isnap)+p_elec(ibasis,isurface_elec)**2
        psit_hole(ibasis,isnap)=psit_hole(ibasis,isnap)+p_hole(ibasis,isurface_hole)**2
        !ksit(ifreem,isnap)=ksit(ifreem,isnap)+0.5d0*mass(ifreem)*v(ifreem)**2
        pes_elec(ibasis,isnap,iaver)=e_elec(ibasis)
        pes_hole(ibasis,isnap,iaver)=e_hole(ibasis)
        !msds(isnap,iaver)=msds(isnap,iaver)+p(ibasis,isurface)**2*(ibasis-icenter_index)**2
        flagd_elec=flagd_elec+p_elec(ibasis,isurface_elec)**4
        flagd_hole=flagd_hole+p_hole(ibasis,isurface_hole)**4
      enddo
      pes_exciton(1:na1site*na2site*Num_occupied,isnap,iaver)=&
          pes_hole(1:na1site*na2site*Num_occupied,isnap,iaver)
      pes_exciton(na1site*na2site*Num_occupied+1:nbasis,isnap,iaver)=&
          pes_elec(na1site*na2site*Num_occupied+1:nbasis,isnap,iaver)
      
      
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



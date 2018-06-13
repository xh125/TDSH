module sh_hamiltonian
  use sh_parameters
  use sh_constants
  use sh_io,only : io_error,io_file_unit,open_file,close_file
  
  implicit none
  
  integer :: nrpts,irpts!n_wann,m_wann
  integer :: n, m
  integer :: ir1, ir2, ir3 
  real(kind=dp) :: ReH, ImH
  character(len=maxlen) :: Hr_name
  integer               :: Hr_unit  
  integer,allocatable   :: ndegen(:)
  ! real and imag of HmnR  
  
  contains
  
  subroutine read_TB_parameters(nnum_wann,nnfreem,HHmnR_Tij_0,HHmnR_Tij_ep)
    implicit none
    !read parameter of HmnR_Tij and HmnR_Tij_ep
    integer,intent(in)::nnum_wann,nnfreem
    real(kind=dp),intent(inout)::HHmnR_Tij_0(nnum_wann,nnum_wann,-1:1,-1:1)
    real(kind=dp),intent(inout)::HHmnR_Tij_ep(nnum_wann,nnum_wann,-1:1,-1:1,nnfreem)
    
    !!!!parameter in set HmnR_Tij_ep
    character(len=maxlen) :: epcfile_name,HmnRdQ_noma_name,eph_noma_name,&
                            TijdQ_noma_name,TiidQ_noma_name,HmnR_Tij_name
    integer               :: epcfile_unit,HmnRdQ_noma_unit,eph_noma_unit,&
                             TijdQ_noma_unit,TiidQ_noma_unit,HmnR_Tij_unit
    integer               :: imode,nmode,idQ,ndQ,NHmn,i,j
    integer,allocatable   :: irvec(:,:,:,:)
    real(kind=dp)         :: ldQ
    character(len=maxlen) :: ctmpmode,ctmpldQ
    complex(kind=dpc),allocatable::HmnR_mode_dQ(:,:,:,:,:),deta_HmnR(:)
    real(kind=dp)         ::dHmn_dQ
    character(len=16),allocatable::charwf(:),charTij(:,:)   
    allocate(charwf(num_wann),charTij(num_wann,num_wann))
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    HmnR_Tij_0  = 0.0d0
    HmnR_Tij_ep = 0.0d0
    adj_Tij     = .False.
    inquire(file="./Tij_parameter/Tij_0",exist=lexist)
    if( .Not. lexist ) then
      !set HmnR_Tij
      Hr_name  ="./wannier/wannier90_hr.dat"
      Hr_unit  = io_file_unit()    
      call open_file(Hr_name,Hr_unit)    
    
      read(Hr_unit, *) ctmp
      read(Hr_unit, *) num_wann
      read(Hr_unit, *) nrpts
      allocate(ndegen(nrpts))    
  
      read(Hr_unit, *) (ndegen(irpts), irpts=1, nrpts)
      !read <m0|H|nR>
      !in wannier90_hr.dat R,m,n,HR,HI
      do irpts=1, nrpts       !R
        do n_wann=1, num_wann   !n
          do m_wann=1, num_wann  !m
            read(Hr_unit,*) ir1, ir2, ir3, m, n, ReH, ImH
            if( abs(ir1)<=1 .and. abs(ir2)<=1 ) then
              HmnR_Tij_0(m,n,ir1,ir2) = ReH
              !<m0|H|nR>
              if(abs(ReH) >= 0.04) then
                adj_Tij(m,n,ir1,ir2) = .TRUE.
              endif
            end if
          enddo
        enddo
      enddo        
      call close_file(Hr_name,Hr_unit)
      !HmnR_Tij_0=HmnR_Tij_0/AU2EV     
      
      HmnR_Tij_name = "./Tij_parameter/Tij_0"
      HmnR_Tij_unit = io_file_unit()
      call open_file(HmnR_Tij_name,HmnR_Tij_unit)
      write(HmnR_Tij_unit,"(4A5,A12)") " ir1 "," ir2 ","  m  ","  n  "," ReH "     
      do ir2=-1,1
        do ir1=-1,1
          do n=1,num_wann
            do m=1,num_wann
              write(HmnR_Tij_unit,"(4I5,F12.6)") ir1,ir2,m,n,HmnR_Tij_0(m,n,ir1,ir2)
            enddo
          enddo
        enddo
      enddo
      call close_file(HmnR_Tij_name,HmnR_Tij_unit)
    
      nmode = nfreem
      ndQ   = 2*nshiftstep+1
      allocate(deta_HmnR(ndQ))
      !!deta_HmnR(:)=HmnR_mode_dQ(m_wann,n_wann,irpts,:,imode)-  &
      !!HmnR_mode_dQ(m_wann,n_wann,irpts,nshiftstep+1,imode)
    
      allocate(HmnR_mode_dQ(num_wann,num_wann,nrpts,ndQ,nmode))
      allocate(irvec(3,nrpts,ndQ,nmode))
      !对于每一个normal mode的处理
      do imode=1,nmode
        write(ctmpmode,*) imode
        deta_HmnR = 0.0d0
        !读取每一个mode下的各个shift下的wannier90_hr.dat文件的H的实部，
        !得到TB参数随该mode上偏移的变化,并写出该mode的phonon与电子之间的电声耦合强度，
        !包括局域与非局域
        do idQ=1,ndQ
          ldQ=(idQ-nshiftstep-1)*dtadQ
          write(ctmpldQ,"(F8.4)") ldQ
        
          Hr_name= "./nomashift/noma_"//trim(adjustl(ctmpmode))//&
                  "/shift_"//trim(adjustl(ctmpldQ))//"/wannier90_hr.dat"
          Hr_unit=io_file_unit()
          call open_file(Hr_name,Hr_unit)        
          read(Hr_unit, *) ctmp
          read(Hr_unit, *) num_wann
          read(Hr_unit, *) nrpts      
          read(Hr_unit, *) (ndegen(irpts), irpts=1, nrpts)
          !read <m0|H|nR>
          !in wannier90_hr.dat R,m,n,HR,HI
          do irpts=1, nrpts       !R
            do n_wann=1, num_wann   !n
              do m_wann=1, num_wann  !m
                read(Hr_unit,*) ir1, ir2, ir3, m, n, ReH, ImH
                !write(*,'(5i5,2f10.5)')i1, i2, i3, i4, i5, r1, r2
                HmnR_mode_dQ(m_wann,n_wann,irpts,idQ,imode)= dcmplx(ReH, ImH)
                if(idQ == ndQ .and. abs(ir1)<=1 .and. abs(ir2)<=1 ) then                
                  HmnR_Tij_ep(m,n,ir1,ir2,imode) = (ReH- &
                  real(HmnR_mode_dQ(m_wann,n_wann,irpts,1,imode)))/(2.0*real(nshiftstep)*dtadQ)
                endif
              enddo
            enddo
            irvec(1, irpts,idQ,imode)=ir1
            irvec(2, irpts,idQ,imode)=ir2
            irvec(3, irpts,idQ,imode)=ir3
          enddo
          call close_file(Hr_name,Hr_unit)        
        enddo
      
        HmnR_Tij_name = "./Tij_parameter/Tij_"//trim(adjustl(ctmpmode))
        HmnR_Tij_unit = io_file_unit()
        call open_file(HmnR_Tij_name,HmnR_Tij_unit)
        write(HmnR_Tij_unit,"(4A5,A12)") " ir1 "," ir2 ","  m  ","  n  "," ReH "
          do ir2=-1,1
            do ir1=-1,1
              do n=1,num_wann
                do m=1,num_wann
                  write(HmnR_Tij_unit,"(4I5,F12.6)") ir1,ir2,m,n,HmnR_Tij_ep(m,n,ir1,ir2,imode)
                enddo
              enddo
            enddo
          enddo    
        call close_file(HmnR_Tij_name,HmnR_Tij_unit)
        
        
      do n_wann=1,num_wann
        write(ctmp,*) n_wann
        charwf(n_wann) = "      WF_"//trim(adjustl(ctmp))
      enddo
      
      do n_wann=1,num_wann
        write(ctmp,*) n_wann
        do m_wann=1,num_wann
          write(ctmp1,*) m_wann
          charTij(m_wann,n_wann) ="   WF"//trim(adjustl(ctmp1))//"-WF"//trim(adjustl(ctmp))
        enddo
      enddo
      
      HmnRdQ_noma_name = "./nomashift/noma_"//trim(adjustl(ctmpmode))//&
                        "/HmnRdQ_noma_"//trim(adjustl(ctmpmode))//".dat"
      HmnRdQ_noma_unit = io_file_unit()
      call open_file(HmnRdQ_noma_name,HmnRdQ_noma_unit)
      
      eph_noma_name    = "./nomashift/noma_"//trim(adjustl(ctmpmode))//&
                        "/eph_noma_"//trim(adjustl(ctmpmode))//".dat"
      eph_noma_unit    =  io_file_unit()
      call open_file(eph_noma_name,eph_noma_unit)
      
      TiidQ_noma_name = "./nomashift/noma_"//trim(adjustl(ctmpmode))//&
                        "/TiidQ_noma_"//trim(adjustl(ctmpmode))//".dat"
      TiidQ_noma_unit = io_file_unit()
      call open_file(TiidQ_noma_name,TiidQ_noma_unit)
      
      TijdQ_noma_name = "./nomashift/noma_"//trim(adjustl(ctmpmode))//&
                        "/TijdQ_noma_"//trim(adjustl(ctmpmode))//".dat"    
      TijdQ_noma_unit = io_file_unit()
      call open_file(TijdQ_noma_name,TijdQ_noma_unit)      

      write(ctmp,*) '(6A5,',ndQ,'(f15.8,1X),','A16)'
      write(HmnRdQ_noma_unit,ctmp) 'NHmn','Rx','Ry','Rz','WF0m','WFRn',(((idQ-nshiftstep-1)*dtadQ),idQ=1,ndQ),'dHmn/dQ'
      write(eph_noma_unit,"(6A5,A16)") 'NHmn','Rx','Ry','Rz','WF0m','WFRn','dHmn/dQ'
      
      write(ctmp,*) '(6I5,',(ndQ+1),'(f16.8))'  
      do irpts=1,nrpts
        do n_wann=1,num_wann
          do m_wann=1,num_wann
            NHmn = (irpts-1)*num_wann*num_wann+(n_wann-1)*num_wann+m_wann
            deta_HmnR(:)=HmnR_mode_dQ(m_wann,n_wann,irpts,:,imode)-HmnR_mode_dQ(m_wann,n_wann,irpts,nshiftstep+1,imode)
            dHmn_dQ=REAL(HmnR_mode_dQ(m_wann,n_wann,irpts,ndQ,imode)-HmnR_mode_dQ(m_wann,n_wann,irpts,1,imode))/(2.0*real(nshiftstep)*dtadQ)
            write(HmnRdQ_noma_unit,ctmp) NHmn,irvec(:,irpts,nshiftstep+1,imode),m_wann,n_wann,Real(deta_HmnR(:)),dHmn_dQ
            write(eph_noma_unit,'(6I5,f16.8)') NHmn,irvec(:,irpts,nshiftstep+1,imode),m_wann,n_wann,dHmn_dQ     
          enddo
        enddo
      enddo
      call close_file(HmnRdQ_noma_name,HmnRdQ_noma_unit)
      call close_file(eph_noma_name,eph_noma_unit)
      
      write(ctmp,*) "(A8,1X,",num_wann,'(A16,1X))'
      write(TiidQ_noma_unit,ctmp) "ldQ",((charwf(n_wann)),n_wann=1,num_wann)
      write(ctmp,*) "(A8,1X,",num_wann*num_wann,'(A16,1X))'
      write(TijdQ_noma_unit,ctmp) "ldQ",((charTij(n_wann,m_wann),n_wann=1,num_wann),m_wann=1,num_wann)
      
      do idQ=1,ndQ
        ldQ= (idQ-nshiftstep-1)*dtadQ
        write(TiidQ_noma_unit,'(F8.4,1X\)') ldQ
        write(TijdQ_noma_unit,'(F8.4,1X\)') ldQ
        do irpts=1,nrpts
          if(irvec(1,irpts,idQ,imode)==0 .and. irvec(2,irpts,idQ,imode)==0 .and. irvec(3,irpts,idQ,imode)==0) then
            write(ctmp,*) '(',num_wann,'(f16.8,1X))'
            write(TiidQ_noma_unit,ctmp) (REAL(HmnR_mode_dQ(i,i,irpts,idQ,imode)-HmnR_mode_dQ(i,i,irpts,nshiftstep+1,imode)),i=1,num_wann)
            write(ctmp,*) "(",num_wann*num_wann,'(A16,1X))'
            write(TijdQ_noma_unit,ctmp) ((REAL(HmnR_mode_dQ(m_wann,n_wann,irpts,idQ,imode)-&
            HmnR_mode_dQ(m_wann,n_wann,irpts,nshiftstep+1,imode)),m_wann=1,num_wann),n_wann=1,num_wann)
            
          end if
        end do
      end do
      call close_file(TijdQ_noma_name,TijdQ_noma_unit)
      call close_file(TiidQ_noma_name,TiidQ_noma_unit)
      
    enddo
    
    else       
      do imode=0,nmode
        write(ctmpmode,*) imode
        HmnR_Tij_name= "./Tij_parameter/Tij_"//trim(adjustl(ctmpmode))
        HmnR_Tij_unit=  io_file_unit()
        read(HmnR_Tij_unit,*)
        call open_file(HmnR_Tij_name,HmnR_Tij_unit)
        do ir2=-1,1
          do ir1=-1,1
            do n=1,num_wann
              do m=1,num_wann
                if(imode==0) then
                  read(HmnR_Tij_unit,"(T20,F12.6)") HmnR_Tij_0(m,n,ir1,ir2)
                else
                  read(HmnR_Tij_unit,"(T20,F12.6)") HmnR_Tij_ep(m,n,ir1,ir2,imode)                  
                endif
              enddo
            enddo
          enddo
        enddo
        call close_file(HmnR_Tij_name,HmnR_Tij_unit)
      end do
    
    endif
  
    HmnR_Tij_0  = HmnR_Tij_0/AU2EV   !!转换 为原子单位
    HmnR_Tij_ep = HmnR_Tij_ep/AU2EV*Au2ang*dsqrt(au2amu)
    
  end subroutine read_TB_parameters
  
  !!set the Tij parameter of elec and hole with out Coulomb interaction with elec and hole
  
  subroutine set_H_without_Coulomb(xx)
    implicit none
    
    real(kind=dp) :: xx(nfreem)
    integer :: initial1,initial2,m,n
    HmnR_Tij_elec=HmnR_Tij_0
    do ifreem=1,nfreem
      HmnR_Tij_elec=HmnR_Tij_elec+xx(ifreem)*HmnR_Tij_ep(:,:,:,:,ifreem)
    enddo
    HmnR_Tij_hole=HmnR_Tij_elec
    do iwann=1,num_wann
      HmnR_Tij_hole(iwann,iwann,0,0)=-HmnR_Tij_hole(iwann,iwann,0,0)
    enddo
    
    !<m0|H|nR>
    !<m0|
    do ia2site=0,na2site-1
      do ia1site=0,na1site-1
        !|nR>    R=(n-ia1site,m-ia2site) 使用了周期性边界条件          
        do m=ia2site-1,ia2site+1
          do n=ia1site-1,ia1site+1
            if ( m>=0 .and. m<na2site-1 .and. n>=0 .and. m<na1site-1 ) then
              initial1=(ia2site)*na1site*num_wann+(ia1site)*num_wann
              initial2=m*na1site*num_wann+n*num_wann
              ir1=n-ia1site
              ir2=m-ia2site
              H_elec(initial1+1:initial1+num_wann,initial2+1:initial2+num_wann)=HmnR_Tij_elec(:,:,ir1,ir2)
              H_hole(initial1+1:initial1+num_wann,initial2+1:initial2+num_wann)=HmnR_Tij_hole(:,:,ir1,ir2)
            endif
          enddo
        enddo 
      enddo
    enddo
    
  end subroutine set_H_without_Coulomb
  
  subroutine set_H_with_Coulomb(nn_elec,nn_hole)
    implicit none
    real(kind=dp)::nn_elec(num_wann,na1site,na2site),nn_hole(num_wann,na1site,na2site)
    integer:: i_WF,j_WF,m,n,initial1,initial2
    do ia2site=1,na2site
      do ia1site=1,na1site
        initial1=(ia2site-1)*na1site*num_wann+(ia1site-1)*num_wann                
        do i_WF=1,num_wann            
          do m=1,na2site
            do n=1,na1site
              initial2=(m-1)*na1site*num_wann+(n-1)*num_wann
              do j_WF=1,num_wann
                H_elec(initial1+i_WF,initial1+i_WF)=H_elec(initial1+i_WF,initial1+i_WF)-&
                  nn_hole(j_WF,n,m)*coulomb(i_WF,ia1site,ia2site,j_WF,n,m)
                H_hole(initial1+i_WF,initial1+i_WF)=H_hole(initial1+i_WF,initial1+i_WF)+&
                  nn_elec(j_WF,n,m)*coulomb(j_WF,n,m,i_WF,ia1site,ia2site)  
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
    
  end subroutine set_H_with_Coulomb
  
  function coulomb(elec_WF,nn,mm,hole_WF,iia1site,iia2site)
    use f95_precision
    use blas95
    implicit none
    real(kind=dp)::coulomb
    integer::elec_WF,nn,mm,hole_WF,iia1site,iia2site
    integer::Relec(3),Rhole(3)
    
    !real(kind=dp)::ee(1,nbasis)
    real(kind=dp)::Relec_hole(3)
    !real(kind=dp)::Reh
    !real(kind=dp)::E_ehExtion
    real(kind=dp)::a(3,3)
    real(kind=dp)::x(3)
    real(kind=dp)::y(3)
    Relec(1)=nn
    Relec(2)=mm
    Relec(3)=1
    Rhole(1)=iia1site
    Rhole(2)=iia2site
    Rhole(3)=1
    x = Relec-Rhole
    a = real_lattice
    call gemv(a,x,y,trans='T')    
    Relec_hole(:) = y(:)+Rwann(:,elec_WF)-Rwann(:,hole_WF)
    Reh = sqrt(Sum(Relec_hole**2))
    if(nn==iia1site .and. mm==iia2site) then
      !coulomb = (sqrt_elem_charge_SI)/(fopieps0*epsr*a_lattice)
      coulomb = 1.0/epsr*a_lattice
    else
      !coulomb = (sqrt_elem_charge_SI)/(fopieps0*epsr*Reh)
      coulomb = 1.0/epsr*Reh
    endif
  
  end function coulomb
  
  function en_addit(nn_elec,nn_hole)
    implicit none
    real(kind=dp) ::nn_elec(num_wann,na1site,na2site),nn_hole(num_wann,na1site,na2site)
    real(kind=dp) ::  en_addit
    integer::l,m
    do ia2site=1,na2site
      do ia1site=1,na1site
        do m_wann=1,num_wann
          do l=1,na2site
            do m=1,na1site
              do n_wann=1,num_wann
              en_addit=en_addit+coulomb(m_wann,ia1site,ia2site,n_wann,m,l)*&
              nn_elec(m_wann,ia1site,ia2site)*nn_hole(n_wann,m,l)
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
    
  end function
  
  
  subroutine calculate_eigen_energy_state(xx,nn_elec,nn_hole) 
    implicit none
    real(kind=dp) ::  xx(nfreem)
    real(kind=dp) ::  nn_elec(num_wann,na1site,na2site),nn_hole(num_wann,na1site,na2site)
    real(kind=dp) ::  e_coulomb
    call set_H_without_Coulomb(xx)
    call set_H_with_Coulomb(nn_elec,nn_hole)
    call dia_H(h_elec,e_elec,p_elec)
    call dia_H(h_hole,e_hole,p_hole)
    e_coulomb = en_addit(nn_elec,nn_hole)
    e_elec=e_elec+e_coulomb/2.0
    e_hole=e_hole-e_coulomb/2.0
    
    
  end subroutine
  !========================================!
  != calculate eigenenergy and eigenstate =!
  !========================================!
  subroutine dia_H(hh,ee,pp)
    use f95_precision
    !use blas95
    use lapack95
    implicit none
    !!use LAPACK with Fortran f95 interface
    !include "mkl_lapack.fi"
  
    integer:: ierror,info
    real(kind=dp) hh(1:nbasis,1:nbasis),ee(1:nbasis),pp(1:nbasis,1:nbasis)
    real(kind=dp),allocatable::fv1(:),fv2(:)
    !!!USED in MKL geev
    real(kind=dp),allocatable::vl(:,:),eei(:,:)
  
    allocate(fv1(1:nbasis))
    allocate(fv2(1:nbasis))
    allocate(vl(nbasis,nbasis))
    allocate(eei(nbasis,nbasis))
    !!!!!!!!!!!!!!!!!!!!!!!
    vl = 0.0
    eei = 0.0
    pp=hh
    !call rs(nbasis,nbasis,hh,ee,1,pp,fv1,fv2,ierror)
    !!     USE MKL lib could have a high speed in dgeev , sgeev   !in page 1131 and 1241
    !call geev(hh,ee,eei,vl,pp,info)
    call syev(pp,ee,'V','U')
    !call geev(hh,ee,eei)
    deallocate(vl,eei)
    deallocate(fv1,fv2)
    !!On exit, hh array is overwritten
  end subroutine dia_H
  
  !===================================!
  != calculate nonadiabatic coupling =!
  !===================================!
  != ref: notebook page 630          =!
  !===================================!

  subroutine calculate_nonadiabatic_coupling(ee,pp,dd)
    implicit none
    integer:: jbasis,ik1site,ik2site,initial1,initial2
    real(kind=dp)::ee(1:nbasis),pp(1:nbasis,1:nbasis),dd(1:nbasis,1:nbasis,1:nfreem)
    dd=0.0d0
    do ibasis=1,nbasis
      do jbasis=1,nbasis        
          do ifreem=1,nfreem
            !HHmnR_Tij_ep(nnum_wann,nnum_wann,-1:1,-1:1,nnfreem)
            do ia2site=0,na2site-1
              do ia1site=0,na1site-1
                initial1=(ia2site)*na1site*num_wann+(ia1site)*num_wann
                do n_wann=1,num_wann
                  ik1site=initial1+n_wann
                  do m=ia2site-1,ia2site+1
                    do n=ia1site-1,ia1site+1
                      if(m>=0 .and. m<=na2site-1 .and. n>=0 .and. n<=na1site-1 ) then
                      initial2=m*na1site*num_wann+n*num_wann
                      do m_wann=1,num_wann
                      ik2site=initial2+m_wann
                      !read <m0|H|nR>
                      if(adj_Tij(n_wann,m_wann,n-ia1site,m-ia2site)) then
                      dd(ibasis,jbasis,ifreem) = dd(ibasis,jbasis,ifreem)+&
                      pp(ik2site,ibasis)*pp(ik1site,jbasis)*HmnR_Tij_ep(n_wann,m_wann,n-ia1site,m-ia2site,ifreem)
                      endif
                      enddo
                      endif
                    enddo
                  enddo
                enddo
              enddo
            enddo 
            if(jbasis /= ibasis) then
              dd(ibasis,jbasis,ifreem)=dd(ibasis,jbasis,ifreem)/(ee(jbasis)-ee(ibasis))
              !dij_dQq
            endif
          end do
        !endif
      enddo
    enddo
    
  end subroutine calculate_nonadiabatic_coupling
  
  !=========================================================!
  != convert wavefunction from diabatix to adiabatic basis =!
  !=========================================================!

  subroutine convert_diabatic_adiabatic(pp,cc,ww)
  
  implicit none

    integer           :: jbasis
    real(kind=dp)     :: pp(1:nbasis,1:nbasis)
    complex(kind=dpc) :: cc(1:nbasis),ww(1:nbasis)

    ww=0.0d0
    do ibasis=1,nbasis
      do jbasis=1,nbasis
        ww(ibasis)=ww(ibasis)+pp(jbasis,ibasis)*cc(jbasis) !
      enddo
    enddo
    
  endsubroutine convert_diabatic_adiabatic
  
end module sh_hamiltonian
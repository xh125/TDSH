module sh_hamiltonian
  use sh_parameters
  use sh_constants
  use sh_io,only : io_error,io_file_unit,open_file,close_file
  
  implicit none
  
  integer :: nrpts,irpts,n_wann,m_wann
  integer :: n, m
  integer :: ir1, ir2, ir3 
  real(kind=dp) :: ReH, ImH
  character(len=maxlen) :: Hr_name
  integer               :: Hr_unit  
  integer,allocatable   :: ndegen(:)
  ! real and imag of HmnR  
  
  contains
  
  subroutine set_HH0()
    implicit none		
    
    !Hr_name  =DIR_APP(1:DIR_LEN)//"wannier/wannier90_hr.dat"
    Hr_name  ="./wannier/wannier90_hr.dat"
    Hr_unit  = io_file_unit()
    
    call open_file(Hr_name,Hr_unit)
    
    read(Hr_unit, *) ctmp
    read(Hr_unit, *) num_wann
    read(Hr_unit, *) nrpts
    allocate(ndegen(nrpts))		
  
    HH0=0.0d0
    
    read(Hr_unit, *) (ndegen(irpts), irpts=1, nrpts)
    !read <m0|H|nR>
    !in wannier90_hr.dat R,m,n,HR,HI
		do irpts=1, nrpts       !R
			do n_wann=1, num_wann   !n
				do m_wann=1, num_wann	!m
					read(Hr_unit,*) ir1, ir2, ir3, m, n, ReH, ImH
					do ia3site=1,na3site
						do ia2site=1,na2site
							do ia1site=1,na1site
								if (((ia1site+ir1) <= na1site .and. (ia1site+ir1) >=1) .and. ((ia2site+ir2) <= na2site .and. (ia2site+ir2) >=1) &
										.and. ((ia3site+ir3) <= na3site .and. (ia3site+ir3) >=1)) then
									HH0((ia3site-1)*na2site*na1site*num_wann+(ia2site-1)*na1site*num_wann+(ia1site-1)*num_wann+m,&
									(ia3site+ir3-1)*na2site*na1site*num_wann+(ia2site+ir2-1)*na1site*num_wann+(ia1site+ir1-1)*num_wann+n)=ReH ! /AU2eV
								end if
							end do
						end do
					end do
					!write(*,'(5i5,2f10.5)')i1, i2, i3, i4, i5, r1, r2
				enddo
			enddo
		enddo
    
    call close_file(Hr_name,Hr_unit)
  
  end subroutine set_HH0
  
  subroutine set_Hep()
    implicit none  	
    character(len=maxlen) :: epcfile_name,HmnRdQ_noma_name,eph_noma_name,&
                            TijdQ_noma_name,TiidQ_noma_name
    integer               :: epcfile_unit,HmnRdQ_noma_unit,eph_noma_unit,&
                             TijdQ_noma_unit,TiidQ_noma_unit
    integer               :: imode,nmode,idQ,ndQ,NHmn,i,j
    integer,allocatable   :: irvec(:,:,:,:)
    real(kind=dp)         :: ldQ
    character(len=maxlen) :: ctmpmode,ctmpldQ
    complex(kind=dpc),allocatable::HmnR_mode_dQ(:,:,:,:,:)
    character(len=16),allocatable::charwf(:),charTij(:,:)
    allocate(Hep(1:nbasis,1:nbasis,1:nfreem))
    allocate(charwf(num_wann),charTij(num_wann,num_wann))
    
    Hep=0.0d0	
    nmode = nfreem
    ndQ   = 2*nshiftstep+1
    dtadQ = real(nshiftdQ)/real(10000)
  
  if ( .Not. Lephfile)  then
    allocate(HmnR_mode_dQ(num_wann,num_wann,nrpts,ndQ,nmode))
    allocate(irvec(3, nrpts,ndQ,nmode))
    !对于每一个normal mode的处理
    do imode=1,nmode
      write(ctmpmode,*) imode
      !读取每一个mode下的各个shift下的wannier90_hr.dat文件的H的实部，
      !得到TB参数随该mode上偏移的变化,并写出该mode的phonon与电子之间的电声耦合强度，
      !包括局域与非局域
      do idQ=1,ndQ
        ldQ=(idQ-nshiftdQ-1)*dtadQ
        write(ctmpldQ,*) ldQ
        
        Hr_name= "./nomashift/noma_"//trim(adjustl(ctmpmode))//&
                  "shift_"//trim(adjustl(ctmpldQ))//"wannier90_hr.dat"
        Hr_unit=io_file_unit()
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
            do m_wann=1, num_wann	!m
              read(Hr_unit,*) ir1, ir2, ir3, m, n, ReH, ImH
              !write(*,'(5i5,2f10.5)')i1, i2, i3, i4, i5, r1, r2
              HmnR_mode_dQ(m_wann,n_wann,irpts,idQ,imode)= dcmplx(ReH, ImH)
            enddo
          enddo
          irvec(1, irpts,idQ,imode)=ir1
          irvec(2, irpts,idQ,imode)=ir2
          irvec(3, irpts,idQ,imode)=ir3
        enddo
        call close_file(Hr_name,Hr_unit)
        
      enddo
      
      do n_wann=1,num_wann
        write(ctmp,*) n_wann
        charwf(n_wann) = "WF_"//trim(adjustl(ctmp))
      enddo
      do n_wann=1,num_wann
        write(ctmp,*) n_wann
        do m_wann=1,num_wann
          write(ctmp1,*) m_wann
          charTij(m_wann,n_wann) ="WF"//trim(adjustl(ctmp1))//"-WF"//trim(adjustl(ctmp))
        enddo
      enddo
      
      HmnRdQ_noma_name = "./nomashift/noma_"//trim(adjustl(ctmpmode))//&
                        "HmnRdQ_noma_"//trim(adjustl(ctmpmode))//".dat"
      HmnRdQ_noma_unit = io_file_unit()
      eph_noma_name    = "./nomashift/noma_"//trim(adjustl(ctmpmode))//&
                        "eph_noma_"//trim(adjustl(ctmpmode))//".dat"
      eph_noma_unit    =  io_file_unit()
      TiidQ_noma_name = "./nomashift/noma_"//trim(adjustl(ctmpmode))//&
                        "TiidQ_noma_"//trim(adjustl(ctmpmode))//".dat"
      TijdQ_noma_name = "./nomashift/noma_"//trim(adjustl(ctmpmode))//&
                        "TijdQ_noma_"//trim(adjustl(ctmpmode))//".dat"
      TiidQ_noma_unit = io_file_unit()
      TijdQ_noma_unit = io_file_unit()
      
      call open_file(HmnRdQ_noma_name,HmnRdQ_noma_unit)
      call open_file(eph_noma_name,eph_noma_unit)
      write(ctmp,*) '(6A5,',ndQ,'(f15.8,1X),','A16)'
      write(HmnRdQ_noma_unit,ctmp) 'NHmn','Rx','Ry','Rz','WF0m','WFRn',(((idQ-nstep-1)*dtadQ),idQ=1,ndQ),'dHmn/dQ'
      write(ctmp,*) '(6A5,A16)'
      write(eph_noma_unit,ctmp) 'NHmn','Rx','Ry','Rz','WF0m','WFRn','dHmn/dQ'
      
      write(ctmp,*) '(6I5,',(ndQ+1),'f16.8)'  
      do irpts=1,nrpts
				do n_wann=1,num_wann
					do m_wann=1,num_wann
            NHmn = (irpts-1)*num_wann*num_wann+(n_wann-1)*num_wann+m_wann
						write(HmnRdQ_noma_unit,ctmp) NHmn,irvec(:,irpts,nshiftdQ+1,imode),m_wann,n_wann,&
							REAL(HmnR_mode_dQ(m_wann,n_wann,irpts,:,imode)-HmnR_mode_dQ(m_wann,n_wann,irpts,nstep+1,imode)),&
							REAL(HmnR_mode_dQ(m_wann,n_wann,irpts,ndQ,imode)-HmnR_mode_dQ(m_wann,n_wann,irpts,1,imode))/(2.0*real(nshiftdQ)*dtadQ)
						write(eph_noma_unit,'(6I5,f16.8)') NHmn,irvec(:,irpts,nshiftdQ+1,imode),m_wann,n_wann,&
							REAL(HmnR_mode_dQ(m_wann,n_wann,irpts,ndQ,imode)-HmnR_mode_dQ(m_wann,n_wann,irpts,1,imode))/(2.0*real(nshiftdQ)*dtadQ)            
          enddo
				enddo
			enddo
      call close_file(HmnRdQ_noma_name,HmnRdQ_noma_unit)
      call close_file(eph_noma_name,eph_noma_unit)
      
      call open_file(TiidQ_noma_name,TiidQ_noma_unit)
      call open_file(TijdQ_noma_name,TijdQ_noma_unit)
      write(ctmp,*) "(A8,1X,",num_wann,'(A16,1X))'
      write(TiidQ_noma_unit,ctmp) "ldQ",((charwf(n_wann)),n_wann=1,num_wann)
      write(ctmp,*) "(A8,1X,",num_wann*num_wann,'(A16,1X))'
      write(TijdQ_noma_unit,ctmp) "ldQ",((charTij(n_wann,m_wann),n_wann=1,num_wann),m_wann=1,num_wann)
      do idQ=1,ndQ
				ldQ= (idQ-nshiftdQ-1)*dtadQ
				write(TiidQ_noma_unit,'(F8.4,1X\)') ldQ
				write(TijdQ_noma_unit,'(F8.4,1X\)') ldQ
				do irpts=1,nrpts
					if(irvec(1,irpts,idQ,imode)==0 .and. irvec(2,irpts,idQ,imode)==0 .and. irvec(3,irpts,idQ,imode)==0) then
						write(ctmp,*) '(',num_wann,'(f16.8,1X))'
						write(TiidQ_noma_unit,ctmp) (REAL(HmnR_mode_dQ(i,i,irpts,idQ,imode)-HmnR_mode_dQ(i,i,irpts,nstep+1,imode)),i=1,num_wann)
						write(ctmp,*) "(",num_wann*num_wann,'(A16,1X))'
            write(TijdQ_noma_unit,ctmp) ((REAL(HmnR_mode_dQ(m_wann,n_wann,irpts,idQ,imode)-&
            HmnR_mode_dQ(m_wann,n_wann,irpts,nstep+1,imode)),m_wann=1,num_wann),n_wann=1,num_wann)
            
					end if
				end do
			end do
      call close_file(TijdQ_noma_name,TijdQ_noma_unit)
      call close_file(TiidQ_noma_name,TiidQ_noma_unit)
      
    enddo
  endif
  
    do imode=1,nmode
      write(ctmpmode,*) imode
      eph_noma_name= "./nomashift/noma_"//trim(adjustl(ctmpmode))//&
                        "eph_noma_"//trim(adjustl(ctmpmode))//".dat"
      eph_noma_unit=  io_file_unit()
      call open_file(eph_noma_name,eph_noma_unit)
      read(eph_noma_unit,*) ctmp
        do irpts=1, nrpts   !R
          do n_wann=1, num_wann   !n
            do m_wann=1, num_wann	!m
              read(epcfile_unit,'(6I5,f16.8)') NHmn,ir1, ir2, ir3, m, n, ReH
              do ia3site=1,na3site
                do ia2site=1,na2site
                  do ia1site=1,na1site
                    if ((ia1site+ir1 <= na1site .and. ia1site+ir1 >=1) .and. (ia2site+ir2 <= na2site .and. ia2site+ir2 >=1) &
                      .and. (ia3site+ir3 <= na3site .and. ia3site+ir3 >=1)) then
                      hep((ia3site-1)*na2site*na1site*num_wann+(ia2site-1)*na1site*num_wann+(ia1site-1)*num_wann+m,&
                      (ia3site+ir3-1)*na2site*na1site*num_wann+(ia2site+ir2-1)*na1site*num_wann+(ia1site+ir1-1)*num_wann+n,imode)=ReH ! /AU2eV
                    end if
                  end do
                end do
              end do
              !write(*,'(5i5,2f10.5)')i1, i2, i3, i4, i5, r1, r2
            enddo
          enddo
        enddo
      call close_file(eph_noma_name,eph_noma_unit)
    end do			
  end subroutine set_Hep
  
  subroutine set_HHt(xx)
  implicit none
    real(kind=dp) xx(1:nfreem)
    
    HHt=HH0
    do ifreem=1,nfreem
      !hh(ibasis,ibasis)=hh(ibasis,ibasis)+alpha*xx(ibasis)
      HHt(:,:)=HHt(:,:)+xx(ifreem)*Hep(:,:,ifreem)
    enddo
  
  end subroutine set_HHt

  !========================================!
  != calculate eigenenergy and eigenstate =!
  !========================================!
  subroutine calculate_eigen_energy_state(xx,ee,pp)
    use f95_precision
    !use blas95
    use lapack95
    implicit none
    !!use LAPACK with Fortran f95 interface
    !include "mkl_lapack.fi"
  
    integer:: ierror,info
    real(kind=dp) xx(1:nfreem),ee(1:nbasis),pp(1:nbasis,1:nbasis)
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
    call set_HHt(xx)
    !call rs(nbasis,nbasis,hh,ee,1,pp,fv1,fv2,ierror)
    !!     USE MKL lib could have a high speed in dgeev , sgeev   !in page 1131 and 1241
    !call geev(hh,ee,eei,vl,pp,info)
    call syev(HHt,ee,'V','U')
    pp=HHt
    !call geev(hh,ee,eei)
    deallocate(vl,eei)
    deallocate(fv1,fv2)
    !!On exit, hh array is overwritten
endsubroutine
  
  !===================================!
  != calculate nonadiabatic coupling =!
  !===================================!
  != ref: notebook page 630          =!
  !===================================!

  subroutine calculate_nonadiabatic_coupling(ee,pp,dd)
    implicit none

    integer:: jbasis,ik1site,ik2site
    real(kind=dp)::ee(1:nbasis),pp(1:nbasis,1:nbasis),dd(1:nbasis,1:nbasis,1:nfreem)

    dd=0.0d0
    do ibasis=1,nbasis
      do jbasis=1,nbasis
        if(jbasis /= ibasis) then
          do ifreem=1,nfreem
            do ik1site=1,nbasis
              do ik2site=1,nbasis
                dd(ibasis,jbasis,ifreem) = dd(ibasis,jbasis,ifreem)+&
                pp(ik2site,ibasis)*pp(ik1site,jbasis)*hep(ik1site,ik2site,ifreem)
              end do
            end do
            dd(ibasis,jbasis,ifreem)=dd(ibasis,jbasis,ifreem)/(ee(jbasis)-ee(ibasis))
          end do
        endif
      enddo
    enddo
    
  end subroutine calculate_nonadiabatic_coupling
  
  !=========================================================!
  != convert wavefunction from diabatix to adiabatic basis =!
  !=========================================================!

  subroutine convert_diabatic_adiabatic(pp,cc,ww)
  
  implicit none

    integer ::jbasis
    real(kind=dp) pp(1:nbasis,1:nbasis)
    complex(kind=dpc) cc(1:nbasis),ww(1:nbasis)

    ww=0.0d0
    do ibasis=1,nbasis
      do jbasis=1,nbasis
        ww(ibasis)=ww(ibasis)+pp(jbasis,ibasis)*cc(jbasis) !
      enddo
    enddo
    
  endsubroutine convert_diabatic_adiabatic
  
end module sh_hamiltonian
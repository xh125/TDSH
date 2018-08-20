module sh_parameters
  !! This module contains parameters to control the actions of SCSH.
  !! Also routines to read the parameters and write them out again.
  use mkl_service
  use sh_constants
  use sh_io       ,only : stdout,maxlen  
  implicit none
  
  integer,public  ::  mkl_threads,max_threads
  integer,public  ::  na1site,ia1site,na2site,ia2site
  integer,public  ::  ia1site_r,ia1site_l,ia2site_r,ia2site_l
  !! Number of cell in xy plan
  integer,public  ::  num_wann,iwann,n_wann,m_wann
  !! Number of wannier function in home-unit cell
  integer,public  ::  Num_occupied,Num_unoccupied,nband_hole,nband_elec
  !! Number of wannier Wavefunction in one cell
  integer,public  ::  nbasis,ibasis
  !! Number of basis in Hamitonian( nbasis = na1site*na2site*num_wann)
  integer,public  ::  nfreem,ifreem
  !! Number of phonon zhishu (nfreem = natoms*3-3)
  integer         ::  naver,iaver
  !! Number of trajects in Sample
  integer         ::  nsnap,isnap
  !! Number of snap in one trajects
  integer,public  ::  nstep,istep
  !! Number of electro MD step in one snap  
  real(kind=dp)   ::  temp
  !!The temperature of the system
  real(kind=dp)   ::  gamma
  !! The fraction coefficient characterizing system-bath coupling strength danwei ps(1-20ps)
  real(kind=dp)   ::  dt
  !! The time interval(fs)
  real(kind=dp),allocatable :: womiga(:)  
  !the w of nomal viboration
  real(kind=dp),allocatable,public :: HmnR_Tij_0(:,:,:,:),HmnR_Tij_ep(:,:,:,:,:)
  logical,allocatable,public       :: adj_Tij(:,:,:,:)
  !平衡位置TB参数，以及TB参数对简正坐标的导数( 电声耦合常数 )
  real(kind=dp),allocatable,public :: HmnR_Tij_total(:,:,:,:)
  real(kind=dp),allocatable,public :: HmnR_Tij_hole(:,:,:,:),HmnR_Tij_elec(:,:,:,:)
  ! 考虑最邻近cell之间的转移积分的TB的参数以及电声相互作用对转移积分的影响，
  ! 包含局域与非局域电声相互作用后的TB参数
  real(kind=dp),allocatable,public :: H_hole(:,:),H_elec(:,:)
  real(kind=dp),allocatable,public :: RE_H_hole(:,:,:,:,:,:),RE_H_elec(:,:,:,:,:,:)
  ! 电子和空穴的TB-Hamitonian
  integer,public  ::  isurface_elec,isurface_hole    
  !! 电子和空穴在势能面上的指标
  !! The elecrto or hole state index
  integer         ::  R_elec(0:2),R_hole(0:2)
  ! 电子和空穴所在的site包含晶胞R，和wannier函数指标
  integer ::  H_shape(2),H_RE_shape(6),P_shape(2),P_RE_shape(4)

  !real(kind=dp),allocatable,public :: q(:),v(:),e(:),p(:,:),d(:,:,:),g(:)
  !real(kind=dp),allocatable,public :: q0(:),v0(:),e0(:),p0(:,:),d0(:,:,:),g1(:)
  real(kind=dp),allocatable,public :: Q(:),Vq(:),Q0(:),Vq0(:)
  real(kind=dp),allocatable,public :: e_elec(:),p_elec(:,:),d_elec(:,:,:),g_elec(:)
  real(kind=dp),allocatable,public :: RE_p_elec(:,:,:,:)
  real(kind=dp),allocatable,public :: e_hole(:),p_hole(:,:),d_hole(:,:,:),g_hole(:)
  real(kind=dp),allocatable,public :: RE_p_hole(:,:,:,:)
  !电子和空穴的本征值本征态，非绝热耦合项，跃迁几率元
  real(kind=dp),allocatable,public :: e0_elec(:),p0_elec(:,:),d0_elec(:,:,:),g1_elec(:)
  real(kind=dp),allocatable,public :: e0_hole(:),p0_hole(:,:),d0_hole(:,:,:),g1_hole(:)
  real(kind=dp),allocatable,public :: RE_p0_elec(:,:,:,:),RE_p0_hole(:,:,:,:)
  
  real(kind=dp),allocatable,public :: pes_elec(:,:,:),inf_elec(:,:,:),csit_elec(:,:),&
                                      wsit_elec(:,:),psit_elec(:,:),ipr_elec(:)                         
  real(kind=dp),allocatable,public :: pes_hole(:,:,:),inf_hole(:,:,:),csit_hole(:,:),&
                                      wsit_hole(:,:),psit_hole(:,:),ipr_hole(:)
  real(kind=dp),allocatable,public :: pes_exciton(:,:,:)
  real(kind=dp),allocatable,public :: RE_csit_elec(:,:,:,:),RE_csit_hole(:,:,:,:)
  
  complex(kind=dpc),allocatable,public :: c_elec(:) , w_elec(:) , w0_elec(:)
  complex(kind=dpc),allocatable,public :: c_hole(:) , w_hole(:) , w0_hole(:)
  real(kind=dp),allocatable,public     :: n_elec(:) , n_hole(:),n0_elec(:),n0_hole(:)
  real(kind=dp),allocatable,public     :: RE_n_elec(:,:,:) , RE_n_hole(:,:,:)
  real(kind=dp),allocatable,public     :: RE_n0_elec(:,:,:), RE_n0_hole(:,:,:)
  !电子和空穴的占据分布
  real(kind=dp),allocatable,public :: xsit(:,:),ksit(:,:),msd(:),msds(:,:) 
  
  
  !initial hot phonon
  logical       :: L_hotphonon
  integer       :: hot_mode
  real(kind=dp) :: hot_scal
  ! write hamitonian parameter and test run one step time
  logical       :: L_test
  
  !initial state elec and hole
  character(len=2):: EHinitial   !BK or EN or WF or ES
  !if EHinitial == "BK" then
  integer ::  init_elec_K,elecK
  integer ::  init_elec_band,elecB
  integer ::  init_hole_K,holeK
  integer ::  init_hole_band,holeB
  !if EHinitial == "WF" then
  integer ::  init_elec_WF
  integer ::  init_hole_WF
  integer ::  index_elec
  integer ::  index_hole
  !if EHinitial == "EN" then
  real(kind=dp) :: init_elec_en
  real(kind=dp) :: init_hole_en
  !if EHinitial == "ES" then (initial the elec or hole energy surface index)
  integer :: init_elec_es
  integer :: init_hole_es
  ! does count of the effect of Coulub intaction of elec and hole
  logical ::  L_exciton
  !! Dielectric Constant ->  $$\epsilon_r$$
  real(kind=dp):: epsr
  !!!band project on WFs
  real(kind=dp),allocatable,public :: bands_projs(:,:,:)
  real(kind=dp)     ::t0,t1,t2,t3,time0,time1,time2,time3
  character(len=9)  ::cdate,ctime 
  
  integer         ::  Rcenter(3)
  ! 体系的中心
 
  ! Atom sites POSCAR information
  !character(len=maxlen),save :: seedname           
  real(kind=dp)        ,allocatable,public, save :: atoms_pos_frac(:,:,:)
  real(kind=dp)        ,allocatable,public, save :: atoms_pos_cart(:,:,:)
  integer              ,allocatable,public, save :: atoms_species_num(:)  
  character(len=maxlen),allocatable,public, save :: atoms_label(:)
  character(len=2)     ,allocatable,public, save :: atoms_symbol(:)
  integer                          ,public, save :: num_atoms,iatom
  !! Number of atoms in one cell
  integer                          ,public, save :: num_species
  !! Number of species of atoms in one cell    
  real(kind=dp),     public, save :: real_lattice(3,3) 
  real(kind=dp),     public, save :: a_lattice
  !晶格常数，用于HUbbard-like electron-hole interaction
  integer,           public, save :: num_kpts
  real(kind=dp),     public, save :: recip_lattice(3,3)
  real(kind=dp),     public, save :: cell_volume
  
  real(kind=dp),allocatable :: Rwann(:,:) 
  !the center of wannier wavefunction
  
  real(kind=dp)     ::sumg0_elec,sumg0_hole,sumg1_elec,sumg1_hole,&
                      minde_elec,minde_hole,flagd_elec,flagd_hole
  ! Private data
  integer                            :: num_lines,line_counter
  character(len=maxlen), allocatable :: in_data(:)
  character(len=maxlen)              :: ctmp,ctmp1,ctmp2
  logical                            :: ltmp,lexist,Lhole,Lephfile
  logical                            :: lbolziman
  integer                            :: itmp
  real(kind=dp)                      :: rtmp
  integer                            :: ierr
  character(len=maxlen)              :: msg
  
  ! Normal shift data
  integer::nshiftstep,ndQ
  real(kind=dp)::dtadQ,ldQ
  !!!!
  
	namelist / shinput / na1site,na2site,num_wann,temp,gamma,dt,&
                       nstep,nsnap,naver,elecb,eleck,holeb,holek,&
                       epsr,nshiftstep,dtadq,lephfile,Num_occupied,&
                       L_hotphonon,hot_mode,hot_scal,lbolziman,L_test,&
                       L_exciton,mkl_threads,&
                       EHinitial,init_elec_K,init_elec_band,init_hole_K,init_hole_band,&
                       init_elec_WF,init_hole_WF,index_elec,index_hole,&
                       init_elec_en,init_hole_en,&
                       init_elec_es,init_hole_es
  
  contains  
  
  subroutine soft_information()
  use sh_io
  implicit none
    time0   = io_time()
    call cpu_time(t0)
    stdout  = io_file_unit()
    open(unit=stdout,file="SCSH.out")
    call io_date(cdate,ctime)
    call mkl_get_version_string( ctmp )
    max_threads = mkl_get_max_threads()
    write(stdout,*) "SCSH complied with using IMKL:"
    write(stdout,*) trim(adjustl(ctmp))
    write(stdout,*) 'SCSH :Execution started on ',cdate,' at ',ctime
    write(ctmp,*) max_threads
    write(stdout,*) "By default, Intel MKL uses "//trim(adjustl(ctmp))//" threads"
    write(stdout,*) "where "//trim(adjustl(ctmp))//" is the number of physical cores on the system"
    
  end subroutine
    
  !==================================================================!
  subroutine read_parameters( )
  !==================================================================!
  !                                                                  !
  !! Read parameters and calculate derived values                    
  !                                                                  !
  !===================================================================  
    use sh_io,        only : io_error,io_file_unit,open_file,close_file
    implicit none
    call param_in_file()
    !用于将输入文件中的每一条写入字符串文件，并且改为小写，去除注释
    call read_param_infile()
    !将输入参数写入namelist文件并读取相关参数
    call set_num_wann()
    ! read num_wann from file of wannier90_hr.dat and set nbasis
    call readPOSCAR()
    
    call readWomiga()
    call readWFcentre()
    call readBandProjs()
    call treat_parameters()
    write(stdout,*) "Read parameter Successful!"
  end subroutine read_parameters
  
  !=======================================!
  subroutine param_in_file()  
  !用于将输入文件中的每一条写入字符串文件，并且改为小写，去除注释          
  !=======================================!
  !! Load the shin file into a character  
  !! array in_file, ignoring comments and  
  !! blank lines and converting everything 
  !! to lowercase characters               
  !=======================================!

    use sh_io,        only : io_file_unit,io_error
    use sh_utility,   only : utility_lowercase
    implicit none

    integer               :: in_unit,tot_num_lines,ierr,loop,in1,in2
    integer               :: ia,iz
    character(len=maxlen) :: dummy
    integer               :: pos
    character, parameter :: TABCHAR = char(9) !char(9)为制表符TAB
    character(len=80) ::msg
    ia = ichar('a')
    iz = ichar('z')
    
    in_unit=io_file_unit( )
    open (unit=in_unit, file='SHIN',form='formatted',status='old',iostat=ierr)
    if(ierr /= 0) then
      call io_error('Error: Problem opening input file SHIN')
    endif
    
    num_lines=0;tot_num_lines=0;ierr=0
    do while( ierr == 0 )
      read(in_unit, '(a)', iostat = ierr ,iomsg=msg) dummy   !参考p.177 P.540
      if(ierr > 0 ) then
        call io_error('Error: Problem reading input file SHIN')
        call io_error(msg)
      elseif(ierr == 0 )then
        
        ! convert all tabulation characters to spaces
        pos = index(dummy,TABCHAR) !查询字符串在字符串中出现的位置,并将制表符改为空格
        do while (pos /= 0)
          dummy(pos:pos) = ' '
          pos = index(dummy,TABCHAR)
        end do
        ! 
        dummy=adjustl(dummy)
        
        tot_num_lines=tot_num_lines+1
        if( dummy(1:1)/='!'  .and. dummy(1:1)/='#' ) then
          if(len_trim(adjustl(dummy)) > 0 ) num_lines=num_lines+1
        endif
      endif
    end do
    !得到SHIN文件中总的行数tot_num_lines以及非注释和空行 num_lines

    rewind(in_unit)

    allocate(in_data(num_lines),stat=ierr)  !字符串数组，内部文件 line=449
    if (ierr/=0) call io_error('Error allocating in_data in param_in_file')

    line_counter=0
    do loop=1,tot_num_lines
      read(in_unit, '(a)', iostat = ierr ,iomsg=msg) dummy
        if(ierr /= 0) then
          call io_error('Error: Problem opening input file SHIN')
          call io_error(msg)
        endif
      !I convert all tabulation characters to spaces
      pos = index(dummy,TABCHAR)
      do while (pos /= 0)
        dummy(pos:pos) = ' '
        pos = index(dummy,TABCHAR)
      end do
      !
      dummy=utility_lowercase(dummy) !将字符串中大写字母全部改为小写
      dummy=trim(adjustl(dummy))     !
      if( dummy(1:1)=='!' .or.  dummy(1:1)=='#' ) cycle
      if(len(trim(dummy)) == 0 ) cycle
      if(index(dummy,'=') <=1 )  cycle  !当该行中没有‘=’ 或‘=’前没有内容则跳过该行
      line_counter=line_counter+1
      
      !去除有效行信息中的注释部分，注释可以采用 ！或者 #
      in1=index(dummy,'!')
      in2=index(dummy,'#')
      if(in1==0 .and. in2==0)  in_data(line_counter)=dummy
      !不存在'!'与'#'
      if(in1==0 .and. in2>0 )  in_data(line_counter)=dummy(:in2-1)
      if(in2==0 .and. in1>0 )  in_data(line_counter)=dummy(:in1-1)
      if(in2> 0 .and. in1>0 )  in_data(line_counter)=dummy(:min(in1,in2)-1)
      
      !如果输入参数为字符串，则给字符串加上"*"
      !itmp = index(dummy,'=')
      !ctmp = dummy(:itmp-1)
      !ctmp1= dummy(itmp+1:)
      !ctmp1= trim(adjustl(ctmp1))
      !itmp = ichar(ctmp1(1:1))
      !if(itmp>=ia .and. itmp<=iz) then
      !  dummy = ctmp//"='"//ctmp1//"'"
      !endif
      
    end do
    !得到包含有效信息的行数line_counter,和相应的数据in_data(line_counter)

    close(in_unit)

  end subroutine param_in_file  
  
  subroutine read_param_infile()
  
    use sh_io
    implicit none
    integer::incar_unit,i
    incar_unit = io_file_unit()
    open(unit=incar_unit,status='SCRATCH',iostat=ierr,iomsg=msg)
    if(ierr > 0 ) then
      call io_error('Error: Problem reading SCRATCH input namelist file')
      call io_error(msg)
    elseif(ierr == 0 )then    
      write(incar_unit,*)"&shinput" 
      do i=1,line_counter
        write(incar_unit,*) trim(adjustl(in_data(i)))
      enddo
      write(incar_unit,"(A1)") "/"
    endif
    rewind(incar_unit)        
    write(stdout,*)   "======================================================"
    write(stdout,"(1X,10X,A)") "The namelist file as follows"
    write(stdout,*)   "======================================================"
    do i=1,line_counter+2
      read(incar_unit,"(A80)") ctmp
      write(stdout,"(A80)") ctmp
    enddo
    write(stdout,*)   "======================================================="
    rewind(incar_unit)
      !initial parameters
      mkl_threads = max_threads
      na1site = 100
      na2site = 100
      num_wann= 22
      temp    = 300
      gamma   = 0.05
      dt      = 0.01
      nstep   = 10000
      nsnap   = 100
      naver   = 500
      elecB   = 18
      elecK   = 100
      holeB   = 17
      holeK   = 100
      Lephfile= .False.
      lbolziman=.TRUE.
      L_test  = .TRUE.
      Num_occupied= 1
      L_hotphonon= .true.
      hot_mode = 1
      hot_scal = 1.5
      !end initial
      read(UNIT=incar_unit,nml=shinput,iostat=ierr,iomsg=msg)
      if(ierr /= 0) then
        call io_error('Error: Problem reading namelist file SHIN')
        call io_error(msg)
      endif
     
      init_elec_K     = elecK
      init_elec_band  = elecB
      init_hole_K     = holeK
      init_hole_band  = holeB
      Num_unoccupied  = num_wann - Num_occupied
      nband_elec      = Num_unoccupied
      nband_hole      = Num_occupied
      
      close(incar_unit)
      
  end subroutine read_param_infile
  
  subroutine set_num_wann()
    use sh_io,        only : io_file_unit,open_file,close_file
    implicit none
  	
    character(len=255) :: Hr_name
    integer            :: Hr_unit   		
    Hr_name  ="./wannier/wannier90_hr.dat"
    inquire(directory = './wannier',exist=lexist)
    if(lexist) then
      Hr_unit = io_file_unit()
      call open_file(Hr_name,Hr_unit)
      read(Hr_unit, *) ctmp
      read(Hr_unit, *) num_wann
      rewind(Hr_unit)
      call close_file(Hr_name,Hr_unit)
    endif
    nbasis = na1site*na2site*num_wann
    write(stdout,*) "NUM_wann=",num_wann
    
  end subroutine set_num_wann
  
  subroutine readPOSCAR()
    
  !use sh_constants, only : dp
  use sh_io,        only : io_file_unit,io_error,open_file,close_file
  use sh_utility,   only : utility_recip_lattice  
  implicit none      
   
    real(kind=dp)::scaling
    integer::i,j    
    integer::poscar_unit
    character(len=maxlen)::poscar_name
    
    poscar_name  = './wannier/POSCAR'
    poscar_unit  = io_file_unit()
    call open_file(poscar_name,poscar_unit)
    read(poscar_unit,*) ctmp
    read(poscar_unit,"(F16.8)") scaling
    read(poscar_unit,'(3F20.12)') ((real_lattice(I,J),J=1,3),I=1,3)
    real_lattice = real_lattice*scaling
    a_lattice = sqrt(Sum(real_lattice(1,:)**2))
    call utility_recip_lattice (real_lattice,recip_lattice,cell_volume)
    
    num_atoms    = 0
    num_species  = 1
    read(poscar_unit,*)
    do while(ierr <= 0)
      allocate(atoms_species_num(num_species))
      read(poscar_unit,FMT=*,iostat=ierr) (atoms_species_num(i),i=1,num_species)
      num_species = num_species + 1
      deallocate (atoms_species_num)
      backspace(poscar_unit)
    end do
    num_species = num_species - 2
    
    allocate(atoms_symbol(num_species))
    allocate(atoms_species_num(num_species))
    backspace(poscar_unit)
    backspace(poscar_unit)
    read(poscar_unit,FMT=*,iostat=ierr) (atoms_symbol(i),i=1,num_species)
    read(poscar_unit,FMT=*,iostat=ierr) (atoms_species_num(i),i=1,num_species)
    num_atoms = sum(atoms_species_num)
    nfreem = 3 * (num_atoms-1)
    call close_file(poscar_name,poscar_unit)
    write(stdout,*) "NUM_atoms=",num_atoms
    write(stdout,*) "Atoms_symbol:",(atoms_symbol(i),i=1,num_species)
    write(stdout,*) "Atoms of each symbol:",(atoms_species_num(i),i=1,num_species)
  end subroutine readPOSCAR
 
  !====================================================================!
  !  read phonon womiga(mev) !use the OUTCAR of phonon calculate       !
  !====================================================================!
  subroutine readwomiga()
    use sh_constants
    use sh_io,only:io_file_unit,io_error,open_file,close_file
    
    implicit none
    integer              ::wvecter_unit
    character(len=maxlen)::wvecter_name
    
    allocate(womiga(nfreem))
    
    wvecter_unit = io_file_unit()
    wvecter_name = "./phonon-gamma/wvecter.txt"
    inquire(file="./phonon-gamma/wvecter.txt",exist=lexist)
    if(.Not. lexist) then
      write(stdout,*) "Wrong: the file './phonon-gamma/wvecter.txt' is not found!"
      stop
    endif
    
    call open_file(wvecter_name,wvecter_unit)
    read(wvecter_unit,'(//,A)') ctmp
    do ifreem=1,nfreem
      read(wvecter_unit,'(A)') ctmp
      read(wvecter_unit,'(T63,F12.7)') womiga(ifreem) !mev
      do iatom=1,num_atoms+1
        read(wvecter_unit,*) ctmp
      enddo
    enddo
    
    call close_file(wvecter_name,wvecter_unit)
    write(stdout,*) "===========Womiga=============="
    do ifreem=1,nfreem
      write(stdout,*) "womiga(",ifreem,")=",womiga(ifreem),"mev"
    enddo
    write(stdout,*) "=======END-Womiga=============="
    
  end subroutine readwomiga
  
  !================================================!
  != read "wannier90.wout" setting Rwann(3,nwann) =!
  !================================================!
  subroutine readWFcentre()
    use sh_constants,only:dp
    use sh_io 
    implicit none
  
    integer::i,wannier_centres_unit
    character(len=maxlen)::wannier_centres_name
    wannier_centres_name = './wannier/wannier90_centres.xyz'    
    wannier_centres_unit = io_file_unit()
    call open_file(wannier_centres_name,wannier_centres_unit)
    allocate(Rwann(3,num_wann))    
    read(wannier_centres_unit,"(/,A)") ctmp
    do iwann=1,num_wann
      read(wannier_centres_unit,*) ctmp,(Rwann(i,iwann),i=1,3)
    enddo
    call close_file(wannier_centres_name,wannier_centres_unit)
    
  end subroutine readWFcentre
  
  !================================================!
  !=get the project of line band on wannier orbital!
  !================================================!
  subroutine readBandProjs()
    use sh_constants
    use sh_io
    implicit none
    integer::total_pts,ipt   !num of kpointer in band line
    real(kind=dp),allocatable ::  xval(:)   !kpointer lable in x value
    real(kind=dp),allocatable ::  eig_int(:,:) !line band
    integer:: i , j
    integer:: bandproj_unit
    character(len=maxlen):: bandproj_name
    
    bandproj_unit = io_file_unit()
    bandproj_name = "./wannier/wannier90_bandproj.dat"
    call open_file(bandproj_name,bandproj_unit)
    read(bandproj_unit,*) ctmp,itmp,ctmp,total_pts
    
    allocate(xval(total_pts))
    allocate(eig_int(total_pts,num_wann))
    allocate(bands_projs(total_pts,num_wann,num_wann))
    ! benzheng tai zai wannier orbital shang de tou ying.
    
    write(ctmp,*) num_wann + 2
    ctmp = '('//trim(adjustl(ctmp))//'E16.8)'
    
    do i=1,num_wann
      do ipt=1,total_pts
        read(bandproj_unit,ctmp) xval(ipt),eig_int(ipt,i),(bands_projs(ipt,i,j),j=1,num_wann)
      enddo
    enddo
    
    call close_file(bandproj_name,bandproj_unit)
    deallocate(xval,eig_int)
  
  end subroutine readBandProjs
  
  subroutine treat_parameters()
    use sh_constants
    use sh_io
    implicit none	
    
    H_shape       = nbasis
    H_RE_shape(1) = num_wann
    H_RE_shape(2) = na1site
    H_RE_shape(3) = na2site
    H_RE_shape(4) = num_wann
    H_RE_shape(5) = na1site
    H_RE_shape(6) = na2site    
    P_shape       = nbasis
    P_RE_shape(1) = num_wann
    P_RE_shape(2) = na1site
    P_RE_shape(3) = na2site
    P_RE_shape(4) = nbasis
    Rcenter(1) =anint(na1site/2.0d0)
    Rcenter(2) =anint(na2site/2.0d0)
    !HmnR_Tij(m,n,ir1, ir2)
    allocate(HmnR_Tij_0(num_wann,num_wann,-1:1,-1:1))
    allocate(adj_Tij(num_wann,num_wann,-1:1,-1:1))
    !HmnR_Tij_ep(m,n,ir1, ir2,nfreem)
    allocate(HmnR_Tij_ep(num_wann,num_wann,-1:1,-1:1,nfreem))
    allocate(HmnR_Tij_total(num_wann,num_wann,-1:1,-1:1))
    allocate(HmnR_Tij_hole(num_wann,num_wann,-1:1,-1:1))
    allocate(HmnR_Tij_elec(num_wann,num_wann,-1:1,-1:1))
    !allocate(HH0(1:nbasis,1:nbasis),HHt(1:nbasis,1:nbasis))
    !allocate(Hep(1:nbasis,1:nbasis,1:nfreem))
    allocate(H_hole(nbasis,nbasis),H_elec(nbasis,nbasis))
    allocate(RE_H_hole(num_wann,na1site,na2site,num_wann,na1site,na2site))
    allocate(RE_H_elec(num_wann,na1site,na2site,num_wann,na1site,na2site))
    allocate(Q(1:nfreem),Vq(1:nfreem),Q0(1:nfreem),Vq0(1:nfreem))
    allocate(e_elec(1:nbasis),p_elec(1:nbasis,1:nbasis),d_elec(1:nbasis,1:nbasis,1:nfreem),g_elec(1:nbasis))
    allocate(e_hole(1:nbasis),p_hole(1:nbasis,1:nbasis),d_hole(1:nbasis,1:nbasis,1:nfreem),g_hole(1:nbasis))
    allocate(RE_p_elec(num_wann,na1site,na2site,nbasis),RE_p_hole(num_wann,na1site,na2site,nbasis))
    allocate(e0_elec(1:nbasis),p0_elec(1:nbasis,1:nbasis),d0_elec(1:nbasis,1:nbasis,1:nfreem),g1_elec(1:nbasis))
    allocate(e0_hole(1:nbasis),p0_hole(1:nbasis,1:nbasis),d0_hole(1:nbasis,1:nbasis,1:nfreem),g1_hole(1:nbasis))
    allocate(RE_p0_elec(num_wann,na1site,na2site,nbasis),RE_p0_hole(num_wann,na1site,na2site,nbasis))
    allocate(pes_elec(0:nbasis,1:nsnap,1:naver),inf_elec(1:3,1:nsnap,1:naver),csit_elec(1:nbasis,1:nsnap),&
           wsit_elec(1:nbasis,1:nsnap),psit_elec(1:nbasis,1:nsnap),ipr_elec(1:nsnap) )           
    allocate(pes_hole(0:nbasis,1:nsnap,1:naver),inf_hole(1:3,1:nsnap,1:naver),csit_hole(1:nbasis,1:nsnap),&
           wsit_hole(1:nbasis,1:nsnap),psit_hole(1:nbasis,1:nsnap),ipr_hole(1:nsnap) )
    allocate(pes_exciton(-1:nbasis,1:nsnap,1:naver))
    allocate(c_elec(1:nbasis),w_elec(1:nbasis),w0_elec(1:nbasis))
    allocate(c_hole(1:nbasis),w_hole(1:nbasis),w0_hole(1:nbasis))
    allocate(n_elec(nbasis),n_hole(nbasis),n0_elec(nbasis),n0_hole(nbasis))
    allocate(RE_n_elec(num_wann,na1site,na2site),RE_n_hole(num_wann,na1site,na2site))
    allocate(RE_n0_elec(num_wann,na1site,na2site),RE_n0_hole(num_wann,na1site,na2site))
    allocate(msds(1:nsnap,1:naver),xsit(1:nbasis,1:nsnap),ksit(1:nbasis,1:nsnap),msd(1:nsnap))
    
    !! Change to Hartree atomic units
    temp          = temp/Au2k
    gamma         = gamma*Au2ps
    kb            = k_B_SI/Au2J*Au2k !in Hartree atomic units kb = 1
    dt            = dt/Au2fs
    real_lattice  = real_lattice/Au2ang
    a_lattice     = a_lattice/Au2ang
    cell_volume   = cell_volume/(Au2ang)**3
    recip_lattice = recip_lattice*Au2ang
    womiga        = womiga/Au2mev
    Rwann         = Rwann/Au2ang
   
    
    !!initiali
    pes_elec  =  0.0d0
    pes_hole  =  0.0d0
    inf_elec  =  0.0d0
    inf_hole  =  0.0d0
    csit_elec =  0.0d0
    csit_hole =  0.0d0
    wsit_elec =  0.0d0
    wsit_hole =  0.0d0
    psit_elec =  0.0d0
    psit_hole =  0.0d0
    xsit      =  0.0d0
    ksit      =  0.0d0
    msd       =  0.0d0
    ipr_elec  =  0.0d0
    ipr_hole  =  0.0d0
    msds      =  0.0d0
    HmnR_Tij_0    =   0.0d0
    HmnR_Tij_elec =   0.0d0
    HmnR_Tij_ep   =   0.0d0
    HmnR_Tij_hole =   0.0d0
    HmnR_Tij_total=   0.0d0
  
  end subroutine  treat_parameters
 
end module sh_parameters
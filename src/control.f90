module control
  !-----------------------------------------------------------------------
  !!
  !! Common variables for the epw program
  !!  
  !-----------------------------------------------------------------------
  implicit none
  use kinds, only :dp
  integer,parameter,public :: maxlen = 120
  
  !!system dimention of real and kpoint parameter 
  integer :: dimention
  integer :: na1
  integer :: na2
  integer :: na3
  integer :: nk1
  integer :: nk2
  integer :: nk3
  integer :: nq1
  integer :: nq2
  integer :: nq3
  !! representation of the system -- "adiabatic","blochstat","wfstat","atomob","molerob"
  character(len=maxlen):: representation
  !!when use wfstat representation ,use method of nomal shift to calculate elec-phonon coupling
  !!the parameter must == in normal shift input
  integer :: nshiftstep
  real(kind=dp) :: dtadq
  
  !!
  real(kind=dp) :: temp
  real(kind=dp) :: gamma
  real(kind=dp) :: dt
  integer :: nstep
  integer :: nsnap
  integer :: naver
  
  !! initial core state init_normal_stat--"class","quantum" 
  character(len=maxlen) ::  initnmstat
  logical       :: L_hotphonon
  integer       :: hot_mode
  real(kind=dp) :: hot_scal
  
  !! initial elec-hole state--"BK","WF","EN","ES"
  !! 分别对应 blochstat 下的kpoint-band,
  !! wfstat -- wannier-function
  !! wfstat -- energy
  !! wfstat -- energy-surface
  integer ::  Num_occupied
  character(len=maxlen) ::  initehstat
  integer :: initek            !init_elec_K
  integer :: initeb            !init_elec_band
  integer :: inithk            !init_hole_K
  integer :: inithb            !init_hole_band
  integer :: initeWF           !init_elec_WF
  integer :: inithWF           !init_hole_WF
  real(kind=dp) :: initeEN     !init_elec_en
  real(kind=dp) :: inithEN     !init_hole_en
  integer :: initeES           !init_elec_esurface
  integer :: inithES           !init_hole_esurface
  
  !! use the exciton effects
  logical       :: L_exciton
  real(kind=dp) :: epsr
  
  !! which SH method do we use
  character(len=maxlen) ::  MSH  !"FSSH","LD-FSSH","SC-FSSH","mSC-FSSH","CC_FSSH"
  !! do we use decoherence method to account for decoherence correction
  logical::Ldecoherece
  !! which methoe of decoherence do we use
  character(len=maxlen) :: Tdecoherence !!"enbased"
  
  !! do we account for feedback of elec to core
  logical:: Lfeedback
  
  
  !! MKL parallel
  integer :: mkl_threads
  
end module control
                       
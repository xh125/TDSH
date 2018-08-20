module io_global
  implicit none
  integer,public,parameter  ::  stdin   = 5 !! unit connected to standard input
  integer,public,parameter  ::  stdout  = 6 ! unit connected to input file (xml or text)
  integer,public,parameter  ::  qestdin = 9 ! unit connected to standard output
  !
  ! For parallel execution: I/O within an image
  ! These are set at startup by calling mp_world_start
  !
  integer :: ionode_id= 0       ! index of the i/o node for this image
  Logical :: Lionode  = .True.  ! true if this processor is a i/o node
                                   ! for this image 
end module io_global
  
  
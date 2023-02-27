!#######################################################################
!
! Data Input/Output System
!
! NAME
!   DIOS - a generic interface to input/output data and files
!
!
! BACKGROUND
!
!   This system refers to the MDF used in TM5.
!
!   Creation of file follows the steps similar to writing a NetCDF file:
!     * open a file
!     * (global attributes)
!     * definition of dimensions (plus attributes)
!     * definition of variables (plus attributes)
!     * end of definition phase
!     * write one or more time records
!     * close file
!
!
! PROCEDURES
!
!   !
!   ! module initialisation
!   !
!
!   subroutine dios_init( stat )
!     integer, intent(out)                ::  stat
!
!   !
!   ! write data
!   !
!
!   ! create a new file for output:
!   subroutine dios_create( file_name, file_type, open_mode, hid, stat )
!     character(len=*), intent(in)        ::  file_name
!     integer, intent(in)                 ::  file_type
!     integer, intent(in)                 ::  open_mode
!     integer, intent(out)                ::  hid
!     integer, intent(out)                ::  stat
!
!   ! ... or create more than one with different formats;
!   ! specify a single base name and an equal number of extensions and type:
!   subroutine DIOS_Create( basename, exts, ftypes, cmode, hid, stat )
!     character(len=*), intent(in)        ::  basename
!     character(len=*), intent(in)        ::  exts(:)
!     integer, intent(in)                 ::  ftypes(:)
!     integer, intent(in)                 ::  cmode
!     integer, intent(out)                ::  hid
!     integer, intent(out)                ::  stat
!
!   subroutine DIOS_Def_Dim( hid, name, length, dimid, stat )
!     integer, intent(in)                 ::  hid
!     character(len=*), intent(in)        ::  name
!     integer, intent(in)                 ::  length
!     integer, intent(out)                ::  dimid
!     integer, intent(out)                ::  stat
!
!   subroutine DIOS_Def_Var( hid, name, xtype, dimids, varid, stat, &
!                                 compression, deflate_level )
!     integer, intent(in)                 ::  hid
!     character(len=*), intent(in)        ::  name
!     integer, intent(in)                 ::  xtype
!     integer, intent(in)                 ::  dimids(:)
!     integer, intent(out)                ::  varid
!     integer, intent(out)                ::  stat
!     integer, intent(in), optional       ::  compression
!     integer, intent(in), optional       ::  deflate_level  ! 0-9
!
!   subroutine DIOS_Put_Att( hid, varid, name, values, stat )
!     integer, intent(in)                 ::  hid
!     integer, intent(in)                 ::  varid
!     character(len=*), intent(in)        ::  name
!     <type>, intent(in)                  ::  value(s)
!     integer, intent(out)                ::  stat
!
!   subroutine DIOS_EndDef( hid, stat )
!     integer, intent(in)                 ::  hid
!     integer, intent(out)                ::  stat
!
!   ! put variable:
!   subroutine DIOS_Put_Var( hid, varid, values, stat, &
!                             start, count, stride, map )
!     integer, intent(in)                 ::  hid
!     integer, intent(in)                 ::  varid
!     <type>, intent(in)                  ::  values(<shape>)
!     integer, intent(out)                ::  stat    
!     integer, intent(in), optional       ::  start(:)        ! (/1,1,..,1/)
!     integer, intent(in), optional       ::  count(:)
!     integer, intent(in), optional       ::  stride(:)
!     integer, intent(in), optional       ::  map(:)
!
!   ! close file(s):
!   subroutine DIOS_Close( hid, stat )
!     integer, intent(out)                ::  hid
!     integer, intent(out)                ::  stat
!
!   !
!   ! read file
!   !
!
!   ! open single file:
!   subroutine DIOS_Open( filename, ftype, mode, hid, stat )
!     character(len=*), intent(in)        ::  filename
!     integer, intent(in)                 ::  ftype
!     integer, intent(in)                 ::  mode
!     integer, intent(out)                ::  hid
!     integer, intent(out)                ::  stat
!
!   subroutine DIOS_Inquire( hid, stat, &
!                             nDimensions, nVariables, nAttributes )
!     integer, intent(in)                 ::  hid
!     integer, intent(out)                ::  stat
!     integer, intent(out), optional      ::  nDimensions
!     integer, intent(out), optional      ::  nVariables
!     integer, intent(out), optional      ::  nAttributes
!
!   subroutine DIOS_Inq_DimID( hid, name, dimid, stat )
!     integer, intent(in)                           ::  hid
!     character(len=*), intent(in)                  ::  name
!     integer, intent(out)                          ::  dimid
!     integer, intent(out)                          ::  stat
!
!   subroutine DIOS_Inquire_Dimension( hid, dimid, stat, name, length, unlimited )
!     integer, intent(in)                           ::  hid
!     integer, intent(in)                           ::  dimid
!     integer, intent(out)                          ::  stat
!     character(len=*), intent(out), optional       ::  name
!     integer, intent(out), optional                ::  length
!     logical, intent(out), optional                ::  unlimited
!
!   subroutine DIOS_Inq_VarID( hid, name, varid, stat )
!     integer, intent(in)                           ::  hid
!     character(len=*), intent(in)                  ::  name
!     integer, intent(out)                          ::  varid
!     integer, intent(out)                          ::  stat
!
!   subroutine DIOS_Inquire_Variable( hid, varid, stat, &
!                                      name, xtype, ndims, dimids, natts )
!     integer, intent(in)                           ::  hid
!     integer, intent(in)                           ::  varid
!     integer, intent(out)                          ::  stat
!     character(len=*), intent(out), optional       ::  name
!     integer, intent(out), optional                ::  xtype
!     integer, intent(out), optional                ::  ndims
!     integer, intent(out), optional                ::  dimids(:)
!     integer, intent(out), optional                ::  natts
!
!   subroutine DIOS_Get_Var( hid, varid, values, stat, &
!                                        start, count, stride, map )
!     integer, intent(in)                 ::  hid
!     integer, intent(in)                 ::  varid
!     <type>, intent(out)                 ::  values(<shape>)
!     integer, intent(out)                ::  stat 
!     integer, intent(in), optional       ::  start (:)
!     integer, intent(in), optional       ::  count (:)
!     integer, intent(in), optional       ::  stride(:)
!     integer, intent(in), optional       ::  map   (:)
!
!   subroutine DIOS_Inq_AttName( hid, varid, attnum, name, stat )
!     integer, intent(in)                           ::  hid
!     integer, intent(in)                           ::  varid
!     integer, intent(in)                           ::  attnum
!     character(len=*), intent(out)                 ::  name
!     integer, intent(out)                          ::  stat
!
!   subroutine DIOS_Inquire_Attribute( hid, varid, name, stat, xtype, length )
!     integer, intent(in)                           ::  hid
!     integer, intent(in)                           ::  varid
!     character(len=*), intent(out)                 ::  name
!     integer, intent(out)                          ::  stat
!     integer, intent(out), optional                ::  xtype
!     integer, intent(out), optional                ::  length
!
!   subroutine DIOS_Get_Att( hid, varid, name, values, stat )
!     integer, intent(in)                 ::  hid
!     integer, intent(in)                 ::  varid
!     character(len=*), intent(in)        ::  name
!     <type>, intent(out)                 ::  value(s)
!     integer, intent(out)                ::  stat
!
!   ! close file:
!   subroutine DIOS_Close( hid, stat )
!     integer, intent(out)                ::  hid
!     integer, intent(out)                ::  stat
!
!   !
!   ! show file content
!   !
!
!   ! show file headers similar to 'ncdump -h' ;
!   ! file type is guessed from extension if not specified directly:
!   subroutine DIOS_Show( filename, stat [,filetype=DIOS_NETCDF4|DIOS_HDF|...] )
!     character(len=*), intent(in)        ::  filename
!     integer, intent(out)                ::  stat
!     integer, intent(in), optional       ::  filetype
!     integer, intent(out)                ::  stat
!
!   !
!   ! end module access
!   !
!
!   ! done with module:
!   subroutine DIOS_Done( stat )
!     integer, intent(out)                ::  stat
!
!
!
! GLOBAL ATTRIBUTES
!
!   To write global attributes, use the constant 'DIOS_GLOBAL' as variable id.
!
!
! UNLIMITED DIMENSION
!
!   To define an unlimited dimension, use the constant 'DIOS_UNLIMITED' 
!   as dimension length.
!
!
!#######################################################################

module dios_mod

  implicit none
  
  ! --- in/out ---------------------------------------

  private

  public  ::  dios_init, dios_done
  
  public  ::  dios_create, dios_open, dios_close
  public  ::  dios_enddef
  public  ::  dios_inquire

  public  ::  dios_def_dim
  public  ::  dios_inq_dimid
  public  ::  dios_inquire_dimension

  public  ::  dios_def_var
  public  ::  dios_inq_varid
  public  ::  dios_inquire_variable
  public  ::  dios_put_var
  public  ::  dios_get_var

  public  ::  dios_put_att
  public  ::  dios_get_att
  
  public  ::  dios_show

  public  ::  DIOS_NONE

  public  ::  DIOS_NEW
  public  ::  DIOS_REPLACE
  public  ::  DIOS_READ
  public  ::  DIOS_WRITE
  
  public  ::  DIOS_NETCDF
  
  public  ::  DIOS_CHAR  
  public  ::  DIOS_BYTE  
  public  ::  DIOS_SHORT 
  public  ::  DIOS_INT   
  public  ::  DIOS_INT64
  public  ::  DIOS_FLOAT 
  public  ::  DIOS_DOUBLE
  public  ::  DIOS_DATATYPE_NAME
  
  public  ::  DIOS_GLOBAL
  public  ::  DIOS_UNLIMITED
  

  !
  ! Creation and open modes
  !
  !   DIOS_NEW     : new file, error if already present
  !   DIOS_REPLACE : new file, overwrite older file if necessary
  !   DIOS_READ    : open existing file for reading
  !   DIOS_WRITE   : open existing file for writing
  !
  
  integer, parameter          ::  DIOS_NEW     = 1
  integer, parameter          ::  DIOS_REPLACE = 2
  integer, parameter          ::  DIOS_READ    = 3
  integer, parameter          ::  DIOS_WRITE   = 4
  !
  integer, parameter          ::  DIOS_CMODE_MAX = DIOS_WRITE
  character(len=*), parameter ::  DIOS_CMODE_NAME(1:DIOS_CMODE_MAX) = &
                                    (/ 'new    ', 'replace', 'read   ', 'write  ' /)

  !
  ! File types
  !
  !   DIOS_ASCII     :  text file
  !   DIOS_NETCDF    :  NetCDF  (clasical format ; via NetCDF-3 or NetCDF-4 library)
  !
  
  integer, parameter          ::  DIOS_ASCII     = 1
  integer, parameter          ::  DIOS_NETCDF    = 2
  !
  integer, parameter          ::  DIOS_FILETYPE_MAX = DIOS_NETCDF
  character(len=*), parameter ::  DIOS_FILETYPE_NAME(1:DIOS_FILETYPE_MAX) = &
                                    (/ 'ASCII  ', 'NetCDF '/)

  !
  ! Data types
  !
  
  integer, parameter          ::  DIOS_CHAR   = 1  ! character
  integer, parameter          ::  DIOS_BYTE   = 2  ! integer(1)
  integer, parameter          ::  DIOS_SHORT  = 3  ! integer(2)
  integer, parameter          ::  DIOS_INT    = 4  ! integer(4)
  integer, parameter          ::  DIOS_FLOAT  = 5  ! real(4)
  integer, parameter          ::  DIOS_DOUBLE = 6  ! real(8)
  integer, parameter          ::  DIOS_INT64  = 7  ! integer(8)
  !
  integer, parameter          ::  DIOS_DATATYPE_MAX = DIOS_INT64
  character(len=*), parameter ::  DIOS_DATATYPE_NAME(1:DIOS_DATATYPE_MAX) = &
                                    (/ 'char  ','byte  ','short ','int   ','float ', 'double', 'int64 ' /)

  !
  ! special parameters
  !
  
  ! dummy ...
  integer, parameter          ::  DIOS_NONE  = -100

  ! special 'variable id' to add global attributes:
  integer, parameter          ::  DIOS_GLOBAL     = -101
  
  ! special dimension 'length' to denote unlimited dimension:
  integer, parameter          ::  DIOS_UNLIMITED  = -102
  
  
  !
  ! internal
  !
  
  ! maximum rank of Fortran arrays:
  integer, parameter    ::  MAX_RANK = 7
  
  ! maximum length for variable names etc:
  integer, parameter    ::  LEN_NAME = 64
  integer, parameter    ::  LEN_FILE = 512
  integer, parameter    ::  LEN_LINE = 4000


  ! --- types ----------------------------------------

  ! interface to DIOS dimension
  
  type dios_dim
    character(len=LEN_NAME)     ::  name
    integer                     ::  length
    logical                     ::  unlimited
    logical                     ::  named
    integer                     ::  netcdf_dimid
  end type dios_dim

  ! Define a structure with a pointer to the type,
  ! this is necessary to create a list of pointers.
  ! Because you can not only deallocate one element of a pointed array, e.g.,
  ! deallocate(dios_dim_list%item(3)) does not work.
  type p_dios_dim
    type(dios_dim), pointer      ::  p
  end type p_dios_dim

  ! define storage type for list with pointers:
  type dios_dim_list
    ! array of pointers; flexible size, increased if necessary
    type(p_dios_dim), pointer ::  item(:)
    ! maximum number of filled items:
    integer                   ::  maxitem
    ! actual number of filled items:
    integer                   ::  nitem
  end type dios_dim_list


  ! interface to DIOS variable
  
  type dios_var
    ! standard fields:
    character(len=LEN_NAME)     ::  name
    integer                     ::  xtype
    integer                     ::  xkind
    integer                     ::  ndim
    integer                     ::  dimids(MAX_RANK)
    integer                     ::  shp(MAX_RANK)
    integer                     ::  natt
    integer                     ::  netcdf_varid
  end type dios_var

  ! Define a structure with a pointer to the type;
  ! this is necessary to create a list of pointers:
  type p_dios_var
    type(dios_var), pointer      ::  p
  end type p_dios_var

  ! define storage type for list with pointers:
  type dios_var_list
    ! array of pointers; flexible size, increased if necessary
    type(p_dios_var), pointer    ::  item(:)
    ! maximum number of filled items:
    integer                         ::  maxitem
    ! actual number of filled items:
    integer                         ::  nitem
  end type dios_var_list


  ! interface to io file
  
  type dios_file
    ! name of the file
    character(len=LEN_FILE)     ::  filename
    ! creation mode:
    integer                     ::  cmode
    ! dimensions:
    type(dios_dim_List)         ::  dim_list
    ! variables:
    type(dios_var_List)         ::  var_list
    ! number of global attributes:
    integer                     ::  natt
    ! file type:
    integer                     ::  ftype
    ! access to file types:
    integer                     ::  netcdf_id
  end type dios_file

  ! Define a structure with a pointer to the type;
  ! this is necessary to create a list of pointers:
  type p_dios_file
    type(dios_file), pointer      ::  p
  end type p_dios_file

  ! define storage type for list with pointers:
  type dios_file_list
    ! array of pointers; flexible size, increased if necessary
    type(p_dios_file), pointer    ::  item(:)
    ! maximum number of filled items:
    integer                         ::  maxitem
    ! actual number of filled items:
    integer                         ::  nitem
  end type dios_file_list



  ! --- interfaces -----------------------------------
  
  interface DIOS_Put_Var
    module procedure DIOS_Put_Var_c1_1d
    module procedure DIOS_Put_Var_c1_2d
    module procedure DIOS_Put_Var_c1_3d
    module procedure DIOS_Put_Var_c1_4d
    module procedure DIOS_Put_Var_c1_5d
    module procedure DIOS_Put_Var_c1_6d
    module procedure DIOS_Put_Var_c1_7d
    !
    module procedure DIOS_Put_Var_i1_1d
    module procedure DIOS_Put_Var_i1_2d
    module procedure DIOS_Put_Var_i1_3d
    module procedure DIOS_Put_Var_i1_4d
    module procedure DIOS_Put_Var_i1_5d
    module procedure DIOS_Put_Var_i1_6d
    module procedure DIOS_Put_Var_i1_7d
    !
    module procedure DIOS_Put_Var_i2_1d
    module procedure DIOS_Put_Var_i2_2d
    module procedure DIOS_Put_Var_i2_3d
    module procedure DIOS_Put_Var_i2_4d
    module procedure DIOS_Put_Var_i2_5d
    module procedure DIOS_Put_Var_i2_6d
    module procedure DIOS_Put_Var_i2_7d
    !
    module procedure DIOS_Put_Var_i4_1d
    module procedure DIOS_Put_Var_i4_2d
    module procedure DIOS_Put_Var_i4_3d
    module procedure DIOS_Put_Var_i4_4d
    module procedure DIOS_Put_Var_i4_5d
    module procedure DIOS_Put_Var_i4_6d
    module procedure DIOS_Put_Var_i4_7d
    !
    module procedure DIOS_Put_Var_r4_1d
    module procedure DIOS_Put_Var_r4_2d
    module procedure DIOS_Put_Var_r4_3d
    module procedure DIOS_Put_Var_r4_4d
    module procedure DIOS_Put_Var_r4_5d
    module procedure DIOS_Put_Var_r4_6d
    module procedure DIOS_Put_Var_r4_7d
    !
    module procedure DIOS_Put_Var_r8_1d
    module procedure DIOS_Put_Var_r8_2d
    module procedure DIOS_Put_Var_r8_3d
    module procedure DIOS_Put_Var_r8_4d
    module procedure DIOS_Put_Var_r8_5d
    module procedure DIOS_Put_Var_r8_6d
    module procedure DIOS_Put_Var_r8_7d
  end interface
  
  interface DIOS_Get_Var
    module procedure DIOS_Get_Var_c1_1d
    module procedure DIOS_Get_Var_c1_2d
    module procedure DIOS_Get_Var_c1_3d
    module procedure DIOS_Get_Var_c1_4d
    module procedure DIOS_Get_Var_c1_5d
    module procedure DIOS_Get_Var_c1_6d
    module procedure DIOS_Get_Var_c1_7d
    !
    module procedure DIOS_Get_Var_i1_1d
    module procedure DIOS_Get_Var_i1_2d
    module procedure DIOS_Get_Var_i1_3d
    module procedure DIOS_Get_Var_i1_4d
    module procedure DIOS_Get_Var_i1_5d
    module procedure DIOS_Get_Var_i1_6d
    module procedure DIOS_Get_Var_i1_7d
    !
    module procedure DIOS_Get_Var_i2_1d
    module procedure DIOS_Get_Var_i2_2d
    module procedure DIOS_Get_Var_i2_3d
    module procedure DIOS_Get_Var_i2_4d
    module procedure DIOS_Get_Var_i2_5d
    module procedure DIOS_Get_Var_i2_6d
    module procedure DIOS_Get_Var_i2_7d
    !
    module procedure DIOS_Get_Var_i4_1d
    module procedure DIOS_Get_Var_i4_2d
    module procedure DIOS_Get_Var_i4_3d
    module procedure DIOS_Get_Var_i4_4d
    module procedure DIOS_Get_Var_i4_5d
    module procedure DIOS_Get_Var_i4_6d
    module procedure DIOS_Get_Var_i4_7d
    !
    module procedure DIOS_Get_Var_r4_1d
    module procedure DIOS_Get_Var_r4_2d
    module procedure DIOS_Get_Var_r4_3d
    module procedure DIOS_Get_Var_r4_4d
    module procedure DIOS_Get_Var_r4_5d
    module procedure DIOS_Get_Var_r4_6d
    module procedure DIOS_Get_Var_r4_7d
    !
    module procedure DIOS_Get_Var_r8_1d
    module procedure DIOS_Get_Var_r8_2d
    module procedure DIOS_Get_Var_r8_3d
    module procedure DIOS_Get_Var_r8_4d
    module procedure DIOS_Get_Var_r8_5d
    module procedure DIOS_Get_Var_r8_6d
    module procedure DIOS_Get_Var_r8_7d
  end interface

  interface DIOS_Put_Att
    module procedure DIOS_Put_Att_c1_0d
    module procedure DIOS_Put_Att_i1_0d
    module procedure DIOS_Put_Att_i1_1d
    module procedure DIOS_Put_Att_i2_0d
    module procedure DIOS_Put_Att_i2_1d
    module procedure DIOS_Put_Att_i4_0d
    module procedure DIOS_Put_Att_i4_1d
    module procedure DIOS_Put_Att_r4_0d
    module procedure DIOS_Put_Att_r4_1d
    module procedure DIOS_Put_Att_r8_0d
    module procedure DIOS_Put_Att_r8_1d
  end interface

  interface DIOS_Get_Att
    module procedure DIOS_Get_Att_c1_0d
    module procedure DIOS_Get_Att_i1_0d
    module procedure DIOS_Get_Att_i1_1d
    module procedure DIOS_Get_Att_i2_0d
    module procedure DIOS_Get_Att_i2_1d
    module procedure DIOS_Get_Att_i4_0d
    module procedure DIOS_Get_Att_i4_1d
    module procedure DIOS_Get_Att_r4_0d
    module procedure DIOS_Get_Att_r4_1d
    module procedure DIOS_Get_Att_r8_0d
    module procedure DIOS_Get_Att_r8_1d
  end interface

  
  ! --- var ------------------------------------------

  ! define lists:
  type(DIOS_File_List)     ::  File_List
  

contains


  ! ********************************************************************
  ! ***
  ! *** DIOS_Dim procedures
  ! ***
  ! ********************************************************************


  !
  ! Initialise a list.
  !

  subroutine dios_dim_list_init( list, stat )

    ! --- in/out -------------------------------------

    type(dios_dim_list), intent(out)  ::  list
    integer, intent(out)              ::  stat

    ! --- const --------------------------------------


    ! --- begin --------------------------------------

    ! empty list:
    nullify( list%item )

    ! set counters:
    list%maxitem = 0
    list%nitem = 0

    ! ok
    stat = 0

  end subroutine dios_dim_list_init


  ! ***


  !
  ! Clear list, deallocate content.
  !

  subroutine dios_dim_list_done( list, stat )

    ! --- in/out -------------------------------------

    type(dios_dim_list), intent(inout)  ::  list
    integer, intent(out)                ::  stat

    ! --- const --------------------------------------


    ! --- local --------------------------------------

    integer           ::  i

    ! --- begin --------------------------------------

    ! list defined ?
    if ( associated(list%item) ) then
      ! loop over all possible indices:
      do i = 1, list%maxitem
        ! filled ?
        if ( associated(list%item(i)%p) ) then
          ! remove structure, reset to save value:
          deallocate( list%item(i)%p )
          nullify( list%item(i)%p )
        end if
      end do
      ! clear, reset to save value:
      deallocate( list%item )
      nullify( list%item )
    end if

    ! set counters:
    list%maxitem = 0
    list%nitem = 0

    ! ok
    stat = 0

  end subroutine dios_dim_list_done


  ! ***


  !
  ! Add new item to list, return id number.
  !

  subroutine dios_dim_list_new_item( list, hid, stat )

    ! --- in/out -------------------------------------

    type(dios_dim_list), intent(inout)  ::  list
    integer, intent(out)                ::  hid
    integer, intent(out)                ::  stat

    ! --- const --------------------------------------


    ! --- local --------------------------------------

    integer                          ::  i
    type(p_dios_dim), pointer        ::  item_new(:)

    ! --- begin --------------------------------------

    ! free item available ?
    if ( list%nitem < list%maxitem ) then
      ! search first empty item:
      hid = -1
      do i = 1, list%maxitem
        if ( .not. associated(list%item(i)%p) ) then
          hid = i
          exit
        end if
      end do
      ! not found ?
      if ( hid < 0 ) then
        write(*, '("all items seem to be associated while counters suggest something else ...")')
        write(*, '("  maxitem : ",i6)') list%maxitem
        write(*, '("  nitem   : ",i6)') list%nitem
        stat=1; return
      end if
    else
      ! allocate extra space:
      allocate( item_new(1:list%maxitem+100) )
      ! copy old pointers:
      do i = 1, list%maxitem
        item_new(i)%p => list%item(i)%p
      end do
      ! init new pointers:
      do i = list%maxitem+1, size(item_new)
        nullify(item_new(i)%p)
      end do
      ! first empty item:
      hid = list%maxitem+1
      ! clear old list if necessary:
      if ( associated(list%item) ) deallocate( list%item )
      ! point to new list:
      list%item => item_new
      ! reset size counter:
      list%maxitem = size(list%item)
      ! clear:
      nullify( item_new )
    end if

    ! allocate structure:
    allocate( list%item(hid)%p )

    ! increase counter:
    list%nitem = list%nitem + 1

    ! ok
    stat = 0

  end subroutine dios_dim_list_new_item


  ! ***


  !
  ! Remove item with given id from list.
  !

  subroutine DIOS_Dim_List_Clear_Item( list, hid, stat )

    ! --- in/out -------------------------------------

    type(DIOS_Dim_List), intent(inout)     ::  list
    integer, intent(inout)                    ::  hid
    integer, intent(out)                      ::  stat

    ! --- const --------------------------------------


    ! --- local --------------------------------------

    ! --- begin --------------------------------------

    ! check index in list ...
    if ( (hid < 0) .or. (hid > list%maxitem) ) then
      write(*, '("handle outside range:")')
      write(*, '("  handle   : ",i6)') hid
      write(*, '("  range    : ",2i6)') 1, list%maxitem
      stat=1; return
    end if

    ! check ...
    if ( .not. associated(list%item(hid)%p) ) then
      write(*, '("handle not in use: ",i6)') hid
      stat=1; return
    end if

    ! clear structure:
    deallocate( list%item(hid)%p )
    ! reset pointer to save value:
    nullify( list%item(hid)%p )

    ! reset counter:
    list%nitem = list%nitem - 1

    ! ok
    stat = 0

  end subroutine DIOS_Dim_List_Clear_Item


  ! ***


  !
  ! Return pointer to user type given id.
  ! stat -1 if id is not in use.
  !

  subroutine DIOS_Dim_List_Get_Pointer( list, hid, p, stat, silent )

    ! --- in/out -------------------------------------

    type(DIOS_Dim_List), intent(inout)     ::  list
    integer, intent(in)                       ::  hid
    type(DIOS_Dim), pointer                ::  p
    integer, intent(out)                      ::  stat

    logical, intent(in), optional             ::  silent

    ! --- const --------------------------------------


    ! --- local --------------------------------------

    logical     ::  shout

    ! --- begin --------------------------------------

    ! messages ?
    shout = .true.
    if ( present(silent) ) shout = .not. silent

    ! check index in list ...
    if ( (hid < 1) .or. (hid > list%maxitem) ) then
      write(*, '("handle outside range:")')
      write(*, '("  handle   : ",i6)') hid
      write(*, '("  range    : ",2i6)') 1, list%maxitem
      stat=1; return
    end if

    ! check if handle is in use ...
    if ( .not. associated(list%item(hid)%p) ) then
      ! error or warning ?
      if ( shout ) then
        ! error stat:
        write(*, '("handle not in use: ",i6)') hid
        stat=1; return
      else
        ! warning stat; this routine is used to test if a handle is in use:
        nullify( p )
        stat = -1 ; return
      end if
    end if

    ! set shorthand:
    p => list%item(hid)%p

    ! ok
    stat = 0

  end subroutine DIOS_Dim_List_Get_Pointer


  ! ***


  !
  ! Return information:
  !   n
  !           Number of elements in use.
  !   maxid
  !           Current possible upper value for id's.
  !           Not all id's in {1,..,maxid} are in use.
  !           Usefull to implement a loop over all possible items.
  !

  subroutine DIOS_Dim_List_Inquire( list, stat, &
                                                  n, maxid )

    ! --- in/out -------------------------------------

    type(DIOS_Dim_List), intent(inout)     ::  list
    integer, intent(out)                      ::  stat

    integer, intent(out), optional            ::  n
    integer, intent(out), optional            ::  maxid

    ! --- const --------------------------------------


    ! --- begin --------------------------------------

    ! set values ?
    if ( present(n    ) ) n     = list%nitem
    if ( present(maxid) ) maxid = list%maxitem

    ! ok
    stat = 0

  end subroutine DIOS_Dim_List_Inquire




  ! ********************************************************************
  ! ***
  ! *** DIOS_Var procedures
  ! ***
  ! ********************************************************************
  


  !
  ! Initialise a list.
  !

  subroutine dios_var_list_init( list, stat )

    ! --- in/out -------------------------------------

    type(dios_var_list), intent(out)  ::  list
    integer, intent(out)              ::  stat

    ! --- const --------------------------------------


    ! --- begin --------------------------------------

    ! empty list:
    nullify( list%item )

    ! set counters:
    list%maxitem = 0
    list%nitem = 0

    ! ok
    stat = 0

  end subroutine dios_var_list_init


  ! ***


  !
  ! Clear list, deallocate content.
  !

  subroutine dios_var_list_done( list, stat )

    ! --- in/out -------------------------------------

    type(dios_var_list), intent(inout)  ::  list
    integer, intent(out)                ::  stat

    ! --- const --------------------------------------


    ! --- local --------------------------------------

    integer           ::  i

    ! --- begin --------------------------------------

    ! list defined ?
    if ( associated(list%item) ) then
      ! loop over all possible indices:
      do i = 1, list%maxitem
        ! filled ?
        if ( associated(list%item(i)%p) ) then
          ! remove structure, reset to save value:
          deallocate( list%item(i)%p )
          nullify( list%item(i)%p )
        end if
      end do
      ! clear, reset to save value:
      deallocate( list%item )
      nullify( list%item )
    end if

    ! set counters:
    list%maxitem = 0
    list%nitem = 0

    ! ok
    stat = 0

  end subroutine dios_var_list_done


  ! ***


  !
  ! Add new item to list, return id number.
  !

  subroutine DIOS_Var_List_New_Item( list, hid, stat )

    ! --- in/out -------------------------------------

    type(DIOS_Var_List), intent(inout)     ::  list
    integer, intent(out)                      ::  hid
    integer, intent(out)                      ::  stat

    ! --- const --------------------------------------


    ! --- local --------------------------------------

    integer                             ::  i
    type(P_DIOS_Var), pointer        ::  item_new(:)

    ! --- begin --------------------------------------

    ! free item available ?
    if ( list%nitem < list%maxitem ) then
      ! search first empty item:
      hid = -1
      do i = 1, list%maxitem
        if ( .not. associated(list%item(i)%p) ) then
          hid = i
          exit
        end if
      end do
      ! not found ?
      if ( hid < 0 ) then
        write(*, '("all items seem to be associated while counters suggest something else ...")')
        write(*, '("  maxitem : ",i6)') list%maxitem
        write(*, '("  nitem   : ",i6)') list%nitem
        stat=1; return
      end if
    else
      ! allocate extra space:
      allocate( item_new(1:list%maxitem+100) )
      ! copy old pointers:
      do i = 1, list%maxitem
        item_new(i)%p => list%item(i)%p
      end do
      ! init new pointers:
      do i = list%maxitem+1, size(item_new)
        nullify(item_new(i)%p)
      end do
      ! first empty item:
      hid = list%maxitem+1
      ! clear old list if necessary:
      if ( associated(list%item) ) deallocate( list%item )
      ! point to new list:
      list%item => item_new
      ! reset size counter:
      list%maxitem = size(list%item)
      ! clear:
      nullify( item_new )
    end if

    ! allocate structure:
    allocate( list%item(hid)%p )

    ! increase counter:
    list%nitem = list%nitem + 1

    ! ok
    stat = 0

  end subroutine DIOS_Var_List_New_Item


  ! ***


  !
  ! Remove item with given id from list.
  !

  subroutine DIOS_Var_List_Clear_Item( list, hid, stat )

    ! --- in/out -------------------------------------

    type(DIOS_Var_List), intent(inout)     ::  list
    integer, intent(inout)                    ::  hid
    integer, intent(out)                      ::  stat

    ! --- const --------------------------------------


    ! --- local --------------------------------------

    ! --- begin --------------------------------------

    ! check index in list ...
    if ( (hid < 0) .or. (hid > list%maxitem) ) then
      write(*, '("handle outside range:")')
      write(*, '("  handle   : ",i6)') hid
      write(*, '("  range    : ",2i6)') 1, list%maxitem
      stat=1; return
    end if

    ! check ...
    if ( .not. associated(list%item(hid)%p) ) then
      write(*, '("handle not in use: ",i6)') hid
      stat=1; return
    end if

    ! clear structure:
    deallocate( list%item(hid)%p )
    ! reset pointer to save value:
    nullify( list%item(hid)%p )

    ! reset counter:
    list%nitem = list%nitem - 1

    ! ok
    stat = 0

  end subroutine DIOS_Var_List_Clear_Item


  ! ***


  !
  ! Return pointer to user type given id.
  ! stat -1 if id is not in use.
  !

  subroutine DIOS_Var_List_Get_Pointer( list, hid, p, stat, silent )

    ! --- in/out -------------------------------------

    type(DIOS_Var_List), intent(inout)     ::  list
    integer, intent(in)                       ::  hid
    type(DIOS_Var), pointer                ::  p
    integer, intent(out)                      ::  stat

    logical, intent(in), optional             ::  silent

    ! --- const --------------------------------------


    ! --- local --------------------------------------

    logical     ::  shout

    ! --- begin --------------------------------------

    ! messages ?
    shout = .true.
    if ( present(silent) ) shout = .not. silent

    ! check index in list ...
    if ( (hid < 1) .or. (hid > list%maxitem) ) then
      write(*, '("handle outside range:")')
      write(*, '("  handle   : ",i6)') hid
      write(*, '("  range    : ",2i6)') 1, list%maxitem
      stat=1; return
    end if

    ! check if handle is in use ...
    if ( .not. associated(list%item(hid)%p) ) then
      ! error or warning ?
      if ( shout ) then
        ! error stat:
        write(*, '("handle not in use: ",i6)') hid
        stat=1; return
      else
        ! warning stat; this routine is used to test if a handle is in use:
        nullify( p )
        stat = -1 ; return
      end if
    end if

    ! set shorthand:
    p => list%item(hid)%p

    ! ok
    stat = 0

  end subroutine DIOS_Var_List_Get_Pointer


  ! ***


  !
  ! Return information:
  !   n
  !           Number of elements in use.
  !   maxid
  !           Current possible upper value for id's.
  !           Not all id's in {1,..,maxid} are in use.
  !           Usefull to implement a loop over all possible items.
  !

  subroutine DIOS_Var_List_Inquire( list, stat, &
                                                  n, maxid )

    ! --- in/out -------------------------------------

    type(DIOS_Var_List), intent(inout)     ::  list
    integer, intent(out)                      ::  stat

    integer, intent(out), optional            ::  n
    integer, intent(out), optional            ::  maxid

    ! --- const --------------------------------------


    ! --- begin --------------------------------------

    ! set values ?
    if ( present(n    ) ) n     = list%nitem
    if ( present(maxid) ) maxid = list%maxitem

    ! ok
    stat = 0

  end subroutine DIOS_Var_List_Inquire




  ! ********************************************************************
  ! ***
  ! *** DIOS procedures
  ! ***
  ! ********************************************************************
  


  !
  ! Initialise a list.
  !

  subroutine dios_file_list_init( list, stat )

    ! --- in/out -------------------------------------

    type(dios_file_list), intent(out)  ::  list
    integer, intent(out)               ::  stat

    ! --- const --------------------------------------


    ! --- begin --------------------------------------

    ! empty list:
    nullify( list%item )

    ! set counters:
    list%maxitem = 0
    list%nitem = 0

    ! ok
    stat = 0

  end subroutine dios_file_list_init


  ! ***


  !
  ! Clear list, deallocate content.
  !

  subroutine dios_file_list_done( list, stat )

    ! --- in/out -------------------------------------

    type(dios_file_list), intent(inout)  ::  list
    integer, intent(out)                 ::  stat

    ! --- const --------------------------------------


    ! --- local --------------------------------------

    integer           ::  i

    ! --- begin --------------------------------------

    ! list defined ?
    if ( associated(list%item) ) then
      ! loop over all possible indices:
      do i = 1, list%maxitem
        ! filled ?
        if ( associated(list%item(i)%p) ) then
          ! remove structure, reset to save value:
          deallocate( list%item(i)%p )
          nullify( list%item(i)%p )
        end if
      end do
      ! clear, reset to save value:
      deallocate( list%item )
      nullify( list%item )
    end if

    ! set counters:
    list%maxitem = 0
    list%nitem = 0

    ! ok
    stat = 0

  end subroutine dios_file_list_done


  ! ***


  !
  ! Add new item to list, return id number.
  !

  subroutine DIOS_File_List_New_Item( list, hid, stat )

    ! --- in/out -------------------------------------

    type(DIOS_File_List), intent(inout)     ::  list
    integer, intent(out)                      ::  hid
    integer, intent(out)                      ::  stat

    ! --- const --------------------------------------


    ! --- local --------------------------------------

    integer                             ::  i
    type(P_DIOS_File), pointer        ::  item_new(:)

    ! --- begin --------------------------------------

    ! free item available ?
    if ( list%nitem < list%maxitem ) then
      ! search first empty item:
      hid = -1
      do i = 1, list%maxitem
        if ( .not. associated(list%item(i)%p) ) then
          hid = i
          exit
        end if
      end do
      ! not found ?
      if ( hid < 0 ) then
        write(*, '("all items seem to be associated while counters suggest something else ...")')
        write(*, '("  maxitem : ",i6)') list%maxitem
        write(*, '("  nitem   : ",i6)') list%nitem
        stat=1; return
      end if
    else
      ! allocate extra space:
      allocate( item_new(1:list%maxitem+100) )
      ! copy old pointers:
      do i = 1, list%maxitem
        item_new(i)%p => list%item(i)%p
      end do
      ! init new pointers:
      do i = list%maxitem+1, size(item_new)
        nullify(item_new(i)%p)
      end do
      ! first empty item:
      hid = list%maxitem+1
      ! clear old list if necessary:
      if ( associated(list%item) ) deallocate( list%item )
      ! point to new list:
      list%item => item_new
      ! reset size counter:
      list%maxitem = size(list%item)
      ! clear:
      nullify( item_new )
    end if

    ! allocate structure:
    allocate( list%item(hid)%p )

    ! increase counter:
    list%nitem = list%nitem + 1

    ! ok
    stat = 0

  end subroutine DIOS_File_List_New_Item


  ! ***


  !
  ! Remove item with given id from list.
  !

  subroutine dios_file_list_clear_item( list, hid, stat )

    ! --- in/out -------------------------------------

    type(dios_file_list), intent(inout)  ::  list
    integer, intent(inout)               ::  hid
    integer, intent(out)                 ::  stat

    ! --- const --------------------------------------


    ! --- local --------------------------------------

    ! --- begin --------------------------------------

    ! check index in list ...
    if ( (hid < 0) .or. (hid > list%maxitem) ) then
      write(*, '("handle outside range:")')
      write(*, '("  handle   : ",i6)') hid
      write(*, '("  range    : ",2i6)') 1, list%maxitem
      stat=1; return
    end if

    ! check ...
    if ( .not. associated(list%item(hid)%p) ) then
      write(*, '("handle not in use: ",i6)') hid
      stat=1; return
    end if

    ! clear structure:
    deallocate( list%item(hid)%p )
    ! reset pointer to save value:
    nullify( list%item(hid)%p )

    ! reset counter:
    list%nitem = list%nitem - 1

    ! ok
    stat = 0

  end subroutine dios_file_list_clear_item


  ! ***


  !
  ! Return pointer to user type given id.
  ! stat -1 if id is not in use.
  !

  subroutine dios_file_list_get_pointer( list, hid, p, stat, silent )

    ! --- in/out -------------------------------------

    type(dios_file_list), intent(inout) :: list
    integer, intent(in)                 :: hid
    type(dios_file), pointer            :: p
    integer, intent(out)                :: stat

    logical, intent(in), optional       :: silent

    ! --- const --------------------------------------


    ! --- local --------------------------------------

    logical     ::  shout

    ! --- begin --------------------------------------

    ! messages ?
    shout = .true.
    if ( present(silent) ) shout = .not. silent

    ! check index in list ...
    if ( (hid < 1) .or. (hid > list%maxitem) ) then
      write(*, '("handle outside range:")')
      write(*, '("  handle   : ",i6)') hid
      write(*, '("  range    : ",2i6)') 1, list%maxitem
      stat=1; return
    end if

    ! check if handle is in use ...
    if ( .not. associated(list%item(hid)%p) ) then
      ! error or warning ?
      if ( shout ) then
        ! error stat:
        write(*, '("handle not in use: ",i6)') hid
        stat=1; return
      else
        ! warning stat; this routine is used to test if a handle is in use:
        nullify( p )
        stat = -1 ; return
      end if
    end if

    ! set shorthand:
    p => list%item(hid)%p

    ! ok
    stat = 0

  end subroutine dios_file_list_get_pointer


  ! ***


  !
  ! Return information:
  !   n
  !           Number of elements in use.
  !   maxid
  !           Current possible upper value for id's.
  !           Not all id's in {1,..,maxid} are in use.
  !           Usefull to implement a loop over all possible items.
  !

  subroutine dios_file_list_inquire( list, stat, n, maxid )

    ! --- in/out -------------------------------------

    type(dios_file_list), intent(inout)     ::  list
    integer, intent(out)                    ::  stat

    integer, intent(out), optional          ::  n
    integer, intent(out), optional          ::  maxid

    ! --- const --------------------------------------


    ! --- begin --------------------------------------

    ! set values ?
    if ( present(n    ) ) n     = list%nitem
    if ( present(maxid) ) maxid = list%maxitem

    ! ok
    stat = 0

  end subroutine dios_file_list_inquire



  ! ********************************************************************
  ! ***
  ! *** tools
  ! ***
  ! ********************************************************************


  subroutine dios_get_kind( xtype, xkind, stat )

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  xtype
    integer, intent(out)                ::  xkind
    integer, intent(out)                ::  stat

    ! --- const --------------------------------------


    ! --- begin --------------------------------------

    ! set kind value given type:
    select case ( xtype )
      case ( DIOS_CHAR   ) ;  xkind = 1
      case ( DIOS_BYTE   ) ;  xkind = 1
      case ( DIOS_SHORT  ) ;  xkind = 2
      case ( DIOS_INT    ) ;  xkind = 4
      case ( DIOS_INT64  ) ;  xkind = 8
      case ( DIOS_FLOAT  ) ;  xkind = 4
      case ( DIOS_DOUBLE ) ;  xkind = 8
      case default
        write(*, '("do not know kind for variable type : ",i6)') xtype
        stat=1; return
    end select
    
    ! ok
    stat = 0
    
  end subroutine dios_get_kind
  
  
  ! ***


  subroutine NetCDF_Get_FileType( fname, ncformat, stat )
  
    use NetCDF, only : NF90_Open, NF90_Close, NF90_Inquire
    use NetCDF, only : NF90_NOWRITE
    use NetCDF, only : NF90_FORMAT_CLASSIC, NF90_FORMAT_64BIT, NF90_FORMAT_NETCDF4, NF90_FORMAT_NETCDF4_CLASSIC
  
    ! --- in/out ---------------------------------
    
    character(len=*), intent(in)        ::  fname
    character(len=*), intent(out)       ::  ncformat
    integer, intent(out)                ::  stat

    ! --- const --------------------------------------

    
    ! --- local ----------------------------------
    
    logical       ::  exist
    integer       ::  ncid
    integer       ::  formatNum
    
    ! --- begin ----------------------------------

    ! check ...
    inquire( file=trim(fname), exist=exist )
    if ( .not. exist ) then
      write(*, '("file to be opened not found : ",a)') trim(fname)
      stat=1; return
    end if

    ! open file for reading:
    stat = NF90_Open( trim(fname), NF90_NOWRITE, ncid )
    
    
    ! get format number:
    stat = NF90_Inquire( ncid, formatNum=formatNum )
    
    
    ! translate ...
    select case ( formatNum )
      case ( NF90_FORMAT_CLASSIC         ) ; ncformat = 'netcdf_classic'
      case ( NF90_FORMAT_64BIT           ) ; ncformat = 'netcdf_64bit'
      case ( NF90_FORMAT_NETCDF4         ) ; ncformat = 'netcdf4'
      case ( NF90_FORMAT_NETCDF4_CLASSIC ) ; ncformat = 'netcdf4_classic'
      case default                         ; ncformat = 'netcdf_unknown'
    end select
    
    ! close file:
    stat = NF90_Close( ncid )
    

    ! ok
    stat = 0
    
  end subroutine NetCDF_Get_FileType


  ! ********************************************************************
  ! ***
  ! *** module init/done
  ! ***
  ! ********************************************************************
  

  subroutine dios_init( stat )
    integer, intent(out)                ::  stat

    ! setup empty list:
    call dios_file_list_init( file_list, stat )

    ! ok
    stat = 0

  end subroutine dios_init


  ! ***

 
  subroutine dios_done( stat )

    ! --- in/out -------------------------------------

    integer, intent(out)                ::  stat

    ! --- const --------------------------------------


    ! --- local --------------------------------------
    
    integer                     ::  maxid
    integer                     ::  id
    type(dios_file), pointer    ::  filep
    integer                     ::  nerror

    ! --- begin --------------------------------------
    
    ! no errors yet ...
    nerror = 0
    
    ! get maximum id number:
    call dios_file_list_inquire( file_list, stat, maxid=maxid )
    
    ! loop over all possible id's:
    do id = 1, maxid
      ! get pointer to file structure; stat -1 if not in use:
      call dios_file_list_get_pointer( file_list, id, filep, stat, silent=.true. )
      if ( stat == -1 ) cycle

      ! error ...
      write(*,'("Called dios_done but file still in use: ",a)') trim(filep%filename)
      nerror = nerror + 1
      !! done with variables:
      !call DIOS_Var_List_Done( filep%Var_List, stat )
      !! done with dimensions:
      !call DIOS_Dim_List_Done( filep%Dim_List, stat )
    end do

    ! clear list:
    call dios_file_list_done( file_list, stat )

    ! ok
    stat = nerror

  end subroutine dios_done


  ! ********************************************************************
  ! ***
  ! *** file create/close
  ! ***
  ! ********************************************************************



  ! ***
  

  subroutine dios_create( filename, ftype, cmode, hid, stat )

    use NetCDF, only : NF90_NOERR, NF90_CLOBBER, NF90_NOCLOBBER
    use NetCDF, only : NF90_Create
    use NetCDF, only : NF90_CLASSIC_MODEL, NF90_NETCDF4
    use NetCDF, only : NF90_Inq_LibVers
    use NetCDF, only : nf90_strerror

    ! --- in/out -------------------------------------

    character(len=*), intent(in)        ::  filename
    integer, intent(in)                 ::  ftype
    integer, intent(in)                 ::  cmode
    integer, intent(out)                ::  hid
    integer, intent(out)                ::  stat
    
    ! --- const --------------------------------------

    
    ! --- external ----------------------------


    ! --- local --------------------------------------

    type(dios_file), pointer          ::  filep
    integer                           ::  netcdf_cmode
    character(len=80)                 ::  netcdf_version

    ! --- begin --------------------------------------

    ! new file:
    call dios_file_list_new_item( file_list, hid, stat )

    ! pointer to file structure:
    call dios_file_list_get_pointer( file_list, hid, filep, stat )

    ! store filename stuff:
    filep%filename = trim(filename)
    
    ! store creation mode:
    filep%cmode = cmode
    
    ! store file types:
    filep%ftype = ftype
    
    ! select appropriate routine for each type:
    select case ( filep%ftype )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_ASCII)
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_NETCDF )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! initial creation mode:
        netcdf_cmode = 0
        ! set creation mode:
        select case ( cmode )
          case ( DIOS_NEW )
            netcdf_cmode = NF90_NOCLOBBER    ! complain if already present
          case ( DIOS_REPLACE )
            netcdf_cmode = NF90_CLOBBER      ! overwrite if necessary
          case default
            write(*,'("unsupported creation mode : ",i6)') cmode
            stat=1; return
        end select

        ! create file:
        stat = NF90_Create( trim(filep%filename), netcdf_cmode, filep%netcdf_id )
        if (stat/=NF90_NOERR) then
          write(*,*) trim(nf90_strerror(stat))
          write(*, '("from creating netcdf4 file :")')
          write(*, '("  ",a)') trim(filep%filename)
          write(*, '("  does directory exist ?")')
          stat=1; return
        end if

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*, '("unsupported filetype : ",i6)') filep%ftype
          write(*, '("  filetype name : ",a)') trim(DIOS_FILETYPE_NAME(ftype))
          write(*, '("  compiled without macro `with_*` defined ?")')
          stat=1; return
    end select
    
    ! init dimension list:
    call dios_dim_list_init( filep%dim_list, stat )
    
    ! init variable list:
    call dios_var_list_init( filep%var_list, stat )
    
    ! no global attributes yet:
    filep%natt = 0
        
    ! ok
    stat = 0

  end subroutine dios_create


  ! ***


  subroutine dios_open( filename, ftype, mode, hid, stat )

    use NetCDF, only : NF90_WRITE, NF90_NOWRITE
    use NetCDF, only : NF90_Open
    use NetCDF, only : NF90_Inquire
    use NetCDF, only : NF90_Inquire_Dimension
    use NetCDF, only : NF90_Inquire_Variable
    use NetCDF, only : NF90_CHAR, NF90_BYTE, NF90_SHORT, NF90_INT, NF90_INT64, NF90_FLOAT, NF90_DOUBLE

    ! --- in/out -------------------------------------

    character(len=*), intent(in)        ::  filename
    integer, intent(in)                 ::  ftype
    integer, intent(in)                 ::  mode
    integer, intent(out)                ::  hid
    integer, intent(out)                ::  stat
    
    ! --- const --------------------------------------

    
    ! --- external ----------------------------


    ! --- local --------------------------------------

    type(DIOS_File), pointer          ::  filep
    type(DIOS_Dim), pointer           ::  dimp
    type(DIOS_Var), pointer           ::  varp
    logical                           ::  exist

    integer                           ::  netcdf_mode
    integer                           ::  netcdf_xtype
    integer                           ::  unlimid

    integer                           ::  ndim, idim, dimid
    integer                           ::  nvar, ivar, varid
    integer                           ::  natt
    character(len=LEN_NAME)           ::  name
    integer                           ::  length
    integer                           ::  dimids(MAX_RANK)
    integer                           ::  shp(MAX_RANK)
    integer                           ::  k, n
    character(len=80)                 ::  netcdf_version

    ! --- begin --------------------------------------

    ! new file:
    call dios_file_list_new_item( file_list, hid, stat )

    ! pointer to file structure:
    call dios_file_list_get_pointer( file_list, hid, filep, stat )

    ! init dimension list:
    call dios_dim_list_init( filep%dim_list, stat )
    
    ! init variable list:
    call dios_var_list_init( filep%var_list, stat )
        
    ! store filename stuff:
    filep%filename = trim(filename)
    
    ! store dummy creation mode:
    filep%cmode = -1

    ! store file type:
    filep%ftype = ftype
    
    ! select appropriate routine for each type:
    select case ( filep%ftype )

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_ASCII )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_NETCDF )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! check ...
        inquire( file=trim(filep%filename), exist=exist )
        if ( .not. exist ) then
          write(*, '("file to be opened not found : ",a)') trim(filep%filename)
          stat=1; return
        end if

        ! set open mode:
        select case ( mode )
          case ( DIOS_READ )
            netcdf_mode = NF90_NOWRITE
          case ( DIOS_WRITE )
            netcdf_mode = NF90_WRITE
          case default
            write(*, '("unsupported creation mode : ",i6)') mode
            stat=1; return
        end select

        ! open file:
        stat = NF90_Open( trim(filep%filename), netcdf_mode, filep%netcdf_id )

        ! get number of global attributes:
        stat = NF90_Inquire( filep%netcdf_id, nAttributes=filep%natt )

        ! get number of dimensions and (dummy) id of unlimitted dimension:
        stat = NF90_Inquire( filep%netcdf_id, nDimensions=ndim, unlimitedDimID=unlimid )

        ! loop over dimensions:
        do idim = 1, ndim
          ! new dimension:
          call dios_dim_list_new_item( filep%dim_list, dimid, stat )
          ! pointer to dimension structure:
          call dios_dim_list_get_pointer( filep%dim_list, dimid, dimp, stat )
          ! netcdf dimension id is number from 1..ndim
          dimp%netcdf_dimid = idim
          ! get info:
          stat = NF90_Inquire_Dimension( filep%netcdf_id, dimp%netcdf_dimid, name=name, len=length )
          ! store:
          dimp%named     = .true.
          dimp%name      = trim(name)
          dimp%length    = length
          dimp%unlimited = dimp%netcdf_dimid == unlimid
        end do

        ! get number of variables:
        stat = NF90_Inquire( filep%netcdf_id, nVariables=nvar )
        ! loop over variables:
        do ivar = 1, nvar
          ! new variable:
          call dios_var_list_new_item( filep%var_list, varid, stat )
          ! pointer to variable structure:
          call dios_var_list_get_pointer( filep%var_list, varid, varp, stat )
          ! netcdf variable id is number from 1..nvar
          varp%netcdf_varid = ivar
          ! get info:
          stat = NF90_Inquire_Variable( filep%netcdf_id, varp%netcdf_varid, &
            name=name, xtype=netcdf_xtype, ndims=ndim )
          ! store name:
          varp%name      = trim(name)
          ! convert type:
          select case ( netcdf_xtype )
            case ( NF90_CHAR   ) ;  varp%xtype = DIOS_CHAR  
            case ( NF90_BYTE   ) ;  varp%xtype = DIOS_BYTE  
            case ( NF90_SHORT  ) ;  varp%xtype = DIOS_SHORT 
            case ( NF90_INT    ) ;  varp%xtype = DIOS_INT   
            case ( NF90_INT64  ) ;  varp%xtype = DIOS_INT64
            case ( NF90_FLOAT  ) ;  varp%xtype = DIOS_FLOAT 
            case ( NF90_DOUBLE ) ;  varp%xtype = DIOS_DOUBLE
            case default
              write(*, '("unsupported data type : ",i6)') netcdf_xtype
              stat=1; return
          end select
          ! set kind given type:
          call dios_get_kind( varp%xtype, varp%xkind, stat )
          ! store number of dimensions:
          varp%ndim = ndim
          ! get netcdf dimension id's now that number is known:
          stat = NF90_Inquire_Variable( filep%netcdf_id, varp%netcdf_varid, dimids=dimids(1:ndim) )
          ! init arrays:
          varp%dimids = -1
          varp%shp    = -1
          ! loop over dimensions:
          do idim = 1, ndim
            ! DIOS dimension id is the same as the netcdf dimension id,
            ! both are numbers 1,..,maxdim :
            dimid = dimids(idim)
            ! pointer to dimension structure:
            call dios_dim_list_get_pointer( filep%dim_list, dimid, dimp, stat )
            ! store:
            varp%dimids(idim) = dimid
            varp%shp   (idim) = dimp%length
          end do
          ! get number of variable attributes:
          stat = nf90_inquire_variable( filep%netcdf_id, varp%netcdf_varid, nAtts=varp%natt )
        end do  ! variables

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        write(*, '("unsupported filetype : ",i6)') filep%ftype
        stat=1; return
    end select
    
    ! ok
    stat = 0

  end subroutine dios_open


  ! ***


  subroutine dios_close( hid, stat )

    use NetCDF, only : NF90_NOERR
    use NetCDF, only : NF90_Close

    ! --- in/out -------------------------------------

    integer, intent(inout)              ::  hid
    integer, intent(out)                ::  stat

    ! --- const --------------------------------------

    
    ! --- external ----------------------------


    ! --- local --------------------------------------

    type(dios_file), pointer           ::  filep
    integer                            ::  iftype
    integer                            ::  ftype
    integer                            ::  ivar, nvar
    type(dios_var), pointer            ::  varp

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call dios_file_list_get_pointer( file_list, hid, filep, stat )

    ! select appropriate routine for each type:
    select case ( filep%ftype )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_ASCII )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_NETCDF )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! close file:
      stat = NF90_Close( filep%netcdf_id )
      if ( stat /= NF90_NOERR ) then
        write(*, '("while closing NetCDF4 file:")')
        write(*, '("  file name : ",a)') trim(filep%filename)
        stat=1; return
      end if

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        write(*, '("unsupported filetype : ",i6)') filep%ftype
        stat=1; return
    end select
    
    ! done with variable list:
    call dios_var_list_done( filep%var_list, stat )

    ! done with dimension list:
    call dios_dim_list_done( filep%dim_list, stat )
        
    ! remove item:
    call dios_file_list_clear_item( file_list, hid, stat )

    ! ok
    stat = 0

  end subroutine dios_close
  
  
  ! ********************************************************************
  ! ***
  ! *** end of definition phase
  ! ***
  ! ********************************************************************


  subroutine dios_enddef( hid, stat )

    use NetCDF, only : NF90_EndDef

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(out)                ::  stat

    ! --- const --------------------------------------


    ! --- local --------------------------------------

    type(dios_file), pointer  ::  filep
    integer                   ::  iftype
    integer                   ::  ftype

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call dios_file_list_get_pointer( file_list, hid, filep, stat )

    ! select appropriate routine for each type:
    select case ( filep%ftype )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_ASCII )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_NETCDF )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! end of definition phase:
        stat = NF90_EndDef( filep%netcdf_id )

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        write(*, '("unsupported filetype : ",i6)') filep%ftype
        stat=1; return
      end select

    ! ok
    stat = 0

  end subroutine dios_enddef


  ! ********************************************************************
  ! ***
  ! *** dimensions
  ! ***
  ! ********************************************************************
  

  subroutine dios_def_dim( hid, name, length, dimid, stat )

    use NetCDF, only : NF90_Def_Dim, NF90_UNLIMITED

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    character(len=*), intent(in)        ::  name
    integer, intent(in)                 ::  length
    integer, intent(out)                ::  dimid
    integer, intent(out)                ::  stat

    ! --- const --------------------------------------


    ! --- local --------------------------------------

    type(dios_file), pointer        ::  filep
    type(dios_dim), pointer         ::  dimp
    integer                         ::  iftype
    integer                         ::  ftype
    integer                         ::  netcdf_length

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call dios_file_list_get_pointer( file_list, hid, filep, stat )

    ! new dimension:
    call dios_dim_list_new_item( filep%dim_list, dimid, stat )

    ! pointer to dimension structure:
    call dios_dim_list_get_pointer( filep%dim_list, dimid, dimp, stat )
    
    ! store:
    dimp%name   = trim(name)
    dimp%length = length

    ! unlimited length ?
    dimp%unlimited = length == DIOS_UNLIMITED

    ! select appropriate routine for each type:
    select case ( filep%ftype )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_ASCII )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_NETCDF )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! set dimension length:
        if ( dimp%unlimited ) then
          netcdf_length = NF90_UNLIMITED
        else
          netcdf_length = length
        end if
        ! define dimension:
        stat = NF90_Def_Dim( filep%netcdf_id, trim(name), netcdf_length, dimp%netcdf_dimid )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        write(*, '("unsupported filetype : ",i6)') filep%ftype
        stat=1; return
    end select
        
    ! ok
    stat = 0

  end subroutine dios_def_dim


  ! ********************************************************************
  ! ***
  ! *** variables
  ! ***
  ! ********************************************************************
  

  subroutine dios_def_var( hid, name, xtype, dimids, varid, stat )

    use NetCDF, only : NF90_CHAR, NF90_BYTE, NF90_SHORT, NF90_INT, NF90_FLOAT, NF90_DOUBLE
    use NetCDF, only : NF90_Def_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    character(len=*), intent(in)        ::  name
    integer, intent(in)                 ::  xtype
    integer, intent(in)                 ::  dimids(:)
    integer, intent(out)                ::  varid
    integer, intent(out)                ::  stat
    
    ! --- const --------------------------------------

    
    ! --- external -----------------------------------


    ! --- local --------------------------------------

    type(dios_file), pointer         ::  filep
    type(dios_dim), pointer          ::  dimp
    type(dios_var), pointer          ::  varp
    integer                         ::  iftype
    integer                         ::  ftype
    integer                         ::  idim

    integer                         ::  netcdf_xtype
    integer                         ::  netcdf_dimids(MAX_RANK)

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call dios_file_list_get_pointer( file_list, hid, filep, stat )
    
    ! new variable:
    call dios_var_list_new_item( filep%var_list, varid, stat )

    ! pointer to variable structure:
    call dios_var_list_get_pointer( filep%var_list, varid, varp, stat )
    
    ! store name;
    varp%name = trim(name)
    
    ! store type:
    varp%xtype = xtype
    
    ! set kind value given type:
    call dios_get_kind( varp%xtype, varp%xkind, stat )
    
    ! number of dimensions:
    varp%ndim = size(dimids)
    
    ! dimension id's :
    varp%dimids(1:varp%ndim) = dimids
    
    ! fill shape:
    do idim = 1, varp%ndim
      ! pointer to  dimension type:
      call dios_dim_list_get_pointer( filep%dim_list, dimids(idim), dimp, stat )
      ! copy dimension id:
      varp%shp(idim) = dimp%length
    end do

    ! select appropriate routine for each type:
    select case ( filep%ftype )

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_ASCII )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_NETCDF )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      
        ! set data type:
        select case ( xtype )
          case ( DIOS_CHAR   ) ;  netcdf_xtype = NF90_CHAR  
          case ( DIOS_BYTE   ) ;  netcdf_xtype = NF90_BYTE  
          case ( DIOS_SHORT  ) ;  netcdf_xtype = NF90_SHORT 
          case ( DIOS_INT    ) ;  netcdf_xtype = NF90_INT   
          case ( DIOS_FLOAT  ) ;  netcdf_xtype = NF90_FLOAT 
          case ( DIOS_DOUBLE ) ;  netcdf_xtype = NF90_DOUBLE
          case default
            write(*, '("unsupported data type : ",i6)') xtype
            stat=1; return
        end select
        
        ! extract dimensions:
        do idim = 1, varp%ndim
          ! pointer to  dimension type:
          call dios_dim_list_get_pointer( filep%dim_list, dimids(idim), dimp, stat )
          ! copy dimension id:
          netcdf_dimids(idim) = dimp%netcdf_dimid
        end do
        
        ! define variable:
        stat = NF90_Def_Var( filep%netcdf_id, trim(name), netcdf_xtype, &
          netcdf_dimids(1:varp%ndim), varp%netcdf_varid )
        
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        write(*, '("unsupported filetype : ",i6)') filep%ftype
        stat=1; return
      end select
    
    ! no attributes yet:
    varp%natt = 0
        
    ! ok
    stat = 0

  end subroutine dios_def_var


  ! ********************************************************************
  ! ***
  ! *** put_var_c1, get_var_c1
  ! ***
  ! ********************************************************************

  
  subroutine dios_put_var_c1_1d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Put_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    character(len=*), intent(in)        ::  values
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(dios_file), pointer        ::  filep
    type(dios_var), pointer         ::  varp

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call dios_file_list_get_pointer( file_list, hid, filep, stat )
    
    ! pointer to variable structure:
    call dios_var_list_get_pointer( filep%var_list, varid, varp, stat )

    ! check ...
    call check_var_values_rank(size(shape(values)), varp%ndim, stat)
    if (stat == 1) return
      
    if ( present(start) ) call check_var_start_size(size(start), varp%ndim, stat)
    if (stat == 1) return

    if ( present(count) ) call check_var_count_size(size(count), varp%ndim, stat)
    if (stat == 1) return

    if ( present(stride) ) call check_var_stride_size(size(stride), varp%ndim, stat)
    if (stat == 1) return

    if ( present(map) ) call check_var_map_size(size(map), varp%ndim, stat)
    if (stat == 1) return

    ! select appropriate routine for each type:
    select case ( filep%ftype )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_ASCII )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_NETCDF )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      
        stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values, start, count, stride, map )
      
        ! just put; let netcdf library convert the right kind:
        !stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values, &
        !                              start, count, stride, map )
        !
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        write(*, '("unsupported filetype : ",i6)') filep%ftype
        stat=1; return
    end select
        
    ! ok
    stat = 0

  end subroutine dios_put_var_c1_1d


  ! ***
  

  subroutine dios_get_var_c1_1d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Get_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    character(len=*), intent(out)       ::  values
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(dios_file), pointer         ::  filep
    type(dios_var), pointer          ::  varp
      

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call dios_file_list_get_pointer( file_list, hid, filep, stat )
    
    ! pointer to variable structure:
    call dios_var_list_get_pointer( filep%var_list, varid, varp, stat )

    ! check ...
    call check_var_values_rank(size(shape(values)), varp%ndim, stat)
    if (stat == 1) return
      
    if ( present(start) ) call check_var_start_size(size(start), varp%ndim, stat)
    if (stat == 1) return

    if ( present(count) ) call check_var_count_size(size(count), varp%ndim, stat)
    if (stat == 1) return

    if ( present(stride) ) call check_var_stride_size(size(stride), varp%ndim, stat)
    if (stat == 1) return

    if ( present(map) ) call check_var_map_size(size(map), varp%ndim, stat)
    if (stat == 1) return
      
    ! select appropriate routine for each type:
    select case ( filep%ftype )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_ASCII )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_NETCDF )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! read values, converted automatically:
        stat = NF90_Get_Var( filep%netcdf_id, varp%netcdf_varid, values, start, count, stride, map )

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        write(*, '("unsupported filetype : ",i6)') filep%ftype
        stat=1; return
    end select
        
    ! ok
    stat = 0

  end subroutine dios_get_var_c1_1d
  
  
  ! ***

  
  subroutine dios_put_var_c1_2d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Put_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    character(len=*), intent(in)        ::  values(:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(dios_file), pointer        ::  filep
    type(dios_var), pointer         ::  varp

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call dios_file_list_get_pointer( file_list, hid, filep, stat )
    
    ! pointer to variable structure:
    call dios_var_list_get_pointer( filep%var_list, varid, varp, stat )

    ! check ...
    call check_var_values_rank(size(shape(values)), varp%ndim, stat)
    if (stat == 1) return
      
    if ( present(start) ) call check_var_start_size(size(start), varp%ndim, stat)
    if (stat == 1) return

    if ( present(count) ) call check_var_count_size(size(count), varp%ndim, stat)
    if (stat == 1) return

    if ( present(stride) ) call check_var_stride_size(size(stride), varp%ndim, stat)
    if (stat == 1) return

    if ( present(map) ) call check_var_map_size(size(map), varp%ndim, stat)
    if (stat == 1) return

    ! select appropriate routine for each type:
    select case ( filep%ftype )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_ASCII )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_NETCDF )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      
        stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values, start, count, stride, map )
      
        ! just put; let netcdf library convert the right kind:
        !stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values, &
        !                              start, count, stride, map )
        !
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        write(*, '("unsupported filetype : ",i6)') filep%ftype
        stat=1; return
    end select
        
    ! ok
    stat = 0

  end subroutine dios_put_var_c1_2d


  ! ***
  

  subroutine dios_get_var_c1_2d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Get_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    character(len=*), intent(out)       ::  values(:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(dios_file), pointer         ::  filep
    type(dios_var), pointer          ::  varp
      

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call dios_file_list_get_pointer( file_list, hid, filep, stat )
    
    ! pointer to variable structure:
    call dios_var_list_get_pointer( filep%var_list, varid, varp, stat )

    ! check ...
    call check_var_values_rank(size(shape(values)), varp%ndim, stat)
    if (stat == 1) return
      
    if ( present(start) ) call check_var_start_size(size(start), varp%ndim, stat)
    if (stat == 1) return

    if ( present(count) ) call check_var_count_size(size(count), varp%ndim, stat)
    if (stat == 1) return

    if ( present(stride) ) call check_var_stride_size(size(stride), varp%ndim, stat)
    if (stat == 1) return

    if ( present(map) ) call check_var_map_size(size(map), varp%ndim, stat)
    if (stat == 1) return
      
    ! select appropriate routine for each type:
    select case ( filep%ftype )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_ASCII )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_NETCDF )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! read values, converted automatically:
        stat = NF90_Get_Var( filep%netcdf_id, varp%netcdf_varid, values, start, count, stride, map )

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        write(*, '("unsupported filetype : ",i6)') filep%ftype
        stat=1; return
    end select
        
    ! ok
    stat = 0

  end subroutine dios_get_var_c1_2d


  ! ***

  
  subroutine dios_put_var_c1_3d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Put_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    character(len=*), intent(in)        ::  values(:, :)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(dios_file), pointer        ::  filep
    type(dios_var), pointer         ::  varp

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call dios_file_list_get_pointer( file_list, hid, filep, stat )
    
    ! pointer to variable structure:
    call dios_var_list_get_pointer( filep%var_list, varid, varp, stat )

    ! check ...
    call check_var_values_rank(size(shape(values)), varp%ndim, stat)
    if (stat == 1) return
      
    if ( present(start) ) call check_var_start_size(size(start), varp%ndim, stat)
    if (stat == 1) return

    if ( present(count) ) call check_var_count_size(size(count), varp%ndim, stat)
    if (stat == 1) return

    if ( present(stride) ) call check_var_stride_size(size(stride), varp%ndim, stat)
    if (stat == 1) return

    if ( present(map) ) call check_var_map_size(size(map), varp%ndim, stat)
    if (stat == 1) return

    ! select appropriate routine for each type:
    select case ( filep%ftype )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_ASCII )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_NETCDF )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      
        stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values, start, count, stride, map )
      
        ! just put; let netcdf library convert the right kind:
        !stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values, &
        !                              start, count, stride, map )
        !
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        write(*, '("unsupported filetype : ",i6)') filep%ftype
        stat=1; return
    end select
        
    ! ok
    stat = 0

  end subroutine dios_put_var_c1_3d


  ! ***
  

  subroutine dios_get_var_c1_3d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Get_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    character(len=*), intent(out)       ::  values(:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(dios_file), pointer         ::  filep
    type(dios_var), pointer          ::  varp
      

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call dios_file_list_get_pointer( file_list, hid, filep, stat )
    
    ! pointer to variable structure:
    call dios_var_list_get_pointer( filep%var_list, varid, varp, stat )

    ! check ...
    call check_var_values_rank(size(shape(values)), varp%ndim, stat)
    if (stat == 1) return
      
    if ( present(start) ) call check_var_start_size(size(start), varp%ndim, stat)
    if (stat == 1) return

    if ( present(count) ) call check_var_count_size(size(count), varp%ndim, stat)
    if (stat == 1) return

    if ( present(stride) ) call check_var_stride_size(size(stride), varp%ndim, stat)
    if (stat == 1) return

    if ( present(map) ) call check_var_map_size(size(map), varp%ndim, stat)
    if (stat == 1) return
      
    ! select appropriate routine for each type:
    select case ( filep%ftype )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_ASCII )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_NETCDF )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! read values, converted automatically:
        stat = NF90_Get_Var( filep%netcdf_id, varp%netcdf_varid, values, start, count, stride, map )

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        write(*, '("unsupported filetype : ",i6)') filep%ftype
        stat=1; return
    end select
        
    ! ok
    stat = 0

  end subroutine dios_get_var_c1_3d
  ! ***

  
  subroutine dios_put_var_c1_4d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Put_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    character(len=*), intent(in)        ::  values(:, :, :)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(dios_file), pointer        ::  filep
    type(dios_var), pointer         ::  varp

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call dios_file_list_get_pointer( file_list, hid, filep, stat )
    
    ! pointer to variable structure:
    call dios_var_list_get_pointer( filep%var_list, varid, varp, stat )

    ! check ...
    call check_var_values_rank(size(shape(values)), varp%ndim, stat)
    if (stat == 1) return
      
    if ( present(start) ) call check_var_start_size(size(start), varp%ndim, stat)
    if (stat == 1) return

    if ( present(count) ) call check_var_count_size(size(count), varp%ndim, stat)
    if (stat == 1) return

    if ( present(stride) ) call check_var_stride_size(size(stride), varp%ndim, stat)
    if (stat == 1) return

    if ( present(map) ) call check_var_map_size(size(map), varp%ndim, stat)
    if (stat == 1) return

    ! select appropriate routine for each type:
    select case ( filep%ftype )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_ASCII )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_NETCDF )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      
        stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values, start, count, stride, map )
      
        ! just put; let netcdf library convert the right kind:
        !stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values, &
        !                              start, count, stride, map )
        !
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        write(*, '("unsupported filetype : ",i6)') filep%ftype
        stat=1; return
    end select
        
    ! ok
    stat = 0

  end subroutine dios_put_var_c1_4d


  ! ***
  

  subroutine dios_get_var_c1_4d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Get_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    character(len=*), intent(out)       ::  values(:,:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(dios_file), pointer         ::  filep
    type(dios_var), pointer          ::  varp
      

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call dios_file_list_get_pointer( file_list, hid, filep, stat )
    
    ! pointer to variable structure:
    call dios_var_list_get_pointer( filep%var_list, varid, varp, stat )

    ! check ...
    call check_var_values_rank(size(shape(values)), varp%ndim, stat)
    if (stat == 1) return
      
    if ( present(start) ) call check_var_start_size(size(start), varp%ndim, stat)
    if (stat == 1) return

    if ( present(count) ) call check_var_count_size(size(count), varp%ndim, stat)
    if (stat == 1) return

    if ( present(stride) ) call check_var_stride_size(size(stride), varp%ndim, stat)
    if (stat == 1) return

    if ( present(map) ) call check_var_map_size(size(map), varp%ndim, stat)
    if (stat == 1) return
      
    ! select appropriate routine for each type:
    select case ( filep%ftype )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_ASCII )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_NETCDF )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! read values, converted automatically:
        stat = NF90_Get_Var( filep%netcdf_id, varp%netcdf_varid, values, start, count, stride, map )

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        write(*, '("unsupported filetype : ",i6)') filep%ftype
        stat=1; return
    end select
        
    ! ok
    stat = 0

  end subroutine dios_get_var_c1_4d


  ! ***

  
  subroutine dios_put_var_c1_5d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Put_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    character(len=*), intent(in)        ::  values(:, :, :, :)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(dios_file), pointer        ::  filep
    type(dios_var), pointer         ::  varp

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call dios_file_list_get_pointer( file_list, hid, filep, stat )
    
    ! pointer to variable structure:
    call dios_var_list_get_pointer( filep%var_list, varid, varp, stat )

    ! check ...
    call check_var_values_rank(size(shape(values)), varp%ndim, stat)
    if (stat == 1) return
      
    if ( present(start) ) call check_var_start_size(size(start), varp%ndim, stat)
    if (stat == 1) return

    if ( present(count) ) call check_var_count_size(size(count), varp%ndim, stat)
    if (stat == 1) return

    if ( present(stride) ) call check_var_stride_size(size(stride), varp%ndim, stat)
    if (stat == 1) return

    if ( present(map) ) call check_var_map_size(size(map), varp%ndim, stat)
    if (stat == 1) return

    ! select appropriate routine for each type:
    select case ( filep%ftype )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_ASCII )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_NETCDF )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      
        stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values, start, count, stride, map )
      
        ! just put; let netcdf library convert the right kind:
        !stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values, &
        !                              start, count, stride, map )
        !
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        write(*, '("unsupported filetype : ",i6)') filep%ftype
        stat=1; return
    end select
        
    ! ok
    stat = 0

  end subroutine dios_put_var_c1_5d


  ! ***
  

  subroutine dios_get_var_c1_5d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Get_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    character(len=*), intent(out)       ::  values(:,:,:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(dios_file), pointer         ::  filep
    type(dios_var), pointer          ::  varp
      

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call dios_file_list_get_pointer( file_list, hid, filep, stat )
    
    ! pointer to variable structure:
    call dios_var_list_get_pointer( filep%var_list, varid, varp, stat )

    ! check ...
    call check_var_values_rank(size(shape(values)), varp%ndim, stat)
    if (stat == 1) return
      
    if ( present(start) ) call check_var_start_size(size(start), varp%ndim, stat)
    if (stat == 1) return

    if ( present(count) ) call check_var_count_size(size(count), varp%ndim, stat)
    if (stat == 1) return

    if ( present(stride) ) call check_var_stride_size(size(stride), varp%ndim, stat)
    if (stat == 1) return

    if ( present(map) ) call check_var_map_size(size(map), varp%ndim, stat)
    if (stat == 1) return
      
    ! select appropriate routine for each type:
    select case ( filep%ftype )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_ASCII )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_NETCDF )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! read values, converted automatically:
        stat = NF90_Get_Var( filep%netcdf_id, varp%netcdf_varid, values, start, count, stride, map )

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        write(*, '("unsupported filetype : ",i6)') filep%ftype
        stat=1; return
    end select
        
    ! ok
    stat = 0

  end subroutine dios_get_var_c1_5d


  ! ***

  
  subroutine dios_put_var_c1_6d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Put_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    character(len=*), intent(in)        ::  values(:, :, :, :, :)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(dios_file), pointer        ::  filep
    type(dios_var), pointer         ::  varp

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call dios_file_list_get_pointer( file_list, hid, filep, stat )
    
    ! pointer to variable structure:
    call dios_var_list_get_pointer( filep%var_list, varid, varp, stat )

    ! check ...
    call check_var_values_rank(size(shape(values)), varp%ndim, stat)
    if (stat == 1) return
      
    if ( present(start) ) call check_var_start_size(size(start), varp%ndim, stat)
    if (stat == 1) return

    if ( present(count) ) call check_var_count_size(size(count), varp%ndim, stat)
    if (stat == 1) return

    if ( present(stride) ) call check_var_stride_size(size(stride), varp%ndim, stat)
    if (stat == 1) return

    if ( present(map) ) call check_var_map_size(size(map), varp%ndim, stat)
    if (stat == 1) return

    ! select appropriate routine for each type:
    select case ( filep%ftype )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_ASCII )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_NETCDF )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      
        stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values, start, count, stride, map )
      
        ! just put; let netcdf library convert the right kind:
        !stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values, &
        !                              start, count, stride, map )
        !
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        write(*, '("unsupported filetype : ",i6)') filep%ftype
        stat=1; return
    end select
        
    ! ok
    stat = 0

  end subroutine dios_put_var_c1_6d


  ! ***
  

  subroutine dios_get_var_c1_6d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Get_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    character(len=*), intent(out)       ::  values(:,:,:,:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(dios_file), pointer         ::  filep
    type(dios_var), pointer          ::  varp
      

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call dios_file_list_get_pointer( file_list, hid, filep, stat )
    
    ! pointer to variable structure:
    call dios_var_list_get_pointer( filep%var_list, varid, varp, stat )

    ! check ...
    call check_var_values_rank(size(shape(values)), varp%ndim, stat)
    if (stat == 1) return
      
    if ( present(start) ) call check_var_start_size(size(start), varp%ndim, stat)
    if (stat == 1) return

    if ( present(count) ) call check_var_count_size(size(count), varp%ndim, stat)
    if (stat == 1) return

    if ( present(stride) ) call check_var_stride_size(size(stride), varp%ndim, stat)
    if (stat == 1) return

    if ( present(map) ) call check_var_map_size(size(map), varp%ndim, stat)
    if (stat == 1) return
      
    ! select appropriate routine for each type:
    select case ( filep%ftype )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_ASCII )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_NETCDF )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! read values, converted automatically:
        stat = NF90_Get_Var( filep%netcdf_id, varp%netcdf_varid, values, start, count, stride, map )

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        write(*, '("unsupported filetype : ",i6)') filep%ftype
        stat=1; return
    end select
        
    ! ok
    stat = 0

  end subroutine dios_get_var_c1_6d


  ! ***

  
  subroutine dios_put_var_c1_7d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Put_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    character(len=*), intent(in)        ::  values(:,:,:,:,:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(dios_file), pointer        ::  filep
    type(dios_var), pointer         ::  varp

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call dios_file_list_get_pointer( file_list, hid, filep, stat )
    
    ! pointer to variable structure:
    call dios_var_list_get_pointer( filep%var_list, varid, varp, stat )

    ! check ...
    call check_var_values_rank(size(shape(values)), varp%ndim, stat)
    if (stat == 1) return
      
    if ( present(start) ) call check_var_start_size(size(start), varp%ndim, stat)
    if (stat == 1) return

    if ( present(count) ) call check_var_count_size(size(count), varp%ndim, stat)
    if (stat == 1) return

    if ( present(stride) ) call check_var_stride_size(size(stride), varp%ndim, stat)
    if (stat == 1) return

    if ( present(map) ) call check_var_map_size(size(map), varp%ndim, stat)
    if (stat == 1) return

    ! select appropriate routine for each type:
    select case ( filep%ftype )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_ASCII )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_NETCDF )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      
        stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values, start, count, stride, map )
      
        ! just put; let netcdf library convert the right kind:
        !stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values, &
        !                              start, count, stride, map )
        !
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        write(*, '("unsupported filetype : ",i6)') filep%ftype
        stat=1; return
    end select
        
    ! ok
    stat = 0

  end subroutine dios_put_var_c1_7d


  ! ***
  

  subroutine dios_get_var_c1_7d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Get_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    character(len=*), intent(out)       ::  values(:,:,:,:,:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(dios_file), pointer         ::  filep
    type(dios_var), pointer          ::  varp
      

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call dios_file_list_get_pointer( file_list, hid, filep, stat )
    
    ! pointer to variable structure:
    call dios_var_list_get_pointer( filep%var_list, varid, varp, stat )

    ! check ...
    call check_var_values_rank(size(shape(values)), varp%ndim, stat)
    if (stat == 1) return
      
    if ( present(start) ) call check_var_start_size(size(start), varp%ndim, stat)
    if (stat == 1) return

    if ( present(count) ) call check_var_count_size(size(count), varp%ndim, stat)
    if (stat == 1) return

    if ( present(stride) ) call check_var_stride_size(size(stride), varp%ndim, stat)
    if (stat == 1) return

    if ( present(map) ) call check_var_map_size(size(map), varp%ndim, stat)
    if (stat == 1) return
      
    ! select appropriate routine for each type:
    select case ( filep%ftype )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_ASCII )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_NETCDF )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! read values, converted automatically:
        stat = NF90_Get_Var( filep%netcdf_id, varp%netcdf_varid, values, start, count, stride, map )

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        write(*, '("unsupported filetype : ",i6)') filep%ftype
        stat=1; return
    end select
        
    ! ok
    stat = 0

  end subroutine dios_get_var_c1_7d


  ! ********************************************************************
  ! ***
  ! *** var_i1
  ! ***
  ! ********************************************************************

  
  subroutine dios_put_var_i1_1d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Put_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    integer(1), intent(in)              ::  values(:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(dios_file), pointer         ::  filep
    type(dios_var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:)
    integer(2), allocatable         ::  values_int2(:)
    integer(4), allocatable         ::  values_int4(:)
    integer(8), allocatable         ::  values_int8(:)
    real(4), allocatable            ::  values_real4(:)
    real(8), allocatable            ::  values_real8(:)
        

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call dios_file_list_get_pointer( file_list, hid, filep, stat )
    
    ! pointer to variable structure:
    call dios_var_list_get_pointer( filep%var_list, varid, varp, stat )

    ! select appropriate routine for each type:
    select case ( filep%ftype )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_ASCII )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_NETCDF )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      
        ! test target type:
        ! convert to required kind wrt var type before entering NF90_Put_Var, 
        ! otherwise segmentation faults on some machines ...
        select case ( varp%xtype )
        
          case ( DIOS_BYTE )
            allocate( values_int1(size(values,1)) )
            values_int1 = int(values,kind=1)
            stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int1, start, count, stride, map )
            deallocate( values_int1 )
        
          case ( DIOS_SHORT )
            allocate( values_int2(size(values,1)) )
            values_int2 = int(values,kind=2)
            stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int2, start, count, stride, map )
            deallocate( values_int2 )
        
          case ( DIOS_INT )
            allocate( values_int4(size(values,1)) )
            values_int4 = int(values,kind=4)
            stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int4, start, count, stride, map )
            deallocate( values_int4 )
        
          case ( DIOS_FLOAT )
            allocate( values_real4(size(values,1)) )
            values_real4 = real(values,kind=4)
            stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real4, start, count, stride, map )
            deallocate( values_real4 )
        
          case ( DIOS_DOUBLE )
            allocate( values_real8(size(values,1)) )
            values_real8 = real(values,kind=8)
            stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real8, start, count, stride, map )
            deallocate( values_real8 )
        
          case default
            write(*, '("not implemented yet for output type : ",i6)') varp%xtype
            stat=1; return
        end select
        
        ! just put; let netcdf library convert the right kind:
        !stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values, start, count, stride, map )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        write(*, '("unsupported filetype : ",i6)') filep%ftype
        stat=1; return
    end select
        
    ! ok
    stat = 0

  end subroutine dios_put_var_i1_1d


  ! ***
  

  subroutine dios_get_var_i1_1d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Get_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    integer(1), intent(out)             ::  values(:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(dios_file), pointer         ::  filep
    type(dios_var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:)
    integer(2), allocatable         ::  values_int2(:)
    integer(4), allocatable         ::  values_int4(:)
    integer(8), allocatable         ::  values_int8(:)
    real(4), allocatable            ::  values_real4(:)
    real(8), allocatable            ::  values_real8(:)

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call dios_file_list_get_pointer( file_list, hid, filep, stat )
    
    ! pointer to variable structure:
    call dios_var_list_get_pointer( filep%var_list, varid, varp, stat )
      
    ! select appropriate routine for each type:
    select case ( filep%ftype )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_ASCII )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_NETCDF )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! read values, converted automatically:
        stat = NF90_Get_Var( filep%netcdf_id, varp%netcdf_varid, values, start, count, stride, map )

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        write(*, '("unsupported filetype : ",i6)') filep%ftype
        stat=1; return
    end select
        
    ! ok
    stat = 0

  end subroutine dios_get_var_i1_1d
  
  
  ! ***

  
  subroutine DIOS_Put_Var_i1_2d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Put_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    integer(1), intent(in)            ::  values(:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:,:)
    integer(2), allocatable         ::  values_int2(:,:)
    integer(4), allocatable         ::  values_int4(:,:)
    integer(8), allocatable         ::  values_int8(:,:)
    real(4), allocatable            ::  values_real4(:,:)
    real(8), allocatable            ::  values_real8(:,:)
        

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call dios_file_list_get_pointer( file_list, hid, filep, stat )
    
    ! pointer to variable structure:
    call dios_var_list_get_pointer( filep%var_list, varid, varp, stat )
      
      ! select appropriate routine for each type:
      select case ( filep%ftype )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_ASCII )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_NETCDF )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      
          ! test target type:
          ! convert to required kind before entering NF90_Put_Var, 
          ! otherwise segmentation faults on some machines ...
          select case ( varp%xtype )
        
            case ( DIOS_BYTE )
              allocate( values_int1(size(values,1),size(values,2)) )
              values_int1 = int(values,kind=1)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int1, &
                                            start, count, stride, map )
              deallocate( values_int1 )
          
            case ( DIOS_SHORT )
              allocate( values_int2(size(values,1),size(values,2)) )
              values_int2 = int(values,kind=2)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int2, &
                                            start, count, stride, map )
              deallocate( values_int2 )
          
            case ( DIOS_INT )
              allocate( values_int4(size(values,1),size(values,2)) )
              values_int4 = int(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int4, &
                                            start, count, stride, map )
              
              deallocate( values_int4 )
          
            case ( DIOS_FLOAT )
              allocate( values_real4(size(values,1),size(values,2)) )
              values_real4 = real(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real4, &
                                            start, count, stride, map )
              
              deallocate( values_real4 )
          
            case ( DIOS_DOUBLE )
              allocate( values_real8(size(values,1),size(values,2)) )
              values_real8 = real(values,kind=8)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real8, &
                                            start, count, stride, map )
              
              deallocate( values_real8 )
          
            case default
              write(*, '("not implemented yet for output type : ",i6)') varp%xtype
              stat=1; return
          end select
        
          ! just put; let netcdf library convert the right kind:
          !stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values, &
          !                              start, count, stride, map )
          !

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*, '("unsupported filetype : ",i6)') filep%ftype
          stat=1; return
      end select
        
    ! ok
    stat = 0

  end subroutine DIOS_Put_Var_i1_2d


  ! ***


  subroutine dios_get_var_i1_2d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Get_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    integer(1), intent(out)             ::  values(:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(dios_file), pointer         ::  filep
    type(dios_var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:,:)
    integer(2), allocatable         ::  values_int2(:,:)
    integer(4), allocatable         ::  values_int4(:,:)
    integer(8), allocatable         ::  values_int8(:,:)
    real(4), allocatable            ::  values_real4(:,:)
    real(8), allocatable            ::  values_real8(:,:)

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call dios_file_list_get_pointer( file_list, hid, filep, stat )
    
    ! pointer to variable structure:
    call dios_var_list_get_pointer( filep%var_list, varid, varp, stat )
      
    ! select appropriate routine for each type:
    select case ( filep%ftype )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_ASCII )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_NETCDF )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! read values, converted automatically:
        stat = NF90_Get_Var( filep%netcdf_id, varp%netcdf_varid, values, start, count, stride, map )

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        write(*, '("unsupported filetype : ",i6)') filep%ftype
        stat=1; return
    end select
        
    ! ok
    stat = 0

  end subroutine dios_get_var_i1_2d
  

  ! ***

  
  subroutine DIOS_Put_Var_i1_3d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Put_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    integer(1), intent(in)            ::  values(:,:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:,:,:)
    integer(2), allocatable         ::  values_int2(:,:,:)
    integer(4), allocatable         ::  values_int4(:,:,:)
    integer(8), allocatable         ::  values_int8(:,:,:)
    real(4), allocatable            ::  values_real4(:,:,:)
    real(8), allocatable            ::  values_real8(:,:,:)
        

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    ! pointer to variable structure:
    call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )

      ! select appropriate routine for each type:
      select case ( filep%ftype )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_ASCII )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_NETCDF )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      
          ! test target type:
          ! convert to required kind before entering NF90_Put_Var, 
          ! otherwise segmentation faults on some machines ...
          select case ( varp%xtype )
        
            case ( DIOS_BYTE )
              allocate( values_int1(size(values,1),size(values,2),size(values,3)) )
              values_int1 = int(values,kind=1)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int1, &
                                            start, count, stride, map )
              
              deallocate( values_int1 )
          
            case ( DIOS_SHORT )
              allocate( values_int2(size(values,1),size(values,2),size(values,3)) )
              values_int2 = int(values,kind=2)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int2, &
                                            start, count, stride, map )
              
              deallocate( values_int2 )
          
            case ( DIOS_INT )
              allocate( values_int4(size(values,1),size(values,2),size(values,3)) )
              values_int4 = int(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int4, &
                                            start, count, stride, map )
              
              deallocate( values_int4 )
          
            case ( DIOS_FLOAT )
              allocate( values_real4(size(values,1),size(values,2),size(values,3)) )
              values_real4 = real(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real4, &
                                            start, count, stride, map )
              
              deallocate( values_real4 )
          
            case ( DIOS_DOUBLE )
              allocate( values_real8(size(values,1),size(values,2),size(values,3)) )
              values_real8 = real(values,kind=8)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real8, &
                                            start, count, stride, map )
              
              deallocate( values_real8 )
          
            case default
              write(*, '("not implemented yet for output type : ",i6)') varp%xtype
              stat=1; return
          end select
        
          ! just put; let netcdf library convert the right kind:
          !stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values, &
          !                              start, count, stride, map )
          !
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*, '("unsupported filetype : ",i6)') filep%ftype
          stat=1; return
      end select
        
    ! ok
    stat = 0

  end subroutine DIOS_Put_Var_i1_3d


  ! ***


  subroutine dios_get_var_i1_3d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Get_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    integer(1), intent(out)             ::  values(:,:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(dios_file), pointer         ::  filep
    type(dios_var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:,:,:)
    integer(2), allocatable         ::  values_int2(:,:,:)
    integer(4), allocatable         ::  values_int4(:,:,:)
    integer(8), allocatable         ::  values_int8(:,:,:)
    real(4), allocatable            ::  values_real4(:,:,:)
    real(8), allocatable            ::  values_real8(:,:,:)

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call dios_file_list_get_pointer( file_list, hid, filep, stat )
    
    ! pointer to variable structure:
    call dios_var_list_get_pointer( filep%var_list, varid, varp, stat )
      
    ! select appropriate routine for each type:
    select case ( filep%ftype )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_ASCII )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_NETCDF )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! read values, converted automatically:
        stat = NF90_Get_Var( filep%netcdf_id, varp%netcdf_varid, values, start, count, stride, map )

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        write(*, '("unsupported filetype : ",i6)') filep%ftype
        stat=1; return
    end select
        
    ! ok
    stat = 0

  end subroutine dios_get_var_i1_3d
  

  ! ***

  
  subroutine DIOS_Put_Var_i1_4d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Put_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    integer(1), intent(in)            ::  values(:,:,:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:,:,:,:)
    integer(2), allocatable         ::  values_int2(:,:,:,:)
    integer(4), allocatable         ::  values_int4(:,:,:,:)
    integer(8), allocatable         ::  values_int8(:,:,:,:)
    real(4), allocatable            ::  values_real4(:,:,:,:)
    real(8), allocatable            ::  values_real8(:,:,:,:)
        

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    
    ! pointer to variable structure:
    call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
    

      
      ! select appropriate routine for each type:
      select case ( filep%ftype )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_ASCII )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_NETCDF )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      
          ! test target type:
          ! convert to required kind before entering NF90_Put_Var, 
          ! otherwise segmentation faults on some machines ...
          select case ( varp%xtype )
        
            case ( DIOS_BYTE )
              allocate( values_int1(size(values,1),size(values,2),size(values,3),size(values,4)) )
              values_int1 = int(values,kind=1)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int1, &
                                            start, count, stride, map )
              
              deallocate( values_int1 )
          
            case ( DIOS_SHORT )
              allocate( values_int2(size(values,1),size(values,2),size(values,3),size(values,4)) )
              values_int2 = int(values,kind=2)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int2, &
                                            start, count, stride, map )
              
              deallocate( values_int2 )
          
            case ( DIOS_INT )
              allocate( values_int4(size(values,1),size(values,2),size(values,3),size(values,4)) )
              values_int4 = int(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int4, &
                                            start, count, stride, map )
              
              deallocate( values_int4 )
          
            case ( DIOS_FLOAT )
              allocate( values_real4(size(values,1),size(values,2),size(values,3),size(values,4)) )
              values_real4 = real(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real4, &
                                            start, count, stride, map )
              
              deallocate( values_real4 )
          
            case ( DIOS_DOUBLE )
              allocate( values_real8(size(values,1),size(values,2),size(values,3),size(values,4)) )
              values_real8 = real(values,kind=8)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real8, &
                                            start, count, stride, map )
              
              deallocate( values_real8 )
          
            case default
              write(*, '("not implemented yet for output type : ",i6)') varp%xtype
              stat=1; return
          end select
        
          ! just put; let netcdf library convert the right kind:
          !stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values, &
          !                              start, count, stride, map )
          !
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*, '("unsupported filetype : ",i6)') filep%ftype
          stat=1; return
      end select
        
    ! ok
    stat = 0

  end subroutine DIOS_Put_Var_i1_4d


  ! ***


  subroutine dios_get_var_i1_4d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Get_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    integer(1), intent(out)             ::  values(:,:,:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(dios_file), pointer         ::  filep
    type(dios_var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:,:,:,:)
    integer(2), allocatable         ::  values_int2(:,:,:,:)
    integer(4), allocatable         ::  values_int4(:,:,:,:)
    integer(8), allocatable         ::  values_int8(:,:,:,:)
    real(4), allocatable            ::  values_real4(:,:,:,:)
    real(8), allocatable            ::  values_real8(:,:,:,:)

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call dios_file_list_get_pointer( file_list, hid, filep, stat )
    
    ! pointer to variable structure:
    call dios_var_list_get_pointer( filep%var_list, varid, varp, stat )
      
    ! select appropriate routine for each type:
    select case ( filep%ftype )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_ASCII )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_NETCDF )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! read values, converted automatically:
        stat = NF90_Get_Var( filep%netcdf_id, varp%netcdf_varid, values, start, count, stride, map )

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        write(*, '("unsupported filetype : ",i6)') filep%ftype
        stat=1; return
    end select
        
    ! ok
    stat = 0

  end subroutine dios_get_var_i1_4d


  ! ***
  
  
  subroutine DIOS_Put_Var_i1_5d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Put_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    integer(1), intent(in)            ::  values(:,:,:,:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:,:,:,:,:)
    integer(2), allocatable         ::  values_int2(:,:,:,:,:)
    integer(4), allocatable         ::  values_int4(:,:,:,:,:)
    integer(8), allocatable         ::  values_int8(:,:,:,:,:)
    real(4), allocatable            ::  values_real4(:,:,:,:,:)
    real(8), allocatable            ::  values_real8(:,:,:,:,:)
        

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    
    ! pointer to variable structure:
    call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
    

      
      select case ( filep%ftype )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_ASCII )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_NETCDF )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      
          ! test target type:
          ! convert to required kind before entering NF90_Put_Var, 
          ! otherwise segmentation faults on some machines ...
          select case ( varp%xtype )
        
            case ( DIOS_BYTE )
              allocate( values_int1(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5)) )
              values_int1 = int(values,kind=1)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int1, &
                                            start, count, stride, map )
              
              deallocate( values_int1 )
          
            case ( DIOS_SHORT )
              allocate( values_int2(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5)) )
              values_int2 = int(values,kind=2)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int2, &
                                            start, count, stride, map )
              
              deallocate( values_int2 )
          
            case ( DIOS_INT )
              allocate( values_int4(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5)) )
              values_int4 = int(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int4, &
                                            start, count, stride, map )
              
              deallocate( values_int4 )
          
            case ( DIOS_FLOAT )
              allocate( values_real4(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5)) )
              values_real4 = real(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real4, &
                                            start, count, stride, map )
              
              deallocate( values_real4 )
          
            case ( DIOS_DOUBLE )
              allocate( values_real8(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5)) )
              values_real8 = real(values,kind=8)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real8, &
                                            start, count, stride, map )
              
              deallocate( values_real8 )
          
            case default
              write(*, '("not implemented yet for output type : ",i6)') varp%xtype
              stat=1; return
          end select
        
          ! just put; let netcdf library convert the right kind:
          !stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values, &
          !                              start, count, stride, map )
          !

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*, '("unsupported filetype : ",i6)') filep%ftype
          stat=1; return
      end select
        
    ! ok
    stat = 0

  end subroutine DIOS_Put_Var_i1_5d


  ! ***


  subroutine dios_get_var_i1_5d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Get_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    integer(1), intent(out)             ::  values(:,:,:,:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(dios_file), pointer         ::  filep
    type(dios_var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:,:,:,:,:)
    integer(2), allocatable         ::  values_int2(:,:,:,:,:)
    integer(4), allocatable         ::  values_int4(:,:,:,:,:)
    integer(8), allocatable         ::  values_int8(:,:,:,:,:)
    real(4), allocatable            ::  values_real4(:,:,:,:,:)
    real(8), allocatable            ::  values_real8(:,:,:,:,:)

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call dios_file_list_get_pointer( file_list, hid, filep, stat )
    
    ! pointer to variable structure:
    call dios_var_list_get_pointer( filep%var_list, varid, varp, stat )
      
    ! select appropriate routine for each type:
    select case ( filep%ftype )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_ASCII )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_NETCDF )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! read values, converted automatically:
        stat = NF90_Get_Var( filep%netcdf_id, varp%netcdf_varid, values, start, count, stride, map )

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        write(*, '("unsupported filetype : ",i6)') filep%ftype
        stat=1; return
    end select
        
    ! ok
    stat = 0

  end subroutine dios_get_var_i1_5d


  ! ***
  
  
  subroutine DIOS_Put_Var_i1_6d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Put_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    integer(1), intent(in)            ::  values(:,:,:,:,:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:,:,:,:,:,:)
    integer(2), allocatable         ::  values_int2(:,:,:,:,:,:)
    integer(4), allocatable         ::  values_int4(:,:,:,:,:,:)
    integer(8), allocatable         ::  values_int8(:,:,:,:,:,:)
    real(4), allocatable            ::  values_real4(:,:,:,:,:,:)
    real(8), allocatable            ::  values_real8(:,:,:,:,:,:)
        

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    
    ! pointer to variable structure:
    call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
    

      select case ( filep%ftype )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_ASCII )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_NETCDF )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      
          ! test target type:
          ! convert to required kind before entering NF90_Put_Var, 
          ! otherwise segmentation faults on some machines ...
          select case ( varp%xtype )
        
            case ( DIOS_BYTE )
              allocate( values_int1(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5),size(values,6)) )
              values_int1 = int(values,kind=1)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int1, &
                                            start, count, stride, map )
              
              deallocate( values_int1 )
          
            case ( DIOS_SHORT )
              allocate( values_int2(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5),size(values,6)) )
              values_int2 = int(values,kind=2)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int2, &
                                            start, count, stride, map )
              
              deallocate( values_int2 )
          
            case ( DIOS_INT )
              allocate( values_int4(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5),size(values,6)) )
              values_int4 = int(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int4, &
                                            start, count, stride, map )
              
              deallocate( values_int4 )
          
            case ( DIOS_FLOAT )
              allocate( values_real4(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5),size(values,6)) )
              values_real4 = real(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real4, &
                                            start, count, stride, map )
              
              deallocate( values_real4 )
          
            case ( DIOS_DOUBLE )
              allocate( values_real8(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5),size(values,6)) )
              values_real8 = real(values,kind=8)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real8, &
                                            start, count, stride, map )
              
              deallocate( values_real8 )
          
            case default
              write(*, '("not implemented yet for output type : ",i6)') varp%xtype
              stat=1; return
          end select
        
          ! just put; let netcdf library convert the right kind:
          !stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values, &
          !                              start, count, stride, map )
          !
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*, '("unsupported filetype : ",i6)') filep%ftype
          stat=1; return
      end select
        
    ! ok
    stat = 0

  end subroutine DIOS_Put_Var_i1_6d


  ! ***


  subroutine dios_get_var_i1_6d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Get_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    integer(1), intent(out)             ::  values(:,:,:,:,:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(dios_file), pointer         ::  filep
    type(dios_var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:,:,:,:,:,:)
    integer(2), allocatable         ::  values_int2(:,:,:,:,:,:)
    integer(4), allocatable         ::  values_int4(:,:,:,:,:,:)
    integer(8), allocatable         ::  values_int8(:,:,:,:,:,:)
    real(4), allocatable            ::  values_real4(:,:,:,:,:,:)
    real(8), allocatable            ::  values_real8(:,:,:,:,:,:)

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call dios_file_list_get_pointer( file_list, hid, filep, stat )
    
    ! pointer to variable structure:
    call dios_var_list_get_pointer( filep%var_list, varid, varp, stat )
      
    ! select appropriate routine for each type:
    select case ( filep%ftype )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_ASCII )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_NETCDF )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! read values, converted automatically:
        stat = NF90_Get_Var( filep%netcdf_id, varp%netcdf_varid, values, start, count, stride, map )

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        write(*, '("unsupported filetype : ",i6)') filep%ftype
        stat=1; return
    end select
        
    ! ok
    stat = 0

  end subroutine dios_get_var_i1_6d


  ! ***
  
  
  subroutine DIOS_Put_Var_i1_7d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Put_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    integer(1), intent(in)            ::  values(:,:,:,:,:,:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:,:,:,:,:,:,:)
    integer(2), allocatable         ::  values_int2(:,:,:,:,:,:,:)
    integer(4), allocatable         ::  values_int4(:,:,:,:,:,:,:)
    integer(8), allocatable         ::  values_int8(:,:,:,:,:,:,:)
    real(4), allocatable            ::  values_real4(:,:,:,:,:,:,:)
    real(8), allocatable            ::  values_real8(:,:,:,:,:,:,:)
        

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    
    ! pointer to variable structure:
    call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
    

      
      ! select appropriate routine for each type:
      select case ( filep%ftype )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_ASCII )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_NETCDF )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      
          ! test target type:
          ! convert to required kind before entering NF90_Put_Var, 
          ! otherwise segmentation faults on some machines ...
          select case ( varp%xtype )
        
            case ( DIOS_BYTE )
              allocate( values_int1(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5),size(values,6),size(values,7)) )
              values_int1 = int(values,kind=1)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int1, &
                                            start, count, stride, map )
              
              deallocate( values_int1 )
          
            case ( DIOS_SHORT )
              allocate( values_int2(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5),size(values,6),size(values,7)) )
              values_int2 = int(values,kind=2)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int2, &
                                            start, count, stride, map )
              
              deallocate( values_int2 )
          
            case ( DIOS_INT )
              allocate( values_int4(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5),size(values,6),size(values,7)) )
              values_int4 = int(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int4, &
                                            start, count, stride, map )
              
              deallocate( values_int4 )
          
            case ( DIOS_FLOAT )
              allocate( values_real4(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5),size(values,6),size(values,7)) )
              values_real4 = real(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real4, &
                                            start, count, stride, map )
              
              deallocate( values_real4 )
          
            case ( DIOS_DOUBLE )
              allocate( values_real8(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5),size(values,6),size(values,7)) )
              values_real8 = real(values,kind=8)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real8, &
                                            start, count, stride, map )
              
              deallocate( values_real8 )
          
            case default
              write(*, '("not implemented yet for output type : ",i6)') varp%xtype
              stat=1; return
          end select
        
          ! just put; let netcdf library convert the right kind:
          !stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values, &
          !                              start, count, stride, map )
          !
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*, '("unsupported filetype : ",i6)') filep%ftype
          stat=1; return
      end select
        
    ! ok
    stat = 0

  end subroutine DIOS_Put_Var_i1_7d


  ! ***


  subroutine dios_get_var_i1_7d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Get_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    integer(1), intent(out)             ::  values(:,:,:,:,:,:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(dios_file), pointer         ::  filep
    type(dios_var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:,:,:,:,:,:,:)
    integer(2), allocatable         ::  values_int2(:,:,:,:,:,:,:)
    integer(4), allocatable         ::  values_int4(:,:,:,:,:,:,:)
    integer(8), allocatable         ::  values_int8(:,:,:,:,:,:,:)
    real(4), allocatable            ::  values_real4(:,:,:,:,:,:,:)
    real(8), allocatable            ::  values_real8(:,:,:,:,:,:,:)

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call dios_file_list_get_pointer( file_list, hid, filep, stat )
    
    ! pointer to variable structure:
    call dios_var_list_get_pointer( filep%var_list, varid, varp, stat )
      
    ! select appropriate routine for each type:
    select case ( filep%ftype )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_ASCII )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_NETCDF )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! read values, converted automatically:
        stat = NF90_Get_Var( filep%netcdf_id, varp%netcdf_varid, values, start, count, stride, map )

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        write(*, '("unsupported filetype : ",i6)') filep%ftype
        stat=1; return
    end select
        
    ! ok
    stat = 0

  end subroutine dios_get_var_i1_7d


  ! ********************************************************************
  ! ***
  ! *** var_i2
  ! ***
  ! ********************************************************************

  
  subroutine dios_put_var_i2_1d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Put_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    integer(2), intent(in)              ::  values(:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(dios_file), pointer         ::  filep
    type(dios_var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:)
    integer(2), allocatable         ::  values_int2(:)
    integer(4), allocatable         ::  values_int4(:)
    integer(8), allocatable         ::  values_int8(:)
    real(4), allocatable            ::  values_real4(:)
    real(8), allocatable            ::  values_real8(:)
        

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call dios_file_list_get_pointer( file_list, hid, filep, stat )
    
    ! pointer to variable structure:
    call dios_var_list_get_pointer( filep%var_list, varid, varp, stat )

    ! select appropriate routine for each type:
    select case ( filep%ftype )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_ASCII )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_NETCDF )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      
        ! test target type:
        ! convert to required kind wrt var type before entering NF90_Put_Var, 
        ! otherwise segmentation faults on some machines ...
        select case ( varp%xtype )
        
          case ( DIOS_BYTE )
            allocate( values_int1(size(values,1)) )
            values_int1 = int(values,kind=1)
            stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int1, start, count, stride, map )
            deallocate( values_int1 )
        
          case ( DIOS_SHORT )
            allocate( values_int2(size(values,1)) )
            values_int2 = int(values,kind=2)
            stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int2, start, count, stride, map )
            deallocate( values_int2 )
        
          case ( DIOS_INT )
            allocate( values_int4(size(values,1)) )
            values_int4 = int(values,kind=4)
            stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int4, start, count, stride, map )
            deallocate( values_int4 )
        
          case ( DIOS_FLOAT )
            allocate( values_real4(size(values,1)) )
            values_real4 = real(values,kind=4)
            stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real4, start, count, stride, map )
            deallocate( values_real4 )
        
          case ( DIOS_DOUBLE )
            allocate( values_real8(size(values,1)) )
            values_real8 = real(values,kind=8)
            stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real8, start, count, stride, map )
            deallocate( values_real8 )
        
          case default
            write(*, '("not implemented yet for output type : ",i6)') varp%xtype
            stat=1; return
        end select
        
        ! just put; let netcdf library convert the right kind:
        !stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values, start, count, stride, map )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        write(*, '("unsupported filetype : ",i6)') filep%ftype
        stat=1; return
    end select
        
    ! ok
    stat = 0

  end subroutine dios_put_var_i2_1d


  ! ***
  

  subroutine dios_get_var_i2_1d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Get_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    integer(2), intent(out)             ::  values(:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(dios_file), pointer         ::  filep
    type(dios_var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:)
    integer(2), allocatable         ::  values_int2(:)
    integer(4), allocatable         ::  values_int4(:)
    integer(8), allocatable         ::  values_int8(:)
    real(4), allocatable            ::  values_real4(:)
    real(8), allocatable            ::  values_real8(:)

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call dios_file_list_get_pointer( file_list, hid, filep, stat )
    
    ! pointer to variable structure:
    call dios_var_list_get_pointer( filep%var_list, varid, varp, stat )
      
    ! select appropriate routine for each type:
    select case ( filep%ftype )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_ASCII )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_NETCDF )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! read values, converted automatically:
        stat = NF90_Get_Var( filep%netcdf_id, varp%netcdf_varid, values, start, count, stride, map )

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        write(*, '("unsupported filetype : ",i6)') filep%ftype
        stat=1; return
    end select
        
    ! ok
    stat = 0

  end subroutine dios_get_var_i2_1d
  
  
  ! ***

  
  subroutine DIOS_Put_var_i2_2d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Put_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    integer(2), intent(in)            ::  values(:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:,:)
    integer(2), allocatable         ::  values_int2(:,:)
    integer(4), allocatable         ::  values_int4(:,:)
    integer(8), allocatable         ::  values_int8(:,:)
    real(4), allocatable            ::  values_real4(:,:)
    real(8), allocatable            ::  values_real8(:,:)
        

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call dios_file_list_get_pointer( file_list, hid, filep, stat )
    
    ! pointer to variable structure:
    call dios_var_list_get_pointer( filep%var_list, varid, varp, stat )
      
      ! select appropriate routine for each type:
      select case ( filep%ftype )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_ASCII )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_NETCDF )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      
          ! test target type:
          ! convert to required kind before entering NF90_Put_Var, 
          ! otherwise segmentation faults on some machines ...
          select case ( varp%xtype )
        
            case ( DIOS_BYTE )
              allocate( values_int1(size(values,1),size(values,2)) )
              values_int1 = int(values,kind=1)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int1, &
                                            start, count, stride, map )
              deallocate( values_int1 )
          
            case ( DIOS_SHORT )
              allocate( values_int2(size(values,1),size(values,2)) )
              values_int2 = int(values,kind=2)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int2, &
                                            start, count, stride, map )
              deallocate( values_int2 )
          
            case ( DIOS_INT )
              allocate( values_int4(size(values,1),size(values,2)) )
              values_int4 = int(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int4, &
                                            start, count, stride, map )
              
              deallocate( values_int4 )
          
            case ( DIOS_FLOAT )
              allocate( values_real4(size(values,1),size(values,2)) )
              values_real4 = real(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real4, &
                                            start, count, stride, map )
              
              deallocate( values_real4 )
          
            case ( DIOS_DOUBLE )
              allocate( values_real8(size(values,1),size(values,2)) )
              values_real8 = real(values,kind=8)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real8, &
                                            start, count, stride, map )
              
              deallocate( values_real8 )
          
            case default
              write(*, '("not implemented yet for output type : ",i6)') varp%xtype
              stat=1; return
          end select
        
          ! just put; let netcdf library convert the right kind:
          !stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values, &
          !                              start, count, stride, map )
          !

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*, '("unsupported filetype : ",i6)') filep%ftype
          stat=1; return
      end select
        
    ! ok
    stat = 0

  end subroutine DIOS_Put_var_i2_2d


  ! ***


  subroutine dios_get_var_i2_2d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Get_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    integer(2), intent(out)             ::  values(:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(dios_file), pointer         ::  filep
    type(dios_var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:,:)
    integer(2), allocatable         ::  values_int2(:,:)
    integer(4), allocatable         ::  values_int4(:,:)
    integer(8), allocatable         ::  values_int8(:,:)
    real(4), allocatable            ::  values_real4(:,:)
    real(8), allocatable            ::  values_real8(:,:)

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call dios_file_list_get_pointer( file_list, hid, filep, stat )
    
    ! pointer to variable structure:
    call dios_var_list_get_pointer( filep%var_list, varid, varp, stat )
      
    ! select appropriate routine for each type:
    select case ( filep%ftype )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_ASCII )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_NETCDF )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! read values, converted automatically:
        stat = NF90_Get_Var( filep%netcdf_id, varp%netcdf_varid, values, start, count, stride, map )

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        write(*, '("unsupported filetype : ",i6)') filep%ftype
        stat=1; return
    end select
        
    ! ok
    stat = 0

  end subroutine dios_get_var_i2_2d
  

  ! ***

  
  subroutine DIOS_Put_var_i2_3d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Put_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    integer(2), intent(in)            ::  values(:,:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:,:,:)
    integer(2), allocatable         ::  values_int2(:,:,:)
    integer(4), allocatable         ::  values_int4(:,:,:)
    integer(8), allocatable         ::  values_int8(:,:,:)
    real(4), allocatable            ::  values_real4(:,:,:)
    real(8), allocatable            ::  values_real8(:,:,:)
        

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    ! pointer to variable structure:
    call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )

      ! select appropriate routine for each type:
      select case ( filep%ftype )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_ASCII )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_NETCDF )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      
          ! test target type:
          ! convert to required kind before entering NF90_Put_Var, 
          ! otherwise segmentation faults on some machines ...
          select case ( varp%xtype )
        
            case ( DIOS_BYTE )
              allocate( values_int1(size(values,1),size(values,2),size(values,3)) )
              values_int1 = int(values,kind=1)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int1, &
                                            start, count, stride, map )
              
              deallocate( values_int1 )
          
            case ( DIOS_SHORT )
              allocate( values_int2(size(values,1),size(values,2),size(values,3)) )
              values_int2 = int(values,kind=2)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int2, &
                                            start, count, stride, map )
              
              deallocate( values_int2 )
          
            case ( DIOS_INT )
              allocate( values_int4(size(values,1),size(values,2),size(values,3)) )
              values_int4 = int(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int4, &
                                            start, count, stride, map )
              
              deallocate( values_int4 )
          
            case ( DIOS_FLOAT )
              allocate( values_real4(size(values,1),size(values,2),size(values,3)) )
              values_real4 = real(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real4, &
                                            start, count, stride, map )
              
              deallocate( values_real4 )
          
            case ( DIOS_DOUBLE )
              allocate( values_real8(size(values,1),size(values,2),size(values,3)) )
              values_real8 = real(values,kind=8)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real8, &
                                            start, count, stride, map )
              
              deallocate( values_real8 )
          
            case default
              write(*, '("not implemented yet for output type : ",i6)') varp%xtype
              stat=1; return
          end select
        
          ! just put; let netcdf library convert the right kind:
          !stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values, &
          !                              start, count, stride, map )
          !
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*, '("unsupported filetype : ",i6)') filep%ftype
          stat=1; return
      end select
        
    ! ok
    stat = 0

  end subroutine DIOS_Put_var_i2_3d


  ! ***


  subroutine dios_get_var_i2_3d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Get_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    integer(2), intent(out)             ::  values(:,:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(dios_file), pointer         ::  filep
    type(dios_var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:,:,:)
    integer(2), allocatable         ::  values_int2(:,:,:)
    integer(4), allocatable         ::  values_int4(:,:,:)
    integer(8), allocatable         ::  values_int8(:,:,:)
    real(4), allocatable            ::  values_real4(:,:,:)
    real(8), allocatable            ::  values_real8(:,:,:)

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call dios_file_list_get_pointer( file_list, hid, filep, stat )
    
    ! pointer to variable structure:
    call dios_var_list_get_pointer( filep%var_list, varid, varp, stat )
      
    ! select appropriate routine for each type:
    select case ( filep%ftype )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_ASCII )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_NETCDF )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! read values, converted automatically:
        stat = NF90_Get_Var( filep%netcdf_id, varp%netcdf_varid, values, start, count, stride, map )

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        write(*, '("unsupported filetype : ",i6)') filep%ftype
        stat=1; return
    end select
        
    ! ok
    stat = 0

  end subroutine dios_get_var_i2_3d
  

  ! ***

  
  subroutine DIOS_Put_var_i2_4d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Put_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    integer(2), intent(in)            ::  values(:,:,:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:,:,:,:)
    integer(2), allocatable         ::  values_int2(:,:,:,:)
    integer(4), allocatable         ::  values_int4(:,:,:,:)
    integer(8), allocatable         ::  values_int8(:,:,:,:)
    real(4), allocatable            ::  values_real4(:,:,:,:)
    real(8), allocatable            ::  values_real8(:,:,:,:)
        

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    
    ! pointer to variable structure:
    call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
    

      
      ! select appropriate routine for each type:
      select case ( filep%ftype )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_ASCII )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_NETCDF )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      
          ! test target type:
          ! convert to required kind before entering NF90_Put_Var, 
          ! otherwise segmentation faults on some machines ...
          select case ( varp%xtype )
        
            case ( DIOS_BYTE )
              allocate( values_int1(size(values,1),size(values,2),size(values,3),size(values,4)) )
              values_int1 = int(values,kind=1)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int1, &
                                            start, count, stride, map )
              
              deallocate( values_int1 )
          
            case ( DIOS_SHORT )
              allocate( values_int2(size(values,1),size(values,2),size(values,3),size(values,4)) )
              values_int2 = int(values,kind=2)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int2, &
                                            start, count, stride, map )
              
              deallocate( values_int2 )
          
            case ( DIOS_INT )
              allocate( values_int4(size(values,1),size(values,2),size(values,3),size(values,4)) )
              values_int4 = int(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int4, &
                                            start, count, stride, map )
              
              deallocate( values_int4 )
          
            case ( DIOS_FLOAT )
              allocate( values_real4(size(values,1),size(values,2),size(values,3),size(values,4)) )
              values_real4 = real(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real4, &
                                            start, count, stride, map )
              
              deallocate( values_real4 )
          
            case ( DIOS_DOUBLE )
              allocate( values_real8(size(values,1),size(values,2),size(values,3),size(values,4)) )
              values_real8 = real(values,kind=8)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real8, &
                                            start, count, stride, map )
              
              deallocate( values_real8 )
          
            case default
              write(*, '("not implemented yet for output type : ",i6)') varp%xtype
              stat=1; return
          end select
        
          ! just put; let netcdf library convert the right kind:
          !stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values, &
          !                              start, count, stride, map )
          !
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*, '("unsupported filetype : ",i6)') filep%ftype
          stat=1; return
      end select
        
    ! ok
    stat = 0

  end subroutine DIOS_Put_var_i2_4d


  ! ***


  subroutine dios_get_var_i2_4d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Get_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    integer(2), intent(out)             ::  values(:,:,:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(dios_file), pointer         ::  filep
    type(dios_var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:,:,:,:)
    integer(2), allocatable         ::  values_int2(:,:,:,:)
    integer(4), allocatable         ::  values_int4(:,:,:,:)
    integer(8), allocatable         ::  values_int8(:,:,:,:)
    real(4), allocatable            ::  values_real4(:,:,:,:)
    real(8), allocatable            ::  values_real8(:,:,:,:)

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call dios_file_list_get_pointer( file_list, hid, filep, stat )
    
    ! pointer to variable structure:
    call dios_var_list_get_pointer( filep%var_list, varid, varp, stat )
      
    ! select appropriate routine for each type:
    select case ( filep%ftype )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_ASCII )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_NETCDF )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! read values, converted automatically:
        stat = NF90_Get_Var( filep%netcdf_id, varp%netcdf_varid, values, start, count, stride, map )

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        write(*, '("unsupported filetype : ",i6)') filep%ftype
        stat=1; return
    end select
        
    ! ok
    stat = 0

  end subroutine dios_get_var_i2_4d


  ! ***
  
  
  subroutine DIOS_Put_var_i2_5d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Put_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    integer(2), intent(in)            ::  values(:,:,:,:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:,:,:,:,:)
    integer(2), allocatable         ::  values_int2(:,:,:,:,:)
    integer(4), allocatable         ::  values_int4(:,:,:,:,:)
    integer(8), allocatable         ::  values_int8(:,:,:,:,:)
    real(4), allocatable            ::  values_real4(:,:,:,:,:)
    real(8), allocatable            ::  values_real8(:,:,:,:,:)
        

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    
    ! pointer to variable structure:
    call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
    

      
      select case ( filep%ftype )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_ASCII )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_NETCDF )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      
          ! test target type:
          ! convert to required kind before entering NF90_Put_Var, 
          ! otherwise segmentation faults on some machines ...
          select case ( varp%xtype )
        
            case ( DIOS_BYTE )
              allocate( values_int1(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5)) )
              values_int1 = int(values,kind=1)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int1, &
                                            start, count, stride, map )
              
              deallocate( values_int1 )
          
            case ( DIOS_SHORT )
              allocate( values_int2(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5)) )
              values_int2 = int(values,kind=2)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int2, &
                                            start, count, stride, map )
              
              deallocate( values_int2 )
          
            case ( DIOS_INT )
              allocate( values_int4(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5)) )
              values_int4 = int(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int4, &
                                            start, count, stride, map )
              
              deallocate( values_int4 )
          
            case ( DIOS_FLOAT )
              allocate( values_real4(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5)) )
              values_real4 = real(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real4, &
                                            start, count, stride, map )
              
              deallocate( values_real4 )
          
            case ( DIOS_DOUBLE )
              allocate( values_real8(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5)) )
              values_real8 = real(values,kind=8)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real8, &
                                            start, count, stride, map )
              
              deallocate( values_real8 )
          
            case default
              write(*, '("not implemented yet for output type : ",i6)') varp%xtype
              stat=1; return
          end select
        
          ! just put; let netcdf library convert the right kind:
          !stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values, &
          !                              start, count, stride, map )
          !

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*, '("unsupported filetype : ",i6)') filep%ftype
          stat=1; return
      end select
        
    ! ok
    stat = 0

  end subroutine DIOS_Put_var_i2_5d


  ! ***


  subroutine dios_get_var_i2_5d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Get_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    integer(2), intent(out)             ::  values(:,:,:,:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(dios_file), pointer         ::  filep
    type(dios_var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:,:,:,:,:)
    integer(2), allocatable         ::  values_int2(:,:,:,:,:)
    integer(4), allocatable         ::  values_int4(:,:,:,:,:)
    integer(8), allocatable         ::  values_int8(:,:,:,:,:)
    real(4), allocatable            ::  values_real4(:,:,:,:,:)
    real(8), allocatable            ::  values_real8(:,:,:,:,:)

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call dios_file_list_get_pointer( file_list, hid, filep, stat )
    
    ! pointer to variable structure:
    call dios_var_list_get_pointer( filep%var_list, varid, varp, stat )
      
    ! select appropriate routine for each type:
    select case ( filep%ftype )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_ASCII )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_NETCDF )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! read values, converted automatically:
        stat = NF90_Get_Var( filep%netcdf_id, varp%netcdf_varid, values, start, count, stride, map )

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        write(*, '("unsupported filetype : ",i6)') filep%ftype
        stat=1; return
    end select
        
    ! ok
    stat = 0

  end subroutine dios_get_var_i2_5d


  ! ***
  
  
  subroutine DIOS_Put_var_i2_6d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Put_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    integer(2), intent(in)            ::  values(:,:,:,:,:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:,:,:,:,:,:)
    integer(2), allocatable         ::  values_int2(:,:,:,:,:,:)
    integer(4), allocatable         ::  values_int4(:,:,:,:,:,:)
    integer(8), allocatable         ::  values_int8(:,:,:,:,:,:)
    real(4), allocatable            ::  values_real4(:,:,:,:,:,:)
    real(8), allocatable            ::  values_real8(:,:,:,:,:,:)
        

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    
    ! pointer to variable structure:
    call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
    

      select case ( filep%ftype )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_ASCII )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_NETCDF )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      
          ! test target type:
          ! convert to required kind before entering NF90_Put_Var, 
          ! otherwise segmentation faults on some machines ...
          select case ( varp%xtype )
        
            case ( DIOS_BYTE )
              allocate( values_int1(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5),size(values,6)) )
              values_int1 = int(values,kind=1)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int1, &
                                            start, count, stride, map )
              
              deallocate( values_int1 )
          
            case ( DIOS_SHORT )
              allocate( values_int2(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5),size(values,6)) )
              values_int2 = int(values,kind=2)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int2, &
                                            start, count, stride, map )
              
              deallocate( values_int2 )
          
            case ( DIOS_INT )
              allocate( values_int4(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5),size(values,6)) )
              values_int4 = int(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int4, &
                                            start, count, stride, map )
              
              deallocate( values_int4 )
          
            case ( DIOS_FLOAT )
              allocate( values_real4(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5),size(values,6)) )
              values_real4 = real(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real4, &
                                            start, count, stride, map )
              
              deallocate( values_real4 )
          
            case ( DIOS_DOUBLE )
              allocate( values_real8(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5),size(values,6)) )
              values_real8 = real(values,kind=8)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real8, &
                                            start, count, stride, map )
              
              deallocate( values_real8 )
          
            case default
              write(*, '("not implemented yet for output type : ",i6)') varp%xtype
              stat=1; return
          end select
        
          ! just put; let netcdf library convert the right kind:
          !stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values, &
          !                              start, count, stride, map )
          !
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*, '("unsupported filetype : ",i6)') filep%ftype
          stat=1; return
      end select
        
    ! ok
    stat = 0

  end subroutine DIOS_Put_var_i2_6d


  ! ***


  subroutine dios_get_var_i2_6d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Get_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    integer(2), intent(out)             ::  values(:,:,:,:,:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(dios_file), pointer         ::  filep
    type(dios_var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:,:,:,:,:,:)
    integer(2), allocatable         ::  values_int2(:,:,:,:,:,:)
    integer(4), allocatable         ::  values_int4(:,:,:,:,:,:)
    integer(8), allocatable         ::  values_int8(:,:,:,:,:,:)
    real(4), allocatable            ::  values_real4(:,:,:,:,:,:)
    real(8), allocatable            ::  values_real8(:,:,:,:,:,:)

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call dios_file_list_get_pointer( file_list, hid, filep, stat )
    
    ! pointer to variable structure:
    call dios_var_list_get_pointer( filep%var_list, varid, varp, stat )
      
    ! select appropriate routine for each type:
    select case ( filep%ftype )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_ASCII )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_NETCDF )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! read values, converted automatically:
        stat = NF90_Get_Var( filep%netcdf_id, varp%netcdf_varid, values, start, count, stride, map )

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        write(*, '("unsupported filetype : ",i6)') filep%ftype
        stat=1; return
    end select
        
    ! ok
    stat = 0

  end subroutine dios_get_var_i2_6d


  ! ***
  
  
  subroutine DIOS_Put_var_i2_7d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Put_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    integer(2), intent(in)            ::  values(:,:,:,:,:,:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:,:,:,:,:,:,:)
    integer(2), allocatable         ::  values_int2(:,:,:,:,:,:,:)
    integer(4), allocatable         ::  values_int4(:,:,:,:,:,:,:)
    integer(8), allocatable         ::  values_int8(:,:,:,:,:,:,:)
    real(4), allocatable            ::  values_real4(:,:,:,:,:,:,:)
    real(8), allocatable            ::  values_real8(:,:,:,:,:,:,:)
        

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    
    ! pointer to variable structure:
    call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
    

      
      ! select appropriate routine for each type:
      select case ( filep%ftype )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_ASCII )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_NETCDF )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      
          ! test target type:
          ! convert to required kind before entering NF90_Put_Var, 
          ! otherwise segmentation faults on some machines ...
          select case ( varp%xtype )
        
            case ( DIOS_BYTE )
              allocate( values_int1(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5),size(values,6),size(values,7)) )
              values_int1 = int(values,kind=1)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int1, &
                                            start, count, stride, map )
              
              deallocate( values_int1 )
          
            case ( DIOS_SHORT )
              allocate( values_int2(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5),size(values,6),size(values,7)) )
              values_int2 = int(values,kind=2)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int2, &
                                            start, count, stride, map )
              
              deallocate( values_int2 )
          
            case ( DIOS_INT )
              allocate( values_int4(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5),size(values,6),size(values,7)) )
              values_int4 = int(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int4, &
                                            start, count, stride, map )
              
              deallocate( values_int4 )
          
            case ( DIOS_FLOAT )
              allocate( values_real4(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5),size(values,6),size(values,7)) )
              values_real4 = real(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real4, &
                                            start, count, stride, map )
              
              deallocate( values_real4 )
          
            case ( DIOS_DOUBLE )
              allocate( values_real8(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5),size(values,6),size(values,7)) )
              values_real8 = real(values,kind=8)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real8, &
                                            start, count, stride, map )
              
              deallocate( values_real8 )
          
            case default
              write(*, '("not implemented yet for output type : ",i6)') varp%xtype
              stat=1; return
          end select
        
          ! just put; let netcdf library convert the right kind:
          !stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values, &
          !                              start, count, stride, map )
          !
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*, '("unsupported filetype : ",i6)') filep%ftype
          stat=1; return
      end select
        
    ! ok
    stat = 0

  end subroutine DIOS_Put_var_i2_7d


  ! ***


  subroutine dios_get_var_i2_7d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Get_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    integer(2), intent(out)             ::  values(:,:,:,:,:,:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(dios_file), pointer         ::  filep
    type(dios_var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:,:,:,:,:,:,:)
    integer(2), allocatable         ::  values_int2(:,:,:,:,:,:,:)
    integer(4), allocatable         ::  values_int4(:,:,:,:,:,:,:)
    integer(8), allocatable         ::  values_int8(:,:,:,:,:,:,:)
    real(4), allocatable            ::  values_real4(:,:,:,:,:,:,:)
    real(8), allocatable            ::  values_real8(:,:,:,:,:,:,:)

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call dios_file_list_get_pointer( file_list, hid, filep, stat )
    
    ! pointer to variable structure:
    call dios_var_list_get_pointer( filep%var_list, varid, varp, stat )
      
    ! select appropriate routine for each type:
    select case ( filep%ftype )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_ASCII )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_NETCDF )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! read values, converted automatically:
        stat = NF90_Get_Var( filep%netcdf_id, varp%netcdf_varid, values, start, count, stride, map )

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        write(*, '("unsupported filetype : ",i6)') filep%ftype
        stat=1; return
    end select
        
    ! ok
    stat = 0

  end subroutine dios_get_var_i2_7d


  ! ********************************************************************
  ! ***
  ! *** var_i4
  ! ***
  ! ********************************************************************


  subroutine dios_put_var_i4_1d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Put_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    integer(4), intent(in)              ::  values(:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(dios_file), pointer         ::  filep
    type(dios_var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:)
    integer(2), allocatable         ::  values_int2(:)
    integer(4), allocatable         ::  values_int4(:)
    integer(8), allocatable         ::  values_int8(:)
    real(4), allocatable            ::  values_real4(:)
    real(8), allocatable            ::  values_real8(:)
        

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call dios_file_list_get_pointer( file_list, hid, filep, stat )
    
    ! pointer to variable structure:
    call dios_var_list_get_pointer( filep%var_list, varid, varp, stat )

    ! select appropriate routine for each type:
    select case ( filep%ftype )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_ASCII )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_NETCDF )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      
        ! test target type:
        ! convert to required kind wrt var type before entering NF90_Put_Var, 
        ! otherwise segmentation faults on some machines ...
        select case ( varp%xtype )
        
          case ( DIOS_BYTE )
            allocate( values_int1(size(values,1)) )
            values_int1 = int(values,kind=1)
            stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int1, start, count, stride, map )
            deallocate( values_int1 )
        
          case ( DIOS_SHORT )
            allocate( values_int2(size(values,1)) )
            values_int2 = int(values,kind=2)
            stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int2, start, count, stride, map )
            deallocate( values_int2 )
        
          case ( DIOS_INT )
            allocate( values_int4(size(values,1)) )
            values_int4 = int(values,kind=4)
            stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int4, start, count, stride, map )
            deallocate( values_int4 )
        
          case ( DIOS_FLOAT )
            allocate( values_real4(size(values,1)) )
            values_real4 = real(values,kind=4)
            stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real4, start, count, stride, map )
            deallocate( values_real4 )
        
          case ( DIOS_DOUBLE )
            allocate( values_real8(size(values,1)) )
            values_real8 = real(values,kind=8)
            stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real8, start, count, stride, map )
            deallocate( values_real8 )
        
          case default
            write(*, '("not implemented yet for output type : ",i6)') varp%xtype
            stat=1; return
        end select
        
        ! just put; let netcdf library convert the right kind:
        !stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values, start, count, stride, map )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        write(*, '("unsupported filetype : ",i6)') filep%ftype
        stat=1; return
    end select
        
    ! ok
    stat = 0

  end subroutine dios_put_var_i4_1d


  ! ***
  

  subroutine dios_get_var_i4_1d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Get_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    integer(4), intent(out)             ::  values(:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(dios_file), pointer         ::  filep
    type(dios_var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:)
    integer(2), allocatable         ::  values_int2(:)
    integer(4), allocatable         ::  values_int4(:)
    integer(8), allocatable         ::  values_int8(:)
    real(4), allocatable            ::  values_real4(:)
    real(8), allocatable            ::  values_real8(:)

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call dios_file_list_get_pointer( file_list, hid, filep, stat )
    
    ! pointer to variable structure:
    call dios_var_list_get_pointer( filep%var_list, varid, varp, stat )
      
    ! select appropriate routine for each type:
    select case ( filep%ftype )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_ASCII )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_NETCDF )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! read values, converted automatically:
        stat = NF90_Get_Var( filep%netcdf_id, varp%netcdf_varid, values, start, count, stride, map )

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        write(*, '("unsupported filetype : ",i6)') filep%ftype
        stat=1; return
    end select
        
    ! ok
    stat = 0

  end subroutine dios_get_var_i4_1d
  
  
  ! ***

  
  subroutine DIOS_Put_var_i4_2d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Put_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    integer(4), intent(in)            ::  values(:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:,:)
    integer(2), allocatable         ::  values_int2(:,:)
    integer(4), allocatable         ::  values_int4(:,:)
    integer(8), allocatable         ::  values_int8(:,:)
    real(4), allocatable            ::  values_real4(:,:)
    real(8), allocatable            ::  values_real8(:,:)
        

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call dios_file_list_get_pointer( file_list, hid, filep, stat )
    
    ! pointer to variable structure:
    call dios_var_list_get_pointer( filep%var_list, varid, varp, stat )
      
      ! select appropriate routine for each type:
      select case ( filep%ftype )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_ASCII )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_NETCDF )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      
          ! test target type:
          ! convert to required kind before entering NF90_Put_Var, 
          ! otherwise segmentation faults on some machines ...
          select case ( varp%xtype )
        
            case ( DIOS_BYTE )
              allocate( values_int1(size(values,1),size(values,2)) )
              values_int1 = int(values,kind=1)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int1, &
                                            start, count, stride, map )
              deallocate( values_int1 )
          
            case ( DIOS_SHORT )
              allocate( values_int2(size(values,1),size(values,2)) )
              values_int2 = int(values,kind=2)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int2, &
                                            start, count, stride, map )
              deallocate( values_int2 )
          
            case ( DIOS_INT )
              allocate( values_int4(size(values,1),size(values,2)) )
              values_int4 = int(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int4, &
                                            start, count, stride, map )
              
              deallocate( values_int4 )
          
            case ( DIOS_FLOAT )
              allocate( values_real4(size(values,1),size(values,2)) )
              values_real4 = real(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real4, &
                                            start, count, stride, map )
              
              deallocate( values_real4 )
          
            case ( DIOS_DOUBLE )
              allocate( values_real8(size(values,1),size(values,2)) )
              values_real8 = real(values,kind=8)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real8, &
                                            start, count, stride, map )
              
              deallocate( values_real8 )
          
            case default
              write(*, '("not implemented yet for output type : ",i6)') varp%xtype
              stat=1; return
          end select
        
          ! just put; let netcdf library convert the right kind:
          !stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values, &
          !                              start, count, stride, map )
          !

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*, '("unsupported filetype : ",i6)') filep%ftype
          stat=1; return
      end select
        
    ! ok
    stat = 0

  end subroutine DIOS_Put_var_i4_2d


  ! ***


  subroutine dios_get_var_i4_2d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Get_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    integer(4), intent(out)             ::  values(:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(dios_file), pointer         ::  filep
    type(dios_var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:,:)
    integer(2), allocatable         ::  values_int2(:,:)
    integer(4), allocatable         ::  values_int4(:,:)
    integer(8), allocatable         ::  values_int8(:,:)
    real(4), allocatable            ::  values_real4(:,:)
    real(8), allocatable            ::  values_real8(:,:)

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call dios_file_list_get_pointer( file_list, hid, filep, stat )
    
    ! pointer to variable structure:
    call dios_var_list_get_pointer( filep%var_list, varid, varp, stat )
      
    ! select appropriate routine for each type:
    select case ( filep%ftype )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_ASCII )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_NETCDF )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! read values, converted automatically:
        stat = NF90_Get_Var( filep%netcdf_id, varp%netcdf_varid, values, start, count, stride, map )

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        write(*, '("unsupported filetype : ",i6)') filep%ftype
        stat=1; return
    end select
        
    ! ok
    stat = 0

  end subroutine dios_get_var_i4_2d
  

  ! ***

  
  subroutine DIOS_Put_var_i4_3d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Put_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    integer(4), intent(in)            ::  values(:,:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:,:,:)
    integer(2), allocatable         ::  values_int2(:,:,:)
    integer(4), allocatable         ::  values_int4(:,:,:)
    integer(8), allocatable         ::  values_int8(:,:,:)
    real(4), allocatable            ::  values_real4(:,:,:)
    real(8), allocatable            ::  values_real8(:,:,:)
        

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    ! pointer to variable structure:
    call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )

      ! select appropriate routine for each type:
      select case ( filep%ftype )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_ASCII )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_NETCDF )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      
          ! test target type:
          ! convert to required kind before entering NF90_Put_Var, 
          ! otherwise segmentation faults on some machines ...
          select case ( varp%xtype )
        
            case ( DIOS_BYTE )
              allocate( values_int1(size(values,1),size(values,2),size(values,3)) )
              values_int1 = int(values,kind=1)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int1, &
                                            start, count, stride, map )
              
              deallocate( values_int1 )
          
            case ( DIOS_SHORT )
              allocate( values_int2(size(values,1),size(values,2),size(values,3)) )
              values_int2 = int(values,kind=2)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int2, &
                                            start, count, stride, map )
              
              deallocate( values_int2 )
          
            case ( DIOS_INT )
              allocate( values_int4(size(values,1),size(values,2),size(values,3)) )
              values_int4 = int(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int4, &
                                            start, count, stride, map )
              
              deallocate( values_int4 )
          
            case ( DIOS_FLOAT )
              allocate( values_real4(size(values,1),size(values,2),size(values,3)) )
              values_real4 = real(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real4, &
                                            start, count, stride, map )
              
              deallocate( values_real4 )
          
            case ( DIOS_DOUBLE )
              allocate( values_real8(size(values,1),size(values,2),size(values,3)) )
              values_real8 = real(values,kind=8)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real8, &
                                            start, count, stride, map )
              
              deallocate( values_real8 )
          
            case default
              write(*, '("not implemented yet for output type : ",i6)') varp%xtype
              stat=1; return
          end select
        
          ! just put; let netcdf library convert the right kind:
          !stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values, &
          !                              start, count, stride, map )
          !
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*, '("unsupported filetype : ",i6)') filep%ftype
          stat=1; return
      end select
        
    ! ok
    stat = 0

  end subroutine DIOS_Put_var_i4_3d


  ! ***


  subroutine dios_get_var_i4_3d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Get_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    integer(4), intent(out)             ::  values(:,:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(dios_file), pointer         ::  filep
    type(dios_var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:,:,:)
    integer(2), allocatable         ::  values_int2(:,:,:)
    integer(4), allocatable         ::  values_int4(:,:,:)
    integer(8), allocatable         ::  values_int8(:,:,:)
    real(4), allocatable            ::  values_real4(:,:,:)
    real(8), allocatable            ::  values_real8(:,:,:)

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call dios_file_list_get_pointer( file_list, hid, filep, stat )
    
    ! pointer to variable structure:
    call dios_var_list_get_pointer( filep%var_list, varid, varp, stat )
      
    ! select appropriate routine for each type:
    select case ( filep%ftype )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_ASCII )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_NETCDF )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! read values, converted automatically:
        stat = NF90_Get_Var( filep%netcdf_id, varp%netcdf_varid, values, start, count, stride, map )

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        write(*, '("unsupported filetype : ",i6)') filep%ftype
        stat=1; return
    end select
        
    ! ok
    stat = 0

  end subroutine dios_get_var_i4_3d
  

  ! ***

  
  subroutine DIOS_Put_var_i4_4d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Put_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    integer(4), intent(in)            ::  values(:,:,:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:,:,:,:)
    integer(2), allocatable         ::  values_int2(:,:,:,:)
    integer(4), allocatable         ::  values_int4(:,:,:,:)
    integer(8), allocatable         ::  values_int8(:,:,:,:)
    real(4), allocatable            ::  values_real4(:,:,:,:)
    real(8), allocatable            ::  values_real8(:,:,:,:)
        

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    
    ! pointer to variable structure:
    call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
    

      
      ! select appropriate routine for each type:
      select case ( filep%ftype )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_ASCII )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_NETCDF )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      
          ! test target type:
          ! convert to required kind before entering NF90_Put_Var, 
          ! otherwise segmentation faults on some machines ...
          select case ( varp%xtype )
        
            case ( DIOS_BYTE )
              allocate( values_int1(size(values,1),size(values,2),size(values,3),size(values,4)) )
              values_int1 = int(values,kind=1)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int1, &
                                            start, count, stride, map )
              
              deallocate( values_int1 )
          
            case ( DIOS_SHORT )
              allocate( values_int2(size(values,1),size(values,2),size(values,3),size(values,4)) )
              values_int2 = int(values,kind=2)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int2, &
                                            start, count, stride, map )
              
              deallocate( values_int2 )
          
            case ( DIOS_INT )
              allocate( values_int4(size(values,1),size(values,2),size(values,3),size(values,4)) )
              values_int4 = int(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int4, &
                                            start, count, stride, map )
              
              deallocate( values_int4 )
          
            case ( DIOS_FLOAT )
              allocate( values_real4(size(values,1),size(values,2),size(values,3),size(values,4)) )
              values_real4 = real(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real4, &
                                            start, count, stride, map )
              
              deallocate( values_real4 )
          
            case ( DIOS_DOUBLE )
              allocate( values_real8(size(values,1),size(values,2),size(values,3),size(values,4)) )
              values_real8 = real(values,kind=8)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real8, &
                                            start, count, stride, map )
              
              deallocate( values_real8 )
          
            case default
              write(*, '("not implemented yet for output type : ",i6)') varp%xtype
              stat=1; return
          end select
        
          ! just put; let netcdf library convert the right kind:
          !stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values, &
          !                              start, count, stride, map )
          !
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*, '("unsupported filetype : ",i6)') filep%ftype
          stat=1; return
      end select
        
    ! ok
    stat = 0

  end subroutine DIOS_Put_var_i4_4d


  ! ***


  subroutine dios_get_var_i4_4d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Get_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    integer(4), intent(out)             ::  values(:,:,:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(dios_file), pointer         ::  filep
    type(dios_var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:,:,:,:)
    integer(2), allocatable         ::  values_int2(:,:,:,:)
    integer(4), allocatable         ::  values_int4(:,:,:,:)
    integer(8), allocatable         ::  values_int8(:,:,:,:)
    real(4), allocatable            ::  values_real4(:,:,:,:)
    real(8), allocatable            ::  values_real8(:,:,:,:)

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call dios_file_list_get_pointer( file_list, hid, filep, stat )
    
    ! pointer to variable structure:
    call dios_var_list_get_pointer( filep%var_list, varid, varp, stat )
      
    ! select appropriate routine for each type:
    select case ( filep%ftype )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_ASCII )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_NETCDF )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! read values, converted automatically:
        stat = NF90_Get_Var( filep%netcdf_id, varp%netcdf_varid, values, start, count, stride, map )

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        write(*, '("unsupported filetype : ",i6)') filep%ftype
        stat=1; return
    end select
        
    ! ok
    stat = 0

  end subroutine dios_get_var_i4_4d


  ! ***
  
  
  subroutine DIOS_Put_var_i4_5d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Put_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    integer(4), intent(in)            ::  values(:,:,:,:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:,:,:,:,:)
    integer(2), allocatable         ::  values_int2(:,:,:,:,:)
    integer(4), allocatable         ::  values_int4(:,:,:,:,:)
    integer(8), allocatable         ::  values_int8(:,:,:,:,:)
    real(4), allocatable            ::  values_real4(:,:,:,:,:)
    real(8), allocatable            ::  values_real8(:,:,:,:,:)
        

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    
    ! pointer to variable structure:
    call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
    

      
      select case ( filep%ftype )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_ASCII )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_NETCDF )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      
          ! test target type:
          ! convert to required kind before entering NF90_Put_Var, 
          ! otherwise segmentation faults on some machines ...
          select case ( varp%xtype )
        
            case ( DIOS_BYTE )
              allocate( values_int1(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5)) )
              values_int1 = int(values,kind=1)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int1, &
                                            start, count, stride, map )
              
              deallocate( values_int1 )
          
            case ( DIOS_SHORT )
              allocate( values_int2(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5)) )
              values_int2 = int(values,kind=2)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int2, &
                                            start, count, stride, map )
              
              deallocate( values_int2 )
          
            case ( DIOS_INT )
              allocate( values_int4(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5)) )
              values_int4 = int(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int4, &
                                            start, count, stride, map )
              
              deallocate( values_int4 )
          
            case ( DIOS_FLOAT )
              allocate( values_real4(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5)) )
              values_real4 = real(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real4, &
                                            start, count, stride, map )
              
              deallocate( values_real4 )
          
            case ( DIOS_DOUBLE )
              allocate( values_real8(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5)) )
              values_real8 = real(values,kind=8)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real8, &
                                            start, count, stride, map )
              
              deallocate( values_real8 )
          
            case default
              write(*, '("not implemented yet for output type : ",i6)') varp%xtype
              stat=1; return
          end select
        
          ! just put; let netcdf library convert the right kind:
          !stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values, &
          !                              start, count, stride, map )
          !

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*, '("unsupported filetype : ",i6)') filep%ftype
          stat=1; return
      end select
        
    ! ok
    stat = 0

  end subroutine DIOS_Put_var_i4_5d


  ! ***


  subroutine dios_get_var_i4_5d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Get_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    integer(4), intent(out)             ::  values(:,:,:,:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(dios_file), pointer         ::  filep
    type(dios_var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:,:,:,:,:)
    integer(2), allocatable         ::  values_int2(:,:,:,:,:)
    integer(4), allocatable         ::  values_int4(:,:,:,:,:)
    integer(8), allocatable         ::  values_int8(:,:,:,:,:)
    real(4), allocatable            ::  values_real4(:,:,:,:,:)
    real(8), allocatable            ::  values_real8(:,:,:,:,:)

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call dios_file_list_get_pointer( file_list, hid, filep, stat )
    
    ! pointer to variable structure:
    call dios_var_list_get_pointer( filep%var_list, varid, varp, stat )
      
    ! select appropriate routine for each type:
    select case ( filep%ftype )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_ASCII )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_NETCDF )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! read values, converted automatically:
        stat = NF90_Get_Var( filep%netcdf_id, varp%netcdf_varid, values, start, count, stride, map )

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        write(*, '("unsupported filetype : ",i6)') filep%ftype
        stat=1; return
    end select
        
    ! ok
    stat = 0

  end subroutine dios_get_var_i4_5d


  ! ***
  
  
  subroutine DIOS_Put_var_i4_6d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Put_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    integer(4), intent(in)            ::  values(:,:,:,:,:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:,:,:,:,:,:)
    integer(2), allocatable         ::  values_int2(:,:,:,:,:,:)
    integer(4), allocatable         ::  values_int4(:,:,:,:,:,:)
    integer(8), allocatable         ::  values_int8(:,:,:,:,:,:)
    real(4), allocatable            ::  values_real4(:,:,:,:,:,:)
    real(8), allocatable            ::  values_real8(:,:,:,:,:,:)
        

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    
    ! pointer to variable structure:
    call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
    

      select case ( filep%ftype )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_ASCII )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_NETCDF )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      
          ! test target type:
          ! convert to required kind before entering NF90_Put_Var, 
          ! otherwise segmentation faults on some machines ...
          select case ( varp%xtype )
        
            case ( DIOS_BYTE )
              allocate( values_int1(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5),size(values,6)) )
              values_int1 = int(values,kind=1)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int1, &
                                            start, count, stride, map )
              
              deallocate( values_int1 )
          
            case ( DIOS_SHORT )
              allocate( values_int2(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5),size(values,6)) )
              values_int2 = int(values,kind=2)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int2, &
                                            start, count, stride, map )
              
              deallocate( values_int2 )
          
            case ( DIOS_INT )
              allocate( values_int4(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5),size(values,6)) )
              values_int4 = int(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int4, &
                                            start, count, stride, map )
              
              deallocate( values_int4 )
          
            case ( DIOS_FLOAT )
              allocate( values_real4(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5),size(values,6)) )
              values_real4 = real(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real4, &
                                            start, count, stride, map )
              
              deallocate( values_real4 )
          
            case ( DIOS_DOUBLE )
              allocate( values_real8(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5),size(values,6)) )
              values_real8 = real(values,kind=8)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real8, &
                                            start, count, stride, map )
              
              deallocate( values_real8 )
          
            case default
              write(*, '("not implemented yet for output type : ",i6)') varp%xtype
              stat=1; return
          end select
        
          ! just put; let netcdf library convert the right kind:
          !stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values, &
          !                              start, count, stride, map )
          !
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*, '("unsupported filetype : ",i6)') filep%ftype
          stat=1; return
      end select
        
    ! ok
    stat = 0

  end subroutine DIOS_Put_var_i4_6d


  ! ***


  subroutine dios_get_var_i4_6d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Get_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    integer(4), intent(out)             ::  values(:,:,:,:,:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(dios_file), pointer         ::  filep
    type(dios_var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:,:,:,:,:,:)
    integer(2), allocatable         ::  values_int2(:,:,:,:,:,:)
    integer(4), allocatable         ::  values_int4(:,:,:,:,:,:)
    integer(8), allocatable         ::  values_int8(:,:,:,:,:,:)
    real(4), allocatable            ::  values_real4(:,:,:,:,:,:)
    real(8), allocatable            ::  values_real8(:,:,:,:,:,:)

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call dios_file_list_get_pointer( file_list, hid, filep, stat )
    
    ! pointer to variable structure:
    call dios_var_list_get_pointer( filep%var_list, varid, varp, stat )
      
    ! select appropriate routine for each type:
    select case ( filep%ftype )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_ASCII )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_NETCDF )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! read values, converted automatically:
        stat = NF90_Get_Var( filep%netcdf_id, varp%netcdf_varid, values, start, count, stride, map )

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        write(*, '("unsupported filetype : ",i6)') filep%ftype
        stat=1; return
    end select
        
    ! ok
    stat = 0

  end subroutine dios_get_var_i4_6d


  ! ***
  
  
  subroutine DIOS_Put_var_i4_7d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Put_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    integer(4), intent(in)            ::  values(:,:,:,:,:,:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:,:,:,:,:,:,:)
    integer(2), allocatable         ::  values_int2(:,:,:,:,:,:,:)
    integer(4), allocatable         ::  values_int4(:,:,:,:,:,:,:)
    integer(8), allocatable         ::  values_int8(:,:,:,:,:,:,:)
    real(4), allocatable            ::  values_real4(:,:,:,:,:,:,:)
    real(8), allocatable            ::  values_real8(:,:,:,:,:,:,:)
        

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    
    ! pointer to variable structure:
    call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
    

      
      ! select appropriate routine for each type:
      select case ( filep%ftype )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_ASCII )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_NETCDF )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      
          ! test target type:
          ! convert to required kind before entering NF90_Put_Var, 
          ! otherwise segmentation faults on some machines ...
          select case ( varp%xtype )
        
            case ( DIOS_BYTE )
              allocate( values_int1(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5),size(values,6),size(values,7)) )
              values_int1 = int(values,kind=1)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int1, &
                                            start, count, stride, map )
              
              deallocate( values_int1 )
          
            case ( DIOS_SHORT )
              allocate( values_int2(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5),size(values,6),size(values,7)) )
              values_int2 = int(values,kind=2)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int2, &
                                            start, count, stride, map )
              
              deallocate( values_int2 )
          
            case ( DIOS_INT )
              allocate( values_int4(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5),size(values,6),size(values,7)) )
              values_int4 = int(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int4, &
                                            start, count, stride, map )
              
              deallocate( values_int4 )
          
            case ( DIOS_FLOAT )
              allocate( values_real4(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5),size(values,6),size(values,7)) )
              values_real4 = real(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real4, &
                                            start, count, stride, map )
              
              deallocate( values_real4 )
          
            case ( DIOS_DOUBLE )
              allocate( values_real8(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5),size(values,6),size(values,7)) )
              values_real8 = real(values,kind=8)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real8, &
                                            start, count, stride, map )
              
              deallocate( values_real8 )
          
            case default
              write(*, '("not implemented yet for output type : ",i6)') varp%xtype
              stat=1; return
          end select
        
          ! just put; let netcdf library convert the right kind:
          !stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values, &
          !                              start, count, stride, map )
          !
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*, '("unsupported filetype : ",i6)') filep%ftype
          stat=1; return
      end select
        
    ! ok
    stat = 0

  end subroutine DIOS_Put_var_i4_7d


  ! ***


  subroutine dios_get_var_i4_7d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Get_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    integer(4), intent(out)             ::  values(:,:,:,:,:,:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(dios_file), pointer         ::  filep
    type(dios_var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:,:,:,:,:,:,:)
    integer(2), allocatable         ::  values_int2(:,:,:,:,:,:,:)
    integer(4), allocatable         ::  values_int4(:,:,:,:,:,:,:)
    integer(8), allocatable         ::  values_int8(:,:,:,:,:,:,:)
    real(4), allocatable            ::  values_real4(:,:,:,:,:,:,:)
    real(8), allocatable            ::  values_real8(:,:,:,:,:,:,:)

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call dios_file_list_get_pointer( file_list, hid, filep, stat )
    
    ! pointer to variable structure:
    call dios_var_list_get_pointer( filep%var_list, varid, varp, stat )
      
    ! select appropriate routine for each type:
    select case ( filep%ftype )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_ASCII )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_NETCDF )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! read values, converted automatically:
        stat = NF90_Get_Var( filep%netcdf_id, varp%netcdf_varid, values, start, count, stride, map )

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        write(*, '("unsupported filetype : ",i6)') filep%ftype
        stat=1; return
    end select
        
    ! ok
    stat = 0

  end subroutine dios_get_var_i4_7d
  

  
  ! ********************************************************************
  ! ***
  ! *** put_var_r4
  ! ***
  ! ********************************************************************

  

  subroutine DIOS_Put_Var_r4_1d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Put_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    real(4), intent(in)            ::  values(:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:)
    integer(2), allocatable         ::  values_int2(:)
    integer(4), allocatable         ::  values_int4(:)
    integer(8), allocatable         ::  values_int8(:)
    real(4), allocatable            ::  values_real4(:)
    real(8), allocatable            ::  values_real8(:)
        

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    
    ! pointer to variable structure:
    call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
    

      
      ! select appropriate routine for each type:
      select case ( filep%ftype )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_ASCII )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_NETCDF )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      
          ! test target type:
          ! convert to required kind before entering NF90_Put_Var, 
          ! otherwise segmentation faults on some machines ...
          select case ( varp%xtype )
        
            case ( DIOS_BYTE )
              allocate( values_int1(size(values,1)) )
              values_int1 = int(values,kind=1)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int1, &
                                            start, count, stride, map )
              
              deallocate( values_int1 )
          
            case ( DIOS_SHORT )
              allocate( values_int2(size(values,1)) )
              values_int2 = int(values,kind=2)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int2, &
                                            start, count, stride, map )
              
              deallocate( values_int2 )
          
            case ( DIOS_INT )
              allocate( values_int4(size(values,1)) )
              values_int4 = int(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int4, &
                                            start, count, stride, map )
              
              deallocate( values_int4 )
          
            case ( DIOS_FLOAT )
              allocate( values_real4(size(values,1)) )
              values_real4 = real(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real4, &
                                            start, count, stride, map )
              
              deallocate( values_real4 )
          
            case ( DIOS_DOUBLE )
              allocate( values_real8(size(values,1)) )
              values_real8 = real(values,kind=8)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real8, &
                                            start, count, stride, map )
              
              deallocate( values_real8 )
          
            case default
              write(*, '("not implemented yet for output type : ",i6)') varp%xtype
              stat=1; return
          end select
        
          ! just put; let netcdf library convert the right kind:
          !stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values, &
          !                              start, count, stride, map )
          !
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*, '("unsupported filetype : ",i6)') filep%ftype
          stat=1; return
      end select
        
    ! ok
    stat = 0

  end subroutine DIOS_Put_Var_r4_1d


  ! ***
  

  subroutine DIOS_Get_Var_r4_1d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Get_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    real(4), intent(out)           ::  values(:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:)
    integer(2), allocatable         ::  values_int2(:)
    integer(4), allocatable         ::  values_int4(:)
    integer(8), allocatable         ::  values_int8(:)
    real(4), allocatable            ::  values_real4(:)
    real(8), allocatable            ::  values_real8(:)
          

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    
    ! pointer to variable structure:
    call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
    

      
      ! select appropriate routine for each type:
      select case ( filep%ftype )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_ASCII )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_NETCDF )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          ! read values, converted automatically:
          stat = NF90_Get_Var( filep%netcdf_id, varp%netcdf_varid, values, start, count, stride, map )
          
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*, '("unsupported filetype : ",i6)') filep%ftype
          stat=1; return
      end select
        
    ! ok
    stat = 0

  end subroutine DIOS_Get_Var_r4_1d
  
  
  ! ***

  
  subroutine DIOS_Put_Var_r4_2d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Put_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    real(4), intent(in)            ::  values(:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:,:)
    integer(2), allocatable         ::  values_int2(:,:)
    integer(4), allocatable         ::  values_int4(:,:)
    integer(8), allocatable         ::  values_int8(:,:)
    real(4), allocatable            ::  values_real4(:,:)
    real(8), allocatable            ::  values_real8(:,:)
        

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    
    ! pointer to variable structure:
    call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
    

      
      ! select appropriate routine for each type:
      select case ( filep%ftype )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_ASCII )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_NETCDF )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      
          ! test target type:
          ! convert to required kind before entering NF90_Put_Var, 
          ! otherwise segmentation faults on some machines ...
          select case ( varp%xtype )
        
            case ( DIOS_BYTE )
              allocate( values_int1(size(values,1),size(values,2)) )
              values_int1 = int(values,kind=1)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int1, &
                                            start, count, stride, map )
              
              deallocate( values_int1 )
          
            case ( DIOS_SHORT )
              allocate( values_int2(size(values,1),size(values,2)) )
              values_int2 = int(values,kind=2)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int2, &
                                            start, count, stride, map )
              
              deallocate( values_int2 )
          
            case ( DIOS_INT )
              allocate( values_int4(size(values,1),size(values,2)) )
              values_int4 = int(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int4, &
                                            start, count, stride, map )
              
              deallocate( values_int4 )
          
            case ( DIOS_FLOAT )
              allocate( values_real4(size(values,1),size(values,2)) )
              values_real4 = real(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real4, &
                                            start, count, stride, map )
              
              deallocate( values_real4 )
          
            case ( DIOS_DOUBLE )
              allocate( values_real8(size(values,1),size(values,2)) )
              values_real8 = real(values,kind=8)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real8, &
                                            start, count, stride, map )
              
              deallocate( values_real8 )
          
            case default
              write(*, '("not implemented yet for output type : ",i6)') varp%xtype
              stat=1; return
          end select
        
          ! just put; let netcdf library convert the right kind:
          !stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values, &
          !                              start, count, stride, map )
          !
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*, '("unsupported filetype : ",i6)') filep%ftype
          stat=1; return
      end select
        
    ! ok
    stat = 0

  end subroutine DIOS_Put_Var_r4_2d


  ! ***
  

  subroutine DIOS_Get_Var_r4_2d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Get_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    real(4), intent(out)           ::  values(:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:,:)
    integer(2), allocatable         ::  values_int2(:,:)
    integer(4), allocatable         ::  values_int4(:,:)
    integer(8), allocatable         ::  values_int8(:,:)
    real(4), allocatable            ::  values_real4(:,:)
    real(8), allocatable            ::  values_real8(:,:)
          

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    
    ! pointer to variable structure:
    call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
    

      
      ! select appropriate routine for each type:
      select case ( filep%ftype )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_ASCII )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_NETCDF )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          ! read values, converted automatically:
          stat = NF90_Get_Var( filep%netcdf_id, varp%netcdf_varid, values, start, count, stride, map )
          
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*, '("unsupported filetype : ",i6)') filep%ftype
          stat=1; return
      end select
        
    ! ok
    stat = 0

  end subroutine DIOS_Get_Var_r4_2d
  
  
  ! ***

  
  subroutine DIOS_Put_Var_r4_3d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Put_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    real(4), intent(in)            ::  values(:,:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:,:,:)
    integer(2), allocatable         ::  values_int2(:,:,:)
    integer(4), allocatable         ::  values_int4(:,:,:)
    integer(8), allocatable         ::  values_int8(:,:,:)
    real(4), allocatable            ::  values_real4(:,:,:)
    real(8), allocatable            ::  values_real8(:,:,:)
        

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    
    ! pointer to variable structure:
    call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
    

      
      ! select appropriate routine for each type:
      select case ( filep%ftype )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_ASCII )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_NETCDF )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      
          ! test target type:
          ! convert to required kind before entering NF90_Put_Var, 
          ! otherwise segmentation faults on some machines ...
          select case ( varp%xtype )
        
            case ( DIOS_BYTE )
              allocate( values_int1(size(values,1),size(values,2),size(values,3)) )
              values_int1 = int(values,kind=1)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int1, &
                                            start, count, stride, map )
              
              deallocate( values_int1 )
          
            case ( DIOS_SHORT )
              allocate( values_int2(size(values,1),size(values,2),size(values,3)) )
              values_int2 = int(values,kind=2)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int2, &
                                            start, count, stride, map )
              
              deallocate( values_int2 )
          
            case ( DIOS_INT )
              allocate( values_int4(size(values,1),size(values,2),size(values,3)) )
              values_int4 = int(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int4, &
                                            start, count, stride, map )
              
              deallocate( values_int4 )
          
            case ( DIOS_FLOAT )
              allocate( values_real4(size(values,1),size(values,2),size(values,3)) )
              values_real4 = real(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real4, &
                                            start, count, stride, map )
              
              deallocate( values_real4 )
          
            case ( DIOS_DOUBLE )
              allocate( values_real8(size(values,1),size(values,2),size(values,3)) )
              values_real8 = real(values,kind=8)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real8, &
                                            start, count, stride, map )
              
              deallocate( values_real8 )
          
            case default
              write(*, '("not implemented yet for output type : ",i6)') varp%xtype
              stat=1; return
          end select
        
          ! just put; let netcdf library convert the right kind:
          !stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values, &
          !                              start, count, stride, map )
          !
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*, '("unsupported filetype : ",i6)') filep%ftype
          stat=1; return
      end select
        
    ! ok
    stat = 0

  end subroutine DIOS_Put_Var_r4_3d


  ! ***
  

  subroutine DIOS_Get_Var_r4_3d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Get_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    real(4), intent(out)           ::  values(:,:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:,:,:)
    integer(2), allocatable         ::  values_int2(:,:,:)
    integer(4), allocatable         ::  values_int4(:,:,:)
    integer(8), allocatable         ::  values_int8(:,:,:)
    real(4), allocatable            ::  values_real4(:,:,:)
    real(8), allocatable            ::  values_real8(:,:,:)
          

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    
    ! pointer to variable structure:
    call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
    

      
      ! select appropriate routine for each type:
      select case ( filep%ftype )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_ASCII )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_NETCDF )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          ! read values, converted automatically:
          stat = NF90_Get_Var( filep%netcdf_id, varp%netcdf_varid, values, start, count, stride, map )
          
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*, '("unsupported filetype : ",i6)') filep%ftype
          stat=1; return
      end select
        
    ! ok
    stat = 0

  end subroutine DIOS_Get_Var_r4_3d
  
  
  ! ***

  
  subroutine DIOS_Put_Var_r4_4d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Put_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    real(4), intent(in)            ::  values(:,:,:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:,:,:,:)
    integer(2), allocatable         ::  values_int2(:,:,:,:)
    integer(4), allocatable         ::  values_int4(:,:,:,:)
    integer(8), allocatable         ::  values_int8(:,:,:,:)
    real(4), allocatable            ::  values_real4(:,:,:,:)
    real(8), allocatable            ::  values_real8(:,:,:,:)
        

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    
    ! pointer to variable structure:
    call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
    

      
      ! select appropriate routine for each type:
      select case ( filep%ftype )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_ASCII )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_NETCDF )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      
          ! test target type:
          ! convert to required kind before entering NF90_Put_Var, 
          ! otherwise segmentation faults on some machines ...
          select case ( varp%xtype )
        
            case ( DIOS_BYTE )
              allocate( values_int1(size(values,1),size(values,2),size(values,3),size(values,4)) )
              values_int1 = int(values,kind=1)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int1, &
                                            start, count, stride, map )
              
              deallocate( values_int1 )
          
            case ( DIOS_SHORT )
              allocate( values_int2(size(values,1),size(values,2),size(values,3),size(values,4)) )
              values_int2 = int(values,kind=2)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int2, &
                                            start, count, stride, map )
              
              deallocate( values_int2 )
          
            case ( DIOS_INT )
              allocate( values_int4(size(values,1),size(values,2),size(values,3),size(values,4)) )
              values_int4 = int(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int4, &
                                            start, count, stride, map )
              
              deallocate( values_int4 )
          
            case ( DIOS_FLOAT )
              allocate( values_real4(size(values,1),size(values,2),size(values,3),size(values,4)) )
              values_real4 = real(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real4, &
                                            start, count, stride, map )
              
              deallocate( values_real4 )
          
            case ( DIOS_DOUBLE )
              allocate( values_real8(size(values,1),size(values,2),size(values,3),size(values,4)) )
              values_real8 = real(values,kind=8)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real8, &
                                            start, count, stride, map )
              
              deallocate( values_real8 )
          
            case default
              write(*, '("not implemented yet for output type : ",i6)') varp%xtype
              stat=1; return
          end select
        
          ! just put; let netcdf library convert the right kind:
          !stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values, &
          !                              start, count, stride, map )
          !
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*, '("unsupported filetype : ",i6)') filep%ftype
          stat=1; return
      end select
        
    ! ok
    stat = 0

  end subroutine DIOS_Put_Var_r4_4d


  ! ***
  

  subroutine DIOS_Get_Var_r4_4d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Get_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    real(4), intent(out)           ::  values(:,:,:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:,:,:,:)
    integer(2), allocatable         ::  values_int2(:,:,:,:)
    integer(4), allocatable         ::  values_int4(:,:,:,:)
    integer(8), allocatable         ::  values_int8(:,:,:,:)
    real(4), allocatable            ::  values_real4(:,:,:,:)
    real(8), allocatable            ::  values_real8(:,:,:,:)
          
    ! --- begin --------------------------------------

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    
    ! pointer to variable structure:
    call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
    

      
      ! select appropriate routine for each type:
      select case ( filep%ftype )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_ASCII )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_NETCDF )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          ! read values, converted automatically:
          stat = NF90_Get_Var( filep%netcdf_id, varp%netcdf_varid, values, start, count, stride, map )
          
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*, '("unsupported filetype : ",i6)') filep%ftype
          stat=1; return
      end select
        
    ! ok
    stat = 0

  end subroutine DIOS_Get_Var_r4_4d
  
  
  ! ***

  
  subroutine DIOS_Put_Var_r4_5d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Put_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    real(4), intent(in)            ::  values(:,:,:,:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp

    integer(1), allocatable         ::  values_int1(:,:,:,:,:)
    integer(2), allocatable         ::  values_int2(:,:,:,:,:)
    integer(4), allocatable         ::  values_int4(:,:,:,:,:)
    integer(8), allocatable         ::  values_int8(:,:,:,:,:)
    real(4), allocatable            ::  values_real4(:,:,:,:,:)
    real(8), allocatable            ::  values_real8(:,:,:,:,:)
        

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    
    ! pointer to variable structure:
    call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
    

      
      ! select appropriate routine for each type:
      select case ( filep%ftype )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_ASCII )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_NETCDF )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      
          ! test target type:
          ! convert to required kind before entering NF90_Put_Var, 
          ! otherwise segmentation faults on some machines ...
          select case ( varp%xtype )
        
            case ( DIOS_BYTE )
              allocate( values_int1(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5)) )
              values_int1 = int(values,kind=1)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int1, &
                                            start, count, stride, map )
              
              deallocate( values_int1 )
          
            case ( DIOS_SHORT )
              allocate( values_int2(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5)) )
              values_int2 = int(values,kind=2)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int2, &
                                            start, count, stride, map )
              
              deallocate( values_int2 )
          
            case ( DIOS_INT )
              allocate( values_int4(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5)) )
              values_int4 = int(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int4, &
                                            start, count, stride, map )
              
              deallocate( values_int4 )
          
            case ( DIOS_FLOAT )
              allocate( values_real4(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5)) )
              values_real4 = real(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real4, &
                                            start, count, stride, map )
              
              deallocate( values_real4 )
          
            case ( DIOS_DOUBLE )
              allocate( values_real8(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5)) )
              values_real8 = real(values,kind=8)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real8, &
                                            start, count, stride, map )
              
              deallocate( values_real8 )
          
            case default
              write(*, '("not implemented yet for output type : ",i6)') varp%xtype
              stat=1; return
          end select
        
          ! just put; let netcdf library convert the right kind:
          !stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values, &
          !                              start, count, stride, map )
          !
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*, '("unsupported filetype : ",i6)') filep%ftype
          stat=1; return
      end select
        
    ! ok
    stat = 0

  end subroutine DIOS_Put_Var_r4_5d


  ! ***
  

  subroutine DIOS_Get_Var_r4_5d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Get_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    real(4), intent(out)           ::  values(:,:,:,:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:,:,:,:,:)
    integer(2), allocatable         ::  values_int2(:,:,:,:,:)
    integer(4), allocatable         ::  values_int4(:,:,:,:,:)
    integer(8), allocatable         ::  values_int8(:,:,:,:,:)
    real(4), allocatable            ::  values_real4(:,:,:,:,:)
    real(8), allocatable            ::  values_real8(:,:,:,:,:)
          

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    
    ! pointer to variable structure:
    call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
    

      
      ! select appropriate routine for each type:
      select case ( filep%ftype )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_ASCII )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_NETCDF )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          ! read values, converted automatically:
          stat = NF90_Get_Var( filep%netcdf_id, varp%netcdf_varid, values, start, count, stride, map )
          
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*, '("unsupported filetype : ",i6)') filep%ftype
          stat=1; return
      end select
        
    ! ok
    stat = 0

  end subroutine DIOS_Get_Var_r4_5d
  
  
  ! ***

  
  subroutine DIOS_Put_Var_r4_6d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Put_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    real(4), intent(in)            ::  values(:,:,:,:,:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:,:,:,:,:,:)
    integer(2), allocatable         ::  values_int2(:,:,:,:,:,:)
    integer(4), allocatable         ::  values_int4(:,:,:,:,:,:)
    integer(8), allocatable         ::  values_int8(:,:,:,:,:,:)
    real(4), allocatable            ::  values_real4(:,:,:,:,:,:)
    real(8), allocatable            ::  values_real8(:,:,:,:,:,:)
        

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    
    ! pointer to variable structure:
    call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
    

      
      ! select appropriate routine for each type:
      select case ( filep%ftype )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_ASCII )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_NETCDF )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      
          ! test target type:
          ! convert to required kind before entering NF90_Put_Var, 
          ! otherwise segmentation faults on some machines ...
          select case ( varp%xtype )
        
            case ( DIOS_BYTE )
              allocate( values_int1(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5),size(values,6)) )
              values_int1 = int(values,kind=1)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int1, &
                                            start, count, stride, map )
              
              deallocate( values_int1 )
          
            case ( DIOS_SHORT )
              allocate( values_int2(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5),size(values,6)) )
              values_int2 = int(values,kind=2)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int2, &
                                            start, count, stride, map )
              
              deallocate( values_int2 )
          
            case ( DIOS_INT )
              allocate( values_int4(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5),size(values,6)) )
              values_int4 = int(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int4, &
                                            start, count, stride, map )
              
              deallocate( values_int4 )
          
            case ( DIOS_FLOAT )
              allocate( values_real4(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5),size(values,6)) )
              values_real4 = real(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real4, &
                                            start, count, stride, map )
              
              deallocate( values_real4 )
          
            case ( DIOS_DOUBLE )
              allocate( values_real8(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5),size(values,6)) )
              values_real8 = real(values,kind=8)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real8, &
                                            start, count, stride, map )
              
              deallocate( values_real8 )
          
            case default
              write(*, '("not implemented yet for output type : ",i6)') varp%xtype
              stat=1; return
          end select
        
          ! just put; let netcdf library convert the right kind:
          !stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values, &
          !                              start, count, stride, map )
          !
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*, '("unsupported filetype : ",i6)') filep%ftype
          stat=1; return
      end select
        
    ! ok
    stat = 0

  end subroutine DIOS_Put_Var_r4_6d


  ! ***
  

  subroutine DIOS_Get_Var_r4_6d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Get_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    real(4), intent(out)           ::  values(:,:,:,:,:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:,:,:,:,:,:)
    integer(2), allocatable         ::  values_int2(:,:,:,:,:,:)
    integer(4), allocatable         ::  values_int4(:,:,:,:,:,:)
    integer(8), allocatable         ::  values_int8(:,:,:,:,:,:)
    real(4), allocatable            ::  values_real4(:,:,:,:,:,:)
    real(8), allocatable            ::  values_real8(:,:,:,:,:,:)
          

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    
    ! pointer to variable structure:
    call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
    

      
      ! select appropriate routine for each type:
      select case ( filep%ftype )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_ASCII )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_NETCDF )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          ! read values, converted automatically:
          stat = NF90_Get_Var( filep%netcdf_id, varp%netcdf_varid, values, start, count, stride, map )
          
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*, '("unsupported filetype : ",i6)') filep%ftype
          stat=1; return
      end select
        
    ! ok
    stat = 0

  end subroutine DIOS_Get_Var_r4_6d
  
  
  ! ***

  
  subroutine DIOS_Put_Var_r4_7d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Put_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    real(4), intent(in)            ::  values(:,:,:,:,:,:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:,:,:,:,:,:,:)
    integer(2), allocatable         ::  values_int2(:,:,:,:,:,:,:)
    integer(4), allocatable         ::  values_int4(:,:,:,:,:,:,:)
    integer(8), allocatable         ::  values_int8(:,:,:,:,:,:,:)
    real(4), allocatable            ::  values_real4(:,:,:,:,:,:,:)
    real(8), allocatable            ::  values_real8(:,:,:,:,:,:,:)
        

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    
    ! pointer to variable structure:
    call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
    

      
      ! select appropriate routine for each type:
      select case ( filep%ftype )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_ASCII )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_NETCDF )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      
          ! test target type:
          ! convert to required kind before entering NF90_Put_Var, 
          ! otherwise segmentation faults on some machines ...
          select case ( varp%xtype )
        
            case ( DIOS_BYTE )
              allocate( values_int1(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5),size(values,6),size(values,7)) )
              values_int1 = int(values,kind=1)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int1, &
                                            start, count, stride, map )
              
              deallocate( values_int1 )
          
            case ( DIOS_SHORT )
              allocate( values_int2(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5),size(values,6),size(values,7)) )
              values_int2 = int(values,kind=2)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int2, &
                                            start, count, stride, map )
              
              deallocate( values_int2 )
          
            case ( DIOS_INT )
              allocate( values_int4(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5),size(values,6),size(values,7)) )
              values_int4 = int(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int4, &
                                            start, count, stride, map )
              
              deallocate( values_int4 )
          
            case ( DIOS_FLOAT )
              allocate( values_real4(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5),size(values,6),size(values,7)) )
              values_real4 = real(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real4, &
                                            start, count, stride, map )
              
              deallocate( values_real4 )
          
            case ( DIOS_DOUBLE )
              allocate( values_real8(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5),size(values,6),size(values,7)) )
              values_real8 = real(values,kind=8)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real8, &
                                            start, count, stride, map )
              
              deallocate( values_real8 )
          
            case default
              write(*, '("not implemented yet for output type : ",i6)') varp%xtype
              stat=1; return
          end select
        
          ! just put; let netcdf library convert the right kind:
          !stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values, &
          !                              start, count, stride, map )
          !
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*, '("unsupported filetype : ",i6)') filep%ftype
          stat=1; return
      end select
        
    ! ok
    stat = 0

  end subroutine DIOS_Put_Var_r4_7d


  ! ***
  

  subroutine DIOS_Get_Var_r4_7d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Get_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    real(4), intent(out)           ::  values(:,:,:,:,:,:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:,:,:,:,:,:,:)
    integer(2), allocatable         ::  values_int2(:,:,:,:,:,:,:)
    integer(4), allocatable         ::  values_int4(:,:,:,:,:,:,:)
    integer(8), allocatable         ::  values_int8(:,:,:,:,:,:,:)
    real(4), allocatable            ::  values_real4(:,:,:,:,:,:,:)
    real(8), allocatable            ::  values_real8(:,:,:,:,:,:,:)
          

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    
    ! pointer to variable structure:
    call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
    

      
      ! select appropriate routine for each type:
      select case ( filep%ftype )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_ASCII )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_NETCDF )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          ! read values, converted automatically:
          stat = NF90_Get_Var( filep%netcdf_id, varp%netcdf_varid, values, start, count, stride, map )
          
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*, '("unsupported filetype : ",i6)') filep%ftype
          stat=1; return
      end select
        
    ! ok
    stat = 0

  end subroutine DIOS_Get_Var_r4_7d
  
  
  ! ********************************************************************
  ! ***
  ! *** put_var_r8
  ! ***
  ! ********************************************************************

  
  subroutine DIOS_Put_var_r8_1d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Put_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    real(8), intent(in)            ::  values(:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:)
    integer(2), allocatable         ::  values_int2(:)
    integer(4), allocatable         ::  values_int4(:)
    integer(8), allocatable         ::  values_int8(:)
    real(4), allocatable            ::  values_real4(:)
    real(8), allocatable            ::  values_real8(:)
        

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    
    ! pointer to variable structure:
    call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
    

      
      ! select appropriate routine for each type:
      select case ( filep%ftype )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_ASCII )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_NETCDF )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      
          ! test target type:
          ! convert to required kind before entering NF90_Put_Var, 
          ! otherwise segmentation faults on some machines ...
          select case ( varp%xtype )
        
            case ( DIOS_BYTE )
              allocate( values_int1(size(values,1)) )
              values_int1 = int(values,kind=1)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int1, &
                                            start, count, stride, map )
              
              deallocate( values_int1 )
          
            case ( DIOS_SHORT )
              allocate( values_int2(size(values,1)) )
              values_int2 = int(values,kind=2)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int2, &
                                            start, count, stride, map )
              
              deallocate( values_int2 )
          
            case ( DIOS_INT )
              allocate( values_int4(size(values,1)) )
              values_int4 = int(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int4, &
                                            start, count, stride, map )
              
              deallocate( values_int4 )
          
            case ( DIOS_FLOAT )
              allocate( values_real4(size(values,1)) )
              values_real4 = real(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real4, &
                                            start, count, stride, map )
              
              deallocate( values_real4 )
          
            case ( DIOS_DOUBLE )
              allocate( values_real8(size(values,1)) )
              values_real8 = real(values,kind=8)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real8, &
                                            start, count, stride, map )
              
              deallocate( values_real8 )
          
            case default
              write(*, '("not implemented yet for output type : ",i6)') varp%xtype
              stat=1; return
          end select
        
          ! just put; let netcdf library convert the right kind:
          !stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values, &
          !                              start, count, stride, map )
          !
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*, '("unsupported filetype : ",i6)') filep%ftype
          stat=1; return
      end select
        
    ! ok
    stat = 0

  end subroutine DIOS_Put_var_r8_1d


  ! ***
  

  subroutine DIOS_Get_var_r8_1d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Get_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    real(8), intent(out)           ::  values(:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:)
    integer(2), allocatable         ::  values_int2(:)
    integer(4), allocatable         ::  values_int4(:)
    integer(8), allocatable         ::  values_int8(:)
    real(4), allocatable            ::  values_real4(:)
    real(8), allocatable            ::  values_real8(:)
          

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    
    ! pointer to variable structure:
    call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
    

      
      ! select appropriate routine for each type:
      select case ( filep%ftype )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_ASCII )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_NETCDF )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          ! read values, converted automatically:
          stat = NF90_Get_Var( filep%netcdf_id, varp%netcdf_varid, values, start, count, stride, map )
          
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*, '("unsupported filetype : ",i6)') filep%ftype
          stat=1; return
      end select
        
    ! ok
    stat = 0

  end subroutine DIOS_Get_var_r8_1d
  
  
  ! ***

  
  subroutine DIOS_Put_var_r8_2d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Put_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    real(8), intent(in)            ::  values(:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:,:)
    integer(2), allocatable         ::  values_int2(:,:)
    integer(4), allocatable         ::  values_int4(:,:)
    integer(8), allocatable         ::  values_int8(:,:)
    real(4), allocatable            ::  values_real4(:,:)
    real(8), allocatable            ::  values_real8(:,:)
        

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    
    ! pointer to variable structure:
    call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
    

      
      ! select appropriate routine for each type:
      select case ( filep%ftype )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_ASCII )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_NETCDF )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      
          ! test target type:
          ! convert to required kind before entering NF90_Put_Var, 
          ! otherwise segmentation faults on some machines ...
          select case ( varp%xtype )
        
            case ( DIOS_BYTE )
              allocate( values_int1(size(values,1),size(values,2)) )
              values_int1 = int(values,kind=1)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int1, &
                                            start, count, stride, map )
              
              deallocate( values_int1 )
          
            case ( DIOS_SHORT )
              allocate( values_int2(size(values,1),size(values,2)) )
              values_int2 = int(values,kind=2)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int2, &
                                            start, count, stride, map )
              
              deallocate( values_int2 )
          
            case ( DIOS_INT )
              allocate( values_int4(size(values,1),size(values,2)) )
              values_int4 = int(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int4, &
                                            start, count, stride, map )
              
              deallocate( values_int4 )
          
            case ( DIOS_FLOAT )
              allocate( values_real4(size(values,1),size(values,2)) )
              values_real4 = real(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real4, &
                                            start, count, stride, map )
              
              deallocate( values_real4 )
          
            case ( DIOS_DOUBLE )
              allocate( values_real8(size(values,1),size(values,2)) )
              values_real8 = real(values,kind=8)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real8, &
                                            start, count, stride, map )
              
              deallocate( values_real8 )
          
            case default
              write(*, '("not implemented yet for output type : ",i6)') varp%xtype
              stat=1; return
          end select
        
          ! just put; let netcdf library convert the right kind:
          !stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values, &
          !                              start, count, stride, map )
          !
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*, '("unsupported filetype : ",i6)') filep%ftype
          stat=1; return
      end select
        
    ! ok
    stat = 0

  end subroutine DIOS_Put_var_r8_2d


  ! ***
  

  subroutine DIOS_Get_var_r8_2d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Get_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    real(8), intent(out)           ::  values(:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:,:)
    integer(2), allocatable         ::  values_int2(:,:)
    integer(4), allocatable         ::  values_int4(:,:)
    integer(8), allocatable         ::  values_int8(:,:)
    real(4), allocatable            ::  values_real4(:,:)
    real(8), allocatable            ::  values_real8(:,:)
          

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    
    ! pointer to variable structure:
    call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
    

      
      ! select appropriate routine for each type:
      select case ( filep%ftype )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_ASCII )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_NETCDF )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          ! read values, converted automatically:
          stat = NF90_Get_Var( filep%netcdf_id, varp%netcdf_varid, values, start, count, stride, map )
          
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*, '("unsupported filetype : ",i6)') filep%ftype
          stat=1; return
      end select
        
    ! ok
    stat = 0

  end subroutine DIOS_Get_var_r8_2d
  
  
  ! ***

  
  subroutine DIOS_Put_var_r8_3d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Put_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    real(8), intent(in)            ::  values(:,:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:,:,:)
    integer(2), allocatable         ::  values_int2(:,:,:)
    integer(4), allocatable         ::  values_int4(:,:,:)
    integer(8), allocatable         ::  values_int8(:,:,:)
    real(4), allocatable            ::  values_real4(:,:,:)
    real(8), allocatable            ::  values_real8(:,:,:)
        

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    
    ! pointer to variable structure:
    call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
    

      
      ! select appropriate routine for each type:
      select case ( filep%ftype )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_ASCII )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_NETCDF )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      
          ! test target type:
          ! convert to required kind before entering NF90_Put_Var, 
          ! otherwise segmentation faults on some machines ...
          select case ( varp%xtype )
        
            case ( DIOS_BYTE )
              allocate( values_int1(size(values,1),size(values,2),size(values,3)) )
              values_int1 = int(values,kind=1)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int1, &
                                            start, count, stride, map )
              
              deallocate( values_int1 )
          
            case ( DIOS_SHORT )
              allocate( values_int2(size(values,1),size(values,2),size(values,3)) )
              values_int2 = int(values,kind=2)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int2, &
                                            start, count, stride, map )
              
              deallocate( values_int2 )
          
            case ( DIOS_INT )
              allocate( values_int4(size(values,1),size(values,2),size(values,3)) )
              values_int4 = int(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int4, &
                                            start, count, stride, map )
              
              deallocate( values_int4 )
          
            case ( DIOS_FLOAT )
              allocate( values_real4(size(values,1),size(values,2),size(values,3)) )
              values_real4 = real(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real4, &
                                            start, count, stride, map )
              
              deallocate( values_real4 )
          
            case ( DIOS_DOUBLE )
              allocate( values_real8(size(values,1),size(values,2),size(values,3)) )
              values_real8 = real(values,kind=8)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real8, &
                                            start, count, stride, map )
              
              deallocate( values_real8 )
          
            case default
              write(*, '("not implemented yet for output type : ",i6)') varp%xtype
              stat=1; return
          end select
        
          ! just put; let netcdf library convert the right kind:
          !stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values, &
          !                              start, count, stride, map )
          !
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*, '("unsupported filetype : ",i6)') filep%ftype
          stat=1; return
      end select
        
    ! ok
    stat = 0

  end subroutine DIOS_Put_var_r8_3d


  ! ***
  

  subroutine DIOS_Get_var_r8_3d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Get_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    real(8), intent(out)           ::  values(:,:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:,:,:)
    integer(2), allocatable         ::  values_int2(:,:,:)
    integer(4), allocatable         ::  values_int4(:,:,:)
    integer(8), allocatable         ::  values_int8(:,:,:)
    real(4), allocatable            ::  values_real4(:,:,:)
    real(8), allocatable            ::  values_real8(:,:,:)
          

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    
    ! pointer to variable structure:
    call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
    

      
      ! select appropriate routine for each type:
      select case ( filep%ftype )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_ASCII )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_NETCDF )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          ! read values, converted automatically:
          stat = NF90_Get_Var( filep%netcdf_id, varp%netcdf_varid, values, start, count, stride, map )
          
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*, '("unsupported filetype : ",i6)') filep%ftype
          stat=1; return
      end select
        
    ! ok
    stat = 0

  end subroutine DIOS_Get_var_r8_3d
  
  
  ! ***

  
  subroutine DIOS_Put_var_r8_4d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Put_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    real(8), intent(in)            ::  values(:,:,:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:,:,:,:)
    integer(2), allocatable         ::  values_int2(:,:,:,:)
    integer(4), allocatable         ::  values_int4(:,:,:,:)
    integer(8), allocatable         ::  values_int8(:,:,:,:)
    real(4), allocatable            ::  values_real4(:,:,:,:)
    real(8), allocatable            ::  values_real8(:,:,:,:)
        

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    
    ! pointer to variable structure:
    call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
    

      
      ! select appropriate routine for each type:
      select case ( filep%ftype )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_ASCII )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_NETCDF )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      
          ! test target type:
          ! convert to required kind before entering NF90_Put_Var, 
          ! otherwise segmentation faults on some machines ...
          select case ( varp%xtype )
        
            case ( DIOS_BYTE )
              allocate( values_int1(size(values,1),size(values,2),size(values,3),size(values,4)) )
              values_int1 = int(values,kind=1)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int1, &
                                            start, count, stride, map )
              
              deallocate( values_int1 )
          
            case ( DIOS_SHORT )
              allocate( values_int2(size(values,1),size(values,2),size(values,3),size(values,4)) )
              values_int2 = int(values,kind=2)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int2, &
                                            start, count, stride, map )
              
              deallocate( values_int2 )
          
            case ( DIOS_INT )
              allocate( values_int4(size(values,1),size(values,2),size(values,3),size(values,4)) )
              values_int4 = int(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int4, &
                                            start, count, stride, map )
              
              deallocate( values_int4 )
          
            case ( DIOS_FLOAT )
              allocate( values_real4(size(values,1),size(values,2),size(values,3),size(values,4)) )
              values_real4 = real(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real4, &
                                            start, count, stride, map )
              
              deallocate( values_real4 )
          
            case ( DIOS_DOUBLE )
              allocate( values_real8(size(values,1),size(values,2),size(values,3),size(values,4)) )
              values_real8 = real(values,kind=8)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real8, &
                                            start, count, stride, map )
              
              deallocate( values_real8 )
          
            case default
              write(*, '("not implemented yet for output type : ",i6)') varp%xtype
              stat=1; return
          end select
        
          ! just put; let netcdf library convert the right kind:
          !stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values, &
          !                              start, count, stride, map )
          !
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*, '("unsupported filetype : ",i6)') filep%ftype
          stat=1; return
      end select
        
    ! ok
    stat = 0

  end subroutine DIOS_Put_var_r8_4d


  ! ***
  

  subroutine DIOS_Get_var_r8_4d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Get_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    real(8), intent(out)           ::  values(:,:,:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:,:,:,:)
    integer(2), allocatable         ::  values_int2(:,:,:,:)
    integer(4), allocatable         ::  values_int4(:,:,:,:)
    integer(8), allocatable         ::  values_int8(:,:,:,:)
    real(4), allocatable            ::  values_real4(:,:,:,:)
    real(8), allocatable            ::  values_real8(:,:,:,:)
          
    ! --- begin --------------------------------------

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    
    ! pointer to variable structure:
    call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
    

      
      ! select appropriate routine for each type:
      select case ( filep%ftype )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_ASCII )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_NETCDF )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          ! read values, converted automatically:
          stat = NF90_Get_Var( filep%netcdf_id, varp%netcdf_varid, values, start, count, stride, map )
          
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*, '("unsupported filetype : ",i6)') filep%ftype
          stat=1; return
      end select
        
    ! ok
    stat = 0

  end subroutine DIOS_Get_var_r8_4d
  
  
  ! ***

  
  subroutine DIOS_Put_var_r8_5d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Put_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    real(8), intent(in)            ::  values(:,:,:,:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp

    integer(1), allocatable         ::  values_int1(:,:,:,:,:)
    integer(2), allocatable         ::  values_int2(:,:,:,:,:)
    integer(4), allocatable         ::  values_int4(:,:,:,:,:)
    integer(8), allocatable         ::  values_int8(:,:,:,:,:)
    real(4), allocatable            ::  values_real4(:,:,:,:,:)
    real(8), allocatable            ::  values_real8(:,:,:,:,:)
        

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    
    ! pointer to variable structure:
    call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
    

      
      ! select appropriate routine for each type:
      select case ( filep%ftype )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_ASCII )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_NETCDF )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      
          ! test target type:
          ! convert to required kind before entering NF90_Put_Var, 
          ! otherwise segmentation faults on some machines ...
          select case ( varp%xtype )
        
            case ( DIOS_BYTE )
              allocate( values_int1(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5)) )
              values_int1 = int(values,kind=1)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int1, &
                                            start, count, stride, map )
              
              deallocate( values_int1 )
          
            case ( DIOS_SHORT )
              allocate( values_int2(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5)) )
              values_int2 = int(values,kind=2)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int2, &
                                            start, count, stride, map )
              
              deallocate( values_int2 )
          
            case ( DIOS_INT )
              allocate( values_int4(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5)) )
              values_int4 = int(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int4, &
                                            start, count, stride, map )
              
              deallocate( values_int4 )
          
            case ( DIOS_FLOAT )
              allocate( values_real4(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5)) )
              values_real4 = real(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real4, &
                                            start, count, stride, map )
              
              deallocate( values_real4 )
          
            case ( DIOS_DOUBLE )
              allocate( values_real8(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5)) )
              values_real8 = real(values,kind=8)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real8, &
                                            start, count, stride, map )
              
              deallocate( values_real8 )
          
            case default
              write(*, '("not implemented yet for output type : ",i6)') varp%xtype
              stat=1; return
          end select
        
          ! just put; let netcdf library convert the right kind:
          !stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values, &
          !                              start, count, stride, map )
          !
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*, '("unsupported filetype : ",i6)') filep%ftype
          stat=1; return
      end select
        
    ! ok
    stat = 0

  end subroutine DIOS_Put_var_r8_5d


  ! ***
  

  subroutine DIOS_Get_var_r8_5d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Get_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    real(8), intent(out)           ::  values(:,:,:,:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:,:,:,:,:)
    integer(2), allocatable         ::  values_int2(:,:,:,:,:)
    integer(4), allocatable         ::  values_int4(:,:,:,:,:)
    integer(8), allocatable         ::  values_int8(:,:,:,:,:)
    real(4), allocatable            ::  values_real4(:,:,:,:,:)
    real(8), allocatable            ::  values_real8(:,:,:,:,:)
          

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    
    ! pointer to variable structure:
    call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
    

      
      ! select appropriate routine for each type:
      select case ( filep%ftype )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_ASCII )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_NETCDF )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          ! read values, converted automatically:
          stat = NF90_Get_Var( filep%netcdf_id, varp%netcdf_varid, values, start, count, stride, map )
          
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*, '("unsupported filetype : ",i6)') filep%ftype
          stat=1; return
      end select
        
    ! ok
    stat = 0

  end subroutine DIOS_Get_var_r8_5d
  
  
  ! ***

  
  subroutine DIOS_Put_var_r8_6d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Put_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    real(8), intent(in)            ::  values(:,:,:,:,:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:,:,:,:,:,:)
    integer(2), allocatable         ::  values_int2(:,:,:,:,:,:)
    integer(4), allocatable         ::  values_int4(:,:,:,:,:,:)
    integer(8), allocatable         ::  values_int8(:,:,:,:,:,:)
    real(4), allocatable            ::  values_real4(:,:,:,:,:,:)
    real(8), allocatable            ::  values_real8(:,:,:,:,:,:)
        

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    
    ! pointer to variable structure:
    call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
    

      
      ! select appropriate routine for each type:
      select case ( filep%ftype )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_ASCII )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_NETCDF )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      
          ! test target type:
          ! convert to required kind before entering NF90_Put_Var, 
          ! otherwise segmentation faults on some machines ...
          select case ( varp%xtype )
        
            case ( DIOS_BYTE )
              allocate( values_int1(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5),size(values,6)) )
              values_int1 = int(values,kind=1)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int1, &
                                            start, count, stride, map )
              
              deallocate( values_int1 )
          
            case ( DIOS_SHORT )
              allocate( values_int2(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5),size(values,6)) )
              values_int2 = int(values,kind=2)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int2, &
                                            start, count, stride, map )
              
              deallocate( values_int2 )
          
            case ( DIOS_INT )
              allocate( values_int4(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5),size(values,6)) )
              values_int4 = int(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int4, &
                                            start, count, stride, map )
              
              deallocate( values_int4 )
          
            case ( DIOS_FLOAT )
              allocate( values_real4(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5),size(values,6)) )
              values_real4 = real(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real4, &
                                            start, count, stride, map )
              
              deallocate( values_real4 )
          
            case ( DIOS_DOUBLE )
              allocate( values_real8(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5),size(values,6)) )
              values_real8 = real(values,kind=8)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real8, &
                                            start, count, stride, map )
              
              deallocate( values_real8 )
          
            case default
              write(*, '("not implemented yet for output type : ",i6)') varp%xtype
              stat=1; return
          end select
        
          ! just put; let netcdf library convert the right kind:
          !stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values, &
          !                              start, count, stride, map )
          !
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*, '("unsupported filetype : ",i6)') filep%ftype
          stat=1; return
      end select
        
    ! ok
    stat = 0

  end subroutine DIOS_Put_var_r8_6d


  ! ***
  

  subroutine DIOS_Get_var_r8_6d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Get_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    real(8), intent(out)           ::  values(:,:,:,:,:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:,:,:,:,:,:)
    integer(2), allocatable         ::  values_int2(:,:,:,:,:,:)
    integer(4), allocatable         ::  values_int4(:,:,:,:,:,:)
    integer(8), allocatable         ::  values_int8(:,:,:,:,:,:)
    real(4), allocatable            ::  values_real4(:,:,:,:,:,:)
    real(8), allocatable            ::  values_real8(:,:,:,:,:,:)
          

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    
    ! pointer to variable structure:
    call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
    

      
      ! select appropriate routine for each type:
      select case ( filep%ftype )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_ASCII )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_NETCDF )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          ! read values, converted automatically:
          stat = NF90_Get_Var( filep%netcdf_id, varp%netcdf_varid, values, start, count, stride, map )
          
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*, '("unsupported filetype : ",i6)') filep%ftype
          stat=1; return
      end select
        
    ! ok
    stat = 0

  end subroutine DIOS_Get_var_r8_6d
  
  
  ! ***

  
  subroutine DIOS_Put_var_r8_7d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Put_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    real(8), intent(in)            ::  values(:,:,:,:,:,:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:,:,:,:,:,:,:)
    integer(2), allocatable         ::  values_int2(:,:,:,:,:,:,:)
    integer(4), allocatable         ::  values_int4(:,:,:,:,:,:,:)
    integer(8), allocatable         ::  values_int8(:,:,:,:,:,:,:)
    real(4), allocatable            ::  values_real4(:,:,:,:,:,:,:)
    real(8), allocatable            ::  values_real8(:,:,:,:,:,:,:)
        

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    
    ! pointer to variable structure:
    call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
    

      
      ! select appropriate routine for each type:
      select case ( filep%ftype )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_ASCII )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_NETCDF )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      
          ! test target type:
          ! convert to required kind before entering NF90_Put_Var, 
          ! otherwise segmentation faults on some machines ...
          select case ( varp%xtype )
        
            case ( DIOS_BYTE )
              allocate( values_int1(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5),size(values,6),size(values,7)) )
              values_int1 = int(values,kind=1)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int1, &
                                            start, count, stride, map )
              
              deallocate( values_int1 )
          
            case ( DIOS_SHORT )
              allocate( values_int2(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5),size(values,6),size(values,7)) )
              values_int2 = int(values,kind=2)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int2, &
                                            start, count, stride, map )
              
              deallocate( values_int2 )
          
            case ( DIOS_INT )
              allocate( values_int4(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5),size(values,6),size(values,7)) )
              values_int4 = int(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_int4, &
                                            start, count, stride, map )
              
              deallocate( values_int4 )
          
            case ( DIOS_FLOAT )
              allocate( values_real4(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5),size(values,6),size(values,7)) )
              values_real4 = real(values,kind=4)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real4, &
                                            start, count, stride, map )
              
              deallocate( values_real4 )
          
            case ( DIOS_DOUBLE )
              allocate( values_real8(size(values,1),size(values,2),size(values,3),size(values,4),size(values,5),size(values,6),size(values,7)) )
              values_real8 = real(values,kind=8)
              stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values_real8, &
                                            start, count, stride, map )
              
              deallocate( values_real8 )
          
            case default
              write(*, '("not implemented yet for output type : ",i6)') varp%xtype
              stat=1; return
          end select
        
          ! just put; let netcdf library convert the right kind:
          !stat = NF90_Put_Var( filep%netcdf_id, varp%netcdf_varid, values, &
          !                              start, count, stride, map )
          !
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*, '("unsupported filetype : ",i6)') filep%ftype
          stat=1; return
      end select
        
    ! ok
    stat = 0

  end subroutine DIOS_Put_var_r8_7d


  ! ***
  

  subroutine DIOS_Get_var_r8_7d( hid, varid, values, stat, start, count, stride, map )

    use NetCDF, only : NF90_Get_Var

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    real(8), intent(out)           ::  values(:,:,:,:,:,:,:)
    integer, intent(out)                ::  stat 
    integer, intent(in), optional       ::  start (:)
    integer, intent(in), optional       ::  count (:)
    integer, intent(in), optional       ::  stride(:)
    integer, intent(in), optional       ::  map   (:)
    
    ! --- const --------------------------------------

      
    ! --- external -----------------------------------
    
    
    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp
      
    integer(1), allocatable         ::  values_int1(:,:,:,:,:,:,:)
    integer(2), allocatable         ::  values_int2(:,:,:,:,:,:,:)
    integer(4), allocatable         ::  values_int4(:,:,:,:,:,:,:)
    integer(8), allocatable         ::  values_int8(:,:,:,:,:,:,:)
    real(4), allocatable            ::  values_real4(:,:,:,:,:,:,:)
    real(8), allocatable            ::  values_real8(:,:,:,:,:,:,:)
          

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    
    ! pointer to variable structure:
    call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
    

      
      ! select appropriate routine for each type:
      select case ( filep%ftype )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_ASCII )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_NETCDF )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          ! read values, converted automatically:
          stat = NF90_Get_Var( filep%netcdf_id, varp%netcdf_varid, values, start, count, stride, map )
          
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*, '("unsupported filetype : ",i6)') filep%ftype
          stat=1; return
      end select
        
    ! ok
    stat = 0

  end subroutine DIOS_Get_var_r8_7d



  ! ********************************************************************
  ! ***
  ! *** attributes
  ! ***
  ! ********************************************************************
  

  
  subroutine DIOS_Put_Att_c1_0d( hid, varid, name, values, stat )

    use NetCDF, only : NF90_Put_Att, NF90_GLOBAL

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    character(len=*), intent(in)        ::  name
    character(len=*), intent(in)            ::  values
    integer, intent(out)                ::  stat

    ! --- const --------------------------------------

    
    ! --- external -------------------------------


    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp

    integer                         ::  netcdf_varid

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    ! global or variable attribute ?
    if ( varid == DIOS_GLOBAL ) then
      ! increase counter:
      filep%natt = filep%natt + 1
    else
      ! pointer to variable structure:
      call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
      ! increase counter:
      varp%natt = varp%natt + 1
    end if
   
      ! select appropriate routine for each type:
      select case ( filep%ftype )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_ASCII )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_NETCDF )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          ! set variable id:
          if ( varid == DIOS_GLOBAL ) then
            netcdf_varid = NF90_GLOBAL
          else
            netcdf_varid = varp%netcdf_varid
          end if
          ! write attribute:
          stat = NF90_Put_Att( filep%netcdf_id, netcdf_varid, trim(name), values )
          
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*, '("unsupported filetype : ",i6)') filep%ftype
          stat=1; return
      end select
    
    ! ok
    stat = 0
      
  end subroutine DIOS_Put_Att_c1_0d
      
  
  ! ***  


  subroutine DIOS_Get_Att_c1_0d( hid, varid, name, values, stat )

    use NetCDF, only : NF90_Get_Att, NF90_GLOBAL

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    character(len=*), intent(in)        ::  name
    character(len=*), intent(out)           ::  values
    integer, intent(out)                ::  stat

    ! --- const --------------------------------------

      
    ! --- external -------------------------------


    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp
    integer                         ::  ftype

    integer                         ::  netcdf_varid

    ! --- begin --------------------------------------

    ! single type:
    call DIOS_Get_Type( hid, ftype, stat )

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    ! pointer to variable structure if possible:
    if ( varid /= DIOS_GLOBAL ) then
      call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
    end if
    
    ! select appropriate routine for each type:
    select case ( filep%ftype )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_ASCII )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_NETCDF )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! set variable id:
        if ( varid == DIOS_GLOBAL ) then
          netcdf_varid = NF90_GLOBAL
        else
          netcdf_varid = varp%netcdf_varid
        end if
        ! read attribute:
        stat = NF90_Get_Att( filep%netcdf_id, netcdf_varid, trim(name), values )
        
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        write(*, '("unsupported filetype : ",i6)') filep%ftype
        stat=1; return
    end select

    ! ok
    stat = 0
      
  end subroutine DIOS_Get_Att_c1_0d
      
  
  
      
  
  subroutine DIOS_Put_Att_i1_0d( hid, varid, name, values, stat )

    use NetCDF, only : NF90_Put_Att, NF90_GLOBAL

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    character(len=*), intent(in)        ::  name
    integer(1), intent(in)            ::  values
    integer, intent(out)                ::  stat

    ! --- const --------------------------------------

    
    ! --- external -------------------------------


    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp

    integer                         ::  netcdf_varid

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    ! global or variable attribute ?
    if ( varid == DIOS_GLOBAL ) then
      ! increase counter:
      filep%natt = filep%natt + 1
    else
      ! pointer to variable structure:
      call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
      ! increase counter:
      varp%natt = varp%natt + 1
    end if
   
      ! select appropriate routine for each type:
      select case ( filep%ftype )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_ASCII )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_NETCDF )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          ! set variable id:
          if ( varid == DIOS_GLOBAL ) then
            netcdf_varid = NF90_GLOBAL
          else
            netcdf_varid = varp%netcdf_varid
          end if
          ! write attribute:
          stat = NF90_Put_Att( filep%netcdf_id, netcdf_varid, trim(name), values )
          
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*, '("unsupported filetype : ",i6)') filep%ftype
          stat=1; return
      end select
    
    ! ok
    stat = 0
      
  end subroutine DIOS_Put_Att_i1_0d
      
  
  ! ***  


  subroutine DIOS_Get_Att_i1_0d( hid, varid, name, values, stat )

    use NetCDF, only : NF90_Get_Att, NF90_GLOBAL

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    character(len=*), intent(in)        ::  name
    integer(1), intent(out)           ::  values
    integer, intent(out)                ::  stat

    ! --- const --------------------------------------

    
    ! --- external -------------------------------


    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp
    integer                         ::  ftype
      
    integer                         ::  netcdf_varid

    ! --- begin --------------------------------------

    ! single type:
    call DIOS_Get_Type( hid, ftype, stat )
    

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    ! pointer to variable structure if possible:
    if ( varid /= DIOS_GLOBAL ) then
      call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
    end if
    
    ! select appropriate routine for each type:
    select case ( filep%ftype )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_ASCII )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_NETCDF )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! set variable id:
        if ( varid == DIOS_GLOBAL ) then
          netcdf_varid = NF90_GLOBAL
        else
          netcdf_varid = varp%netcdf_varid
        end if
        ! read attribute:
        stat = NF90_Get_Att( filep%netcdf_id, netcdf_varid, trim(name), values )
        
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        write(*, '("unsupported filetype : ",i6)') filep%ftype
        stat=1; return
    end select

    ! ok
    stat = 0
      
  end subroutine DIOS_Get_Att_i1_0d
      
  
  
      
  
  subroutine DIOS_Put_Att_i1_1d( hid, varid, name, values, stat )

    use NetCDF, only : NF90_Put_Att, NF90_GLOBAL

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    character(len=*), intent(in)        ::  name
    integer(1), intent(in)            ::  values(:)
    integer, intent(out)                ::  stat

    ! --- const --------------------------------------

    
    ! --- external -------------------------------


    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp

    integer                         ::  netcdf_varid

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    ! global or variable attribute ?
    if ( varid == DIOS_GLOBAL ) then
      ! increase counter:
      filep%natt = filep%natt + 1
    else
      ! pointer to variable structure:
      call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
      ! increase counter:
      varp%natt = varp%natt + 1
    end if
   
      ! select appropriate routine for each type:
      select case ( filep%ftype )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_ASCII )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_NETCDF )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          ! set variable id:
          if ( varid == DIOS_GLOBAL ) then
            netcdf_varid = NF90_GLOBAL
          else
            netcdf_varid = varp%netcdf_varid
          end if
          ! write attribute:
          stat = NF90_Put_Att( filep%netcdf_id, netcdf_varid, trim(name), values )
          
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*, '("unsupported filetype : ",i6)') filep%ftype
          stat=1; return
      end select
    
    ! ok
    stat = 0
      
  end subroutine DIOS_Put_Att_i1_1d
      
  
  ! ***  


  subroutine DIOS_Get_Att_i1_1d( hid, varid, name, values, stat )

    use NetCDF, only : NF90_Get_Att, NF90_GLOBAL

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    character(len=*), intent(in)        ::  name
    integer(1), intent(out)           ::  values(:)
    integer, intent(out)                ::  stat

    ! --- const --------------------------------------

    
    ! --- external -------------------------------


    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp
    integer                         ::  ftype

    integer                         ::  netcdf_varid

    ! --- begin --------------------------------------

    ! single type:
    call DIOS_Get_Type( hid, ftype, stat )
    

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    
    ! pointer to variable structure if possible:
    if ( varid /= DIOS_GLOBAL ) then
      call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
      
    end if
    
    ! select appropriate routine for each type:
    select case ( ftype )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_ASCII )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_NETCDF )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! set variable id:
        if ( varid == DIOS_GLOBAL ) then
          netcdf_varid = NF90_GLOBAL
        else
          netcdf_varid = varp%netcdf_varid
        end if
        ! read attribute:
        stat = NF90_Get_Att( filep%netcdf_id, netcdf_varid, trim(name), values )
        
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        write(*, '("unsupported filetype : ",i6)') filep%ftype
        stat=1; return
    end select

    ! ok
    stat = 0
      
  end subroutine DIOS_Get_Att_i1_1d
      
  
  subroutine DIOS_Put_Att_i2_0d( hid, varid, name, values, stat )

    use NetCDF, only : NF90_Put_Att, NF90_GLOBAL

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    character(len=*), intent(in)        ::  name
    integer(2), intent(in)            ::  values
    integer, intent(out)                ::  stat

    ! --- const --------------------------------------

    
    ! --- external -------------------------------


    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp

    integer                         ::  netcdf_varid

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    ! global or variable attribute ?
    if ( varid == DIOS_GLOBAL ) then
      ! increase counter:
      filep%natt = filep%natt + 1
    else
      ! pointer to variable structure:
      call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
      ! increase counter:
      varp%natt = varp%natt + 1
    end if
   
      ! select appropriate routine for each type:
      select case ( filep%ftype )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_ASCII )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_NETCDF )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          ! set variable id:
          if ( varid == DIOS_GLOBAL ) then
            netcdf_varid = NF90_GLOBAL
          else
            netcdf_varid = varp%netcdf_varid
          end if
          ! write attribute:
          stat = NF90_Put_Att( filep%netcdf_id, netcdf_varid, trim(name), values )
          
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*, '("unsupported filetype : ",i6)') filep%ftype
          stat=1; return
      end select
    
    ! ok
    stat = 0
      
  end subroutine DIOS_Put_Att_i2_0d
      
  
  ! ***  


  subroutine DIOS_Get_Att_i2_0d( hid, varid, name, values, stat )

    use NetCDF, only : NF90_Get_Att, NF90_GLOBAL

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    character(len=*), intent(in)        ::  name
    integer(2), intent(out)           ::  values
    integer, intent(out)                ::  stat

    ! --- const --------------------------------------

    
    ! --- external -------------------------------


    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp
    integer                         ::  ftype

    integer                         ::  netcdf_varid

    ! --- begin --------------------------------------

    ! single type:
    call DIOS_Get_Type( hid, ftype, stat )

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    ! pointer to variable structure if possible:
    if ( varid /= DIOS_GLOBAL ) then
      call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
    end if
    
    ! select appropriate routine for each type:
    select case ( ftype )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_ASCII )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_NETCDF )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! set variable id:
        if ( varid == DIOS_GLOBAL ) then
          netcdf_varid = NF90_GLOBAL
        else
          netcdf_varid = varp%netcdf_varid
        end if
        ! read attribute:
        stat = NF90_Get_Att( filep%netcdf_id, netcdf_varid, trim(name), values )
        
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        write(*, '("unsupported filetype : ",i6)') filep%ftype
        stat=1; return
    end select

    ! ok
    stat = 0
      
  end subroutine DIOS_Get_Att_i2_0d
      
  
  
      
  
  subroutine DIOS_Put_Att_i2_1d( hid, varid, name, values, stat )

    use NetCDF, only : NF90_Put_Att, NF90_GLOBAL

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    character(len=*), intent(in)        ::  name
    integer(2), intent(in)            ::  values(:)
    integer, intent(out)                ::  stat

    ! --- const --------------------------------------

    
    ! --- external -------------------------------


    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp

    integer                         ::  netcdf_varid

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    ! global or variable attribute ?
    if ( varid == DIOS_GLOBAL ) then
      ! increase counter:
      filep%natt = filep%natt + 1
    else
      ! pointer to variable structure:
      call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
      ! increase counter:
      varp%natt = varp%natt + 1
    end if
   
      ! select appropriate routine for each type:
      select case ( filep%ftype )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_ASCII )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_NETCDF )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          ! set variable id:
          if ( varid == DIOS_GLOBAL ) then
            netcdf_varid = NF90_GLOBAL
          else
            netcdf_varid = varp%netcdf_varid
          end if
          ! write attribute:
          stat = NF90_Put_Att( filep%netcdf_id, netcdf_varid, trim(name), values )
          
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*, '("unsupported filetype : ",i6)') filep%ftype
          stat=1; return
      end select
    
    ! ok
    stat = 0
      
  end subroutine DIOS_Put_Att_i2_1d
      
  
  ! ***  


  subroutine DIOS_Get_Att_i2_1d( hid, varid, name, values, stat )

    use NetCDF, only : NF90_Get_Att, NF90_GLOBAL

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    character(len=*), intent(in)        ::  name
    integer(2), intent(out)           ::  values(:)
    integer, intent(out)                ::  stat

    ! --- const --------------------------------------

    
    ! --- external -------------------------------


    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp
    integer                         ::  ftype

    integer                         ::  netcdf_varid

    ! --- begin --------------------------------------

    ! single type:
    call DIOS_Get_Type( hid, ftype, stat )

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    ! pointer to variable structure if possible:
    if ( varid /= DIOS_GLOBAL ) then
      call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
    end if
    
    ! select appropriate routine for each type:
    select case ( ftype )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_ASCII )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_NETCDF )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! set variable id:
        if ( varid == DIOS_GLOBAL ) then
          netcdf_varid = NF90_GLOBAL
        else
          netcdf_varid = varp%netcdf_varid
        end if
        ! read attribute:
        stat = NF90_Get_Att( filep%netcdf_id, netcdf_varid, trim(name), values )
        
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        write(*, '("unsupported filetype : ",i6)') filep%ftype
        stat=1; return
    end select

    ! ok
    stat = 0
      
  end subroutine DIOS_Get_Att_i2_1d
      
  
  
      
  
  subroutine DIOS_Put_Att_i4_0d( hid, varid, name, values, stat )

    use NetCDF, only : NF90_Put_Att, NF90_GLOBAL

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    character(len=*), intent(in)        ::  name
    integer(4), intent(in)            ::  values
    integer, intent(out)                ::  stat

    ! --- const --------------------------------------

    
    ! --- external -------------------------------


    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp

    integer                         ::  netcdf_varid

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    ! global or variable attribute ?
    if ( varid == DIOS_GLOBAL ) then
      ! increase counter:
      filep%natt = filep%natt + 1
    else
      ! pointer to variable structure:
      call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
      ! increase counter:
      varp%natt = varp%natt + 1
    end if
   
      ! select appropriate routine for each type:
      select case ( filep%ftype )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_ASCII )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_NETCDF )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          ! set variable id:
          if ( varid == DIOS_GLOBAL ) then
            netcdf_varid = NF90_GLOBAL
          else
            netcdf_varid = varp%netcdf_varid
          end if
          ! write attribute:
          stat = NF90_Put_Att( filep%netcdf_id, netcdf_varid, trim(name), values )
          
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*, '("unsupported filetype : ",i6)') filep%ftype
          stat=1; return
      end select
    
    ! ok
    stat = 0
      
  end subroutine DIOS_Put_Att_i4_0d
      
  
  ! ***  


  subroutine DIOS_Get_Att_i4_0d( hid, varid, name, values, stat )

    use NetCDF, only : NF90_Get_Att, NF90_GLOBAL

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    character(len=*), intent(in)        ::  name
    integer(4), intent(out)           ::  values
    integer, intent(out)                ::  stat

    ! --- const --------------------------------------

    
    ! --- external -------------------------------


    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp
    integer                         ::  ftype

    integer                         ::  netcdf_varid

    ! --- begin --------------------------------------

    ! single type:
    call DIOS_Get_Type( hid, ftype, stat )

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    ! pointer to variable structure if possible:
    if ( varid /= DIOS_GLOBAL ) then
      call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
    end if
    
    ! select appropriate routine for each type:
    select case ( ftype )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_ASCII )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_NETCDF )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! set variable id:
        if ( varid == DIOS_GLOBAL ) then
          netcdf_varid = NF90_GLOBAL
        else
          netcdf_varid = varp%netcdf_varid
        end if
        ! read attribute:
        stat = NF90_Get_Att( filep%netcdf_id, netcdf_varid, trim(name), values )
        
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        write(*, '("unsupported filetype : ",i6)') filep%ftype
        stat=1; return
    end select

    ! ok
    stat = 0
      
  end subroutine DIOS_Get_Att_i4_0d
      
  
  
      
  
  subroutine DIOS_Put_Att_i4_1d( hid, varid, name, values, stat )

    use NetCDF, only : NF90_Put_Att, NF90_GLOBAL

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    character(len=*), intent(in)        ::  name
    integer(4), intent(in)            ::  values(:)
    integer, intent(out)                ::  stat

    ! --- const --------------------------------------

    
    ! --- external -------------------------------


    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp

    integer                         ::  netcdf_varid

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    ! global or variable attribute ?
    if ( varid == DIOS_GLOBAL ) then
      ! increase counter:
      filep%natt = filep%natt + 1
    else
      ! pointer to variable structure:
      call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
      ! increase counter:
      varp%natt = varp%natt + 1
    end if
   
      ! select appropriate routine for each type:
      select case ( filep%ftype )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_ASCII )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_NETCDF )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          ! set variable id:
          if ( varid == DIOS_GLOBAL ) then
            netcdf_varid = NF90_GLOBAL
          else
            netcdf_varid = varp%netcdf_varid
          end if
          ! write attribute:
          stat = NF90_Put_Att( filep%netcdf_id, netcdf_varid, trim(name), values )
          
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*, '("unsupported filetype : ",i6)') filep%ftype
          stat=1; return
      end select
    
    ! ok
    stat = 0
      
  end subroutine DIOS_Put_Att_i4_1d
      
  
  ! ***  


  subroutine DIOS_Get_Att_i4_1d( hid, varid, name, values, stat )

    use NetCDF, only : NF90_Get_Att, NF90_GLOBAL

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    character(len=*), intent(in)        ::  name
    integer(4), intent(out)           ::  values(:)
    integer, intent(out)                ::  stat

    ! --- const --------------------------------------

    
    ! --- external -------------------------------


    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp
    integer                         ::  ftype

    integer                         ::  netcdf_varid

    ! --- begin --------------------------------------

    ! single type:
    call DIOS_Get_Type( hid, ftype, stat )

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    ! pointer to variable structure if possible:
    if ( varid /= DIOS_GLOBAL ) then
      call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
    end if
    
    ! select appropriate routine for each type:
    select case ( ftype )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_ASCII )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_NETCDF )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! set variable id:
        if ( varid == DIOS_GLOBAL ) then
          netcdf_varid = NF90_GLOBAL
        else
          netcdf_varid = varp%netcdf_varid
        end if
        ! read attribute:
        stat = NF90_Get_Att( filep%netcdf_id, netcdf_varid, trim(name), values )
        
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        write(*, '("unsupported filetype : ",i6)') filep%ftype
        stat=1; return
    end select

    ! ok
    stat = 0
      
  end subroutine DIOS_Get_Att_i4_1d
      
  
  
      
  
  subroutine DIOS_Put_Att_r4_0d( hid, varid, name, values, stat )

    use NetCDF, only : NF90_Put_Att, NF90_GLOBAL

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    character(len=*), intent(in)        ::  name
    real(4), intent(in)            ::  values
    integer, intent(out)                ::  stat

    ! --- const --------------------------------------

    
    ! --- external -------------------------------


    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp

    integer                         ::  netcdf_varid

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    ! global or variable attribute ?
    if ( varid == DIOS_GLOBAL ) then
      ! increase counter:
      filep%natt = filep%natt + 1
    else
      ! pointer to variable structure:
      call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
      ! increase counter:
      varp%natt = varp%natt + 1
    end if
   
      ! select appropriate routine for each type:
      select case ( filep%ftype )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_ASCII )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_NETCDF )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          ! set variable id:
          if ( varid == DIOS_GLOBAL ) then
            netcdf_varid = NF90_GLOBAL
          else
            netcdf_varid = varp%netcdf_varid
          end if
          ! write attribute:
          stat = NF90_Put_Att( filep%netcdf_id, netcdf_varid, trim(name), values )
          
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*, '("unsupported filetype : ",i6)') filep%ftype
          stat=1; return
      end select
    
    ! ok
    stat = 0
      
  end subroutine DIOS_Put_Att_r4_0d
      
  
  ! ***  


  subroutine DIOS_Get_Att_r4_0d( hid, varid, name, values, stat )

    use NetCDF, only : NF90_Get_Att, NF90_GLOBAL

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    character(len=*), intent(in)        ::  name
    real(4), intent(out)           ::  values
    integer, intent(out)                ::  stat

    ! --- const --------------------------------------

    
    ! --- external -------------------------------


    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp
    integer                         ::  ftype

    integer                         ::  netcdf_varid

    ! --- begin --------------------------------------

    ! single type:
    call DIOS_Get_Type( hid, ftype, stat )

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    ! pointer to variable structure if possible:
    if ( varid /= DIOS_GLOBAL ) then
      call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
    end if
    
    ! select appropriate routine for each type:
    select case ( ftype )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_ASCII )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_NETCDF )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! set variable id:
        if ( varid == DIOS_GLOBAL ) then
          netcdf_varid = NF90_GLOBAL
        else
          netcdf_varid = varp%netcdf_varid
        end if
        ! read attribute:
        stat = NF90_Get_Att( filep%netcdf_id, netcdf_varid, trim(name), values )
        
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        write(*, '("unsupported filetype : ",i6)') filep%ftype
        stat=1; return
    end select

    ! ok
    stat = 0
      
  end subroutine DIOS_Get_Att_r4_0d
      
  
  
      
  
  subroutine DIOS_Put_Att_r4_1d( hid, varid, name, values, stat )

    use NetCDF, only : NF90_Put_Att, NF90_GLOBAL

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    character(len=*), intent(in)        ::  name
    real(4), intent(in)            ::  values(:)
    integer, intent(out)                ::  stat

    ! --- const --------------------------------------

    
    ! --- external -------------------------------


    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp

    integer                         ::  netcdf_varid

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    
    ! global or variable attribute ?
    if ( varid == DIOS_GLOBAL ) then
      ! increase counter:
      filep%natt = filep%natt + 1
    else
      ! pointer to variable structure:
      call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
      
      ! increase counter:
      varp%natt = varp%natt + 1
    end if
   
      ! select appropriate routine for each type:
      select case ( filep%ftype )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_ASCII )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_NETCDF )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          ! set variable id:
          if ( varid == DIOS_GLOBAL ) then
            netcdf_varid = NF90_GLOBAL
          else
            netcdf_varid = varp%netcdf_varid
          end if
          ! write attribute:
          stat = NF90_Put_Att( filep%netcdf_id, netcdf_varid, trim(name), values )
          
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*, '("unsupported filetype : ",i6)') filep%ftype
          stat=1; return
      end select
    
    ! ok
    stat = 0
      
  end subroutine DIOS_Put_Att_r4_1d
      
  
  ! ***  


  subroutine DIOS_Get_Att_r4_1d( hid, varid, name, values, stat )

    use NetCDF, only : NF90_Get_Att, NF90_GLOBAL

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    character(len=*), intent(in)        ::  name
    real(4), intent(out)           ::  values(:)
    integer, intent(out)                ::  stat

    ! --- const --------------------------------------

    
    ! --- external -------------------------------


    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp
    integer                         ::  ftype

    integer                         ::  netcdf_varid

    ! --- begin --------------------------------------

    ! single type:
    call DIOS_Get_Type( hid, ftype, stat )

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    ! pointer to variable structure if possible:
    if ( varid /= DIOS_GLOBAL ) then
      call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
    end if
    
    ! select appropriate routine for each type:
    select case ( ftype )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_ASCII )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_NETCDF )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! set variable id:
        if ( varid == DIOS_GLOBAL ) then
          netcdf_varid = NF90_GLOBAL
        else
          netcdf_varid = varp%netcdf_varid
        end if
        ! read attribute:
        stat = NF90_Get_Att( filep%netcdf_id, netcdf_varid, trim(name), values )
        
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        write(*, '("unsupported filetype : ",i6)') filep%ftype
        stat=1; return
    end select

    ! ok
    stat = 0
      
  end subroutine DIOS_Get_Att_r4_1d
      
  
  
      
  
  subroutine DIOS_Put_Att_r8_0d( hid, varid, name, values, stat )

    use NetCDF, only : NF90_Put_Att, NF90_GLOBAL

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    character(len=*), intent(in)        ::  name
    real(8), intent(in)            ::  values
    integer, intent(out)                ::  stat

    ! --- const --------------------------------------

    
    ! --- external -------------------------------


    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp

    integer                         ::  netcdf_varid

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    
    ! global or variable attribute ?
    if ( varid == DIOS_GLOBAL ) then
      ! increase counter:
      filep%natt = filep%natt + 1
    else
      ! pointer to variable structure:
      call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
      
      ! increase counter:
      varp%natt = varp%natt + 1
    end if
   
      ! select appropriate routine for each type:
      select case ( filep%ftype )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_ASCII )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_NETCDF )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          ! set variable id:
          if ( varid == DIOS_GLOBAL ) then
            netcdf_varid = NF90_GLOBAL
          else
            netcdf_varid = varp%netcdf_varid
          end if
          ! write attribute:
          stat = NF90_Put_Att( filep%netcdf_id, netcdf_varid, trim(name), values )
          
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*, '("unsupported filetype : ",i6)') filep%ftype
          stat=1; return
      end select
    
    ! ok
    stat = 0
      
  end subroutine DIOS_Put_Att_r8_0d
      
  
  ! ***  


  subroutine DIOS_Get_Att_r8_0d( hid, varid, name, values, stat )

    use NetCDF, only : NF90_Get_Att, NF90_GLOBAL

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    character(len=*), intent(in)        ::  name
    real(8), intent(out)           ::  values
    integer, intent(out)                ::  stat

    ! --- const --------------------------------------

    
    ! --- external -------------------------------


    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp
    integer                         ::  ftype

    integer                         ::  netcdf_varid

    ! --- begin --------------------------------------

    ! single type:
    call DIOS_Get_Type( hid, ftype, stat )

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    ! pointer to variable structure if possible:
    if ( varid /= DIOS_GLOBAL ) then
      call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
    end if
    
    ! select appropriate routine for each type:
    select case ( ftype )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_ASCII )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_NETCDF )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! set variable id:
        if ( varid == DIOS_GLOBAL ) then
          netcdf_varid = NF90_GLOBAL
        else
          netcdf_varid = varp%netcdf_varid
        end if
        ! read attribute:
        stat = NF90_Get_Att( filep%netcdf_id, netcdf_varid, trim(name), values )
        
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        write(*, '("unsupported filetype : ",i6)') filep%ftype
        stat=1; return
    end select

    ! ok
    stat = 0
      
  end subroutine DIOS_Get_Att_r8_0d
      
  
  
      
  
  subroutine DIOS_Put_Att_r8_1d( hid, varid, name, values, stat )

    use NetCDF, only : NF90_Put_Att, NF90_GLOBAL

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    character(len=*), intent(in)        ::  name
    real(8), intent(in)            ::  values(:)
    integer, intent(out)                ::  stat

    ! --- const --------------------------------------

    
    ! --- external -------------------------------


    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp

    integer                         ::  netcdf_varid

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    
    ! global or variable attribute ?
    if ( varid == DIOS_GLOBAL ) then
      ! increase counter:
      filep%natt = filep%natt + 1
    else
      ! pointer to variable structure:
      call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
      
      ! increase counter:
      varp%natt = varp%natt + 1
    end if
   
      ! select appropriate routine for each type:
      select case ( filep%ftype )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_ASCII )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_NETCDF )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          ! set variable id:
          if ( varid == DIOS_GLOBAL ) then
            netcdf_varid = NF90_GLOBAL
          else
            netcdf_varid = varp%netcdf_varid
          end if
          ! write attribute:
          stat = NF90_Put_Att( filep%netcdf_id, netcdf_varid, trim(name), values )
          
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*, '("unsupported filetype : ",i6)') filep%ftype
          stat=1; return
      end select
    
    ! ok
    stat = 0
      
  end subroutine DIOS_Put_Att_r8_1d
      
  
  ! ***  


  subroutine DIOS_Get_Att_r8_1d( hid, varid, name, values, stat )

    use NetCDF, only : NF90_Get_Att, NF90_GLOBAL

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  varid
    character(len=*), intent(in)        ::  name
    real(8), intent(out)           ::  values(:)
    integer, intent(out)                ::  stat

    ! --- const --------------------------------------

    
    ! --- external -------------------------------


    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp
    integer                         ::  ftype

    integer                         ::  netcdf_varid

    ! --- begin --------------------------------------

    ! single type:
    call DIOS_Get_Type( hid, ftype, stat )
    

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    
    ! pointer to variable structure if possible:
    if ( varid /= DIOS_GLOBAL ) then
      call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
      
    end if
    
    ! select appropriate routine for each type:
    select case ( ftype )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_ASCII )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_NETCDF )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! set variable id:
        if ( varid == DIOS_GLOBAL ) then
          netcdf_varid = NF90_GLOBAL
        else
          netcdf_varid = varp%netcdf_varid
        end if
        ! read attribute:
        stat = NF90_Get_Att( filep%netcdf_id, netcdf_varid, trim(name), values )
        
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        write(*, '("unsupported filetype : ",i6)') filep%ftype
        stat=1; return
    end select

    ! ok
    stat = 0
      
  end subroutine DIOS_Get_Att_r8_1d
      
  
  
      

  ! ********************************************************************
  ! ***
  ! *** inquire
  ! ***
  ! ********************************************************************


  subroutine DIOS_Get_Type( hid, ftype, stat )
  
    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(out)                ::  ftype
    integer, intent(out)                ::  stat

    ! --- const --------------------------------------

    
    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    
    ! return single type:
    ftype = filep%ftype
    
    ! ok
    stat = 0

  end subroutine DIOS_Get_Type
  
  
  ! ***
  
  
  subroutine DIOS_Inquire( hid, stat, nDimensions, nVariables, nAttributes )
  
    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(out)                ::  stat
    
    integer, intent(out), optional      ::  nDimensions
    integer, intent(out), optional      ::  nVariables
    integer, intent(out), optional      ::  nAttributes

    ! --- const --------------------------------------

    
    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    ! return number of dimensions ?
    if ( present(nDimensions) ) then
      ! get number of elements in list:
      call DIOS_Dim_List_Inquire( filep%Dim_List, stat, n=nDimensions )
    end if
    
    ! return number of variables ?
    if ( present(nVariables) ) then
      ! get number of elements in list:
      call DIOS_Var_List_Inquire( filep%Var_List, stat, n=nVariables )
    end if
    
    ! return number of global attributes ?
    if ( present(nAttributes) ) then
      ! copy from structure:
      nAttributes = filep%natt
    end if
    
    ! ok
    stat = 0
    
  end subroutine DIOS_Inquire
  
  
  ! ***
  
  
  subroutine DIOS_Inq_DimID( hid, name, dimid, stat )
  
    ! --- in/out -------------------------------------

    integer, intent(in)                           ::  hid
    character(len=*), intent(in)                  ::  name
    integer, intent(out)                          ::  dimid
    integer, intent(out)                          ::  stat

    ! --- const --------------------------------------

    
    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    integer                         ::  ndim, idim
    character(len=LEN_NAME)         ::  dimname

    ! --- begin --------------------------------------
    
    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )

    ! dummy ...
    dimid = -1

    ! number of variables:
    call DIOS_Inquire( hid, stat, ndimensions=ndim )
    ! list variables ?
    if ( ndim > 0 ) then
      ! loop over variables:
      do idim = 1, ndim
        ! get name:
        call DIOS_Inquire_Dimension( hid, idim, stat, name=dimname )
        ! similar ?
        if ( trim(name) == trim(dimname) ) then
          ! store current id:
          dimid = idim
          ! leave:
          exit
        end if
      end do  ! variables
    end if
    
    ! check ...
    if ( dimid < 0 ) then
      write(*, '("no dimension `",a,"` found in file : ",a)') trim(name), trim(filep%filename)
      stat=1; return
    end if
    
    ! ok
    stat = 0
    
  end subroutine DIOS_Inq_DimID


  ! ***
  
  
  subroutine DIOS_Inquire_Dimension( hid, dimid, stat, name, length, unlimited, named )
  
    ! --- in/out -------------------------------------

    integer, intent(in)                           ::  hid
    integer, intent(in)                           ::  dimid
    integer, intent(out)                          ::  stat
    
    character(len=*), intent(out), optional       ::  name
    integer, intent(out), optional                ::  length
    logical, intent(out), optional                ::  unlimited
    logical, intent(out), optional                ::  named

    ! --- const --------------------------------------

    
    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    type(DIOS_Dim), pointer          ::  dimp

    ! --- begin --------------------------------------

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )

    ! pointer to dimension structure:
    call DIOS_Dim_List_Get_Pointer( filep%Dim_List, dimid, dimp, stat )
    
    ! return value ?
    if ( present(name     ) ) name      = trim(dimp%name)
    if ( present(length   ) ) length    = dimp%length
    if ( present(unlimited) ) unlimited = dimp%unlimited
    if ( present(named    ) ) named     = dimp%named
    
    ! ok
    stat = 0
    
  end subroutine DIOS_Inquire_Dimension
  
  
  ! ***
  
  
  subroutine DIOS_Inq_VarID( hid, name, varid, stat )
  
    ! --- in/out -------------------------------------

    integer, intent(in)                           ::  hid
    character(len=*), intent(in)                  ::  name
    integer, intent(out)                          ::  varid
    integer, intent(out)                          ::  stat

    ! --- const --------------------------------------

    
    ! --- local --------------------------------------

    type(DIOS_File), pointer         ::  filep
    integer                         ::  nvar, ivar
    character(len=LEN_NAME)         ::  varname

    ! --- begin --------------------------------------
    
    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )

    ! dummy ...
    varid = -1

    ! number of variables:
    call DIOS_Inquire( hid, stat, nVariables=nvar )
    ! list variables ?
    if ( nvar > 0 ) then
      ! loop over variables:
      do ivar = 1, nvar
        ! get name:
        call DIOS_Inquire_Variable( hid, ivar, stat, name=varname )
        ! similar ?
        if ( trim(name) == trim(varname) ) then
          ! store current id:
          varid = ivar
          ! leave:
          exit
        end if
      end do  ! variables
    end if
    
    ! check ...
    if ( varid < 0 ) then
      write(*, '("no variable `",a,"` found in file : ",a)') trim(name), trim(filep%filename)
      stat=varid; return
    end if
    
    ! ok
    stat = 0
    
  end subroutine DIOS_Inq_VarID
  
  
  ! ***
  
  
  subroutine DIOS_Inquire_Variable( hid, varid, stat, &
                                     name, xtype, ndims, dimids, natts, &
                                     shp )
  
    ! --- in/out -------------------------------------

    integer, intent(in)                           ::  hid
    integer, intent(in)                           ::  varid
    integer, intent(out)                          ::  stat
    
    character(len=*), intent(out), optional       ::  name
    integer, intent(out), optional                ::  xtype
    integer, intent(out), optional                ::  ndims
    integer, intent(out), optional                ::  dimids(:)
    integer, intent(out), optional                ::  natts
    integer, intent(out), optional                ::  shp(:)

    ! --- const --------------------------------------

    
    ! --- local --------------------------------------

    integer                         ::  ftype
    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp

    ! --- begin --------------------------------------
    
    ! single type:
    call DIOS_Get_Type( hid, ftype, stat )

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )

    ! pointer to variable structure:
    call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
    
    ! return value ?
    if ( present(name     ) ) name  = trim(varp%name)
    if ( present(xtype    ) ) xtype = varp%xtype
    if ( present(ndims    ) ) ndims = varp%ndim
    if ( present(dimids   ) ) then
      if ( size(dimids) /= varp%ndim ) then
        write(*, '("size of dimension id array (",i6,") should equal number of dimensions (",i6,")")') size(dimids), varp%ndim
        stat=1; return
      end if
      dimids = varp%dimids(1:varp%ndim)
    end if
    if ( present(natts) ) then
      natts = varp%natt
    end if
    
    ! special:
    if ( present(shp) ) then
      if ( size(shp) /= varp%ndim ) then
        write(*, '("size of shape array (",i6,") should equal number of dimensions (",i6,")")') size(shp), varp%ndim
        stat=1; return
      end if
      shp = varp%shp(1:varp%ndim)
    end if
    
    ! ok
    stat = 0
    
  end subroutine DIOS_Inquire_Variable
  
  
  ! ***
  
  
  subroutine DIOS_Inquire_Attribute( hid, varid, name, stat, xtype, length )
  
    use NetCDF, only : NF90_Inquire_Attribute
    use NetCDF, only : NF90_GLOBAL
    use NetCDF, only : NF90_CHAR, NF90_BYTE, NF90_SHORT, NF90_INT, NF90_FLOAT, NF90_DOUBLE

    ! --- in/out -------------------------------------

    integer, intent(in)                           ::  hid
    integer, intent(in)                           ::  varid
    character(len=*), intent(in)                  ::  name
    integer, intent(out)                          ::  stat
    
    integer, intent(out), optional                ::  xtype
    integer, intent(out), optional                ::  length

    ! --- const --------------------------------------

    
    ! --- external -------------------------------

    
    ! --- local --------------------------------------

    integer                         ::  ftype
    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp

    integer                         ::  netcdf_id

    ! --- begin --------------------------------------
    
    ! single type:
    call DIOS_Get_Type( hid, ftype, stat )
    

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    

    ! select appropriate routine:
    select case ( ftype )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_ASCII )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case ( DIOS_NETCDF )
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! global or variable attribute ?
        if ( varid == DIOS_GLOBAL ) then
          ! file id:
          netcdf_id = NF90_GLOBAL
        else
          ! pointer to variable structure:
          call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
          ! variable id:
          netcdf_id = varp%netcdf_varid
        end if
        ! get type etc:
        stat = NF90_Inquire_Attribute( filep%netcdf_id, netcdf_id, trim(name), &
                         xtype=xtype, len=length )
        
        ! return type ?
        if ( present(xtype) ) then
          ! convert:
          select case ( xtype )
            case ( NF90_CHAR   ) ;  xtype = DIOS_CHAR  
            case ( NF90_BYTE   ) ;  xtype = DIOS_BYTE  
            case ( NF90_SHORT  ) ;  xtype = DIOS_SHORT 
            case ( NF90_INT    ) ;  xtype = DIOS_INT   
            case ( NF90_FLOAT  ) ;  xtype = DIOS_FLOAT 
            case ( NF90_DOUBLE ) ;  xtype = DIOS_DOUBLE
            case default
              write(*, '("unsupported data type : ",i6)') xtype
              stat=1; return
          end select
        end if

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        write(*, '("unsupported filetype : ",i6)') filep%ftype
        stat=1; return
    end select

    ! ok
    stat = 0
    
  end subroutine DIOS_Inquire_Attribute
  
  
  ! ***
  
  
  subroutine DIOS_Inq_AttName( hid, varid, attnum, name, stat )
  
    use NetCDF, only : NF90_Inq_AttName, NF90_GLOBAL

    ! --- in/out -------------------------------------

    integer, intent(in)                           ::  hid
    integer, intent(in)                           ::  varid
    integer, intent(in)                           ::  attnum  ! 1,..,natt
    character(len=*), intent(out)                 ::  name
    integer, intent(out)                          ::  stat

    ! --- const --------------------------------------

    
    ! --- external -------------------------------

    
    ! --- local --------------------------------------

    integer                         ::  ftype
    type(DIOS_File), pointer         ::  filep
    type(DIOS_Var), pointer          ::  varp

    ! --- begin --------------------------------------
    
    ! single type:
    call DIOS_Get_Type( hid, ftype, stat )
    

    ! pointer to file structure:
    call DIOS_File_List_Get_Pointer( File_List, hid, filep, stat )
    
    
    ! global attribute ?
    if ( varid == DIOS_GLOBAL ) then
      
      ! select appropriate routine:
      select case ( ftype )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_ASCII )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_NETCDF )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          ! get name:
          stat = NF90_Inq_AttName( filep%netcdf_id, NF90_GLOBAL, attnum, name )
          
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*, '("unsupported filetype : ",i6)') filep%ftype
          stat=1; return
      end select
    
    else  ! variable attribute

      ! pointer to variable structure:
      call DIOS_Var_List_Get_Pointer( filep%Var_List, varid, varp, stat )
      
      ! select appropriate routine:
      select case ( ftype )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_ASCII )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case ( DIOS_NETCDF )
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          ! get name:
          stat = NF90_Inq_AttName( filep%netcdf_id, varp%netcdf_varid, attnum, name )
          
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
        case default
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*, '("unsupported filetype : ",i6)') filep%ftype
          stat=1; return
      end select

    end if  ! global or variable attribute
    
    ! ok
    stat = 0
    
  end subroutine DIOS_Inq_AttName


  ! ********************************************************************
  ! ***
  ! *** show content
  ! ***
  ! ********************************************************************
  

  subroutine DIOS_Show( filename, stat, filetype, show_data )

    ! --- in/out -------------------------------------

    character(len=*), intent(in)        ::  filename
    integer, intent(out)                ::  stat
    integer, intent(in), optional       ::  filetype
    logical, intent(in), optional       ::  show_data

    ! --- const --------------------------------------

    
    ! --- local --------------------------------------
    
    integer                           ::  l
    integer                           ::  ftype
    logical                           ::  do_show_data
    integer                           ::  hid
    integer                           ::  ndim, idim
    integer                           ::  nvar, ivar
    integer                           ::  natt
    logical                           ::  named
    character(len=LEN_NAME)           ::  ftype_name
    character(len=LEN_NAME)           ::  name
    integer                           ::  length
    logical                           ::  unlimited
    logical                           ::  isfirst
    integer                           ::  xtype
    integer                           ::  dimids(MAX_RANK)
    integer                           ::  shp(MAX_RANK)
    character(len=LEN_LINE)           ::  line
    character(len=LEN_NAME)           ::  val

    ! --- begin --------------------------------------
    
    ! show data ?
    do_show_data = .false.
    if ( present(show_data) ) do_show_data = show_data

    ! guess file type by default, or set to optional argument:
    ftype = DIOS_NONE
    if ( present(filetype) ) ftype = filetype
    
    ! guess file type ?
    if ( ftype == DIOS_NONE ) then
      ! length of filename:
      l = len_trim(filename)
      ! guess ...
      if ( (l > 3) .and. (filename(l-2:l) == '.nc') ) then
        ftype = DIOS_NETCDF
      else if ( (l > 4) .and. ( (filename(l-3:l) == '.dat') .or. (filename(l-3:l) == '.txt') ) ) then
        ftype = DIOS_ASCII
      else
        write(*, '("could not guess file type from file name:")')
        write(*, '("  ",a)') trim(filename)
        stat=1; return
      end if
    end if
    
    ! filetype name:
    ftype_name = trim(DIOS_FILETYPE_NAME(ftype))
    if ( ftype == DIOS_NETCDF ) then
      call NetCDF_Get_FileType( trim(filename), ftype_name, stat )
    end if

    ! open file:
    call DIOS_Open( trim(filename), ftype, DIOS_READ, hid, stat )

    ! header line:
    !   <filetype> <filename> {
    write(*, '(a," ",a," {")') trim(ftype_name), trim(filename)
    
    ! number of dimensions:
    call DIOS_Inquire( hid, stat, nDimensions=ndim )
    ! list dimensions ?
    if ( ndim > 0 ) then
      ! init flag:
      isfirst = .true.
      ! loop over dimensions:
      do idim = 1, ndim
        ! write lines:
        !   x = 4 ;
        !   t = UNLIMITED ; // (5 currently)
        !   ...
        ! get name and length:
        call DIOS_Inquire_Dimension( hid, idim, stat, name=name, named=named, &
                                          length=length, unlimited=unlimited )
        ! skip ?
        if ( .not. named ) cycle
        ! display header ?
        if ( isfirst ) then
          write(*, '("dimensions:")')
          isfirst = .false.
        end if
        ! display:
        if ( unlimited ) then
          write(val,*) length
          val = adjustl(val)
          write(*, '("  ",a," = UNLIMITED ; // (",a," currently)")') trim(name), trim(val)
        else if ( named ) then
          write(*, '("  ",a," = ",i6," ;")') trim(name), length
        end if
      end do   ! dimensions
    end if   ! ndim > 0
    
    ! number of variables:
    call DIOS_Inquire( hid, stat, nVariables=nvar )
    ! list variables ?
    if ( nvar > 0 ) then
      ! start:
      write(*, '("variables:")')
      ! loop over variables:
      do ivar = 1, nvar
        ! write lines:
        !   float afield(y, x) ;
        !     afield:unit = "m" ;
        !     ...
        ! get name etc:
        call DIOS_Inquire_Variable( hid, ivar, stat, name=name, xtype=xtype, ndims=ndim, natts=natt )
        ! get dimension id's now the number is known:
        call DIOS_Inquire_Variable( hid, ivar, stat, dimids=dimids(1:ndim) )
        ! start line with type and variable name:
        write(line,'("  ",a," ",a,"(")') trim(DIOS_DATATYPE_NAME(xtype)), trim(name)
        ! loop over dimensions:
        shp = 1
        do idim = 1, ndim
          ! get dimension name:
          call DIOS_Inquire_Dimension( hid, dimids(idim), stat, name=name, named=named, &
                                              length=length, unlimited=unlimited )
          ! add to line:
          if ( idim > 1 ) line = trim(line)//','
          ! name or number ...
          if ( named ) then
            line = trim(line)//' '//trim(name)
          else
            write(val,*) length
            line = trim(line)//' '//adjustl(val)
            if (unlimited) line = trim(line)//'/Inf'
          end if
          ! store for show_data:
          shp(idim) = length
        end do
        ! close line:
        line = trim(line)//' ) ;'
        ! display dimension name(dims) line:
        write(*, '(a)') trim(line)
        ! write attributes:
        call DIOS_Show_Attributes( hid, ivar, natt, stat )
        ! show data ?
        if ( do_show_data ) then
          ! display data:
          call DIOS_Show_Data( hid, ivar, xtype, ndim, shp, stat )
        end if
      end do  ! variables
    end if  ! nvar > 0
    
    ! number of global attributes:
    call DIOS_Inquire( hid, stat, nAttributes=natt )
    ! global attributes ?
    if ( natt > 0 ) then
      ! intro:
      write(*, '("")')
      write(*, '("// global attributes:")')
      ! display attributes:
      call DIOS_Show_Attributes( hid, DIOS_GLOBAL, natt, stat )
    end if

    ! closure:
    !   }
    write(*, '("}")');

    ! close file:
    call DIOS_Close( hid, stat )

    ! ok
    stat = 0
    
  end subroutine DIOS_Show
  
  
  ! ***
  
  
  subroutine DIOS_Show_Attributes( hid, ivar, natt, stat )

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  ivar
    integer, intent(in)                 ::  natt
    integer, intent(out)                ::  stat

    ! --- const --------------------------------------

    
    ! --- local --------------------------------------
    
    integer                           ::  l
    integer                           ::  iatt
    character(len=LEN_NAME)           ::  name
    integer                           ::  length
    integer                           ::  xtype
    character(len=LEN_LINE)           ::  line
    character(len=LEN_NAME)           ::  val
    integer                           ::  ival
    integer(4), allocatable           ::  values_i4(:)
    real(8), allocatable              ::  values_r8(:)

    ! --- begin --------------------------------------
  
    ! loop over attributes:
    do iatt = 1, natt
      ! get name:
      call DIOS_Inq_AttName( hid, ivar, iatt, name, stat )
      ! get type and length:
      call DIOS_Inquire_Attribute( hid, ivar, name, stat, xtype=xtype, length=length )
      !! info ...
      !write(*, '("    ",a," <",a,"> ",i4," ;")') trim(name), trim(DIOS_DATATYPE_NAME(xtype)), length
      ! different per type ...
      select case ( xtype )
        !character values
        case ( DIOS_CHAR )
          ! get value:
          call DIOS_Get_Att( hid, ivar, name, line, stat )
          if (stat/=0) then
            ! somthing went wrong (attribute too large ?)
            line = '...'
          else
            ! ok; but not too much to the screen ...
            if ( len_trim(line) > 400 )  line = line(1:400)//'...'
          end if
          ! display ..


          ! variable rank etc:
          write(*, '("    ",a," = """,a,""" ;")') trim(name), trim(line);

        !integer values
        case ( DIOS_BYTE, DIOS_SHORT, DIOS_INT )
          ! storage:
          allocate( values_i4(length) )
          ! fill:
          call DIOS_Get_Att( hid, ivar, name, values_i4, stat )
          if (stat/=0) then
            ! somthing went wrong (attribute too large ?)
            line = '...'
          else
            ! loop over values:
            line = ''
            do ival = 1, min(length,50)
              ! add seperation if necessary:
              if ( ival > 1 ) line = trim(line)//','
              ! dump value:
              write(val,*) values_i4(ival)
              ! shift to left:
              val = adjustl(val)
              ! add to line:
              line = trim(line)//' '//trim(val)
              ! add type indicator if necessary:
              if ( xtype == DIOS_BYTE  ) line = trim(line)//'b'
              if ( xtype == DIOS_SHORT ) line = trim(line)//'s'
            end do
            if ( ival < length ) line=trim(line)//' ...'
          end if
          ! display:
          write(*, '("    ",a," =",a)') trim(name), trim(line)
          ! clear:
          deallocate( values_i4 )
        !floating point values
        case ( DIOS_FLOAT, DIOS_DOUBLE )
          ! storage:
          allocate( values_r8(length) )
          ! fill:
          call DIOS_Get_Att( hid, ivar, name, values_r8, stat )
          if (stat/=0) then
            ! somthing went wrong (attribute too large ?)
            line = '...'
          else
            ! loop over values:
            line = ''
            do ival = 1, min(length,50)
              ! add seperation if necessary:
              if ( ival > 1 ) line = trim(line)//','
              ! dump value:
              write(val,*) values_r8(ival)
              ! remove tailing zeros:
              do l = len_trim(val), 1, -1
                if ( val(l:l) /= '0' ) exit
                val(l:l) = ' '
              end do
              ! shift to left:
              val = adjustl(val)
              ! add to line:
              line = trim(line)//' '//trim(val)
              ! add type indicator if necessary:
              if ( xtype == DIOS_FLOAT  ) line = trim(line)//'f'
            end do
            if ( ival < length ) line=trim(line)//' ...'
          end if
          ! display:
          write(*, '("    ",a," =",a)') trim(name), trim(line)
          ! clear:
          deallocate( values_r8 )
        !other ...
        case default
          write(*, '("INTERNAL - unsupported data type : ",i6)') xtype
          stat=1; return
      end select
    end do
    
    ! ok
    stat = 0
    
  end subroutine DIOS_Show_Attributes
  
  
  ! ***
  
  
  subroutine DIOS_Show_Data( hid, ivar, xtype, rank, shp, stat )

    ! --- in/out -------------------------------------

    integer, intent(in)                 ::  hid
    integer, intent(in)                 ::  ivar
    integer, intent(in)                 ::  xtype
    integer, intent(in)                 ::  rank
    integer, intent(in)                 ::  shp(MAX_RANK)
    integer, intent(out)                ::  stat

    ! --- const --------------------------------------

    
    ! --- local --------------------------------------
    
    character(len=LEN_LINE)                 ::  line
    character(len=LEN_NAME)                 ::  val
    integer                                 ::  l
    integer                                 ::  i1,i2,i3,i4,i5,i6,i7
    integer(4), allocatable                 ::  values_i4(:,:,:,:,:,:,:)
    real(8), allocatable                    ::  values_r8(:,:,:,:,:,:,:)
    character(len=shp(1)), allocatable      ::  values_c   (:,:,:,:,:,:)

    ! --- begin --------------------------------------
    
    ! per type:
    select case ( xtype )

      ! character values:
      case ( DIOS_CHAR )
        ! storage:
        allocate( values_c(shp(2),shp(3),shp(4),shp(5),shp(6),shp(7)), stat=stat )
        
        ! read:
        select case ( rank )
          case ( 1 ) ; call DIOS_Get_Var( hid, ivar, values_c(1,1,1,1,1,1), stat, count=shp(1:rank) )
          case ( 2 ) ; call DIOS_Get_Var( hid, ivar, values_c(:,1,1,1,1,1), stat, count=shp(1:rank) )
          case ( 3 ) ; call DIOS_Get_Var( hid, ivar, values_c(:,:,1,1,1,1), stat, count=shp(1:rank) )
          case ( 4 ) ; call DIOS_Get_Var( hid, ivar, values_c(:,:,:,1,1,1), stat, count=shp(1:rank) )
          case ( 5 ) ; call DIOS_Get_Var( hid, ivar, values_c(:,:,:,:,1,1), stat, count=shp(1:rank) )
          case ( 6 ) ; call DIOS_Get_Var( hid, ivar, values_c(:,:,:,:,:,1), stat, count=shp(1:rank) )
          case ( 7 ) ; call DIOS_Get_Var( hid, ivar, values_c(:,:,:,:,:,:), stat, count=shp(1:rank) )
          case default
            write(*, '("INTERNAL - unsupported rank : ",i6)') rank
            stat=1; return
        end select
        
        ! loop over higher dimensions:
        do i7 = 1, shp(7)
          do i6 = 1, shp(6)
            do i5 = 1, shp(5)
              do i4 = 1, shp(4)
                do i3 = 1, shp(3)
                  ! plot index of higer dimensions ?
                  if ( rank > 2 ) then
                    line = '      (:,:'
                    if ( rank >= 3 ) write(line,'(a,",",i4)') trim(line), i3
                    if ( rank >= 4 ) write(line,'(a,",",i4)') trim(line), i4
                    if ( rank >= 5 ) write(line,'(a,",",i4)') trim(line), i5
                    if ( rank >= 6 ) write(line,'(a,",",i4)') trim(line), i6
                    if ( rank >= 7 ) write(line,'(a,",",i4)') trim(line), i7
                    line = trim(line)//')'
                    write(*, '(a)') trim(line)
                  end if
                  ! display matrix:
                  do i2 = 1, shp(2)
                    ! copy value:
                    line = values_c(i2,i3,i4,i5,i6,i7)
                    ! display:
                    write(*, '("        `",a,"` ;")') trim(line)
                  end do  ! i2
                end do  ! i3
              end do  ! i4
            end do  ! i5
          end do  ! i6
        end do  ! i7
        ! clear:
        deallocate( values_c, stat=stat )
        

      ! integer values:
      case ( DIOS_BYTE, DIOS_SHORT, DIOS_INT )
        ! storage:
        allocate( values_i4(shp(1),shp(2),shp(3),shp(4),shp(5),shp(6),shp(7)), stat=stat )
        
        ! read:
        select case ( rank )
          case ( 1 ) ; call DIOS_Get_Var( hid, ivar, values_i4(:,1,1,1,1,1,1), stat )
          case ( 2 ) ; call DIOS_Get_Var( hid, ivar, values_i4(:,:,1,1,1,1,1), stat )
          case ( 3 ) ; call DIOS_Get_Var( hid, ivar, values_i4(:,:,:,1,1,1,1), stat )
          case ( 4 ) ; call DIOS_Get_Var( hid, ivar, values_i4(:,:,:,:,1,1,1), stat )
          case ( 5 ) ; call DIOS_Get_Var( hid, ivar, values_i4(:,:,:,:,:,1,1), stat )
          case ( 6 ) ; call DIOS_Get_Var( hid, ivar, values_i4(:,:,:,:,:,:,1), stat )
          case ( 7 ) ; call DIOS_Get_Var( hid, ivar, values_i4(:,:,:,:,:,:,:), stat )
          case default
            write(*, '("INTERNAL - unsupported rank : ",i6)') rank
            stat=1; return
        end select
        
        ! loop over higher dimensions:
        do i7 = 1, shp(7)
          do i6 = 1, shp(6)
            do i5 = 1, shp(5)
              do i4 = 1, shp(4)
                do i3 = 1, shp(3)
                  ! plot index of higer dimensions ?
                  if ( rank > 2 ) then
                    line = '      (:,:'
                    if ( rank >= 3 ) write(line,'(a,",",i4)') trim(line), i3
                    if ( rank >= 4 ) write(line,'(a,",",i4)') trim(line), i4
                    if ( rank >= 5 ) write(line,'(a,",",i4)') trim(line), i5
                    if ( rank >= 6 ) write(line,'(a,",",i4)') trim(line), i6
                    if ( rank >= 7 ) write(line,'(a,",",i4)') trim(line), i7
                    line = trim(line)//')'
                    write(*, '(a)') trim(line)
                  end if
                  ! display matrix:
                  do i2 = 1, shp(2)
                    line = ''
                    do i1 = 1, shp(1)
                      !! not all ?
                      !if ( i1 > 10 ) then
                      !  line = trim(line)//' ...'
                      !  exit
                      !end if
                      ! add seperation if necessary:
                      if ( i1 > 1 ) line = trim(line)//','
                      ! dump value:
                      write(val,*) values_i4(i1,i2,i3,i4,i5,i6,i7)
                      ! shift to left:
                      val = adjustl(val)
                      ! add to line:
                      line = trim(line)//' '//trim(val)
                      ! add type indicator if necessary:
                      if ( xtype == DIOS_BYTE  ) line = trim(line)//'b'
                      if ( xtype == DIOS_SHORT ) line = trim(line)//'s'
                      ! line too long already ?
                      if ( len_trim(line) > 72 ) then
                        ! display:
                        write(*, '("        ",a," ;")') trim(line)
                        ! empty:
                        line = ''
                      end if
                    end do  ! i1
                    ! display:
                    if ( len_trim(line) > 0 ) then
                      write(*, '("        ",a," ;")') trim(line)
                    end if
                  end do  ! i2
                end do  ! i3
              end do  ! i4
            end do  ! i5
          end do  ! i6
        end do  ! i7
        ! clear:
        deallocate( values_i4, stat=stat )
        

      ! real values:
      case ( DIOS_FLOAT, DIOS_DOUBLE )
        ! storage:
        allocate( values_r8(shp(1),shp(2),shp(3),shp(4),shp(5),shp(6),shp(7)), stat=stat )
        
        ! read:
        select case ( rank )
          case ( 1 ) ; call DIOS_Get_Var( hid, ivar, values_r8(:,1,1,1,1,1,1), stat )
          case ( 2 ) ; call DIOS_Get_Var( hid, ivar, values_r8(:,:,1,1,1,1,1), stat )
          case ( 3 ) ; call DIOS_Get_Var( hid, ivar, values_r8(:,:,:,1,1,1,1), stat )
          case ( 4 ) ; call DIOS_Get_Var( hid, ivar, values_r8(:,:,:,:,1,1,1), stat )
          case ( 5 ) ; call DIOS_Get_Var( hid, ivar, values_r8(:,:,:,:,:,1,1), stat )
          case ( 6 ) ; call DIOS_Get_Var( hid, ivar, values_r8(:,:,:,:,:,:,1), stat )
          case ( 7 ) ; call DIOS_Get_Var( hid, ivar, values_r8(:,:,:,:,:,:,:), stat )
          case default
            write(*, '("INTERNAL - unsupported rank : ",i6)') rank
            stat=1; return
        end select
        
        ! loop over higher dimensions:
        do i7 = 1, shp(7)
          do i6 = 1, shp(6)
            do i5 = 1, shp(5)
              do i4 = 1, shp(4)
                do i3 = 1, shp(3)
                  ! plot index of higer dimensions ?
                  if ( rank > 2 ) then
                    line = '      (:,:'
                    if ( rank >= 3 ) write(line,'(a,",",i4)') trim(line), i3
                    if ( rank >= 4 ) write(line,'(a,",",i4)') trim(line), i4
                    if ( rank >= 5 ) write(line,'(a,",",i4)') trim(line), i5
                    if ( rank >= 6 ) write(line,'(a,",",i4)') trim(line), i6
                    if ( rank >= 7 ) write(line,'(a,",",i4)') trim(line), i7
                    line = trim(line)//')'
                    write(*, '(a)') trim(line)
                  end if
                  ! display matrix:
                  do i2 = 1, shp(2)
                    line = ''
                    do i1 = 1, shp(1)
                      !! not all ?
                      !if ( i1 > 10 ) then
                      !  line = trim(line)//' ...'
                      !  exit
                      !end if
                      ! add seperation if necessary:
                      if ( i1 > 1 ) line = trim(line)//','
                      ! dump value:
                      write(val,*) values_r8(i1,i2,i3,i4,i5,i6,i7)
                      ! remove tailing zeros:
                      do l = len_trim(val), 1, -1
                        if ( val(l:l) /= '0' ) exit
                        val(l:l) = ' '
                      end do
                      ! shift to left:
                      val = adjustl(val)
                      ! add to line:
                      line = trim(line)//' '//trim(val)
                      ! add type indicator if necessary:
                      if ( xtype == DIOS_FLOAT  ) line = trim(line)//'f'
                      ! line too long already ?
                      if ( len_trim(line) > 72 ) then
                        ! display:
                        write(*, '("        ",a," ;")') trim(line)
                        ! empty:
                        line = ''
                      end if
                    end do  ! i1
                    ! display:
                    if ( len_trim(line) > 0 ) then
                      write(*, '("        ",a," ;")') trim(line)
                    end if
                  end do  ! i2
                end do  ! i3
              end do  ! i4
            end do  ! i5
          end do  ! i6
        end do  ! i7
        ! clear:
        deallocate( values_r8, stat=stat )
        

      case default
        write(*, '("INTERNAL - unsupported data type : ",i6)') xtype
        stat=1; return
    end select
    
    ! ok
    stat = 0
    
  end subroutine DIOS_Show_Data


  subroutine handle_err(stat)
    use NetCDF, only : NF90_NOERR
    use NetCDF, only : nf90_strerror

    integer, intent ( in) :: stat
    if(stat /= NF90_NOERR) then
      write(*,*) trim(nf90_strerror(stat))
      ! stop "Stopped"
    end if
  end subroutine handle_err


  subroutine check_var_values_rank(values_rank, var_ndim, stat)
    integer, intent(in)   ::  values_rank
    integer, intent(in)   ::  var_ndim
    integer, intent(out)  ::  stat 

    ! check ...
    if ( values_rank > var_ndim ) then
      write(*, '("Rank of values to be written exceeds variable dimension:")')
      write(*, '("  rank values          : ",i6)') values_rank
      write(*, '("  variable dimension   : ",i6)') var_ndim
      stat=1; return
    end if

    ! ok
    stat = 0
  end subroutine check_var_values_rank


  subroutine check_var_start_size(start_size, var_ndim, stat)
    integer, intent(in)   ::  start_size
    integer, intent(in)   ::  var_ndim
    integer, intent(out)  ::  stat 
      
    ! check ...
    if ( start_size /= var_ndim ) then
      write(*, '("Size of position argument not equal to variable dimension:")')
      write(*, '("  size start  : ",i6)') start_size
      write(*, '("  var dim     : ",i6)') var_ndim
      stat=1; return
    end if

    ! ok
    stat = 0
  end subroutine check_var_start_size


  subroutine check_var_count_size(count_size, var_ndim, stat)
    integer, intent(in)   ::  count_size
    integer, intent(in)   ::  var_ndim
    integer, intent(out)  ::  stat 

    ! check ...
    if ( count_size /= var_ndim ) then
      write(*, '("Size of position argument not equal to variable dimension:")')
      write(*, '("  size count  : ",i6)') count_size
      write(*, '("  var dim     : ",i6)') var_ndim
      stat=1; return
    end if

    ! ok
    stat = 0
  end subroutine check_var_count_size


  subroutine check_var_stride_size(stride_size, var_ndim, stat)
    integer, intent(in)   ::  stride_size
    integer, intent(in)   ::  var_ndim
    integer, intent(out)  ::  stat 

    ! check ...
    if ( stride_size /= var_ndim ) then
      write(*, '("Size of position argument not equal to variable dimension:")')
      write(*, '("  size stride  : ",i6)') stride_size
      write(*, '("  var dim     : ",i6)') var_ndim
      stat=1; return
    end if

    ! ok
    stat = 0
  end subroutine check_var_stride_size


  subroutine check_var_map_size(map_size, var_ndim, stat)
    integer, intent(in)   ::  map_size
    integer, intent(in)   ::  var_ndim
    integer, intent(out)  ::  stat 

    ! check ...
    if ( map_size /= var_ndim ) then
      write(*, '("Size of position argument not equal to variable dimension:")')
      write(*, '("  size map  : ",i6)') map_size
      write(*, '("  var dim     : ",i6)') var_ndim
      stat=1; return
    end if

    ! ok
    stat = 0
  end subroutine check_var_map_size


end module dios_mod

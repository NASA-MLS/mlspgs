
! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module Hdf

  use MLSCommon, only: R8, R4
  implicit none
  public                                ! Most things are public

  private :: Id                         ! Except Id

  !------------------- RCS Ident Info -----------------------
  character(LEN=130) :: Id = &        
    "$Id$"
  !----------------------------------------------------------

  ! Contents: 

  !     Tag definitions
  !     Error return codes
  !     Logical constants

  ! Remarks: This file can be included with Fortran user programs.  As a
  !          general rule, don't use DFNT constants that don't include a
  !          number in their name.  E.g., don't use DFNT_FLOAT, use
  !          DFNT_FLOAT32 or DFNT_FLOAT64.  The DFNT constants that don't
  !          include numbers are for backward compatibility only.  Also,
  !          there are no current plans to support 128-bit number types.
  !          For more information about constants in this file, see the
  !          equivalent constant declarations in the C include file 'hdf.h'

  ! Now define f90 interfaces for some toolkit routines.

  private :: sfrdata
  integer, external :: sfrdata
  integer, external :: sfwdata

  interface

    integer function sfstart(filename, access_mode)
      character (LEN=*), intent (IN) :: filename
      integer, intent (IN) :: access_mode
    end function sfstart

    integer function sfcreate(sd_id, sds_name, data_type,rank, dim_sizes)
      character (LEN=*), intent (IN) :: sds_name
      integer, intent (IN) :: sd_id, data_type, rank
      integer, intent (IN) :: dim_sizes(*)
    end function sfcreate

    integer function sfdimid(sds_id, dim_index)
      integer, intent (IN) :: sds_id, dim_index
    end function sfdimid

    integer function sfsdmname(dim_id, dim_name)
      character (LEN=*), intent (IN) :: dim_name
      integer, intent (IN) :: dim_id
    end function sfsdmname

    integer function sfsdscale(dim_id, n_values, data_type, data)
      integer, intent(IN) :: dim_id
      integer, intent(IN) :: n_values
      integer, intent(IN) :: data_type
      double precision, intent(IN), dimension(*) :: data
    end function sfsdscale

!    integer function sfwdata(sds_id, start, stride, edges , data)      
!      integer, intent (IN) :: sds_id
!      integer, intent (IN), dimension (*) :: start, stride, edges
!      double precision, intent (IN), dimension (*) :: data
!    end function sfwdata

    integer function sfendacc(sds_id)
      integer, intent (IN) :: sds_id
    end function sfendacc

    integer function sfend(sd_id)
      integer, intent (IN) :: sd_id
    end function sfend

    integer function sfrcdata(sds_id, start, stride, edge, buffer)
      integer, intent (IN) :: sds_id
      integer, intent (IN), dimension (*) :: start, stride, edge
      character(len=*), intent (OUT), dimension (*) :: buffer
    end function sfrcdata

    integer function sfn2index(sd_id, sds_name)      
      integer, intent (IN) :: sd_id
      character (LEN=*), intent (IN) :: sds_name
    end function sfn2index

    integer function sfselect(sd_id, sds_index)      
      integer, intent (IN) :: sd_id
      integer, intent (IN) :: sds_index
    end function sfselect

    integer function sfginfo(sds_id, sds_name, rank, dimsizes, data_type, &
      & num_attrs)      
      integer, intent (IN) :: sds_id
      integer, intent (OUT) :: rank, data_type, num_attrs
      character (LEN=*), intent (OUT) :: sds_name
      integer, intent (IN), dimension (*) :: dimsizes
    end function sfginfo

    integer function sfgdinfo(dim_id, dim_name, size, data_type, &
      & num_attrs)      
      integer, intent (IN) :: dim_id
      integer, intent (OUT) :: size, data_type, num_attrs
      character (LEN=*), intent (OUT) :: dim_name
    end function sfgdinfo

  end interface

  interface sfrdata_f90
    module procedure sfrdata_i1, sfrdata_i2, sfrdata_i3, &
      & sfrdata_r1, sfrdata_r2, sfrdata_r3, &
      & sfrdata_dp1, sfrdata_dp2, sfrdata_dp3
  end interface

  interface sfwdata_f90
    module procedure sfwdata_i1, sfwdata_dp1, sfwdata_r1
    module procedure sfwdata_r3, sfwdata_dp3
  end interface


  ! For routines where one of the arguments can be any valid data type !
  ! an interface block causes trouble. One either builds a generic
  ! interface or plays fast-and-loose by declaring the routine to be
  ! "external" in the old-fashioned Fortran-77 way. We do the latter for 
  ! sfsfill
  integer::sfsfill
  external::sfsfill

  !-------------------
  ! Error Return Codes 
  !-------------------

  integer, parameter :: DFE_NOERROR              =     0
  integer, parameter :: DFE_NONE                 =     0
  integer, parameter :: DFE_FNF                  =    -1
  integer, parameter :: DFE_DENIED               =    -2
  integer, parameter :: DFE_ALROPEN              =    -3
  integer, parameter :: DFE_TOOMANY              =    -4
  integer, parameter :: DFE_BADNAME              =    -5
  integer, parameter :: DFE_BADACC               =    -6
  integer, parameter :: DFE_BADOPEN              =    -7
  integer, parameter :: DFE_NOTOPEN              =    -8
  integer, parameter :: DFE_CANTCLOSE            =    -9
  integer, parameter :: DFE_DFNULL               =   -10
  integer, parameter :: DFE_ILLTYPE              =   -11
  integer, parameter :: DFE_UNSUPPORTED          =   -12
  integer, parameter :: DFE_BADDDLIST            =   -13
  integer, parameter :: DFE_NOTDFFILE            =   -14
  integer, parameter :: DFE_SEEDTWICE            =   -15
  integer, parameter :: DFE_NOSPACE              =   -16
  integer, parameter :: DFE_NOSUCHTAG            =   -17
  integer, parameter :: DFE_READERROR            =   -18

  integer, parameter :: DFE_WRITEERROR           =   -19
  integer, parameter :: DFE_SEEKERROR            =   -20
  integer, parameter :: DFE_NOFREEDD             =   -21
  integer, parameter :: DFE_BADTAG               =   -22
  integer, parameter :: DFE_BADREF               =   -23
  integer, parameter :: DFE_RDONLY               =   -24
  integer, parameter :: DFE_BADCALL              =   -25
  integer, parameter :: DFE_BADPTR               =   -26
  integer, parameter :: DFE_BADLEN               =   -27
  integer, parameter :: DFE_BADSEEK              =   -28
  integer, parameter :: DFE_NOMATCH              =   -29
  integer, parameter :: DFE_NOTINSET             =   -30
  integer, parameter :: DFE_BADDIM               =   -31
  integer, parameter :: DFE_BADOFFSET            =   -32
  integer, parameter :: DFE_BADSCHEME            =   -33
  integer, parameter :: DFE_NODIM                =   -34
  integer, parameter :: DFE_NOTENOUGH            =   -35
  integer, parameter :: DFE_NOVALS               =   -36
  integer, parameter :: DFE_CORRUPT              =   -37
  integer, parameter :: DFE_BADFP                =   -38

  integer, parameter :: DFE_NOREF                =   -39
  integer, parameter :: DFE_BADDATATYPE          =   -40
  integer, parameter :: DFE_BADMCTYPE            =   -41
  integer, parameter :: DFE_BADNUMTYPE           =   -42
  integer, parameter :: DFE_BADORDER             =   -43
  integer, parameter :: DFE_ARGS                 =   -44
  integer, parameter :: DFE_INTERNAL             =   -45
  integer, parameter :: DFE_DUPDD                =   -46
  integer, parameter :: DFE_CANTMOD              =   -47
  integer, parameter :: DFE_RANGE                =   -48
  integer, parameter :: DFE_BADTABLE             =   -49
  integer, parameter :: DFE_BADSDG               =   -50
  integer, parameter :: DFE_BADNDG               =   -51
  integer, parameter :: DFE_BADFIELDS            =   -52
  integer, parameter :: DFE_NORESET              =   -53
  integer, parameter :: DFE_NOVS                 =   -54
  integer, parameter :: DFE_VGSIZE               =   -55
  integer, parameter :: DFE_DIFFFILES            =   -56
  integer, parameter :: DFE_VTAB                 =   -57
  integer, parameter :: DFE_BADAID               =   -58

  integer, parameter :: DFE_OPENAID              =   -59
  integer, parameter :: DFE_BADCONV              =   -60
  integer, parameter :: DFE_GENAPP               =   -61
  integer, parameter :: DFE_CANTFLUSH            =   -62
  integer, parameter :: DFE_BADTYPE              =   -63
  integer, parameter :: DFE_SYMSIZE              =   -64
  integer, parameter :: DFE_BADATTACH            =   -65
  integer, parameter :: DFE_CANTDETACH           =   -66

  !---------------------------
  ! internal file access codes
  !---------------------------

  integer, parameter :: DFACC_READ               =     1
  integer, parameter :: DFACC_WRITE              =     2
  integer, parameter :: DFACC_CREATE             =     4
  integer, parameter :: DFACC_ALL                =     7
  integer, parameter :: DFACC_RDONLY             =     1
  integer, parameter :: DFACC_RDWR               =     3
  integer, parameter :: DFACC_CLOBBER            =     4

  !---------------------------------
  ! Access types for SDsetaccesstype
  !---------------------------------

  integer, parameter :: DFACC_DEFAULT            =     0
  integer, parameter :: DFACC_SERIAL             =     1
  integer, parameter :: DFACC_PARALLEL           =     9

  !---------------------------
  ! Constants for DFSDsetorder
  !---------------------------

  integer, parameter :: DFO_FORTRAN              =     1
  integer, parameter :: DFO_C                    =     2

  !----------------------------------
  ! Definitions of storage convention
  !----------------------------------

  integer, parameter :: DFNTF_IEEE               =     1
  integer, parameter :: DFNTF_VAX                =     2
  integer, parameter :: DFNTF_CRAY               =     3
  integer, parameter :: DFNTF_PC                 =     4
  integer, parameter :: DFNTF_CONVEX             =     5
  integer, parameter :: DFNTF_VP                 =     6

  !----------------
  ! Masks for types
  !----------------

  integer, parameter :: DFNT_HDF                 =     0
  integer, parameter :: DFNT_NATIVE              =  4096
  integer, parameter :: DFNT_CUSTOM              =  8192
  integer, parameter :: DFNT_LITEND              = 16384

  !-----------------------
  ! Number type info codes 
  !-----------------------

  integer, parameter :: DFNT_NONE                =     0
  integer, parameter :: DFNT_QUERY               =     0
  integer, parameter :: DFNT_VERSION             =     1

  integer, parameter :: DFNT_FLOAT32             =     5
  integer, parameter :: DFNT_FLOAT               =     5
  integer, parameter :: DFNT_FLOAT64             =     6
  integer, parameter :: DFNT_DOUBLE              =     6
  integer, parameter :: DFNT_FLOAT128            =     7

  integer, parameter :: DFNT_INT8                =    20
  integer, parameter :: DFNT_UINT8               =    21
  integer, parameter :: DFNT_INT16               =    22
  integer, parameter :: DFNT_UINT16              =    23
  integer, parameter :: DFNT_INT32               =    24
  integer, parameter :: DFNT_UINT32              =    25
  integer, parameter :: DFNT_INT64               =    26
  integer, parameter :: DFNT_UINT64              =    27
  integer, parameter :: DFNT_INT128              =    28
  integer, parameter :: DFNT_UINT128             =    29

  integer, parameter :: DFNT_UCHAR8              =     3
  integer, parameter :: DFNT_UCHAR               =     3
  integer, parameter :: DFNT_CHAR8               =     4
  integer, parameter :: DFNT_CHAR                =     4
  integer, parameter :: DFNT_CHAR16              =    42
  integer, parameter :: DFNT_UCHAR16             =    43

  integer, parameter :: DFNT_NFLOAT32            =  4101
  integer, parameter :: DFNT_NFLOAT              =  4101
  integer, parameter :: DFNT_NFLOAT64            =  4102
  integer, parameter :: DFNT_NDOUBLE             =  4102
  integer, parameter :: DFNT_NFLOAT128           =  4103

  integer, parameter :: DFNT_NINT8               =  4116
  integer, parameter :: DFNT_NUINT8              =  4117
  integer, parameter :: DFNT_NINT16              =  4118
  integer, parameter :: DFNT_NUINT16             =  4119
  integer, parameter :: DFNT_NINT32	          =  4120
  integer, parameter :: DFNT_NUINT32             =  4121
  integer, parameter :: DFNT_NINT64              =  4122
  integer, parameter :: DFNT_NUINT64             =  4123
  integer, parameter :: DFNT_NINT128             =  4124
  integer, parameter :: DFNT_NUINT128            =  4125

  integer, parameter :: DFNT_NUCHAR8             =  4099
  integer, parameter :: DFNT_NUCHAR              =  4099
  integer, parameter :: DFNT_NCHAR8              =  4100
  integer, parameter :: DFNT_NCHAR               =  4100
  integer, parameter :: DFNT_NCHAR16             =  4138
  integer, parameter :: DFNT_NUCHAR16            =  4139

  integer, parameter :: DFNT_LFLOAT32            = 16389
  integer, parameter :: DFNT_LFLOAT              = 16389
  integer, parameter :: DFNT_LFLOAT64            = 16390
  integer, parameter :: DFNT_LDOUBLE             = 16390
  integer, parameter :: DFNT_LFLOAT128           = 16391

  integer, parameter :: DFNT_LINT8               = 16404
  integer, parameter :: DFNT_LUINT8              = 16405
  integer, parameter :: DFNT_LINT16              = 16406
  integer, parameter :: DFNT_LUINT16             = 16407
  integer, parameter :: DFNT_LINT32              = 16408
  integer, parameter :: DFNT_LUINT32             = 16409
  integer, parameter :: DFNT_LINT64              = 16410
  integer, parameter :: DFNT_LUINT64             = 16411
  integer, parameter :: DFNT_LINT128             = 16412
  integer, parameter :: DFNT_LUINT128            = 16413

  integer, parameter :: DFNT_LUCHAR8             = 16387
  integer, parameter :: DFNT_LUCHAR              = 16387
  integer, parameter :: DFNT_LCHAR8              = 16388
  integer, parameter :: DFNT_LCHAR               = 16388
  integer, parameter :: DFNT_LCHAR16             = 16426
  integer, parameter :: DFNT_LUCHAR16            = 16427

  !--------------
  ! tags and refs
  !--------------

  integer, parameter :: DFREF_WILDCARD           =     0
  integer, parameter :: DFTAG_WILDCARD           =     0
  integer, parameter :: DFTAG_NULL               =     1
  integer, parameter :: DFTAG_LINKED             =    20
  integer, parameter :: DFTAG_VERSION            =    30
  integer, parameter :: DFTAG_COMPRESSED         =    40

  !------------
  ! utility set
  !------------

  integer, parameter :: DFTAG_FID                =   100
  integer, parameter :: DFTAG_FD                 =   101
  integer, parameter :: DFTAG_TID                =   102
  integer, parameter :: DFTAG_TD                 =   103
  integer, parameter :: DFTAG_DIL                =   104
  integer, parameter :: DFTAG_DIA                =   105
  integer, parameter :: DFTAG_NT                 =   106
  integer, parameter :: DFTAG_MT                 =   107

  !-------------
  ! raster-8 set 
  !-------------

  integer, parameter :: DFTAG_ID8                =   200
  integer, parameter :: DFTAG_IP8                =   201
  integer, parameter :: DFTAG_RI8                =   202 
  integer, parameter :: DFTAG_CI8                =   203
  integer, parameter :: DFTAG_II8                =   204

  !-----------------
  ! Raster Image set
  !-----------------

  integer, parameter :: DFTAG_ID                 =   300 
  integer, parameter :: DFTAG_LUT                =   301
  integer, parameter :: DFTAG_RI                 =   302 
  integer, parameter :: DFTAG_CI                 =   303

  integer, parameter :: DFTAG_RIG                =   306 
  integer, parameter :: DFTAG_LD                 =   307
  integer, parameter :: DFTAG_MD                 =   308 
  integer, parameter :: DFTAG_MA                 =   309
  integer, parameter :: DFTAG_CCN                =   310 
  integer, parameter :: DFTAG_CFM                =   311
  integer, parameter :: DFTAG_AR                 =   312

  integer, parameter :: DFTAG_DRAW               =   400 
  integer, parameter :: DFTAG_RUN                =   401
  integer, parameter :: DFTAG_XYP                =   500 
  integer, parameter :: DFTAG_MTO                =   501

  !----------
  ! Tektronix 
  !----------

  integer, parameter :: DFTAG_T14                =   602 
  integer, parameter :: DFTAG_T105               =   603

  !--------------------
  ! Scientific Data set 
  !--------------------

  integer, parameter :: DFTAG_SDG                =   700 
  integer, parameter :: DFTAG_SDD                =   701
  integer, parameter :: DFTAG_SD                 =   702 
  integer, parameter :: DFTAG_SDS                =   703
  integer, parameter :: DFTAG_SDL                =   704 
  integer, parameter :: DFTAG_SDU                =   705
  integer, parameter :: DFTAG_SDF                =   706 
  integer, parameter :: DFTAG_SDM                =   707
  integer, parameter :: DFTAG_SDC                =   708 
  integer, parameter :: DFTAG_SDT                =   709
  integer, parameter :: DFTAG_SDLNK              =   710 
  integer, parameter :: DFTAG_NDG                =   720
  integer, parameter :: DFTAG_CAL                =   731 
  integer, parameter :: DFTAG_FV                 =   732
  integer, parameter :: DFTAG_BREQ               =   799
  integer, parameter :: DFTAG_EREQ               =   780

  !------
  ! VSets 
  !------

  integer, parameter :: DFTAG_VG                 =  1965 
  integer, parameter :: DFTAG_VH                 =  1962
  integer, parameter :: DFTAG_VS                 =  1963

  !--------------------
  ! compression schemes 
  !--------------------

  integer, parameter :: DFTAG_RLE                =    11 
  integer, parameter :: DFTAG_IMC                =    12
  integer, parameter :: DFTAG_IMCOMP             =    12 
  integer, parameter :: DFTAG_JPEG               =    13
  integer, parameter :: DFTAG_GREYJPEG           =    14

  !--------------
  ! SPECIAL CODES 
  !--------------

  integer, parameter :: SPECIAL_LINKED 	 =     1
  integer, parameter :: SPECIAL_EXT 	         =     2

  !-----------
  ! PARAMETERS 
  !-----------

  integer, parameter :: DF_MAXFNLEN              =   256 
  integer, parameter :: SD_UNLIMITED             =     0
  integer, parameter :: SD_DIMVAL_BW_COMP        =     1 
  integer, parameter :: SD_DIMVAL_BW_INCOMP      =     0
  integer, parameter :: SD_FILL                  =     0 
  integer, parameter :: SD_NOFILL                =   256

  integer, parameter :: HDF_VDATA                =    -1

  !----------------------
  ! Standard return codes       
  !----------------------

  integer, parameter :: SUCCEED                  =     0 
  integer, parameter :: FAIL     	         =    -1

  !------------------
  ! Compression Types 
  !------------------

  integer, parameter :: COMP_NONE                =     0
  integer, parameter :: COMP_RLE                 =    11
  integer, parameter :: COMP_IMCOMP              =    12 
  integer, parameter :: COMP_JPEG                =     2

  !--------------------------------------------------------------------
  ! Fortran chunking (SD and GR interfaces and compression routines use
  ! the following compression types:
  !--------------------------------------------------------------------

  integer, parameter :: COMP_CODE_NONE           =     0
  integer, parameter :: COMP_CODE_RLE            =     1
  integer, parameter :: COMP_CODE_NBIT           =     2
  integer, parameter :: COMP_CODE_SKPHUFF        =     3
  integer, parameter :: COMP_CODE_DEFLATE        =     4

  !----------------
  ! Interlace Types 
  !----------------

  integer, parameter :: MFGR_INTERLACE_PIXEL     =     0
  integer, parameter :: MFGR_INTERLACE_LINE      =     1
  integer, parameter :: MFGR_INTERLACE_COMPONENT =     2

  integer, parameter :: FULL_INTERLACE	         =     0 
  integer, parameter :: NO_INTERLACE             =     1

  !---------------------------
  ! Vdata fields packing types
  !---------------------------

  integer, parameter :: HDF_VSPACK               =     0
  integer, parameter :: HDF_VSUNPACK             =     1

  !----------------------------
  ! Multi-file Annotation types
  !----------------------------

  integer, parameter :: AN_DATA_LABEL            =     0 
  integer, parameter :: AN_DATA_DESC             =     1
  integer, parameter :: AN_FILE_LABEL            =     2
  integer, parameter :: AN_FILE_DESC             =     3

contains ! ============================= Local wrappers ======================

  integer function sfrdata_i1 (sds_id, start, stride, edges, data)      
    integer, intent (in) :: sds_id
    integer, intent (in), dimension (:) :: start, stride, edges
    integer, intent (out), dimension (:) :: data
    
    sfrdata_i1 = sfrdata ( sds_id, start, stride, edges, data )
  end function sfrdata_i1

  integer function sfrdata_i2 (sds_id, start, stride, edges, data)      
    integer, intent (in) :: sds_id
    integer, intent (in), dimension (:) :: start, stride, edges
    integer, intent (out), dimension (:,:) :: data

    integer, dimension(size(data)) :: tmp
    
    sfrdata_i2 = sfrdata ( sds_id, start, stride, edges, tmp )
    data = reshape ( tmp, shape(data) )
  end function sfrdata_i2

  integer function sfrdata_i3 (sds_id, start, stride, edges, data)      
    integer, intent (in) :: sds_id
    integer, intent (in), dimension (:) :: start, stride, edges
    integer, intent (out), dimension (:,:,:) :: data

    integer, dimension(size(data)) :: tmp
    
    sfrdata_i3 = sfrdata ( sds_id, start, stride, edges, tmp )
    data = reshape ( tmp, shape(data) )
  end function sfrdata_i3

  integer function sfrdata_r1 (sds_id, start, stride, edges, data)      
    integer, intent (in) :: sds_id
    integer, intent (in), dimension (:) :: start, stride, edges
    real(r4), intent (out), dimension (:) :: data
    
    sfrdata_r1 = sfrdata ( sds_id, start, stride, edges, data )
  end function sfrdata_r1

  integer function sfrdata_r2 (sds_id, start, stride, edges, data)      
    integer, intent (in) :: sds_id
    integer, intent (in), dimension (:) :: start, stride, edges
    real(r4), intent (out), dimension (:,:) :: data

    real(r4), dimension(size(data)) :: tmp
    
    sfrdata_r2 = sfrdata ( sds_id, start, stride, edges, tmp )
    data = reshape ( tmp, shape(data) )
  end function sfrdata_r2

  integer function sfrdata_r3 (sds_id, start, stride, edges, data)      
    integer, intent (in) :: sds_id
    integer, intent (in), dimension (:) :: start, stride, edges
    real(r4), intent (out), dimension (:,:,:) :: data

    real(r4), dimension(size(data)) :: tmp
    
    sfrdata_r3 = sfrdata ( sds_id, start, stride, edges, tmp )
    data = reshape ( tmp, shape(data) )
  end function sfrdata_r3

  integer function sfrdata_dp1 (sds_id, start, stride, edges, data)      
    integer, intent (in) :: sds_id
    integer, intent (in), dimension (:) :: start, stride, edges
    real(r8), intent (out), dimension (:) :: data
    
    sfrdata_dp1 = sfrdata ( sds_id, start, stride, edges, data )
  end function sfrdata_dp1

  integer function sfrdata_dp2 (sds_id, start, stride, edges, data)      
    integer, intent (in) :: sds_id
    integer, intent (in), dimension (:) :: start, stride, edges
    real(r8), intent (out), dimension (:,:) :: data

    real(r8), dimension(size(data)) :: tmp
    
    sfrdata_dp2 = sfrdata ( sds_id, start, stride, edges, tmp )
    data = reshape ( tmp, shape(data) )
  end function sfrdata_dp2

  integer function sfrdata_dp3 (sds_id, start, stride, edges, data)      
    integer, intent (in) :: sds_id
    integer, intent (in), dimension (:) :: start, stride, edges
    real(r8), intent (out), dimension (:,:,:) :: data

    real(r8), dimension(size(data)) :: tmp
    
    sfrdata_dp3 = sfrdata ( sds_id, start, stride, edges, tmp )
    data = reshape ( tmp, shape(data) )
  end function sfrdata_dp3

  integer function sfwdata_i1 (sds_id, start, stride, edges, data)      
    integer, intent (in) :: sds_id
    integer, intent (in), dimension (:) :: start, stride, edges
    integer, intent (in), dimension (:) :: data
    
    sfwdata_i1 = sfwdata ( sds_id, start, stride, edges, data )
  end function sfwdata_i1

  integer function sfwdata_r1 (sds_id, start, stride, edges, data)      
    integer, intent (in) :: sds_id
    integer, intent (in), dimension (:) :: start, stride, edges
    real, intent (in), dimension (:) :: data
    
    sfwdata_r1 = sfwdata ( sds_id, start, stride, edges, data )
  end function sfwdata_r1

  integer function sfwdata_dp1 (sds_id, start, stride, edges, data)      
    integer, intent (in) :: sds_id
    integer, intent (in), dimension (:) :: start, stride, edges
    double precision, intent (in), dimension (:) :: data
    
    sfwdata_dp1 = sfwdata ( sds_id, start, stride, edges, data )
  end function sfwdata_dp1

  integer function sfwdata_dp3 (sds_id, start, stride, edges, data)      
    integer, intent (in) :: sds_id
    integer, intent (in), dimension (:) :: start, stride, edges
    double precision, intent (in), dimension (:,:,:) :: data
    
    sfwdata_dp3 = sfwdata ( sds_id, start, stride, edges, data )
  end function sfwdata_dp3

  integer function sfwdata_r3 (sds_id, start, stride, edges, data)      
    integer, intent (in) :: sds_id
    integer, intent (in), dimension (:) :: start, stride, edges
    real, intent (in), dimension (:,:,:) :: data
    
    sfwdata_r3 = sfwdata ( sds_id, start, stride, edges, data )
  end function sfwdata_r3

end module Hdf

! $Log$
! Revision 2.7  2001/11/01 21:00:25  pwagner
! Replaced sfwdata generic with sfwdata_f90
!
! Revision 2.6  2001/05/30 23:51:35  livesey
! New version tidied up an some wrappers written.
!
! Revision 2.5  2001/03/06 22:41:34  livesey
! Changed from assumed size to assumed shape
!
! Revision 2.4  2001/01/03 00:45:27  pwagner
! Changed sfgetinfo to sfginfo
!
! Revision 2.3  2000/12/04 21:50:24  pwagner
! Added sfgdinfo
!
! Revision 2.2  2000/12/02 00:05:22  pwagner
! Added sfginfo, sfrdata, sfn2index, sfselect
!
! Revision 2.1  2000/11/29 19:33:18  pumphrey
! Added EXTERNAL reference for sfsfill
!
! Revision 2.0  2000/09/05 17:41:06  dcuddy
! Change revision to 2.0
!
! Revision 1.10  2000/06/30 01:30:26  lungu
!  Replaced Real (r8) by DOUBLE PRECISION
! to make the NAG compiler happy.
!

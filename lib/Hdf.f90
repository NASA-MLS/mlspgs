
! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===============================================================================
MODULE Hdf
!===============================================================================
   USE MLSCommon
   IMPLICIT NONE
   PUBLIC 		! Everything is public.

   PRIVATE :: Id	! Except Id
 
!------------------- RCS Ident Info -----------------------
   CHARACTER(LEN=130) :: Id = &        
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

!PAGE_BREAK

! Now define f90 interfaces for some toolkit routines.

   INTERFACE
      INTEGER FUNCTION sfstart(filename, access_mode)
        CHARACTER (LEN=*), INTENT (IN) :: filename
        INTEGER, INTENT (IN) :: access_mode
      END FUNCTION sfstart
      INTEGER FUNCTION sfcreate(sd_id, sds_name, data_type,rank, dim_sizes)
        CHARACTER (LEN=*), INTENT (IN) :: sds_name
        INTEGER, INTENT (IN) :: sd_id, data_type, rank
        INTEGER, INTENT (IN) :: dim_sizes(*)
      END FUNCTION sfcreate
      INTEGER FUNCTION sfdimid(sds_id, dim_index)
        INTEGER, INTENT (IN) :: sds_id, dim_index
      END FUNCTION sfdimid
      INTEGER FUNCTION sfsdmname(dim_id, dim_name)
        CHARACTER (LEN=*), INTENT (IN) :: dim_name
        INTEGER, INTENT (IN) :: dim_id
      END FUNCTION sfsdmname   
      INTEGER FUNCTION sfwdata(sds_id, start, stride, edges , data)      
        INTEGER, INTENT (IN) :: sds_id
        INTEGER, INTENT (IN), DIMENSION (:) :: start, stride, edges
        DOUBLE PRECISION, INTENT (IN), DIMENSION (:,:,:) :: data
      END FUNCTION sfwdata
      INTEGER FUNCTION sfendacc(sds_id)
        INTEGER, INTENT (IN) :: sds_id
      END FUNCTION sfendacc     
      INTEGER FUNCTION sfend(sd_id)
        INTEGER, INTENT (IN) :: sd_id
      END FUNCTION sfend 
      INTEGER FUNCTION sfrdata(sds_id, start, stride, edges , data)      
        INTEGER, INTENT (IN) :: sds_id
        INTEGER, INTENT (IN), DIMENSION (:) :: start, stride, edges
        DOUBLE PRECISION, INTENT (OUT), DIMENSION (:,:,:) :: data
      END FUNCTION sfrdata
      INTEGER FUNCTION sfn2index(sd_id, sds_name)      
        INTEGER, INTENT (IN) :: sd_id
        CHARACTER (LEN=*), INTENT (IN) :: sds_name
      END FUNCTION sfn2index
      INTEGER FUNCTION sfselect(sd_id, sds_index)      
        INTEGER, INTENT (IN) :: sd_id
        INTEGER, INTENT (IN) :: sds_index
      END FUNCTION sfselect
      INTEGER FUNCTION sfgetinfo(sds_id, sds_name, rank, dimsizes, data_type, &
      & num_attrs)      
        INTEGER, INTENT (IN) :: sds_id
        INTEGER, INTENT (OUT) :: rank, data_type, num_attrs
        CHARACTER (LEN=*), INTENT (OUT) :: sds_name
        INTEGER, INTENT (IN), DIMENSION (:) :: dimsizes
      END FUNCTION sfgetinfo
      INTEGER FUNCTION sfgdinfo(dim_id, dim_name, size, data_type, &
      & num_attrs)      
        INTEGER, INTENT (IN) :: dim_id
        INTEGER, INTENT (OUT) :: size, data_type, num_attrs
        CHARACTER (LEN=*), INTENT (OUT) :: dim_name
      END FUNCTION sfgdinfo

   END INTERFACE

   ! For routines where one of the arguments can be any valid data type !
   ! an interface block causes trouble. One either builds a generic
   ! interface or plays fast-and-loose by declaring the routine to be
   ! "external" in the old-fashioned Fortran-77 way. We do the latter for 
   ! sfsfill
   INTEGER::sfsfill
   EXTERNAL::sfsfill

!-------------------
! Error Return Codes 
!-------------------

   INTEGER, PARAMETER :: DFE_NOERROR              =     0
   INTEGER, PARAMETER :: DFE_NONE                 =     0
   INTEGER, PARAMETER :: DFE_FNF                  =    -1
   INTEGER, PARAMETER :: DFE_DENIED               =    -2
   INTEGER, PARAMETER :: DFE_ALROPEN              =    -3
   INTEGER, PARAMETER :: DFE_TOOMANY              =    -4
   INTEGER, PARAMETER :: DFE_BADNAME              =    -5
   INTEGER, PARAMETER :: DFE_BADACC               =    -6
   INTEGER, PARAMETER :: DFE_BADOPEN              =    -7
   INTEGER, PARAMETER :: DFE_NOTOPEN              =    -8
   INTEGER, PARAMETER :: DFE_CANTCLOSE            =    -9
   INTEGER, PARAMETER :: DFE_DFNULL               =   -10
   INTEGER, PARAMETER :: DFE_ILLTYPE              =   -11
   INTEGER, PARAMETER :: DFE_UNSUPPORTED          =   -12
   INTEGER, PARAMETER :: DFE_BADDDLIST            =   -13
   INTEGER, PARAMETER :: DFE_NOTDFFILE            =   -14
   INTEGER, PARAMETER :: DFE_SEEDTWICE            =   -15
   INTEGER, PARAMETER :: DFE_NOSPACE              =   -16
   INTEGER, PARAMETER :: DFE_NOSUCHTAG            =   -17
   INTEGER, PARAMETER :: DFE_READERROR            =   -18

   INTEGER, PARAMETER :: DFE_WRITEERROR           =   -19
   INTEGER, PARAMETER :: DFE_SEEKERROR            =   -20
   INTEGER, PARAMETER :: DFE_NOFREEDD             =   -21
   INTEGER, PARAMETER :: DFE_BADTAG               =   -22
   INTEGER, PARAMETER :: DFE_BADREF               =   -23
   INTEGER, PARAMETER :: DFE_RDONLY               =   -24
   INTEGER, PARAMETER :: DFE_BADCALL              =   -25
   INTEGER, PARAMETER :: DFE_BADPTR               =   -26
   INTEGER, PARAMETER :: DFE_BADLEN               =   -27
   INTEGER, PARAMETER :: DFE_BADSEEK              =   -28
   INTEGER, PARAMETER :: DFE_NOMATCH              =   -29
   INTEGER, PARAMETER :: DFE_NOTINSET             =   -30
   INTEGER, PARAMETER :: DFE_BADDIM               =   -31
   INTEGER, PARAMETER :: DFE_BADOFFSET            =   -32
   INTEGER, PARAMETER :: DFE_BADSCHEME            =   -33
   INTEGER, PARAMETER :: DFE_NODIM                =   -34
   INTEGER, PARAMETER :: DFE_NOTENOUGH            =   -35
   INTEGER, PARAMETER :: DFE_NOVALS               =   -36
   INTEGER, PARAMETER :: DFE_CORRUPT              =   -37
   INTEGER, PARAMETER :: DFE_BADFP                =   -38

!PAGE_BREAK
   INTEGER, PARAMETER :: DFE_NOREF                =   -39
   INTEGER, PARAMETER :: DFE_BADDATATYPE          =   -40
   INTEGER, PARAMETER :: DFE_BADMCTYPE            =   -41
   INTEGER, PARAMETER :: DFE_BADNUMTYPE           =   -42
   INTEGER, PARAMETER :: DFE_BADORDER             =   -43
   INTEGER, PARAMETER :: DFE_ARGS                 =   -44
   INTEGER, PARAMETER :: DFE_INTERNAL             =   -45
   INTEGER, PARAMETER :: DFE_DUPDD                =   -46
   INTEGER, PARAMETER :: DFE_CANTMOD              =   -47
   INTEGER, PARAMETER :: DFE_RANGE                =   -48
   INTEGER, PARAMETER :: DFE_BADTABLE             =   -49
   INTEGER, PARAMETER :: DFE_BADSDG               =   -50
   INTEGER, PARAMETER :: DFE_BADNDG               =   -51
   INTEGER, PARAMETER :: DFE_BADFIELDS            =   -52
   INTEGER, PARAMETER :: DFE_NORESET              =   -53
   INTEGER, PARAMETER :: DFE_NOVS                 =   -54
   INTEGER, PARAMETER :: DFE_VGSIZE               =   -55
   INTEGER, PARAMETER :: DFE_DIFFFILES            =   -56
   INTEGER, PARAMETER :: DFE_VTAB                 =   -57
   INTEGER, PARAMETER :: DFE_BADAID               =   -58

   INTEGER, PARAMETER :: DFE_OPENAID              =   -59
   INTEGER, PARAMETER :: DFE_BADCONV              =   -60
   INTEGER, PARAMETER :: DFE_GENAPP               =   -61
   INTEGER, PARAMETER :: DFE_CANTFLUSH            =   -62
   INTEGER, PARAMETER :: DFE_BADTYPE              =   -63
   INTEGER, PARAMETER :: DFE_SYMSIZE              =   -64
   INTEGER, PARAMETER :: DFE_BADATTACH            =   -65
   INTEGER, PARAMETER :: DFE_CANTDETACH           =   -66

!---------------------------
! internal file access codes
!---------------------------

   INTEGER, PARAMETER :: DFACC_READ               =     1
   INTEGER, PARAMETER :: DFACC_WRITE              =     2
   INTEGER, PARAMETER :: DFACC_CREATE             =     4
   INTEGER, PARAMETER :: DFACC_ALL                =     7
   INTEGER, PARAMETER :: DFACC_RDONLY             =     1
   INTEGER, PARAMETER :: DFACC_RDWR               =     3
   INTEGER, PARAMETER :: DFACC_CLOBBER            =     4

!PAGE_BREAK
!---------------------------------
! Access types for SDsetaccesstype
!---------------------------------

   INTEGER, PARAMETER :: DFACC_DEFAULT            =     0
   INTEGER, PARAMETER :: DFACC_SERIAL             =     1
   INTEGER, PARAMETER :: DFACC_PARALLEL           =     9

!---------------------------
! Constants for DFSDsetorder
!---------------------------

   INTEGER, PARAMETER :: DFO_FORTRAN              =     1
   INTEGER, PARAMETER :: DFO_C                    =     2

!----------------------------------
! Definitions of storage convention
!----------------------------------

   INTEGER, PARAMETER :: DFNTF_IEEE               =     1
   INTEGER, PARAMETER :: DFNTF_VAX                =     2
   INTEGER, PARAMETER :: DFNTF_CRAY               =     3
   INTEGER, PARAMETER :: DFNTF_PC                 =     4
   INTEGER, PARAMETER :: DFNTF_CONVEX             =     5
   INTEGER, PARAMETER :: DFNTF_VP                 =     6

!----------------
! Masks for types
!----------------

   INTEGER, PARAMETER :: DFNT_HDF                 =     0
   INTEGER, PARAMETER :: DFNT_NATIVE              =  4096
   INTEGER, PARAMETER :: DFNT_CUSTOM              =  8192
   INTEGER, PARAMETER :: DFNT_LITEND              = 16384

!PAGE_BREAK
!-----------------------
! Number type info codes 
!-----------------------

   INTEGER, PARAMETER :: DFNT_NONE                =     0
   INTEGER, PARAMETER :: DFNT_QUERY               =     0
   INTEGER, PARAMETER :: DFNT_VERSION             =     1
      
   INTEGER, PARAMETER :: DFNT_FLOAT32             =     5
   INTEGER, PARAMETER :: DFNT_FLOAT               =     5
   INTEGER, PARAMETER :: DFNT_FLOAT64             =     6
   INTEGER, PARAMETER :: DFNT_DOUBLE              =     6
   INTEGER, PARAMETER :: DFNT_FLOAT128            =     7

   INTEGER, PARAMETER :: DFNT_INT8                =    20
   INTEGER, PARAMETER :: DFNT_UINT8               =    21
   INTEGER, PARAMETER :: DFNT_INT16               =    22
   INTEGER, PARAMETER :: DFNT_UINT16              =    23
   INTEGER, PARAMETER :: DFNT_INT32               =    24
   INTEGER, PARAMETER :: DFNT_UINT32              =    25
   INTEGER, PARAMETER :: DFNT_INT64               =    26
   INTEGER, PARAMETER :: DFNT_UINT64              =    27
   INTEGER, PARAMETER :: DFNT_INT128              =    28
   INTEGER, PARAMETER :: DFNT_UINT128             =    29

   INTEGER, PARAMETER :: DFNT_UCHAR8              =     3
   INTEGER, PARAMETER :: DFNT_UCHAR               =     3
   INTEGER, PARAMETER :: DFNT_CHAR8               =     4
   INTEGER, PARAMETER :: DFNT_CHAR                =     4
   INTEGER, PARAMETER :: DFNT_CHAR16              =    42
   INTEGER, PARAMETER :: DFNT_UCHAR16             =    43

   INTEGER, PARAMETER :: DFNT_NFLOAT32            =  4101
   INTEGER, PARAMETER :: DFNT_NFLOAT              =  4101
   INTEGER, PARAMETER :: DFNT_NFLOAT64            =  4102
   INTEGER, PARAMETER :: DFNT_NDOUBLE             =  4102
   INTEGER, PARAMETER :: DFNT_NFLOAT128           =  4103

!PAGE_BREAK
   INTEGER, PARAMETER :: DFNT_NINT8               =  4116
   INTEGER, PARAMETER :: DFNT_NUINT8              =  4117
   INTEGER, PARAMETER :: DFNT_NINT16              =  4118
   INTEGER, PARAMETER :: DFNT_NUINT16             =  4119
   INTEGER, PARAMETER :: DFNT_NINT32	          =  4120
   INTEGER, PARAMETER :: DFNT_NUINT32             =  4121
   INTEGER, PARAMETER :: DFNT_NINT64              =  4122
   INTEGER, PARAMETER :: DFNT_NUINT64             =  4123
   INTEGER, PARAMETER :: DFNT_NINT128             =  4124
   INTEGER, PARAMETER :: DFNT_NUINT128            =  4125

   INTEGER, PARAMETER :: DFNT_NUCHAR8             =  4099
   INTEGER, PARAMETER :: DFNT_NUCHAR              =  4099
   INTEGER, PARAMETER :: DFNT_NCHAR8              =  4100
   INTEGER, PARAMETER :: DFNT_NCHAR               =  4100
   INTEGER, PARAMETER :: DFNT_NCHAR16             =  4138
   INTEGER, PARAMETER :: DFNT_NUCHAR16            =  4139

   INTEGER, PARAMETER :: DFNT_LFLOAT32            = 16389
   INTEGER, PARAMETER :: DFNT_LFLOAT              = 16389
   INTEGER, PARAMETER :: DFNT_LFLOAT64            = 16390
   INTEGER, PARAMETER :: DFNT_LDOUBLE             = 16390
   INTEGER, PARAMETER :: DFNT_LFLOAT128           = 16391

   INTEGER, PARAMETER :: DFNT_LINT8               = 16404
   INTEGER, PARAMETER :: DFNT_LUINT8              = 16405
   INTEGER, PARAMETER :: DFNT_LINT16              = 16406
   INTEGER, PARAMETER :: DFNT_LUINT16             = 16407
   INTEGER, PARAMETER :: DFNT_LINT32              = 16408
   INTEGER, PARAMETER :: DFNT_LUINT32             = 16409
   INTEGER, PARAMETER :: DFNT_LINT64              = 16410
   INTEGER, PARAMETER :: DFNT_LUINT64             = 16411
   INTEGER, PARAMETER :: DFNT_LINT128             = 16412
   INTEGER, PARAMETER :: DFNT_LUINT128            = 16413

   INTEGER, PARAMETER :: DFNT_LUCHAR8             = 16387
   INTEGER, PARAMETER :: DFNT_LUCHAR              = 16387
   INTEGER, PARAMETER :: DFNT_LCHAR8              = 16388
   INTEGER, PARAMETER :: DFNT_LCHAR               = 16388
   INTEGER, PARAMETER :: DFNT_LCHAR16             = 16426
   INTEGER, PARAMETER :: DFNT_LUCHAR16            = 16427

!PAGE_BREAK
!--------------
! tags and refs
!--------------

   INTEGER, PARAMETER :: DFREF_WILDCARD           =     0
   INTEGER, PARAMETER :: DFTAG_WILDCARD           =     0
   INTEGER, PARAMETER :: DFTAG_NULL               =     1
   INTEGER, PARAMETER :: DFTAG_LINKED             =    20
   INTEGER, PARAMETER :: DFTAG_VERSION            =    30
   INTEGER, PARAMETER :: DFTAG_COMPRESSED         =    40

!------------
! utility set
!------------

   INTEGER, PARAMETER :: DFTAG_FID                =   100
   INTEGER, PARAMETER :: DFTAG_FD                 =   101
   INTEGER, PARAMETER :: DFTAG_TID                =   102
   INTEGER, PARAMETER :: DFTAG_TD                 =   103
   INTEGER, PARAMETER :: DFTAG_DIL                =   104
   INTEGER, PARAMETER :: DFTAG_DIA                =   105
   INTEGER, PARAMETER :: DFTAG_NT                 =   106
   INTEGER, PARAMETER :: DFTAG_MT                 =   107

!PAGE_BREAK
!-------------
! raster-8 set 
!-------------

   INTEGER, PARAMETER :: DFTAG_ID8                =   200
   INTEGER, PARAMETER :: DFTAG_IP8                =   201
   INTEGER, PARAMETER :: DFTAG_RI8                =   202 
   INTEGER, PARAMETER :: DFTAG_CI8                =   203
   INTEGER, PARAMETER :: DFTAG_II8                =   204

!-----------------
! Raster Image set
!-----------------

   INTEGER, PARAMETER :: DFTAG_ID                 =   300 
   INTEGER, PARAMETER :: DFTAG_LUT                =   301
   INTEGER, PARAMETER :: DFTAG_RI                 =   302 
   INTEGER, PARAMETER :: DFTAG_CI                 =   303

   INTEGER, PARAMETER :: DFTAG_RIG                =   306 
   INTEGER, PARAMETER :: DFTAG_LD                 =   307
   INTEGER, PARAMETER :: DFTAG_MD                 =   308 
   INTEGER, PARAMETER :: DFTAG_MA                 =   309
   INTEGER, PARAMETER :: DFTAG_CCN                =   310 
   INTEGER, PARAMETER :: DFTAG_CFM                =   311
   INTEGER, PARAMETER :: DFTAG_AR                 =   312

   INTEGER, PARAMETER :: DFTAG_DRAW               =   400 
   INTEGER, PARAMETER :: DFTAG_RUN                =   401
   INTEGER, PARAMETER :: DFTAG_XYP                =   500 
   INTEGER, PARAMETER :: DFTAG_MTO                =   501

!----------
! Tektronix 
!----------

   INTEGER, PARAMETER :: DFTAG_T14                =   602 
   INTEGER, PARAMETER :: DFTAG_T105               =   603

!PAGE_BREAK
!--------------------
! Scientific Data set 
!--------------------

   INTEGER, PARAMETER :: DFTAG_SDG                =   700 
   INTEGER, PARAMETER :: DFTAG_SDD                =   701
   INTEGER, PARAMETER :: DFTAG_SD                 =   702 
   INTEGER, PARAMETER :: DFTAG_SDS                =   703
   INTEGER, PARAMETER :: DFTAG_SDL                =   704 
   INTEGER, PARAMETER :: DFTAG_SDU                =   705
   INTEGER, PARAMETER :: DFTAG_SDF                =   706 
   INTEGER, PARAMETER :: DFTAG_SDM                =   707
   INTEGER, PARAMETER :: DFTAG_SDC                =   708 
   INTEGER, PARAMETER :: DFTAG_SDT                =   709
   INTEGER, PARAMETER :: DFTAG_SDLNK              =   710 
   INTEGER, PARAMETER :: DFTAG_NDG                =   720
   INTEGER, PARAMETER :: DFTAG_CAL                =   731 
   INTEGER, PARAMETER :: DFTAG_FV                 =   732
   INTEGER, PARAMETER :: DFTAG_BREQ               =   799
   INTEGER, PARAMETER :: DFTAG_EREQ               =   780

!------
! VSets 
!------

   INTEGER, PARAMETER :: DFTAG_VG                 =  1965 
   INTEGER, PARAMETER :: DFTAG_VH                 =  1962
   INTEGER, PARAMETER :: DFTAG_VS                 =  1963

!--------------------
! compression schemes 
!--------------------

   INTEGER, PARAMETER :: DFTAG_RLE                =    11 
   INTEGER, PARAMETER :: DFTAG_IMC                =    12
   INTEGER, PARAMETER :: DFTAG_IMCOMP             =    12 
   INTEGER, PARAMETER :: DFTAG_JPEG               =    13
   INTEGER, PARAMETER :: DFTAG_GREYJPEG           =    14

!--------------
! SPECIAL CODES 
!--------------

   INTEGER, PARAMETER :: SPECIAL_LINKED 	 =     1
   INTEGER, PARAMETER :: SPECIAL_EXT 	         =     2

!PAGE_BREAK
!-----------
! PARAMETERS 
!-----------

   INTEGER, PARAMETER :: DF_MAXFNLEN              =   256 
   INTEGER, PARAMETER :: SD_UNLIMITED             =     0
   INTEGER, PARAMETER :: SD_DIMVAL_BW_COMP        =     1 
   INTEGER, PARAMETER :: SD_DIMVAL_BW_INCOMP      =     0
   INTEGER, PARAMETER :: SD_FILL                  =     0 
   INTEGER, PARAMETER :: SD_NOFILL                =   256

   INTEGER, PARAMETER :: HDF_VDATA                =    -1

!----------------------
! Standard return codes       
!----------------------
 
   INTEGER, PARAMETER :: SUCCEED                  =     0 
   INTEGER, PARAMETER :: FAIL     	          =    -1

!------------------
! Compression Types 
!------------------

   INTEGER, PARAMETER :: COMP_NONE                =     0
   INTEGER, PARAMETER :: COMP_RLE                 =    11
   INTEGER, PARAMETER :: COMP_IMCOMP              =    12 
   INTEGER, PARAMETER :: COMP_JPEG                =     2

!--------------------------------------------------------------------
! Fortran chunking (SD and GR interfaces and compression routines use
! the following compression types:
!--------------------------------------------------------------------
  
   INTEGER, PARAMETER :: COMP_CODE_NONE           =     0
   INTEGER, PARAMETER :: COMP_CODE_RLE            =     1
   INTEGER, PARAMETER :: COMP_CODE_NBIT           =     2
   INTEGER, PARAMETER :: COMP_CODE_SKPHUFF        =     3
   INTEGER, PARAMETER :: COMP_CODE_DEFLATE        =     4

!PAGE_BREAK
!----------------
! Interlace Types 
!----------------

   INTEGER, PARAMETER :: MFGR_INTERLACE_PIXEL     =     0
   INTEGER, PARAMETER :: MFGR_INTERLACE_LINE      =     1
   INTEGER, PARAMETER :: MFGR_INTERLACE_COMPONENT =     2

   INTEGER, PARAMETER :: FULL_INTERLACE	          =     0 
   INTEGER, PARAMETER :: NO_INTERLACE             =     1

!---------------------------
! Vdata fields packing types
!---------------------------
  
   INTEGER, PARAMETER :: HDF_VSPACK               =     0
   INTEGER, PARAMETER :: HDF_VSUNPACK             =     1

!----------------------------
! Multi-file Annotation types
!----------------------------

   INTEGER, PARAMETER :: AN_DATA_LABEL            =     0 
   INTEGER, PARAMETER :: AN_DATA_DESC             =     1
   INTEGER, PARAMETER :: AN_FILE_LABEL            =     2
   INTEGER, PARAMETER :: AN_FILE_DESC             =     3

!===============================================================================
END MODULE Hdf
!===============================================================================

!# $Log$
!# Revision 2.3  2000/12/04 21:50:24  pwagner
!# Added sfgdinfo
!#
!# Revision 2.2  2000/12/02 00:05:22  pwagner
!# Added sfgetinfo, sfrdata, sfn2index, sfselect
!#
!# Revision 2.1  2000/11/29 19:33:18  pumphrey
!# Added EXTERNAL reference for sfsfill
!#
!# Revision 2.0  2000/09/05 17:41:06  dcuddy
!# Change revision to 2.0
!#
!# Revision 1.10  2000/06/30 01:30:26  lungu
!#  Replaced Real (r8) by DOUBLE PRECISION
!# to make the NAG compiler happy.
!#

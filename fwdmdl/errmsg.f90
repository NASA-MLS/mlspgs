! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!
!-----------------------------------------------------------------------
!=====================  For NAG Fortran 95 only =====================
!
 
SUBROUTINE lu_errmsg(lu,msg,ier)
 
!
!  Lu - unit number to write into (usually 6)
!  msg - the User's message if he/she wishes. It could be empty.
!  ier  the iostat or STAT integer flag (error number)
!
Implicit NONE
!
 
INTEGER, INTENT(IN)            :: lu, ier
CHARACTER (LEN=*), INTENT(IN)  :: msg
 
INTEGER, PARAMETER :: numerr = 215
 
INTEGER  :: i,j,errno(numerr)
CHARACTER (LEN=80) :: a,sysmsg(numerr)
CHARACTER (LEN=8) :: ax
 
DATA (errno(i),i=1,numerr)/ -2, -1,                                 &
     1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13, &
    14,  15,  16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26, &
    27,  28,  29,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39, &
    40,  41,  42,  43,  44,  45,  46,  47,  48,  49,  50,  51,  52, &
    53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,  64,  65, &
    66,  67,  68,  69,  70,  71,  72,  74,  77,  78,  79,  80,  81, &
    82,  83,  84,  85,  86,  87,  88,  89,  90,  91,  92,  93,  94, &
    95,  96,  97,  98,  99, 100, 101, 102, 103, 104, 105, 106, 107, &
   108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, &
   121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, &
   134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, &
   147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, &
   160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, &
   173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, &
   186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, &
   199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, &
   212, 213, 214, 215, 216 /
 
DATA (sysmsg(i),i=1,91)/                                       &
   "Trying to read past End-Of-Record",                            &
   "Trying to read past End-Of-File",                              &
   "Not super-user","No such file or directory","No such process", &
   "interrupted system call","I/O error","No such device or address", &
   "Arg list too long","Exec format error","Bad file number", &
   "No children","Resource temporarily unavailable","Not enough core", &
   "Permission denied","Bad address","Block device required", &
   "Mount device busy","File exists","Cross-device link","No such device", &
   "Not a directory","Is a directory","Invalid argument", &
   "File table overflow","Too many open files", &
   "Inappropriate ioctl for device","Text file busy","File too large", &
   "No space left on device","Illegal seek","Read only file system", &
   "Too many links","Broken pipe","Math arg out of domain of func", &
   "Math result not representable","No message of desired type", &
   "Identifier removed","Channel number out of range", &
   "Level 2 not synchronized","Level 3 halted","Level 3 reset", &
   "Link number out of range","Protocol driver not attached", &
   "No CSI structure available","Level 2 halted","Deadlock condition", &
   "No record locks available","Operation canceled", &
   "Operation not supported","Disc quota exceeded","invalid exchange", &
   "invalid request descriptor","exchange full","no anode", &
   "invalid request code","invalid slot","file locking deadlock error", &
   "bad font file fmt","process died with the lock", &
   "lock is not recoverable","Device not a stream", &
   "no data (for no delay io)","timer expired", &
   "out of streams resources","Machine is not on the network", &
   "Package not installed","The object is remote", &
   "the link has been severed","advertise error","srmount error", &
   "Communication error on send","Protocol error", &
   "locked lock was unmapped","multihop attempted", &
   "trying to read unreadable message","path name is too long", &
   "value too large to be stored in data type","given log. name not unique", &
   "f.d. invalid for this operation","Remote address changed", &
   "Can't access a needed shared lib","Accessing a corrupted shared lib", &
   ".lib section in a.out corrupted","Attempting to link in too many libs", &
   "Attempting to exec a shared library","Illegal byte sequence", &
   "Unsupported file system operation","Symbolic link loop", &
   "Restartable system call","if pipe/FIFO, don't sleep in stream head" /

DATA (sysmsg(i),i=92,136)/                                       &
   "directory not empty","Too many users (for UFS)", &
   "Socket operation on non-socket","Destination address required", &
   "Message too long","Protocol wrong type for socket", &
   "Protocol not available","Buffer overflow on output", &
   "Internal file overflow","Scale factor out of range", &
   "Exponent too large for w.d format","Record too long for input buffer", &
   "Zero repeat factor in list-directed input", &
   "Invalid input for integer editing", &
   "Input value too large for INTEGER(KIND=1)", &
   "Input value too large for INTEGER(KIND=2)", &
   "Repeat factor in list-directed input larger than HUGE(0)", &
   "Input value too large for default INTEGER type", &
   "Invalid input for real editing", &
   "Invalid input for logical editing", &
   "Invalid input for complex editing", &
   "Invalid input for character editing", &
   "Format specification does not begin with a left parenthesis", &
   "Format specification does not end with a right parenthesis", &
   "No data edit descriptor in tail of format specification after reversion", &
   "Sub-format groups nested too deeply", &
   "Unexpected end of format specification", &
   "Expected integer literal constant in format specification", &
   "Field/exponent width or repeat in format specification must be non-zero", &
   "Expected decimal point in format specification", &
   "Expected P following signed integer constant in format specification", &
   "Expected N or Z following B in format specification", &
   "Invalid edit descriptor","No edit descriptor following repeat factor", &
   "Repeat factor given for sign edit descriptor", &
   "Repeat factor given for blank-interpretation edit descriptor", &
   "Missing length of H edit descriptor", &
   "Repeat factor given for character string edit descriptor", &
   "No spacing specified for X edit descriptor", &
   "Repeat factor given for position edit descriptor", &
   "Character string edit descriptor used on input", &
   "Invalid edit descriptor for real i/o-list item", &
   "Invalid edit descriptor for integer i/o-list item", &
   "Invalid edit descriptor for logical i/o-list item", &
   "Character string edit descriptor does not terminate before format end" /

DATA (sysmsg(i),i=137,177)/                                       &
   "Sign in a numeric input field not followed by any digits", &
   "Invalid exponent in real input field", &
   "Invalid character in real input field", &
   "Invalid character in integer input field", &
   "Invalid character in binary integer input field", &
   "Invalid character in octal integer input field", &
   "Invalid character in hexadecimal integer input field", &
   "Invalid edit descriptor for character i/o-list item", &
   "READ after WRITE with no intervening file positioning", &
   "Unit number out of range","Unit is not connected", &
   "File connected to unit is not capable of BACKSPACE", &
   "Unit is not connected for SEQUENTIAL i/o", &
   "Unit is not connected for READ action", &
   "Unit is not connected for FORMATTED i/o", &
   "Unit is not connected for WRITE action", &
   "Unit is not connected for UNFORMATTED i/o", &
   "Unit is not connected and no FILE= specifier on OPEN", &
   "FILE= specifier on OPEN with STATUS='SCRATCH'", &
   "OPEN on connected unit has different STATUS= specifier", &
   "OPEN on connected unit has different ACCESS= specifier", &
   "OPEN on connected unit has different FORM= specifier", &
   "OPEN on connected unit has different RECL= specifier", &
   "OPEN on connected unit has different ACTION= specifier", &
   "OPEN on connected unit has different POSITION= specifier", &
   "Invalid value for STATUS= specifier", &
   "Invalid value for ACCESS= specifier", &
   "Invalid value for FORM= specifier", &
   "Invalid value for BLANKS= specifier", &
   "Invalid value for POSITION= specifier", &
   "INvalid value for ACTION= specifier", &
   "Invalid value for DELIM= specifier", &
   "Invalid value for PAD= specifier", &
   "The RECL= specifier must be given for DIRECT access OPEN", &
   "STATUS='KEEP' is invalid for a SCRATCH file", &
   "ENDFILE applied twice to unit with no intervening file positioning", &
   "File name too long","Cannot find OLD file","NEW file already exists", &
   "File connected to unit is not capable of REWIND", &
   "BACKSPACE failed to find the beginning of the previous record" /

DATA (sysmsg(i),i=178,numerr)/                                       &
   "Unit is not connected for DIRECT i/o", &
   "Record number out of range", &
   "Illegal character in LOGICAL input field", &
   "No value found in LOGICAL input field", &
   "Unknown OPEN failure on unit", &
   "Expected '&' but found s/g else in NAMELIST input", &
   "NAMELIST group name in input is too long", &
   "Expected NAMELIST group but found anther name", &
   "Invalid character in NAMELIST input", &
   "NAMELIST group name in input of NAMELIST is too long", &
   "Expected '=' but found s/g else in NAMELIST input", &
   "Unknown group object in input for NAMELIST", &
   "Unexpected subscript for object of NAMELIST", &
   "Unexpected component specifier for object of NAMELIST", &
   "Component name too long in input for object of NAMELIST", &
   "Unknown component in input for object of NAMELIST", &
   "Array component of array parent in input for object of NAMELIST", &
   "Invalid integer literal in input for object of NAMELIST", &
   "Expected ':' but found s/g else in input for object of NAMELIST", &
   "Substring has zero length in input for object of NAMELIST", &
   "Substring (%d:%d) out of bounds in input for object of NAMELIST", &
   "Expected ',' but found s/g else in input for object of NAMELIST", &
   "Expected ')' but found s/g else in input for object of NAMELIST", &
   "Subscript (%d) out of range in input for object of NAMELIST", &
   "Array section has zero size in input for object of NAMELIST", &
   "Section stride has value zero in input for object of NAMELIST", &
   "Missing namelist group name in input of NAMELIST", &
   "Input list bigger than record length in unformatted READ on unit", &
   "Record too short for format requirement and PAD='NO' on unit", &
   "Unformatted data file open on unit %d is corrupt", &
   "File truncation on unit failed", &
   "READ/WRITE attempted after ENDFILE on unit", &
   "Input value too large for 64-bit integer", &
   "No FILE= specifier with STATUS='REPLACE", &
   "Invalid value for RECL= specifier (must be positive)", &
   "READ beyond end of direct access file on unit", &
   "Floating overflow on real number input", &
   "Direct access is incompatible with the POSITION= specifier" /
 
i = Len_Trim(msg)
IF(i > 0) WRITE(lu,'(1x,A)') msg(1:i)
 
ax(1:)=' '
WRITE(ax,*) ier
ax = AdjustL(ax)
 
  j = 0
  a(1:)=' '
  do i = 1, numerr
    if(errno(i) == ier) THEN
      a = sysmsg(i)
      j = Len_Trim(a)
      EXIT
    endif
  end do
 
  if(j < 1) THEN
    a = '** Unknown system error **'
    j = Len_Trim(a)
  endif
 
  i = Len_Trim(ax)
  WRITE(lu,900) ax(1:i),a(1:j)
 
900  FORMAT(' ** Error code: ',a,' --> ',a)
 
RETURN
END SUBROUTINE lu_errmsg
 
!-----------------------------------------------------------------------
 
SUBROUTINE errmsg(msg,ier)
 
Implicit NONE
 
INTEGER, INTENT(IN)            :: ier
CHARACTER (LEN=*), INTENT(IN)  :: msg
 
INTEGER :: lu
 
lu = 6
CALL lu_errmsg(lu,msg,ier)
 
RETURN
END SUBROUTINE errmsg
 

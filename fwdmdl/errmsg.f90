!
!-----------------------------------------------------------------------
!=====================  For PC AbSoft Fortran 77 only =====================
!

SUBROUTINE lu_errmsg(lu,msg,ier)

  use STRINGS, only: LEFTJ
!

INTEGER, INTENT(IN OUT)                  :: lu
CHARACTER (LEN=*), INTENT(IN OUT)        :: msg
INTEGER, INTENT(IN)                      :: ier
INTEGER, PARAMETER :: numerr = 158


CHARACTER (LEN=60) :: a,sysmsg(numerr)
CHARACTER (LEN=8) :: ax

DATA (sysmsg(i),i=1,18)/ 'End-of-file encountered', &
    'Operation not permitted', 'No such file or directory', &
    'No such process', 'Interrupted system call', &
    'I/O error', 'No such device or address', &
    'Arg list too long', 'Exec format error', &
    'Bad file number', 'No child processes', &
    'Try again', 'Out of memory', &
    'Permission denied', 'Bad address', &
    'Block device required', 'Device or resource busy', &
    'File exists'/

DATA (sysmsg(i),i=19,36)/ 'Cross-device link', &
    'No such device', 'Not a directory', &
    'Is a directory', 'Invalid argument', &
    'File table overflow', 'Too many open files', &
    'Not a typewriter', 'Text file busy', &
    'File too large', 'No space left on device', &
    'Illegal seek', 'Read-only file system', &
    'Too many links', 'Broken pipe', &
    'Math argument out of domain of func', &
    'Math result not representable', &
    'Resource deadlock would occur'/

DATA (sysmsg(i),i=37,54)/ 'File name too long', &
    'No record locks available', 'Function not implemented', &
    'Directory not empty', 'Too many symbolic links encountered', &
    'Operation would block', 'No message of desired type', &
    'Identifier removed', 'Channel number out of range', &
    'Level 2 not synchronized', 'Level 3 halted', &
    'Level 3 reset', 'Link number out of range', &
    'Protocol driver not attached', 'No CSI structure available', &
    'Level 2 halted', 'Invalid exchange', &
    'Invalid request descriptor'/

DATA (sysmsg(i),i=55,72)/ 'Exchange full', &
    'No anode', 'Invalid request code', &
    'Invalid slot', 'Resource deadlock would occur', &
    'Bad font file format', 'Device not a stream', &
    'No data available', 'Timer expired', &
    'Out of streams resources', 'Machine is not on the network', &
    'Package not installed', 'Object is remote', &
    'Link has been severed', 'Advertise error', &
    'Srmount error', 'Communication error on send', &
    'Protocol error'/

DATA (sysmsg(i),i=73,90)/ 'Multihop attempted', &
    'RFS specific error', 'Not a data message', &
    'Value too large for defined data type', &
    'Name not unique on network', &
    'File descriptor in bad state', 'Remote address changed', &
    'Can not access a needed shared library', &
    'Accessing a corrupted shared library', &
    '.lib section in a.out corrupted', &
    'Attempting to link in too many shared libraries', &
    'Cannot exec a shared library directly', &
    'Illegal byte sequence', &
    'Interrupted system call should be restarted', &
    'Streams pipe error', &
    'Too many users', 'Socket operation on non-socket', &
    'Destination address required'/

DATA (sysmsg(i),i=91,108)/ 'Message too long', &
    'Protocol wrong type for socket', 'Protocol not available', &
    'Protocol not supported', 'Socket type not supported', &
    'Operation not supported on transport endpoint', &
    'Protocol family not supported', &
    'Address family not supported by protocol', &
    'Address already in use', &
    'Cannot assign requested address', 'Network is down', &
    'Network is unreachable', &
    'Network dropped connection because of reset', &
    'Software caused connection abort', &
    'Connection reset by peer', &
    'No buffer space available', &
    'Transport endpoint is already connected', &
    'Transport endpoint is not connected'/

DATA (sysmsg(i),i=109,125)/  &
    'Cannot send after transport endpoint shutdown', &
    'Too many references: cannot splice', 'Connection timed out', &
    'Connection refused', 'Host is down', &
    'No route to host', 'Operation already in progress', &
    'Operation now in progress', 'Stale NFS file handle', &
    'Structure needs cleaning', 'Not a XENIX named type file', &
    'No XENIX semaphores available', 'Is a named type file', &
    'Remote I/O error', 'Quota exceeded', &
    'No medium found', 'Wrong medium type'/

! The following are FORTRAN erros, not SYSTEM errors:

DATA (sysmsg(i),i=126,143)/ 'File not open for read', &
    'File not open for write', 'File not found', &
    'Record length negative or zero', 'Buffer allocation failed', &
    'Bad iolist specifier', 'Error in format string', &
    'Illegal repeat count', &
    'Hollerith count exceeded remaining format string', &
    'Format string missing opening "("', &
    'Format string has unmatched parens', &
    'Format string has unmatched quotes', &
    'Non-repeatable format descriptor', &
    'Attempt to read past end-of-file', 'Bad file specification', &
    'Format group table overflow', &
    'Illegal character in numeric input', &
    'No record specified for direct access'/

DATA (sysmsg(i),i=144,158)/ 'Maximum record number exceeded', &
    'Illegal file type for namelist directed I/O', &
    'Illegal input for namelist directed I/O', &
    'Variable not present in current namelist', &
    'Variable type or size does not match edit descriptor', &
    'Illegal direct access record number', &
    'Illegal use of internal file', &
    'RECL= only valid for direct access files', &
    'BLOCK= only valid for unformatted sequetial files', &
    'Unable to truncate file after rewind, backspace, or endfile', &
    'Cannot do formatted I/O on an entire structure', &
    'Illegal (negative) unit specified', &
    'Specification in re-open do not match previous open', &
    'No implicit OPEN for direct access files', &
    'Cannot open an existing file with STATUS="NEW"'/

i = Len(Trim(msg))
IF(i > 0) WRITE(lu,'(1x,A)') msg(1:i)

ax(1:)=' '
WRITE(ax,*) ier
CALL Leftj(ax)

jer = -4
IF(ier == -1) THEN
  jer = 1
ELSE IF(ier < 125) THEN
  jer = ier + 1
ELSE IF(ier > 9999) THEN
  jer = ier - 10000 + 126
ELSE IF(ier > 999) THEN
  jer = ier - 1000 + 126
END IF

a(1:)=' '
IF(jer > 0.AND.jer <= numerr) THEN
  a = sysmsg(jer)
ELSE
  a = '** UnKnown System Error **'
END IF

10  j = Len(a)
i = Len(ax)
WRITE(lu,900) ax(1:i),a(1:j)

900  FORMAT(' ** Error code: ',a,' --> ',a)

RETURN
END SUBROUTINE lu_errmsg

!-----------------------------------------------------------------------
!=====================  For PC AbSoft Fortran 77 only =====================

SUBROUTINE errmsg(msg,ier)

CHARACTER (LEN=*), INTENT(IN OUT)        :: msg
INTEGER, INTENT(IN OUT)                  :: ier

lu = 6
CALL lu_errmsg(lu,msg,ier)

RETURN
END SUBROUTINE errmsg

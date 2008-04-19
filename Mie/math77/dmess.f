      subroutine DMESS (MACT, TEXT, IDAT, FDAT)
c     .  Copyright (C) 1991, California Institute of Technology.
c     .  All rights reserved.  U. S. Government sponsorship under
c     .  NASA contract NAS7-918 is acknowledged.
c>> 1997-06-17 DMESS Krogh  In C code made messxc, static.
c>> 1996-07-12 DMESS Krogh  Changes to use .C. and C%%.
c>> 1996-03-30 DMESS Krogh  Added external statement.
c>> 1994-10-20 DMESS Krogh  Changes to use M77CON
c>> 1994-09-21 DMESS Krogh  Added CHGTYP code.
c>> 1994-09-08 DMESS Krogh  Added new matrix/vector capabilities.
c>> 1994-08-17 DMESS Krogh  Removed duplicate save statement.
c>> 1994-04-19 DMESS Krogh  Removed blank line from DMESS.
c>> 1993-05-14 DMESS Krogh  Changed TEXT to array of character strings.
c>> 1993-04-14 DMESS Krogh  Fixes for conversion to C. (C%% comments.)
c>> 1992-07-12 DMESS Krogh  Fixed so negative KDFDEF works.
c>> 1992-05-27 DMESS Krogh  Initialized LDFDEF in a data statement.
c>> 1992-05-14 DMESS Krogh  Put common blocks in save statement.
c>> 1992-04-28 DMESS Krogh  Corrected minor error in floating pt. format
c>> 1992-02-28 DMESS Krogh  Initial Code.
c
c--D replaces "?": ?MESS
c
c Processes Messages -- Actions are controlled by MACT().  See
c comment is subroutine MESS.  This program is for the extra
c argument of type real.
c
c BUF    In common CMESSC, see MESS.
c DOLS   In common for intitialization, not used here.  See MESS.
c EUNIT  In common for intitialization, not used here.  See MESS.
c FDAT   Formal argument -- gives floating point data to print.
c FMAG   Magnitude of floating point number to output, with negative
c   sign if floating point number is < 0.
c FMIN   Minimum floating point number.
c FMTF   In common CMESSC, format for printing floating point number.
c FMTG   In common CMESSC, user format to use in place of FMTF.
c FNMAX  Maximum negative floating point number.
c FOUT   Floating point number to be output.
c FSMA   Smallest postitive floating point number.
c ICOL   In common CMESSI, see MESS.
c ID     Number of decimal digits for floating point format statement.
c IDAT   Integer data -- passed to MESS.
c IVAR   In common CMESSI, see MESS.
c IWF    In common CMESSI, see MESS.
c IWG    In common CMESSI, see MESS.
c J      Temporary index.
c K      Temporary index.
c KAP    Number of extra 0's after the decimal point in "F" format.  If
c        < 0, -KAP gives the number of extra digits to the left of the
c        decimal point.  KAP depends on abs(smallest number printed).
c KBP    Number of extra digits before the decimal point required by the
c        largest number to be printed.
c KDF    In common CMESSI, see MESS.
c KDFDEF In common CMESSI, see MESS.
c KDIAG  In common CMESSI, not used here, see MESS.
c KEXE   Extra space required for E format.
c KF     In common CMESSI, see MESS.
c KLINE  In common CMESSI, see MESS.
c KSCRN  In common CMESSI, see MESS.
c KRES1  In common CMESSI, see MESS.
c KSPEC  In common CMESSI, see MESS.
c LASTER In common CMESSI, not used here, see MESS.
c LASTI  In common CMESSI, see MESS.
c LBUF   In common CMESSI, see MESS.
c LDFDEF Value of KDFDEF for this routine.  (Saved)
c LENBUF In common CMESSI, see MESS.
c LENLIN In common CMESSI, not used here, see MESS.
c LENTRY In common CMESSI, see MESS.
c LHEAD  In common CMESSI, not used here, see MESS.
c LINERR In common CMESSI, not used here, see MESS.
c LINMSG In common CMESSI, not used here, see MESS.
c LOCBEG In common CMESSI, see MESS.
c LPRINT In common CMESSI, not used here, see MESS.
c LSTOP  In common CMESSI, not used here, see MESS.
c LSTRT  In common CMESSI, see MESS.
c MACT   Formal argument, see MESS.
c MDAT   In common CMESSI, not used here, see MESS.
c MEMDA5 In common CMESSI, see MESS.
c MESS   Program called for most of the message processing.
c MPT    In common CMESSI, see MESS.
c MUNIT  In common CMESSI, not used here, see MESS.
c NCOL   In common CMESSI, see MESS.
c NDIM   In common CMESSI, see MESS.
c NFDAT  In common CMESSI, see MESS.
c NIDAT  In common CMESSI, not used here, see MESS.
c NMDAT  In common CMESSI, not used here, see MESS.
c NTEXT  In common CMESSI, not used here, see MESS.
c OUNIT  In common CMESSI, not used here, see MESS.
c D1MACH External func. giving floating pt. info. about the environment.
c SUNIT  In common CMESSI, not used here, see MESS.
c TEXT   Formal argument, passed to MESS, see there.
c XARGOK In common CMESSI, see MESS.
c
C%% static void messxc();
c
      external         D1MACH
      integer          MACT(*), IDAT(*)
      double precision FDAT(*)
      character        TEXT(*)*(*)
      integer          ICOL, ID, J, K, KAP, KBP, KEXE, LDFDEF
      double precision FMAG, FMIN, FOUT, FNMAX, FSMA, D1MACH
      save LDFDEF
      save /CMESSI/, /CMESSC/
c++ CODE for .C. is inactive
c      integer  kciwid, kccwid, kcrwid, lbeg, lend, lfprec, lgprec
c      common /MESSCC/ kciwid,kccwid,kcrwid,lbeg,lend,lfprec,lgprec
c++ END
c
c ************************** Data from common block ********************
c
c For comments on these variables, see the listing for MESS.
c
      integer   LENBUF, MEVBAS, MEVLAS
      parameter (LENBUF = 250)
      parameter (MEVBAS = 10)
      parameter (MEVLAS = 32)
      logical          XARG, GOTFMT, XARGOK
      integer          EUNIT, ICHAR0, IRC, IVAR(MEVBAS:MEVLAS), IMAG,
     1   INC, ITEXT, IWF, IWG, KDF, KDFDEF, KDI, KDIAG, KDJ, KLINE,
     2   KSCRN, KSHIFT, KSPEC, KT, MAXERR, LASTI, LBUF, LENLIN, LENOUT,
     3   LENTRY, LENTXT, LHEAD, LINERR, LINMSG, LOCBEG, LPRINT, LSTOP,
     4   LSTRT, LTEXT, MAXWID(2), MDAT(5), MPT, MUNIT, NCOL, NDIM,
     5   NFDAT, NIDAT, NMDAT, NROW, NTEXT, OUNIT, SUNIT, TABSPA
c
      character BUF*(LENBUF), DOLS*72, FMTC*7, FMTF*15, FMTG*15,
     1  FMTI*7, FMTIM(2)*7, FMTJ*7, FMTR*7, FMTT*15
      common /CMESSI/ SUNIT, LHEAD, KDFDEF, LINMSG, LINERR, MUNIT,
     1   EUNIT, KSCRN, KDIAG, MAXERR, LSTOP, LPRINT, KDF, NTEXT, NIDAT,
     2   NFDAT, NMDAT, MDAT, TABSPA, ICHAR0, IMAG, INC, IRC, ITEXT,
     3   IWF, IWG, KDI, KDJ, KLINE, KSHIFT, KSPEC, KT, LASTI, LBUF,
     4   LENLIN, LENOUT, LENTRY, LENTXT, LOCBEG, LSTRT, LTEXT, MAXWID,
     5   MPT, NROW, NCOL, NDIM, OUNIT, GOTFMT, XARG, XARGOK
      common /CMESSC / BUF, DOLS, FMTF, FMTG, FMTI, FMTJ, FMTT, FMTIM
      equivalence (IVAR(MEVBAS), SUNIT)
      equivalence (FMTIM(1), FMTR), (FMTIM(2), FMTC)
c
      data LDFDEF / 0 /
c
c ************************* Start of Executable Code *******************
c
      XARGOK = .true.
      if (LDFDEF .eq. 0) then
         LDFDEF = 1 - int(log10(d1mach(3)))
      end if
      KDFDEF = LDFDEF
      KDF = KDFDEF
   10 call MESS (MACT, TEXT, IDAT)
      go to (20, 100, 100, 200, 300, 400), LENTRY-3
      XARGOK = .false.
      LDFDEF = KDFDEF
      return
c                                      Print from FDAT
   20 J = LBUF + 1
      FOUT = FDAT(NFDAT)
      NFDAT = NFDAT + 1
      if (KSPEC .ge. 8) then
         LBUF = LBUF + IWG
C%%      messcc.lend = cmessi.lbuf;
C%%      cmessc.buf[messcc.lend] = ' ';
C%%      if ((j > 1) && (cmessc.buf[j-2] >= '0') &&
C%%         (cmessc.buf[j-2] <= '9')) cmessc.buf[j++ - 1] = ' ';
C%%      sprintf(&cmessc.buf[j-1], cmessc.fmtg, cmessi.iwg,
C%%         messcc.lgprec, fout);
C%%      if (cmessc.buf[messcc.lend] != 0) {messcc.lbeg = j; messxc( );}
         write (BUF(J:LBUF), FMTG) FOUT
         go to 10
      end if
      FMAG = FOUT
      FSMA = abs(FOUT)
      IWF = 1
c                                      Get the format.
   40 continue
      KEXE = 3
      KBP = 0
      KAP = 0
      if (FMAG .ne. 0.D0) then
         if (FMAG .lt. 0.D0) then
            FMAG = -FMAG
            IWF = IWF + 1
         end if
         if (FSMA .ge. 1.D0) then
            KAP = -log10(FSMA) - 1.000001D0
         else
            KAP = -log10(FSMA) + 1.D-6
         end if
         if (FMAG .ne. FSMA) then
            if (FMAG .ge. 1.D0) KBP = 1 + int(log10(FMAG) + 1.D-6)
         else
            KBP = max(-KAP, 0)
         end if
         K = max(KAP+1, KBP)
         if (K .ge. 10) KEXE = 3 + int(log10(dble(K)) + 1.D-5)
      end if
      if (KDF .eq. 0) KDF = KDFDEF
      if (KDF .gt. 0) then
         if (KBP + max(KAP, 0) .ge. KEXE) go to 50
         ID = max(KDF + KAP, 0)
         IWF = IWF + ID + KBP
c If LENTRY is 5 or 6 (Vector or Matrix output), and KBP is 0, need 1
c extra since with the extra space for a blank an extra 0 is printed.
         if (LENTRY - 2*KBP .gt. 4) IWF = IWF + 1
      else
         if (KBP .gt. KEXE) go to 50
         ID = -KDF
         IWF = IWF + KBP + ID
      end if
c++ CODE for ~.C. is active
      write (FMTF, '(6H(0P99F,I2,1H.,I2,1H))') IWF,ID
c++ CODE for .C. is inactive
c%%    strcpy(cmessc.fmtf, "%*.*f\0");
c      lfprec = id
c++ END
      go to 60
   50 ID = KDF - 1
      if (ID .lt. 0) ID = -ID - 1
c++ CODE for ~.C. is active
      IWF = IWF + ID + KEXE + 1
      write (FMTF, '(6H(1P99E,I2,1H.,I2,1HE,I1,1H))') IWF,ID,KEXE-2
c++ CODE for .C. is inactive
cc "+ 6" needed for WATCOM C (and others?), some can use "+ 5"
c      iwf = iwf + id + 6
c%%    strcpy(cmessc.fmtf, "%*.*E\0");
c      lfprec = id
c++ END
   60 if (LENTRY .ne. 4) go to 10
c
      LBUF = LBUF + IWF
C%%   messcc.lend = cmessi.lbuf;
C%%   cmessc.buf[messcc.lend] = ' ';
C%%   if ((j > 1) && (cmessc.buf[j-2] >= '0') &&
C%%      (cmessc.buf[j-2] <= '9')) cmessc.buf[j++ - 1] = ' ';
C%%   sprintf(&cmessc.buf[j-1], cmessc.fmtf, cmessi.iwf,
C%%      messcc.lfprec, fout);
C%%   if (cmessc.buf[messcc.lend] != 0) {messcc.lbeg = j; messxc( );}
      write (BUF(J:LBUF),FMTF) FOUT
      go to 10
c                                     Get format for a vector or matrix
  100 if (KDF .eq. 0) KDF = KDFDEF
      ICOL = 1
      FMAG = FDAT(LOCBEG)
      FMIN = FMAG
      FSMA = 1.D0
      FNMAX = -1.D0
  110 do 120 J = LOCBEG, LASTI, INC
         if (FDAT(J) .le. 0.D0) then
            FMIN = min(FMIN, FDAT(J))
            if (FDAT(J) .ne. 0.D0) FNMAX = max(FNMAX, FDAT(J))
         else
            FMAG = max(FMAG, FDAT(J))
            FSMA = min(FSMA, FDAT(J))
         end if
  120 continue
      if (NCOL .ne. 0) then
         ICOL = ICOL + 1
         LOCBEG = LOCBEG + NDIM
         LASTI = LASTI + NDIM
         if (ICOL .le. NCOL) go to 110
      end if
      if (FMIN .lt. 0.D0) then
         FMAG = min(FMIN, -0.1D0*FMAG)
         FSMA = min(-FNMAX, FSMA)
      end if
      if (FSMA .ge. 1.D0) then
         FSMA = abs(FMAG)
      else if (abs(FMAG) .lt. 1.D0) then
         FMAG = sign(FSMA, FMAG)
      end if
      IWF = 2
      go to 40
c                                    Floating point vector output
  200 continue
C%%   messcc.lend = cmessi.lstrt-1;
C%%   for (k=cmessi.mpt; k<cmessi.mpt+cmessi.kline; k++){
C%%      messcc.lend = messcc.lend + cmessi.iwf;
C%%      cmessc.buf[messcc.lend] = ' ';
C%%      messcc.lbeg = messcc.lend - cmessi.iwf;
C%%      if ((messcc.lbeg > 1) && (cmessc.buf[messcc.lbeg-1] >= '0') &&
C%%         (cmessc.buf[messcc.lbeg-1] <= '9'))
C%%         cmessc.buf[messcc.lbeg++] = ' ';
C%%      sprintf(&cmessc.buf[messcc.lbeg], cmessc.fmtf,
C%%         cmessi.iwf, messcc.lfprec, fdat[cmessi.inc*k-1]);
C%%      if (cmessc.buf[messcc.lend] != 0) messxc( ); }
      write(BUF(LSTRT:LBUF),FMTF) (FDAT(K),K=MPT,MPT+INC*(KLINE-1),INC)
      MPT = MPT + KLINE * INC
      go to 10
c                                    Floating point matrix output
  300 continue
C%%   j = 0;
C%%   messcc.lend = cmessi.lstrt-1;
C%%   for (k=cmessi.mpt; k<=cmessi.lasti; k+=cmessi.ndim){
C%%      messcc.lend = messcc.lend + cmessi.iwf;
C%%      cmessc.buf[messcc.lend] = ' ';
C%%      messcc.lbeg = messcc.lend - cmessi.iwf;
C%%      if ((messcc.lbeg > 1) && (cmessc.buf[messcc.lbeg-1] >= '0') &&
C%%         (cmessc.buf[messcc.lbeg-1] <= '9'))
C%%         cmessc.buf[messcc.lbeg++] = ' ';
C%%      sprintf(&cmessc.buf[messcc.lbeg],
C%%         cmessc.fmtf, cmessi.iwf, messcc.lfprec, fdat[k-1]);
C%%      if (cmessc.buf[messcc.lend] != 0) messxc( ); }
      write (BUF(LSTRT:LBUF), FMTF) (FDAT(K), K = MPT, LASTI, NDIM)
      go to 10
c                                    Table output
  400 continue
C%%   messcc.lend = cmessi.lstrt-1;
C%%   for (k=cmessi.mpt; k<cmessi.mpt+cmessi.kline; k++){
C%%      messcc.lend = messcc.lend + cmessi.iwf;
C%%      cmessc.buf[messcc.lend] = ' ';
C%%      messcc.lbeg = messcc.lend - cmessi.iwf;
C%%      if ((messcc.lbeg > 1) && (cmessc.buf[messcc.lbeg-1] >= '0') &&
C%%         (cmessc.buf[messcc.lbeg-1] <= '9'))
C%%         cmessc.buf[messcc.lbeg++] = ' ';
C%%      sprintf(&cmessc.buf[messcc.lbeg], cmessc.fmtf,
C%%         cmessi.iwf, messcc.lfprec, fdat[k-1]);
C%%      if (cmessc.buf[messcc.lend] != 0) messxc( ); }
      write (BUF(LSTRT:LBUF), FMTF) (FDAT(K), K = MPT, MPT+KLINE-1)
      go to 10
      end

C%%  static void /*FUNCTION*/ messxc()
C%%{
C%%  /* Adjust for lack of control on digits in exponent */
C%%   long int k, k0;
C%%   if (cmessc.fmtf[4] == 'f') return;
C%%   k0 = messcc.lend-2;
C%%   if (cmessc.buf[k0] != '0') k0++;
C%%   while (cmessc.buf[messcc.lend] != 0) {
C%%      if (cmessc.buf[k0] != '0') {
C%%         if ((cmessc.buf[messcc.lbeg-1] == ' ') &&
C%%            (cmessc.buf[messcc.lbeg] < '0') ||
C%%            (cmessc.buf[messcc.lbeg] > '9')) {
C%%               messcc.lbeg = messcc.lbeg - 1;
C%%               for (k = messcc.lbeg; cmessc.buf[k] != 0; k++)
C%%                  cmessc.buf[k] = cmessc.buf[k+1];
C%%                  continue; }
C%%         else cmessc.buf[k0] = '?';}
C%%      else cmessc.buf[k0] = cmessc.buf[k0+1];
C%%      for (k = k0+1; cmessc.buf[k] != 0; k++)
C%%         cmessc.buf[k] = cmessc.buf[k+1]; }
C%%      return;
C%%} /* end of function */









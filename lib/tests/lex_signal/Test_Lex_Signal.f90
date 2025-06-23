program Test_Lex_Signal

  use Lexer_Core, only: Init_Lexer, Token
  use Lexer_M, only: Lex_Signal
  use Symbol_Types, only: T_End_of_input
  use Toggles, only: Lex, Toggle

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

  character(len=80) :: Buf    ! Input buffer
  type(token) :: The_Token
  integer :: Where            ! Last examined character in the input buffer

  call init_lexer ( n_chars=10000, n_symbols=1000, hash_table_size=1003 )
  toggle(lex) = .true.
  do
    read ( *, '(a)', end=99 ) buf
    where = 0
    do
      call lex_signal ( buf, where, the_token )
      if ( the_token%class == t_end_of_input ) exit
    end do
  end do
99 stop
end program Test_Lex_Signal

! $Log$

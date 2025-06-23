  program test_streq
    use MLSStrings, ONLY: streq
    character(len=256) :: str1, str2,str3
    str1="*rad*"
    str2="RadiancesFile.dat"
    str3="*snoop*"
    str3="snoopradiances"
    print *,"wc, *rad*, RadiancesFile.dat ",streq(trim(str1),trim(str2),"wc")
    print *,"wc, *snoop*, snoopradiances ",streq(trim(str3),trim(str4),"wc")

    
  end program test_streq

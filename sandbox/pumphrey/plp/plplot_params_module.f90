module plplot_params_module
  !Contains some parameters that plplot was compiled with
  integer,parameter,public::PLFLT=selected_real_kind(13)
  integer,public,parameter::PLINT=selected_int_kind(7)
end module plplot_params_module

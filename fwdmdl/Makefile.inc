PROG =	driver0

SRCS =	MLSCommon.f90 abs_cs_air_cont_m.f90 abs_cs_h2o_213g_cont_m.f90 \
	abs_cs_liq_h2o_m.f90 abs_cs_o2_cont_m.f90 aitken_int_m.f90 \
	complete_i_k_m.f90 create_beta_coeff_m.f90 \
	create_pfa_beta_coeff_m.f90 cs_intrp_m.f90 d_brkline_m.f90 \
	d_csinterp_m.f90 d_cspline_m.f90 d_delta_dt_m.f90 d_fft.f90 \
	d_get_one_eta_m.f90 d_hunt_m.f90 d_lintrp_m.f90 d_pcspl_m.f90 \
	d_q_log_m.f90 d_solve_quad_m.f90 d_stat_temp_m.f90 \
	d_t_script_dtnp_m.f90 daitken_module.f90 dcspline_der_m.f90 \
	do_t_script_m.f90 driver0.f90 dsimpson_module.f90 \
	dump_k_records_m.f90 earth_intersection_m.f90 ellipse.f90 eos_mdb.f90 \
	ermsg_m.f90 format_key_m.f90 fov_convolve_m.f90 fwd_mdl_m.f90 \
	gauss_legendre_m.f90 get_delta_m.f90 get_drad_m.f90 get_eta_m.f90 \
	get_pfa_delta_m.f90 gh06.f90 gl6p.f90 hydrostatic_intrp.f90 \
	i_and_k_star_m.f90 i_k_star_no_conv_m.f90 l2pc_file.f90 \
	l2pc_file_parameters.f90 l2pc_file_structures.f90 \
	l2pc_pfa_structures.f90 l2pcdim.f90 mdbeta.f90 pfa_beta_intrp_m.f90 \
	pfa_db_delta_m.f90 pfa_dbeta_intrp_m.f90 pfa_fwd_model_m.f90 \
	physics.f90 read_uif_m.f90 s_brkline_m.f90 s_csinterp_m.f90 \
	s_cspline_m.f90 s_get_one_eta_m.f90 s_hunt_m.f90 s_lintrp_m.f90 \
	s_pcspl_m.f90 s_solve_quad_m.f90 s_stat_temp_m.f90 saitken_module.f90 \
	scrt_dn_m.f90 scspline_der_m.f90 slabs.f90 slabs_prep_m.f90 \
	spectro_derivative_m.f90 ssimpson_module.f90 strings.f90 \
	temperature_deriv_m.f90 two_d_polate.f90 units.f90 \
	write_k_records_m.f90 write_one_record_m.f90 write_x_i_records_m.f90 \
	zatmos_deriv_m.f90 zopacity_m.f90

OBJS =	MLSCommon.o abs_cs_air_cont_m.o abs_cs_h2o_213g_cont_m.o \
	abs_cs_liq_h2o_m.o abs_cs_o2_cont_m.o aitken_int_m.o complete_i_k_m.o \
	create_beta_coeff_m.o create_pfa_beta_coeff_m.o cs_intrp_m.o \
	d_brkline_m.o d_csinterp_m.o d_cspline_m.o d_delta_dt_m.o d_fft.o \
	d_get_one_eta_m.o d_hunt_m.o d_lintrp_m.o d_pcspl_m.o d_q_log_m.o \
	d_solve_quad_m.o d_stat_temp_m.o d_t_script_dtnp_m.o daitken_module.o \
	dcspline_der_m.o do_t_script_m.o driver0.o dsimpson_module.o \
	dump_k_records_m.o earth_intersection_m.o ellipse.o eos_mdb.o \
	ermsg_m.o format_key_m.o fov_convolve_m.o fwd_mdl_m.o \
	gauss_legendre_m.o get_delta_m.o get_drad_m.o get_eta_m.o \
	get_pfa_delta_m.o gh06.o gl6p.o hydrostatic_intrp.o i_and_k_star_m.o \
	i_k_star_no_conv_m.o l2pc_file.o l2pc_file_parameters.o \
	l2pc_file_structures.o l2pc_pfa_structures.o l2pcdim.o mdbeta.o \
	pfa_beta_intrp_m.o pfa_db_delta_m.o pfa_dbeta_intrp_m.o \
	pfa_fwd_model_m.o physics.o read_uif_m.o s_brkline_m.o s_csinterp_m.o \
	s_cspline_m.o s_get_one_eta_m.o s_hunt_m.o s_lintrp_m.o s_pcspl_m.o \
	s_solve_quad_m.o s_stat_temp_m.o saitken_module.o scrt_dn_m.o \
	scspline_der_m.o slabs.o slabs_prep_m.o spectro_derivative_m.o \
	ssimpson_module.o strings.o temperature_deriv_m.o two_d_polate.o \
	units.o write_k_records_m.o write_one_record_m.o \
	write_x_i_records_m.o zatmos_deriv_m.o zopacity_m.o


ALLOBJS = $(OBJS) machine.o
LIBS =	

all: $(PROG)
objs: $(ALLOBJS)

$(PROG): $(ALLOBJS)
	$(F90) $(LDFLAGS) -o $@ $(ALLOBJS) $(LIBS)

clean:
	rm -f $(PROG) $(ALLOBJS) *.mod

#.f90.o:
%.o %.mod: %.f90
	$(F90) $(F90FLAGS) -c $(MODULEPATH) $(INCLUDE) $<

#.f90.mod:
#	/bin/rm -f $@

#.mod.o:
#	/bin/rm -f $<

%.mod::

.SUFFIXES: 
.SUFFIXES: $(SUFFIXES) .f90

driver0: driver0.o machine.o 

aitken_int_m.o: MLSCommon.o d_cspline_m.o daitken_module.o dsimpson_module.o
complete_i_k_m.o: l2pc_file_parameters.o l2pc_pfa_structures.o l2pcdim.o
create_beta_coeff_m.o: MLSCommon.o abs_cs_air_cont_m.o \
	abs_cs_h2o_213g_cont_m.o abs_cs_liq_h2o_m.o d_hunt_m.o eos_mdb.o \
	l2pc_pfa_structures.o s_hunt_m.o slabs.o
create_pfa_beta_coeff_m.o: MLSCommon.o create_beta_coeff_m.o eos_mdb.o gl6p.o \
	l2pc_pfa_structures.o l2pcdim.o
cs_intrp_m.o: MLSCommon.o l2pcdim.o mdbeta.o s_hunt_m.o
d_brkline_m.o: brkline.f9h
d_csinterp_m.o: MLSCommon.o csinterp.f9h
d_cspline_m.o: d_hunt_m.o d_pcspl_m.o cspline.f9h
d_delta_dt_m.o: MLSCommon.o cs_intrp_m.o ellipse.o gl6p.o mdbeta.o \
	s_get_one_eta_m.o two_d_polate.o
d_fft.o: ermsg_m.o cft.f9h fft.f9h rft.f9h rft1.f9h tcst.f9h
d_get_one_eta_m.o: MLSCommon.o get_one_eta.f9h
d_hunt_m.o: MLSCommon.o hunt.f9h
d_lintrp_m.o: d_hunt_m.o lintrp.f9h
d_pcspl_m.o: pcspl.f9h
d_solve_quad_m.o: solve_quad.f9h
d_stat_temp_m.o: MLSCommon.o physics.o stat_temp.f9h
d_t_script_dtnp_m.o: MLSCommon.o d_stat_temp_m.o l2pcdim.o physics.o \
	s_get_one_eta_m.o
daitken_module.o: aitken.f9h
dcspline_der_m.o: d_hunt_m.o d_pcspl_m.o cspline_der.f9h
do_t_script_m.o: MLSCommon.o l2pcdim.o s_stat_temp_m.o
driver0.o: MLSCommon.o eos_mdb.o fwd_mdl_m.o l2pc_file_parameters.o \
	l2pc_file_structures.o l2pc_pfa_structures.o l2pcdim.o mdbeta.o
dsimpson_module.o: MLSCommon.o d_cspline_m.o simpson.inc
dump_k_records_m.o: MLSCommon.o format_key_m.o l2pc_file_parameters.o \
	l2pc_file_structures.o l2pcdim.o strings.o
earth_intersection_m.o: MLSCommon.o ellipse.o
ellipse.o: MLSCommon.o
eos_mdb.o: MLSCommon.o
format_key_m.o: MLSCommon.o machine.o strings.o
fov_convolve_m.o: MLSCommon.o d_fft.o l2pcdim.o machine.o s_cspline_m.o \
	units.o
fwd_mdl_m.o: MLSCommon.o complete_i_k_m.o d_t_script_dtnp_m.o do_t_script_m.o \
	dump_k_records_m.o earth_intersection_m.o ellipse.o eos_mdb.o \
	get_delta_m.o get_eta_m.o gl6p.o i_and_k_star_m.o \
	i_k_star_no_conv_m.o l2pc_file_parameters.o l2pc_file_structures.o \
	l2pc_pfa_structures.o l2pcdim.o mdbeta.o pfa_fwd_model_m.o \
	s_lintrp_m.o scrt_dn_m.o temperature_deriv_m.o write_k_records_m.o \
	write_x_i_records_m.o zatmos_deriv_m.o zopacity_m.o
gauss_legendre_m.o: MLSCommon.o
get_delta_m.o: MLSCommon.o cs_intrp_m.o ellipse.o gl6p.o mdbeta.o \
	s_get_one_eta_m.o
get_drad_m.o: MLSCommon.o
get_eta_m.o: MLSCommon.o
get_pfa_delta_m.o: MLSCommon.o ellipse.o gl6p.o pfa_beta_intrp_m.o \
	s_get_one_eta_m.o
gh06.o: MLSCommon.o
gl6p.o: MLSCommon.o
hydrostatic_intrp.o: MLSCommon.o d_solve_quad_m.o s_hunt_m.o
i_and_k_star_m.o: MLSCommon.o fov_convolve_m.o hydrostatic_intrp.o \
	l2pc_file_parameters.o l2pc_file_structures.o l2pc_pfa_structures.o \
	l2pcdim.o s_csinterp_m.o s_cspline_m.o s_lintrp_m.o scspline_der_m.o
i_k_star_no_conv_m.o: MLSCommon.o hydrostatic_intrp.o l2pc_file_parameters.o \
	l2pc_file_structures.o l2pc_pfa_structures.o l2pcdim.o s_csinterp_m.o \
	s_cspline_m.o s_lintrp_m.o scspline_der_m.o
l2pc_file.o: MLSCommon.o eos_mdb.o fwd_mdl_m.o gl6p.o l2pc_file_parameters.o \
	l2pc_file_structures.o l2pc_pfa_structures.o l2pcdim.o mdbeta.o \
	read_uif_m.o time_mod.o
l2pc_file_parameters.o: MLSCommon.o l2pcdim.o
l2pc_file_structures.o: MLSCommon.o l2pc_file_parameters.o
l2pc_pfa_structures.o: MLSCommon.o l2pc_file_parameters.o l2pcdim.o
l2pcdim.o: MLSCommon.o
mdbeta.o: MLSCommon.o
pfa_beta_intrp_m.o: MLSCommon.o s_hunt_m.o
pfa_db_delta_m.o: MLSCommon.o ellipse.o gl6p.o pfa_dbeta_intrp_m.o \
	s_get_one_eta_m.o two_d_polate.o
pfa_dbeta_intrp_m.o: MLSCommon.o s_hunt_m.o
pfa_fwd_model_m.o: MLSCommon.o aitken_int_m.o create_beta_coeff_m.o \
	create_pfa_beta_coeff_m.o do_t_script_m.o earth_intersection_m.o \
	ellipse.o eos_mdb.o fov_convolve_m.o get_pfa_delta_m.o gl6p.o \
	hydrostatic_intrp.o l2pc_file_parameters.o l2pc_file_structures.o \
	l2pc_pfa_structures.o l2pcdim.o s_hunt_m.o s_lintrp_m.o scrt_dn_m.o \
	scspline_der_m.o spectro_derivative_m.o strings.o zopacity_m.o
physics.o: MLSCommon.o
read_uif_m.o: MLSCommon.o l2pc_file_parameters.o l2pc_file_structures.o \
	l2pc_pfa_structures.o l2pcdim.o machine.o s_brkline_m.o strings.o
s_brkline_m.o: brkline.f9h
s_csinterp_m.o: MLSCommon.o csinterp.f9h
s_cspline_m.o: s_hunt_m.o s_pcspl_m.o cspline.f9h
s_get_one_eta_m.o: MLSCommon.o get_one_eta.f9h
s_hunt_m.o: MLSCommon.o hunt.f9h
s_lintrp_m.o: s_hunt_m.o lintrp.f9h
s_pcspl_m.o: pcspl.f9h
s_solve_quad_m.o: solve_quad.f9h
s_stat_temp_m.o: MLSCommon.o physics.o stat_temp.f9h
saitken_module.o: aitken.f9h
scrt_dn_m.o: MLSCommon.o
scspline_der_m.o: s_hunt_m.o s_pcspl_m.o cspline_der.f9h
slabs.o: MLSCommon.o gh06.o
slabs_prep_m.o: MLSCommon.o d_q_log_m.o
spectro_derivative_m.o: MLSCommon.o get_drad_m.o l2pc_pfa_structures.o \
	l2pcdim.o pfa_db_delta_m.o
ssimpson_module.o: MLSCommon.o s_cspline_m.o simpson.inc
temperature_deriv_m.o: MLSCommon.o d_delta_dt_m.o get_drad_m.o \
	l2pc_file_parameters.o l2pc_pfa_structures.o l2pcdim.o
two_d_polate.o: MLSCommon.o s_get_one_eta_m.o
write_k_records_m.o: MLSCommon.o l2pc_file_structures.o machine.o \
	write_one_record_m.o
write_one_record_m.o: MLSCommon.o format_key_m.o l2pc_file_parameters.o \
	l2pc_file_structures.o machine.o
write_x_i_records_m.o: MLSCommon.o l2pc_file_parameters.o \
	l2pc_file_structures.o l2pc_pfa_structures.o l2pcdim.o \
	write_one_record_m.o
zatmos_deriv_m.o: MLSCommon.o get_drad_m.o l2pc_file_parameters.o \
	l2pc_pfa_structures.o l2pcdim.o
zopacity_m.o: MLSCommon.o

# $Log$
# Revision 1.3  2000/05/05 17:18:14  vsnyder
# Replaced because CVS lost v1.1
#
# Revision 1.1  2000/05/04 18:12:04  vsnyder
# Initial conversion to Fortran 90
#

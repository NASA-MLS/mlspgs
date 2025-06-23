% LR output converted to TeX by TeX_Grammar program
  \item
    \begin{tabbing}
      $<$GOAL$>$ \= $\rightarrow$ SOG cf EOG\\
    \end{tabbing}
  \item
    \begin{tabbing}
      cf \= $\rightarrow$ cfs\\
    \end{tabbing}
  \item
    \begin{tabbing}
      cfs \= $\rightarrow$ one_cf EOS\\
      \> $\rightarrow$ include\\
      \> $\rightarrow$ cfs one_cf EOS\\
    \end{tabbing}
  \item
    \begin{tabbing}
      one_cf \= $\rightarrow$\\
      \> $\rightarrow$ construct_label select case '(' expr ')' 'EOS+' cases_outer end select $\Rightarrow$ n_select\\
      \> $\rightarrow$ construct_label do do_header EOS cfs end do $\Rightarrow$ n_do\\
      \> $\rightarrow$ construct_label do while '(' expr ')' EOS cfs end do $\Rightarrow$ n_while\\
      \> $\rightarrow$ if_block_outer more_if_outer end if $\Rightarrow$ n_if\\
      \> $\rightarrow$ if_label if_block_outer more_if_outer end if $\Rightarrow$ n_if\\
      \> $\rightarrow$ if_stmt_outer $\Rightarrow$ n_if\\
      \> $\rightarrow$ if_consequent_outer\\
      \> $\rightarrow$ begin name EOS specs end name $\Rightarrow$ n_cf\\
    \end{tabbing}
  \item
    \begin{tabbing}
      cases_outer \= $\rightarrow$ case_test_outer\\
      \> $\rightarrow$ case default EOS cfs $\Rightarrow$ n_default\\
      \> $\rightarrow$ case_test_outer cases_outer\\
    \end{tabbing}
  \item
    \begin{tabbing}
      case_test_outer \= $\rightarrow$ case '(' expr ')' EOS cfs $\Rightarrow$ n_test\\
    \end{tabbing}
  \item
    \begin{tabbing}
      'EOS+' \= $\rightarrow$ EOS\\
      \> $\rightarrow$ 'EOS+' EOS\\
    \end{tabbing}
  \item
    \begin{tabbing}
      construct_label \= $\rightarrow$\\
      \> $\rightarrow$ if_label\\
    \end{tabbing}
  \item
    \begin{tabbing}
      if_label \= $\rightarrow$ name ':' $\Rightarrow$ n_named\\
    \end{tabbing}
  \item
    \begin{tabbing}
      do_header \= $\rightarrow$ name ':=' expr_array $\Rightarrow$ n_do_head\\
      \> $\rightarrow$ name ':=' expr ',' expr $\Rightarrow$ n_do_head\\
      \> $\rightarrow$ name ':=' expr ',' expr ',' expr $\Rightarrow$ n_do_head\\
    \end{tabbing}
  \item
    \begin{tabbing}
      expr_array \= $\rightarrow$ expr $\Rightarrow$ n_array\\
      \> $\rightarrow$ '[' exprs ']' $\Rightarrow$ n_array\\
    \end{tabbing}
  \item
    \begin{tabbing}
      if_block_outer \= $\rightarrow$ if '(' expr ')' then EOS cfs $\Rightarrow$ n_test\\
    \end{tabbing}
  \item
    \begin{tabbing}
      more_if_outer \= $\rightarrow$\\
      \> $\rightarrow$ else if_block_outer more_if_outer\\
      \> $\rightarrow$ else EOS cfs $\Rightarrow$ n_else\\
    \end{tabbing}
  \item
    \begin{tabbing}
      if_stmt_outer \= $\rightarrow$ if '(' expr ')' if_consequent_outer $\Rightarrow$ n_test\\
    \end{tabbing}
  \item
    \begin{tabbing}
      if_consequent_outer \= $\rightarrow$ cycle_stmt\\
      \> $\rightarrow$ exit_stmt\\
      \> $\rightarrow$ variable_def\\
    \end{tabbing}
  \item
    \begin{tabbing}
      cycle_stmt \= $\rightarrow$ cycle $\Rightarrow$ n_cycle\\
      \> $\rightarrow$ cycle name $\Rightarrow$ n_cycle\\
    \end{tabbing}
  \item
    \begin{tabbing}
      exit_stmt \= $\rightarrow$ exit $\Rightarrow$ n_exit\\
      \> $\rightarrow$ exit name $\Rightarrow$ n_exit\\
    \end{tabbing}
  \item
    \begin{tabbing}
      variable_def \= $\rightarrow$ name ':=' value $\Rightarrow$ n_variable\\
      \> $\rightarrow$ name ':=' $\Rightarrow$ n_variable\\
    \end{tabbing}
  \item
    \begin{tabbing}
      specs \= $\rightarrow$\\
      \> $\rightarrow$ specs spec\\
    \end{tabbing}
  \item
    \begin{tabbing}
      spec \= $\rightarrow$ one_spec EOS\\
      \> $\rightarrow$ include\\
    \end{tabbing}
  \item
    \begin{tabbing}
      one_spec \= $\rightarrow$\\
      \> $\rightarrow$ construct_label select case '(' expr ')' 'EOS+' cases end select $\Rightarrow$ n_select\\
      \> $\rightarrow$ construct_label do while '(' expr ')' EOS specs end do $\Rightarrow$ n_while\\
      \> $\rightarrow$ construct_label do do_header EOS specs end do $\Rightarrow$ n_do\\
      \> $\rightarrow$ if_label if_block more_if end if $\Rightarrow$ n_if\\
      \> $\rightarrow$ if_block more_if end if $\Rightarrow$ n_if\\
      \> $\rightarrow$ if_stmt $\Rightarrow$ n_if\\
      \> $\rightarrow$ if_consequent\\
    \end{tabbing}
  \item
    \begin{tabbing}
      include \= $\rightarrow$ '\#include' string EOS $\Rightarrow$ 9\\
    \end{tabbing}
  \item
    \begin{tabbing}
      cases \= $\rightarrow$ case_test\\
      \> $\rightarrow$ case default EOS specs $\Rightarrow$ n_default\\
      \> $\rightarrow$ case_test cases\\
    \end{tabbing}
  \item
    \begin{tabbing}
      case_test \= $\rightarrow$ case '(' expr ')' EOS specs $\Rightarrow$ n_test\\
    \end{tabbing}
  \item
    \begin{tabbing}
      if_block \= $\rightarrow$ if '(' expr ')' then EOS specs $\Rightarrow$ n_test\\
    \end{tabbing}
  \item
    \begin{tabbing}
      more_if \= $\rightarrow$\\
      \> $\rightarrow$ else if_block more_if\\
      \> $\rightarrow$ else EOS specs $\Rightarrow$ n_else\\
    \end{tabbing}
  \item
    \begin{tabbing}
      if_stmt \= $\rightarrow$ if '(' expr ')' if_consequent $\Rightarrow$ n_test\\
    \end{tabbing}
  \item
    \begin{tabbing}
      if_consequent \= $\rightarrow$ cycle_stmt\\
      \> $\rightarrow$ exit_stmt\\
      \> $\rightarrow$ variable_def\\
      \> $\rightarrow$ name ':' spec_rest $\Rightarrow$ n_named\\
      \> $\rightarrow$ spec_rest\\
      \> $\rightarrow$ name '=' value $\Rightarrow$ n_equal\\
    \end{tabbing}
  \item
    \begin{tabbing}
      spec_rest \= $\rightarrow$ name fields $\Rightarrow$ n_spec_args\\
    \end{tabbing}
  \item
    \begin{tabbing}
      fields \= $\rightarrow$\\
      \> $\rightarrow$ ',' field fields\\
    \end{tabbing}
  \item
    \begin{tabbing}
      field \= $\rightarrow$ name '=' value $\Rightarrow$ n_asg\\
      \> $\rightarrow$ '/' name $\Rightarrow$ n_set_one\\
    \end{tabbing}
  \item
    \begin{tabbing}
      value \= $\rightarrow$ expr\\
      \> $\rightarrow$ '[' value2_list ']'\\
    \end{tabbing}
  \item
    \begin{tabbing}
      value2_list \= $\rightarrow$ value2\\
      \> $\rightarrow$ value2_list ',' value2\\
    \end{tabbing}
  \item
    \begin{tabbing}
      value2 \= $\rightarrow$ expr\\
      \> $\rightarrow$ '[' exprs ']' $\Rightarrow$ n_array\\
    \end{tabbing}
  \item
    \begin{tabbing}
      exprs \= $\rightarrow$ expr\\
      \> $\rightarrow$ exprs ',' expr\\
    \end{tabbing}
  \item
    \begin{tabbing}
      expr \= $\rightarrow$ cond\\
      \> $\rightarrow$ test '?' expr '!' expr $\Rightarrow$ n_cond\\
    \end{tabbing}
  \item
    \begin{tabbing}
      cond \= $\rightarrow$ limit\\
      \> $\rightarrow$ limit ':' limit $\Rightarrow$ n_colon\\
      \> $\rightarrow$ limit ':$<$' limit $\Rightarrow$ n_less\\
      \> $\rightarrow$ limit '$<$:' limit $\Rightarrow$ n_less_colon\\
      \> $\rightarrow$ limit '$<$:$<$' limit $\Rightarrow$ n_less\\
    \end{tabbing}
  \item
    \begin{tabbing}
      limit \= $\rightarrow$ lterm\\
      \> $\rightarrow$ limit or lterm $\Rightarrow$ n_or\\
    \end{tabbing}
  \item
    \begin{tabbing}
      lterm \= $\rightarrow$ lnot\\
      \> $\rightarrow$ lterm and lnot $\Rightarrow$ n_and\\
    \end{tabbing}
  \item
    \begin{tabbing}
      lnot \= $\rightarrow$ test\\
      \> $\rightarrow$ not test $\Rightarrow$ n_not\\
    \end{tabbing}
  \item
    \begin{tabbing}
      test \= $\rightarrow$ lfactor\\
      \> $\rightarrow$ lfactor '$<$' lfactor $\Rightarrow$ n_less\\
      \> $\rightarrow$ lfactor '$<$=' lfactor $\Rightarrow$ n_less_eq\\
      \> $\rightarrow$ lfactor '$>$' lfactor $\Rightarrow$ n_greater\\
      \> $\rightarrow$ lfactor '$>$=' lfactor $\Rightarrow$ n_greater_eq\\
      \> $\rightarrow$ lfactor '==' lfactor $\Rightarrow$ n_equal_equal\\
      \> $\rightarrow$ lfactor '/=' lfactor $\Rightarrow$ n_not_equal\\
    \end{tabbing}
  \item
    \begin{tabbing}
      lfactor \= $\rightarrow$ term\\
      \> $\rightarrow$ '+' term $\Rightarrow$ n_plus\\
      \> $\rightarrow$ '-' term $\Rightarrow$ n_minus\\
      \> $\rightarrow$ lfactor '+' term $\Rightarrow$ n_plus\\
      \> $\rightarrow$ lfactor '-' term $\Rightarrow$ n_minus\\
    \end{tabbing}
  \item
    \begin{tabbing}
      term \= $\rightarrow$ factor\\
      \> $\rightarrow$ term '*' factor $\Rightarrow$ n_mult\\
      \> $\rightarrow$ term '/' factor $\Rightarrow$ n_div\\
      \> $\rightarrow$ term '$\backslash$' factor $\Rightarrow$ n_into\\
    \end{tabbing}
  \item
    \begin{tabbing}
      factor \= $\rightarrow$ primary\\
      \> $\rightarrow$ primary '$\wedge$' factor $\Rightarrow$ n_pow\\
    \end{tabbing}
  \item
    \begin{tabbing}
      primary \= $\rightarrow$ name dots $\Rightarrow$ n_dot ?\\
      \> $\rightarrow$ number\\
      \> $\rightarrow$ number name $\Rightarrow$ n_unit\\
      \> $\rightarrow$ string\\
      \> $\rightarrow$ '(' expr ')'\\
      \> $\rightarrow$ '(' expr ')' name $\Rightarrow$ n_unit\\
      \> $\rightarrow$ name '[' expr ']' $\Rightarrow$ n_subscript\\
      \> $\rightarrow$ func_ref\\
      \> $\rightarrow$ func_ref name $\Rightarrow$ n_unit\\
    \end{tabbing}
  \item
    \begin{tabbing}
      func_ref \= $\rightarrow$ name '(' value2_list ')' $\Rightarrow$ n_func_ref\\
    \end{tabbing}
  \item
    \begin{tabbing}
      dots \= $\rightarrow$\\
      \> $\rightarrow$ dot name dots\\
    \end{tabbing}

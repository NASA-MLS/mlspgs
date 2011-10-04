pro icfm_unpack_matrixblock, block
    ;; matrix block kind
    m_absent = 0
    m_banded = 1
    m_column_sparse = 2
    m_full = 3
    m_unknown = 4

    flags = pvm_unpack_quantity()
    kind = flags[0]
    nrows = flags[1]
    ncols = flags[2]

    if kind ne m_absent and kind ne m_full then begin
        MyMessage, /error, "Matrix block is neither absent nor full."
    endif
    if kind eq m_absent then begin
        kind = 'A'
    endif else if kind eq m_full then begin
        kind = 'F'
    endif

    if kind eq 'F' then begin
        values = pvm_unpack_quantity()
        block = {kind: kind, values: values}
    endif else block = {kind: kind}
end

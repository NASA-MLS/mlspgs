; This subroutine utilize some subroutines in AtGod
; Use pvm to send input, receive output, and contact
; the Fortran library's Forward Model server (addressable by its tid)
pro ForwardModel, tid, info, mafNo, state, stateExtra, outputTemplate, output, jacobian=jacobian, doDerivative=doDerivative
    if keyword_set(doDerivative) then doDerivative = 1L else doDerivative = 0L

    bufid = pvm_initsend()

    pvm_pack_idltypeinfo, 2L ;sig_fwdmdl
    pvm_pack_idlvariable, 2L

    icfm_sendvector, state, /packonly
    icfm_sendvector, stateExtra, /packonly
    icfm_sendvector, outputTemplate, /packonly

    pvm_pack_idltypeinfo, doDerivative
    pvm_pack_idlvariable, doDerivative

    mafNo = long(mafNo)
    pvm_pack_idltypeinfo, mafNo
    pvm_pack_idlvariable, mafNo

    msgtag = 200L

    if n_elements(tid) eq 1 then info = pvm_send(tid=tid(0), msgtag=msgtag) $
    else info = pvm_mcast(tids=tid, msgtag=msgtag)

    ; using AtMatrix requires using AtGod
    icfm_receivequantities, output, tid=tid, msgtag=msgtag
    ;createvector, output, measurement.name
    ;for i=0, n_tags(quantities)-1 do begin
    ;    addquantity2vector, output, quantities.(i), /overwrite
    ;endfor

    msgtag = 202L ; different for receiving matrix
    if doDerivative then icfm_receivematrix, jacobian, tid, msgtag
end
; $Log$

; This subroutine utilize some subroutines in AtGod
pro ForwardModel, tid, info, mafNo, state, stateExtra, measurement, output
    bufid = pvm_initsend()

    pvm_pack_idltypeinfo, 2L
    pvm_pack_idlvariable, 2L

    icfm_sendvector, state, /packonly
    icfm_sendvector, stateExtra, /packonly
    icfm_sendvector, measurement, /packonly

    mafNo = long(mafNo)
    pvm_pack_idltypeinfo, mafNo
    pvm_pack_idlvariable, mafNo

    msgtag = 200L

    if n_elements(tid) eq 1 then info = pvm_send(tid=tid(0), msgtag=msgtag) $
    else info = pvm_mcast(tids=tid, msgtag=msgtag)

    icfm_receivequantities, output, tid, msgtag
    ;createvector, output, measurement.name
    ;for i=0, n_tags(quantities)-1 do begin
    ;    addquantity2vector, output, quantities.(i), /overwrite
    ;endfor
end

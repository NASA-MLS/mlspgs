; This subroutine commands the server to die
; using a signal sent via pvm
pro ICFM_die, tid, info
    bufid = pvm_initsend()

    pvm_pack_idltypeinfo, 4L ;sig_die
    pvm_pack_idlvariable, 4L

    msgtag=200L

    if n_elements(tid) eq 1 then info = pvm_send(tid=tid(0), msgtag=msgtag) $
    else info = pvm_mcast(tids=tid, msgtag=msgtag)
end
; $Log$

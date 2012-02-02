; This subroutine utilize some subroutines in AtGod
pro ICFM_Cleanup, tid, info
    bufid = pvm_initsend()

    pvm_pack_idltypeinfo, 1L ;sig_cleanup
    pvm_pack_idlvariable, 1L

    msgtag=200L

    if n_elements(tid) eq 1 then info = pvm_send(tid=tid(0), msgtag=msgtag) $
    else info = pvm_mcast(tids=tid, msgtag=msgtag)
end

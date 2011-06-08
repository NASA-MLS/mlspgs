; This subroutine utilize some subroutines in AtGod
pro icfm_receivequantity, quantity, tid=tid, msgtag=msgtag, justunpack=justunpack
    if n_elements(justunpack) eq 0 or (not justunpack) then begin
        if n_elements(tid) eq 0 then MyMessage, /error, "Missing tid"
        if n_elements(msgtag) eq 0 then MyMessage, /error, "Missing msgtag"
        bufid = pvm_recv(tid=tid, msgtag=msgtag)
    endif

    ;l28 = pvm_unpack_quantity()

    ;print, l28

    quantity = unpackvectorquantity()
end

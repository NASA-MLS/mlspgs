; This subroutine utilize some subroutines in AtGod
pro icfm_receivequantities, quantities, tid=tid, msgtag=msgtag, justunpack=justunpack
    if n_elements(justunpack) eq 0 then bufid = pvm_recv(tid=tid, msgtag=msgtag)

    numQty = pvm_unpack_quantity()
    if numQty eq 0 then return

    icfm_receivequantity, qty, /justunpack
    quantities = create_struct('q0', qty)

    for i=1, numQty-1 do begin
        icfm_receivequantity, qty, /justunpack
        quantities = create_struct(quantities, 'q' + strtrim(i,1), qty)
    endfor
end

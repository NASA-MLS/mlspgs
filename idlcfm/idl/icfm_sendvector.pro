; This subroutine utilize some subroutines in AtGod
pro icfm_sendvector, vector, tid=tid, info=info, packonly=packonly

    if n_elements(packonly) eq 0 then begin
        bufid = pvm_initsend()
        pvm_pack_idltypeinfo, 3L ;sig_vector
        pvm_pack_idlvariable, 3L
    endif 

    tags = tag_names(vector)
    numQty = n_elements(tags) - 1 ; NAME is not a quantity

    pvm_pack_idltypeinfo, vector.name
    pvm_pack_idlvariable, vector.name

    pvm_pack_idltypeinfo, numQty
    pvm_pack_idlvariable, numQty

    for i=0,numQty do begin
        if (tags[i] eq 'NAME') then continue 
        icfm_sendquantity, vector.(i), /packonly
    endfor

    msgtag = 200L

    if n_elements(packonly) eq 0 then begin
        ; send data
        if n_elements(tid) eq 1 then info = pvm_send(tid=tid(0), msgtag=msgtag) $
        else info = pvm_mcast(tids=tid, msgtag=msgtag)
    endif
end
; $Log$

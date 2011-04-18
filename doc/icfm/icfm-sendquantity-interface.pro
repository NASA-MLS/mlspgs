pro ICFM_SendQuantity, qty, tid=tid, info=info, msgtag=msgtag, packonly=packonly
    ; qty: quantity to send
    ; tid (optional): server's tid, only needed if packonly is not set
    ; info (optional): output status info
    ; msgtag(optional): should be 200L, but only needed if packonly is not set
    ; packonly (optional): if set, this procedure will only pack the quantity
    ; into the send buffer without sending it
end

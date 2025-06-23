pro icfm_sendvector, vector, tid=tid, info=info, msgtag=msgtag, packonly=packonly
    ; vector: vector to send or pack
    ; tid (optional): server's tid, needed only if packonly is not set
    ; info (optional): output, status report
    ; msgtag (optional): should be 200L, only needed if packonly is not set
    ; packonly (optional): if set, this procedure will only pack the vector
    ; into the send buffer without sending it
end

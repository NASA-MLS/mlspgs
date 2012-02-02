; $Id$

pro icfm_receivematrix, matrix, tid, msgtag

    if n_elements(tid) eq 0 then MyMessage, /error, "Missing tid"
    if n_elements(msgtag) eq 0 then MyMessage, /error, "Missing msgtag"
    bufid = pvm_recv(tid=tid, msgtag=msgtag)

    icfm_unpack_rc, col
    icfm_unpack_rc, row

    numCol = col.numBlocks
    numRow = row.numBlocks

    matrix = creatematrix(row.vector, col.vector, colinstfirst=col.instfirst, rowinstfirst=row.instfirst)
    for i=0,numCol-1 do begin
        for j=0, numRow-1 do begin
            icfm_unpack_matrixblock, block
            a = createblock(matrix, j, i, type=block.kind)
            if block.kind ne 'A' then a.values = block.values
            matrix.blocks(j,i) = ptr_new(a)
        endfor
    endfor
end

; $Log$
; Revision 1.1  2011/09/07 06:33:35  honghanh
; Add 3 more procedures to receive matrix
; Add 'stacked' concept to createquantity, and icfm_sendquantity
; and related procedures.
;

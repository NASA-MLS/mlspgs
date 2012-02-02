pro icfm_unpack_rc, rc
    nb = pvm_unpack_quantity()
    instFirst = pvm_unpack_quantity()
    nelts = pvm_unpack_quantity()
    inst = pvm_unpack_quantity()
    quant = pvm_unpack_quantity()
    icfm_receivequantities, quantities, /justunpack

    rc = {numBlocks:nb, instfirst: instfirst, vector:quantities}
end

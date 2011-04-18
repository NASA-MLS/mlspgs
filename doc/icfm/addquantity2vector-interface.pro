pro AddQuantity2Vector, vector, qty, overwrite=overwrite
    ; vector: the vector to which the quantity is added
    ; qty: the quantity to add, the quantity's name must not be 'name'
    ; overwrite (optional): if false, then the procedure will generate
    ; an error if there is already a quantity of the same name
    ; as qty's in vector. If true, then the procedure will replace
    ; the quantity in vector with qty.
end

pro AddQuantity2Vector, vector, qty, overwrite=overwrite

    if strlowcase(qty.name) eq 'name' then MyMessage, /error, "qty's name must not be 'name'"

    tags = tag_names(vector)
    if n_elements(overwrite) eq 0 then overwrite = 0
    w = where(tags eq strupcase(qty.name))
    if w gt -1 and (not overwrite) then MyMessage, /error, "a quantity of that name has already been added and overwrite is not set"
    if w eq -1 then vector = create_struct([qty.name], qty, vector) else vector.(w) = qty
    
end

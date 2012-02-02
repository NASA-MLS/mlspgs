pro CreateVector, vector, name

    if size(name, /type) ne 7 then MyMessage, /error, "name must be a string"

    vector = {name: name}
end 

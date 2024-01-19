function rename_tags, struct, oldtags, newtags
    tags = tag_names(struct)
    FOR i=0L, n_elements(tags)-1 DO BEGIN
       w=where(oldtags EQ tags[i], nw)
       IF nw NE 0 THEN taguse = newtags[w[0]] ELSE taguse = tags[i]
       IF i EQ 0 THEN newst = create_struct(taguse, struct[0].(i)) $
	       ELSE newst = create_struct(newst, taguse, struct[0].(i))
    ENDFOR
    newstruct = replicate(newst, n_elements(struct))
    FOR i=0L, n_elements(tags)-1 DO newstruct.(i) = struct.(i)
    return, newstruct
END

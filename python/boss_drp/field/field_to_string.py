
def field_to_string(field):
    if field == '*':
        field = str(0).zfill(6)
        field = field.replace('0','?')
        return(field)
    return(str(field).zfill(6))

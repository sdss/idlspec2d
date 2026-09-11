import re

def parse_format_to_regex(fmt):
    """
    Convert a simple parse-style format string to a compiled regex.

    Example: 'spAll-{run2d}-{mjd}_{obs}.parquet'
              becomes a regex with named capture groups:
                    (?P<run2d>...)
                    (?P<mjd>...)
                    (?P<obs>...)
    """
    parts = re.split(r'(\{[^{}]+\})', fmt)

    regex_parts = []

    for part in parts:
        if part.startswith("{") and part.endswith("}"):
            field = part[1:-1]
            regex_parts.append(f"(?P<{field}>.+?)")
        else:
            regex_parts.append(re.escape(part))

    return re.compile("^" + "".join(regex_parts) + "$")

def parse(fmt, value):
    """
    Parse `value` according to a parse-style format.
    Returns a dict on success, None on failure.
    """
    pattern = parse_format_to_regex(fmt)
    match = pattern.fullmatch(value)

    return match.groupdict() if match else None
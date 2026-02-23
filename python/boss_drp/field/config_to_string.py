def config_to_string(config, legacy = False):
    # If `config` is a list or array, handle each element recursively
    if isinstance(config, (list, tuple)):
        return [config_to_string(c) for c in config]

    # Ensure `config` is non-negative
    if config < 0:
        config = 0

    # Format the value
    if not legacy:
        if config < 1000000:
            return f"{int(config):06d}"  # Format as 6-digit integer with leading zeros
        else:
            return str(config)  # Convert to string for values >= 1,000,000
    else:
        return str(config).zfill(4)

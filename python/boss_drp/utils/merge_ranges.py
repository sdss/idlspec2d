import math

def merge_ranges(ranges):
    if not ranges:
        return []

    def left(x):
        return -math.inf if math.isnan(x) else x

    def right(x):
        return math.inf if math.isnan(x) else x

    ranges = sorted(ranges, key=lambda x: left(x[0]))
    merged = [ranges[0].copy()]

    for start, end in ranges[1:]:
        last_start, last_end = merged[-1]

        if left(start) <= right(last_end):
            # Extend right boundary if necessary
            if right(end) > right(last_end):
                merged[-1][1] = end
        else:
            merged.append([start, end])

    return merged
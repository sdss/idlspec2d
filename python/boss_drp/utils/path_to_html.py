import os
import urllib.parse


def path_to_html(path, dir=False):
    if dir and not path.endswith(os.sep):
        path = path+os.sep
    # Normalize the path to use the correct slashes
    normalized_path = os.path.normpath(path)
    # Replace backslashes with forward slashes
    web_friendly_path = normalized_path.replace(os.sep, '/')
    # URL-encode the path to escape special characters
    html_path = urllib.parse.quote(web_friendly_path)
    return html_path

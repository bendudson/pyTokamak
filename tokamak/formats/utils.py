# Useful utilities when reading files

try:
    import re
except:
    print "Regular expression module 're' not available"
    raise

def file_tokens(fp):
    """ A generator to split a file into tokens
    """
    toklist = []
    while True:
        line = fp.readline()
        if not line: break
        toklist = line.split()
        for tok in toklist:
            yield tok

def file_numbers(fp):
    """Generator to get numbers from a text file"""
    toklist = []
    while True:
        line = fp.readline()
        if not line: break
        # Match numbers in the line using regular expression
        pattern = r'[+-]?\d*[\.]?\d+(?:[Ee][+-]?\d+)?'
        toklist = re.findall(pattern, line)
        for tok in toklist:
            yield tok

from __future__ import print_function

import numpy as np
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
import sys

if len(sys.argv) != 2:
    raise RuntimeError('Need exactly one argument')

# Read from argument
filename = sys.argv[1]

tree = ET.parse(filename)
root = tree.getroot()

# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 15:11:20 2021
@author: norab
"""

import sys
import os
sys.path


# Make sure system modules are contained in sys.path
# so modules can be found upon import
abspath = os.path.abspath('.') + '\\Willow1.0'
add_paths = [x[0] for x in os.walk(abspath)]
sys.path = list(set(sys.path + add_paths))

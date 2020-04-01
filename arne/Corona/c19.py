from lxml import html
#!/usr/bin/env python3

import os
from pathlib import Path
import re

#os.system('./c19.bash')
inp_dir=home = str(Path.home())+"/Desktop/Corona/data/"

file_list=[]
for f in os.listdir(inp_dir):
    if (f.find(".txt")!=-1):
        region=re.sub(r'.txt','',f)
        print (region)
        

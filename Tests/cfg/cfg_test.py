#import HTPolyNet.configuration as c
import logging
import yaml
import sys
#logging.basicConfig(filename='cfg_test.log',encoding='utf-8',filemode='w',format='%(asctime)s %(message)s',level=logging.DEBUG)

#C=c.Configuration.read('FDE-DFDA.yaml')
#C.calculate_maximum_conversion()
filename='test.yaml'
with open(filename,'r') as f:
    mydict=yaml.safe_load(f)

with open('reyamled.yaml','w') as f:
    yaml.dump(mydict,f)

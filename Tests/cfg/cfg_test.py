import HTPolyNet.configuration as c
import logging

logging.basicConfig(filename='cfg_test.log',encoding='utf-8',filemode='w',format='%(asctime)s %(message)s',level=logging.DEBUG)

C=c.Configuration.read('FDE-DFDA.yaml')
C.calculate_maximum_conversion()

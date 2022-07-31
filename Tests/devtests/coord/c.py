from  HTPolyNet.coordinates import Coordinates
import logging
import argparse as ap

parser=ap.ArgumentParser()
parser.add_argument('f',type=str,default='FDE-p.gro')
parser.add_argument('-i',type=int,default=0)
parser.add_argument('-j',type=int,default=0)

args=parser.parse_args()
logging.basicConfig(filename='out.log',encoding='utf-8',filemode='w',format='%(asctime)s %(message)s',level=logging.DEBUG)
c=Coordinates.read_gro(args.f)
c.calc_distance_matrix()
print(c.distance_matrix)

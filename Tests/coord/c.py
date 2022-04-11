from  HTPolyNet.coordinates import Coordinates
import logging
import argparse as ap

parser=ap.ArgumentParser()
parser.add_argument('f',type=str,default='npt-1.gro')
args=parser.parse_args()
logging.basicConfig(filename='out.log',encoding='utf-8',filemode='w',format='%(asctime)s %(message)s',level=logging.DEBUG)
c=Coordinates.read_gro(args.f)

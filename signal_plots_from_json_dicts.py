import matplotlib
import matplotlib.pyplot as plt
import h5py
import glob
import json
import argparse
import numpy as np
import twoBytwo_defs
import auxiliary
import signal_characterization as sig_char
from plot_signal_muons import plot_muons
from plot_signal_hadrons import plot_hadrons

def main(scale_factor, muon_json_file, hadron_json_file, wrong_sign):

    if wrong_sign == False:
        sb = 0
    if wrong_sign == True:
        sb = 3

    muon_file = open(muon_json_file)
    muon_dict=json.load(muon_file)

    hadron_file = open(hadron_json_file)
    hadron_dict=json.load(hadron_file)

    plot_muons(muon_dict, scale_factor, sig_bkg = sb)
    plot_hadrons(hadron_dict, scale_factor, sig_bkg = sb)


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-sf', '--scale_factor', default=1, required=True, type=float, \
                        help='''Scale factor related to input dictionary''')
    parser.add_argument('-mu', '--muon_json_file', default=None, required=True, type=str, \
                        help='''string corresponding to the path of the muon info JSON file''')
    parser.add_argument('-had', '--hadron_json_file', default=None, required=True, type=str, \
                        help='''string corresponding to the path of the hadron info JSON file''')
    parser.add_argument('-ws', '--wrong_sign', default=False, required=True, type=bool, \
                        help='''bool corresponding to true for wrong sign background and false for signal''')
    args = parser.parse_args()
    main(**vars(args))

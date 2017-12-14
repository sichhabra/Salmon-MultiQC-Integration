import pandas as pd
import numpy as np


class SeqModel:
    def __init__(self):
        self.obs3_ = None
        self.obs5_ = None
        self.exp3_ = None
        self.exp5_ = None
        self.dims_ = None
        self.valid_ = False

    def populate_model_(self, data_):
        import struct
        from numpy.linalg import norm

        model = None
        offset = 0
        int_struct = struct.Struct('@i')
        long_struct = struct.Struct('@q')

        mspace1 = int_struct.unpack_from(data_[offset:])[0]
        offset += int_struct.size

        mspace2 = int_struct.unpack_from(data_[offset:])[0]
        offset += int_struct.size

        mspace3 = int_struct.unpack_from(data_[offset:])[0]
        offset += int_struct.size

        model_struct1 = struct.Struct('@' + mspace1 * 'i')
        model1 = model_struct1.unpack_from(data_[offset:])
        offset += model_struct1.size

        model_struct2 = struct.Struct('@' + mspace1 * 'i')
        model2 = model_struct2.unpack_from(data_[offset:])
        offset += model_struct2.size

        model_struct3 = struct.Struct('@' + mspace1 * 'i')
        model3 = model_struct3.unpack_from(data_[offset:])
        offset += model_struct3.size

        nrow = long_struct.unpack_from(data_[offset:])[0]
        offset += long_struct.size

        ncol = long_struct.unpack_from(data_[offset:])[0]
        offset += long_struct.size

        model_struct = struct.Struct('@' + nrow * ncol * 'd')
        model4 = model_struct.unpack_from(data_[offset:])
        offset += model_struct.size

        nrow = int_struct.unpack_from(data_[offset:])[0]
        offset += int_struct.size

        ncol = int_struct.unpack_from(data_[offset:])[0]
        offset += int_struct.size

        model_struct4 = struct.Struct('@' + nrow * mspace1 * 'd')
        model = model_struct4.unpack_from(data_[offset:])
        model = np.array(model)
        model = model.reshape(mspace1, nrow).T
        model = (model.T / model.sum(axis=1)).T
        return model

    # dname is the root directory of salmon output
    def from_file(self, dname):
        import os
        import gzip
        obs3_name = os.path.sep.join([dname, 'aux_info', 'obs3_seq.gz'])
        obs5_name = os.path.sep.join([dname, 'aux_info', 'obs5_seq.gz'])
        exp3_name = os.path.sep.join([dname, 'aux_info', 'exp3_seq.gz'])
        exp5_name = os.path.sep.join([dname, 'aux_info', 'exp5_seq.gz'])

        obs3_dat, obs5_dat, exp3_dat, exp5_dat = None, None, None, None
        try:
            with gzip.open(obs3_name) as obs_file:
                obs3_dat = obs_file.read()
            self.obs3_ = self.populate_model_(obs3_dat)
        except IOError:
            print("Could not open file {}".format(obs3_name))
            return False

        try:
            with gzip.open(obs5_name) as obs_file:
                obs5_dat = obs_file.read()
            self.obs5_ = self.populate_model_(obs5_dat)
        except IOError:
            print("Could not open file {}".format(obs5_name))
            return False

        try:
            with gzip.open(exp3_name) as exp_file:
                exp3_dat = exp_file.read()
            self.exp3_ = self.populate_model_(exp3_dat)
        except IOError:
            print("Could not open file {}".format(exp3_name))
            return False

        try:
            with gzip.open(exp5_name) as exp_file:
                exp5_dat = exp_file.read()
            self.exp5_ = self.populate_model_(exp5_dat)
        except IOError:
            print("Could not open file {}".format(exp5_name))
            return False

        self.valid_ = True
        return True

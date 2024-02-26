import numpy as np
import pandas as pd
import os
import uuid

class predictor:
    def __init__(self, download_folder = "downloads"):
        # set const
        self.bpstep_features = {"Shift", "Slide", "Rise", "Tilt", "Roll", "HelT", "Shift-FL", "Slide-FL", "Rise-FL", "Tilt-FL", "Roll-FL", "HelT-FL"}
        self.bp_features = {"MGW", "EP", "Opening", "ProT", "Buckle", "Stagger", "Stretch", "Shear", 
                         "MGW-FL", "Opening-FL", "ProT-FL", "Buckle-FL", "Stagger-FL", "Stretch-FL", "Shear-FL"}
        self.allfeatures = ["MGW", "EP", "Opening", "ProT", "Buckle", "Stagger", "Stretch", "Shear", 
                         "MGW-FL", "Opening-FL", "ProT-FL", "Buckle-FL", "Stagger-FL", "Stretch-FL", "Shear-FL", 
                         "Shift", "Slide", "Rise", "Tilt", "Roll", "HelT", 
                         "Shift-FL", "Slide-FL", "Rise-FL", "Tilt-FL", "Roll-FL", "HelT-FL"]
        self.seqlimit = 1000000
        self.download_folder = download_folder
        bpfile = "/srv/www/deepdnashape/querytables/bp.parquet"
        bpstepfile = "/srv/www/deepdnashape/querytables/bpstep.parquet"
        if os.path.exists(bpfile):
            self.bp_data = pd.read_parquet(bpfile)
        else:
            # load k-mer sequences
            self.bp_data = pd.DataFrame({"Sequence": self.loadkmer("odd")})
            # load feature tables
            for bp_feature in self.bp_features:
                self.bp_data[bp_feature] = np.array(self.load(bp_feature), dtype = "float32")
            self.bp_data = self.bp_data.set_index("Sequence")
            print(self.bp_data)
            self.bp_data.to_parquet(bpfile, compression='snappy')
        
        if os.path.exists(bpstepfile):
            self.bpstep_data = pd.read_parquet(bpstepfile)
        else:
            # load k-mer sequences
            self.bpstep_data = pd.DataFrame({"Sequence": self.loadkmer("even")})
            # load feature tables
            for bpstep_feature in self.bpstep_features:
                self.bpstep_data[bpstep_feature] = np.array(self.load(bpstep_feature), dtype = "float32")
            self.bpstep_data = self.bpstep_data.set_index("Sequence")
            print(self.bpstep_data)
            self.bpstep_data.to_parquet(bpstepfile, compression='snappy')

        print("Finished Loading Database.")

    def query(self, seq: str, feature: str):
        if feature in self.bp_data.columns:
            return self.bp_data[feature].get(seq, np.nan)
        if feature in self.bpstep_data.columns:
            return self.bpstep_data[feature].get(seq, np.nan)
        return np.nan

    def predictSeq(self, seq, feature, layer):
        if feature in self.bpstep_features:
            if layer >= 5:
                return []
            k = layer * 2 + 2
        else:
            k = layer * 2 + 1
        seq = "N" * layer + seq + "N" * layer
        n = len(seq)
        subseqs = [seq[i:i+k] for i in range(n - k + 1)]
        return self.query(subseqs, feature).tolist()

    def loadkmer(self, kmer):
        kmer_file = "/srv/www/deepdnashape/querytables/" + kmer + ".txt"
        seqs = []
        with open(kmer_file) as fin:
            for line in fin:
                seqs.append(line.strip())
        return seqs
    
    def load(self, feature):
        print("Loading:", feature)
        d = []
        table_name = feature + ".txt"
        with open("/srv/www/deepdnashape/querytables/" + table_name) as fin:
            for line in fin:
                d.append(float(line.strip()))
        return d
    
    def is_valid_dna(self, s):
        valid_chars = {'A', 'C', 'G', 'T', 'N'}
        return set(s).difference(valid_chars) == set()

    def predictFile(self, file, feature, layer):
        seqCount = 0
        lines = []
        with open(file) as fin:
            for line in fin:
                seq = line.strip()
                if not self.is_valid_dna(seq):
                    return None
                predicted = list(map(str, self.predictSeq(seq, feature, layer)))
                lines.append(" ".join(predicted))
                seqCount += 1
                if seqCount >= self.seqlimit:
                    return None
        #Write file to download
        outfilename = str(uuid.uuid4()) + ".txt"
        outfilepath = os.path.join(self.download_folder, outfilename)
        with open(outfilepath, "w") as fout:
            for line in lines:
                fout.write(line)
                fout.write("\n")
        return outfilename


if __name__ == "__main__":
    #run test
    pass
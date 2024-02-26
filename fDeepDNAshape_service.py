import numpy as np
import pandas as pd
import os
import uuid
import Pyro4
import logging
from Bio import SeqIO

logging.basicConfig(level=logging.DEBUG)
Pyro4.config.DETAILED_TRACEBACK = True

class predictor:
    def __init__(self, download_folder = "/srv/www/deepdnashape/downloads"):
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

    @Pyro4.expose
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
    
    # @Pyro4.expose
    # def predictFile(self, file, feature, layer):
    #     seqCount = 0
    #     lines = []
    #     with open(file) as fin:
    #         for line in fin:
    #             seq = line.strip()
    #             if not self.is_valid_dna(seq):
    #                 return None
    #             predicted = list(map(str, self.predictSeq(seq, feature, layer)))
    #             lines.append(" ".join(predicted))
    #             seqCount += 1
    #             if seqCount >= self.seqlimit:
    #                 return None
    #     #Write file to download
    #     outfilename = str(uuid.uuid4()) + ".txt"
    #     outfilepath = os.path.join(self.download_folder, outfilename)
    #     with open(outfilepath, "w") as fout:
    #         for line in lines:
    #             fout.write(line)
    #             fout.write("\n")
    #     return outfilename
    
    @Pyro4.expose
    def predictFile(self, file, feature, layer):
        seqCount = 0
        sequences = []
        
        # Determine file format based on extension
        if file.endswith( ".txt" ):
            # Read sequences from a text file
            with open(file, "r") as fin:
                for line in fin:
                    seq = line.strip()
                    sequences.append(("sequence_" + str(seqCount), seq))
                    seqCount += 1
                    if seqCount >= self.seqlimit:
                        break
            output_extension = ".txt"
        elif file.endswith( ".fasta" ) or file.endswith( ".fa" ):
            # Read sequences from a FASTA file using Biopython
            with open(file, "r") as fin:
                for record in SeqIO.parse(fin, "fasta"):
                    sequences.append((record.id, str(record.seq)))
                    seqCount += 1
                    if seqCount >= self.seqlimit:
                        break
            output_extension = ".fasta"
        else:
            return None
        
        predicted_sequences = []
        for identifier, seq in sequences:
            if not self.is_valid_dna(seq):
                return None
            predicted = list(map(str, self.predictSeq(seq, feature, layer)))
            predicted_sequences.append((identifier, predicted))
        
        # Write predictions to a new file with the same format as the input file
        outfilename = str(uuid.uuid4()) + output_extension
        outfilepath = os.path.join(self.download_folder, outfilename)
        with open(outfilepath, "w") as fout:
            if output_extension == ".txt":
                for identifier, predicted in predicted_sequences:
                    fout.write(" ".join(predicted) + "\n")
            elif output_extension == ".fasta":
                for identifier, predicted in predicted_sequences:
                    fout.write(">" + identifier + "\n")
                    fout.write(",".join(predicted) + "\n")
        
        return outfilename

def start_server():
    # Start the Pyro4 server with your exposed class
    daemon = Pyro4.Daemon()                # Create a Pyro daemon
    my_predictor = predictor()
    my_predictor.query("A", "MGW")
    my_predictor.query("A", "Roll")
    uri = daemon.register(my_predictor)         # Register your class as a Pyro object
    ns = Pyro4.locateNS()                  # Find the name server
    ns.register("deepdnashape.db", uri)    # Register your class with a name in the name server

    print("Server is ready.")
    print(f"Object URI: {uri}")
    daemon.requestLoop()                   # Start the event loop of the server to wait for calls

if __name__ == "__main__":
    start_server()
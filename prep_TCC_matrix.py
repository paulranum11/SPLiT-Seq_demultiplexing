import os
import sys, gc

#json_path=os.path.abspath(sys.argv[1])
#if not os.path.isfile(json_path):
#    print("ERROR: Please provide path to a valid config.json file...")
#    print(sys.argv[1])
#    exit(1)
    
#import json
#with open(json_path) as json_file:
#    parameter = json.load(json_file)
#
#print("Number of threads", parameter["NUM_THREADS"])

# Load dataset
import numpy as np
from scipy.sparse import coo_matrix
from sklearn.preprocessing import normalize
import sys
import argparse

# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument('-E', '--inputMatrixec', required=True, help='Input Matrix.ec')
parser.add_argument('-T', '--inputMatrixtsv', required=True, help='Input Matrix.tsv')
parser.add_argument('-O', '--outputDir', required=True, help='Output Dir')
parser.add_argument('-I', '--indexFasta', required=True, help='Provide the index fasta file used to generate the kallisto .idx')
parser.add_argument('-G', '--GenesOrTranscripts', required=False, help='enter "geneIDs" or "transcriptIDs" to return a rownames as gene or transcript IDs')
#/mnt/isilon/davidson_lab/ranum/Tools/Kallisto_Index/Mus_musculus.GRCm38.cdna.all.fa
args = parser.parse_args()

#matrix.ec file
#ecfile_dir = parameter["kallisto"]["TCC_output"]+"matrix.ec"
#tsvfile_dir = parameter["iskallisto"]["TCC_output"]+"matrix.tsv"
exfile_dir = args.inputMatrixec
tsvfile_dir = args.inputMatrixtsv
output_dir = args.outputDir

#exfile_dir = /mnt/isilon/davidson_lab/ranum/Data/SPLiT-Seq/AlignTest/SPLiT-Seq_demultiplexing/kallisto_output/matrix.ec
#tsvfile_dir = /mnt/isilon/davidson_lab/ranum/Data/SPLiT-Seq/AlignTest/SPLiT-Seq_demultiplexing/kallisto_output/matrix.tsv

print("Loading index fasta..")

if (args.GenesOrTranscripts == "transcriptIDs"):
    print("Extracting transcript IDs..")
    transcriptID_Dict = {}
    with open(args.indexFasta, "r") as infile:
        line_ct = 0
        id_ct = 0
        for line in infile:
            if (">" in line):
                split1 = line.split(" ")
                split2 = split1[0].split(">")
                transcriptID = split2[1]
                #print(str(id_ct) + " " + transcriptID)
                transcriptID_Dict[id_ct] = transcriptID
                id_ct = id_ct + 1
            line_ct = line_ct + 1

if (args.GenesOrTranscripts == "geneIDs"):
    print("Extracting gene IDs..")
    geneID_Dict = {}
    with open(args.indexFasta, "r") as infile:
        line_ct = 0
        id_ct = 0
        for line in infile:
            if (">" in line):
                split1 = line.split(" ")
                split2 = split1[3].split(":")
                geneID = split2[1]
                #print(str(id_ct) + " " + transcriptID)
                geneID_Dict[id_ct] = geneID
                id_ct = id_ct + 1
            line_ct = line_ct + 1

print("Loading input matrix..")
matrixTSV_List = []
with open(args.inputMatrixtsv, "r") as infile2:
    for line in infile2:
        split1 = line.split("\t")
        transcriptNum = split1[0]
        matrixTSV_List.append(int(transcriptNum))

print("Loading in equivalence classes..")
equivalenceClass_Dict = {}
with open(args.inputMatrixec, "r") as infile3:
    for line in infile3:
        split1 = line.split("\t")
        split2 = split1[1].rstrip()
        split3 = split2.split(",")
        equivalenceClass_Dict[int(split1[0])] = split3

print("Loading TCCs..")

COOinput = np.loadtxt( tsvfile_dir, delimiter='\t' , dtype=int)
rows,cols,data = COOinput.T
nonzero_ec = np.unique(rows)
map_rows = { val:ind for ind,val in enumerate( nonzero_ec ) }
map_cols = { val:ind for ind,val in enumerate( np.unique(cols) ) }
TCCmatrix   = coo_matrix( (data.astype(float),( [map_rows[r] for r in rows], [map_cols[c] for c in cols]) ) ) 

NUM_OF_CELLS = TCCmatrix.shape[1]
print("NUM_OF_CELLS =", NUM_OF_CELLS)
      
T = TCCmatrix.tocsr()
T_norm = normalize(T, norm='l1', axis=0) 
T_normT = T_norm.transpose()
del TCCmatrix;
_ = gc.collect()


# Pairwise_distances
from sklearn.metrics.pairwise import pairwise_distances
from scipy.spatial.distance import *
from scipy.stats import entropy

def L1_distance(p,q):
    return cityblock(p,q).sum()

# def jensen_shannon(p, q):
#     m=0.5*p+0.5*q
#     p = np.transpose(p[p > 0])
#     q = np.transpose(q[q > 0])
#     m = np.transpose(m[m > 0])
#     return np.sqrt(entropy(m)-0.5*entropy(q)-0.5*entropy(p))


#num_of_threads = parameter["NUM_THREADS"]
num_of_threads = 5
print("Calculating pairwise L1 distances... ( num_threads =",num_of_threads,")")

# D_js = pairwise_distances(T_normT,metric=jensen_shannon,n_jobs=num_of_threads)
D_l1 = pairwise_distances(T_normT,metric=L1_distance,n_jobs=num_of_threads)

print("writing data...")

#Save data
T = T.todense()
#D_l1 = D_l1.todense()
#nonzero_ex = nonzero_ex.todense()

with open(output_dir+"2_rowNames.txt", 'w') as f:
    for row in map_rows:
        print(row, file=f)

#with open(output_dir+"3_colNames.txt", 'w') as f:
#    for row in map_cols:
#        print(row, file=f)

if (args.GenesOrTranscripts == "transcriptIDs"):
    with open(output_dir+"3_transcriptIDs.txt", 'w') as f:
        for row in map_rows:
            EQgenelist = equivalenceClass_Dict[row]
            Plist = []
            for EQkey in EQgenelist:
                Plist.append(transcriptID_Dict[int(EQkey)]) 
            print('[%s]' % ', '.join(map(str, Plist)), file=f)  
            #print('[%s]' % ', '.join(map(str, EQgenelist)))

if (args.GenesOrTranscripts == "geneIDs"):
    with open(output_dir+"3_geneIDs.txt", 'w') as f:
        for row in map_rows:
            EQgenelist = equivalenceClass_Dict[row]
            Plist = []
            for EQkey in EQgenelist:
                Plist.append(geneID_Dict[int(EQkey)])
            print('[%s]' % ', '.join(map(str, Plist)), file=f)
            #print('[%s]' % ', '.join(map(str, EQgenelist)))

with open(output_dir+"1_expressionMatrix.txt", 'wb') as f:
    np.savetxt(f,T, delimiter="\t")
#with open(output_dir+"2_colNames.txt", 'wb') as f:
    
#with open(output_dir+"pwise_dist_L1.dat", 'wb') as f:
#    np.savetxt(f,D_l1, delimiter="\t")
#with open(output_dir+"nonzero_ec.dat", 'wb') as f:
#    np.savetxt(f,nonzero_ec, delimiter="\t")  
print("DONE.")


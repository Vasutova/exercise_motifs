import itertools
import numpy as np
from numpy.matlib import zeros
from Bio import motifs
from Bio.Seq import Seq
import random
from Bio import SeqIO


lecture_dna = [
    "TGACGTATAAGTTGCGATGGACGAGATAGCAGAGAATAGGCAACGAGAGATAAGCAG",
    "GACGGTAGCAGATAGACAGATGAAGAGTATGAATTGCACAGATAGCAGATAGCAGAT",
    "GGAGTGTGACGTAGCAGAGACGAAAGACGTAGAGTAGCAGTAGCAGATAGAGGGAGT",
    "TAGACAGTATAGAGACAGCGAGTCGGATAGCACCCAGTATGACGATAGCAATGACAG",
    "GCAGTAGAGCAGATTAGCATTGACAGATAGACGATTGGAGAGATGTGTGGATGACGA",
    "GGCAGGTAGCACACTGGGTCGATAAAGAGTAGCATAGAGACATAGACATATTTTAGC",
]


def count_matrix(motifs):
    a_list = [0] * len(motifs[0])
    t_list = [0] * len(motifs[0])
    c_list = [0] * len(motifs[0])
    g_list = [0] * len(motifs[0])

    for i in range(len(motifs)):
        mot = motifs[i]
        for j in range(len(motifs[0])):
            if mot[j] == "A":
                a_list[j] += 1
            elif mot[j] == "G":
                g_list[j] += 1
            elif mot[j] == "C":
                c_list[j] += 1
            else:
                t_list[j] += 1
    return {
        "A": a_list,
        "C": c_list,
        "G": g_list,
        "T": t_list
    }





def score(motifs):
    counts = count_matrix(motifs)
    score = 0

    for j in range(len(motifs[0])):
        max_count = 0

        for base in "ACGT":
            if counts[base][j] > max_count:
                max_count = counts[base][j]

        score += max_count

    return score


def consensus(motifs):
    counts = count_matrix(motifs)
    result = ""

    for j in range(len(motifs[0])):
        best_base = "A"
        best_count = counts["A"][j]

        for base in "ACGT":
            if counts[base][j] > best_count:
                best_base = base
                best_count = counts[base][j]

        result += best_base

    return result



def hamming_distance(a, b):
    count = 0
    for i in range(len(a)):
        if a[i] != b[i]:
            count += 1
    return count

def total_distance(pattern, sequences):
    total = 0
    l = len(pattern)

    for sequence in sequences:
        min_distance = l

        for i in range(len(sequence) - l + 1):
            lmer = sequence[i:i + l]
            distance = hamming_distance(pattern, lmer)

            if distance < min_distance:
                min_distance = distance

        total += min_distance

    return total


red = ["TAAGTT", "TGAATT", "GGAGTG", "CGAGTC", "TGTGTG", "TGGGTC"]  # slide 19
best = ["AGATAG", "AGATAG", "AGATAG", "AGACAG", "AGATAG", "AGGTAG"]

#print(count_matrix(red))


class MotifProfile:
    def __init__(self, motifs, pseudocount=1):
        self.l = len(motifs[0])
        motif_number = len(motifs)

        counts = count_matrix(motifs)

        new_A = []
        new_C = []
        new_G = []
        new_T = []

        for i in range(self.l):
            new_A.append((counts["A"][i] + pseudocount) / (motif_number + 4 * pseudocount))
            new_C.append((counts["C"][i] + pseudocount) / (motif_number + 4 * pseudocount))
            new_G.append((counts["G"][i] + pseudocount) / (motif_number + 4 * pseudocount))
            new_T.append((counts["T"][i] + pseudocount) /(motif_number + 4 * pseudocount))

        self.ppm = {
            "A": new_A,
            "C": new_C,
            "G": new_G,
            "T": new_T}


    def lmer_probability(self, lmer):
        probability = 1

        for i in range(self.l):
            probability *= self.ppm[lmer[i]][i]

        return probability

    def most_probable_lmer(self, sequence):
        best_lmer = sequence[:self.l]
        best_probability = self.lmer_probability(best_lmer)

        for i in range(1, len(sequence) - self.l + 1):
            lmer = sequence[i:i + self.l]
            probability = self.lmer_probability(lmer)

            if probability > best_probability:
                best_lmer = lmer
                best_probability = probability

        return best_lmer

    def consensus(self):
        result = ""

        for i in range(self.l):
            best_base = "A"
            best_probability = self.ppm["A"][i]

            for base in "ACGT":
                if self.ppm[base][i] > best_probability:
                    best_base = base
                    best_probability = self.ppm[base][i]

            result += best_base

        return result


profile = MotifProfile(["ATCCGTA", "GTGCATA", "AAGCGTA", "ATGCGTG"])
print(profile.consensus())                       # ATGCGTA
print(round(profile.lmer_probability("ATGCGTA"), 4))  # 0.0122

two = MotifProfile(["GTAC", "TTAA"])
print(two.most_probable_lmer("ACTGGATGACCC"))    # TGAC
print(round(two.lmer_probability("TGAC"), 4))         # 0.0093






bio = motifs.create([Seq(site) for site in ["ATCCGTA", "GTGCATA", "AAGCGTA", "ATGCGTG"]])
bio.pseudocounts = 1
print(bio.consensus)     # ATGCGTA
print(bio.pwm["A"])      # the same numbers as your profile.ppm["A"]



rng = random.Random(1)                 # a random number generator with seed 1
i = rng.randint(0, 50)                 # random integer, 0 <= i <= 50 (both ends included!)
lmer = rng.choice(["ACG", "CGT", "GTA"])   # one random item of a list


# TASK 3
sequences = [str(record.seq) for record in SeqIO.parse("planted_motif.fasta", "fasta")]


class MotifFinder:
    def __init__(self, sequences, l, seed=None):
        self.l = l
        self.sequences = sequences
        self.rng = random.Random(seed)
        self.windows = []
        
        for sequence in sequences:
            sequence_windows = []

            for i in range(len(sequence) - l + 1):
                sequence_windows.append(sequence[i:i + l])

            self.windows.append(sequence_windows)


    def total_distance(self, pattern):
        total = 0
        l = len(pattern)
        min_distance = l
        for seq in self.windows:
            distance = hamming_distance(pattern, seq)
            if distance < min_distance:
                 min_distance = distance

            total += min_distance

        return total
  
    def median_string(self):
        itertools.product("ACGT", repeat=1)
   #def randomized_search(self):

   #def best_of(self, runs):

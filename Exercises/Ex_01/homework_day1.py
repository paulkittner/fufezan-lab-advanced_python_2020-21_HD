import collections
import csv
import seaborn as sns
import matplotlib.pyplot as plt
import warnings
import re
import pandas as pd

human = "data/human.fasta"
bacteria = "data/bacteria.fasta" #(Prochlorococcus marinus)
archea = "data/archea.fasta" #(N(atronomonas salsuginis)
plantae = "data/plantae.fasta" #(Arabidopsis thaliana)
animalia  = "data/animalia.fasta" #(Mus musculus)


def aa_counter(data):
    '''
    Function to go over every input line, decide wether its a name or not,
    count the aas in the sequence and its length for each protein.

    :param data: fasta data sequence
    :return:
        collections.Counter(counted_values) : total counted aas in the genome
        names : names of all the ind. proteins
        counts : aas counts of all the ind. proteins
    '''
    counted_values = ""
    names = []
    counts = []
    with open(data) as file:
        file_dict = csv.reader(file)
        for line in file_dict:
            if ">" in line[0]:
                names.append(line)
                length = ''
            else:
                counted_values += line[0]
                length += line[0]
            counts.append(collections.Counter(length))

    return collections.Counter(counted_values), names, counts


def plot_aminoacid_histogram(as_counter):
    '''
    Function to plot the amounts of aas discovered in the given genome

    :param as_counter: total aas count
    :return: histogram
    '''
    aminoacids = ['G', 'P', 'A', 'V', 'L', 'I', 'M', 'C', 'F', 'Y', 'W', 'H', 'K', 'R', 'Q', 'N', 'E', 'D', 'S', 'T']
    aminoacid_count = []
    for aminoacid in aminoacids:
        aminoacid_count.append(as_counter[aminoacid])

    plt.figure(figsize=(10, 6))
    plt.grid(color='w', linestyle='--', linewidth=1, zorder=0)
    ax = plt.axes()
    ax.set_facecolor("#ADDFFF")
    ax.set_axisbelow(True)

    plt.bar(aminoacids, aminoacid_count, width=0.7, color='#0000A0', zorder=1)

    plt.xlabel('1-Letter name for the amino acids', fontsize=14)
    plt.ylabel('Count of the amino acids (#)', fontsize=14)
    plt.title('Amino Acid Distribution', fontsize=20)
    plt.tight_layout();

def create_csv(aas):
    with open('Hausaufgabe1.csv', 'w') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(["aa, count"])
        for key, value in aas.items():
            writer.writerow([key, value])


def calc_protein_attributs(full_names, as_counts, organism):
    '''

    :param full_names: the list with the full names taken from the fasta set
    :param as_counts: The amino acid counts for each protein
    :param organism: the name of the inspected organism
    :return: Dataframe with all the protein names, their aas count and total length
    '''
    protein_names = []
    amino_acid_count = []
    empty_counter = collections.Counter('')
    protein_length = []

    organism_names = {'human': 'HUMAN(.+?)OS=', 'plantae': 'ARATH(.+?)OS=', 'archea': '9EURY(.+?)OS=',
                      'animalia': 'MOUSE(.+?)OS=', 'bacteria': 'PROMT(.+?)OS='}

    keys = list(organism_names.keys())
    values = list(organism_names.values())
    str_to_search = (values[keys.index(organism)])

    for i in range(0, len(full_names)):
        protein_name = re.search(str_to_search, str(full_names[i]))
        if protein_name:
            protein_names.append(protein_name.group(1))

    for i in range(1, len(as_counts)):
        if as_counts[i] == empty_counter:
            amino_acid_count.append(as_counts[i - 1])
            protein_length.append(sum(as_counts[i - 1].values()))

    amino_acid_count.append(as_counts[int(len(as_counts)) - 1])
    protein_length.append(sum(as_counts[int(len(as_counts)) - 1].values()))

    protein_info = {'name': protein_names, 'AS count': amino_acid_count, 'length': protein_length}
    protein_attributs_df = pd.DataFrame(protein_info)

    return protein_attributs_df

if __name__ == '__main__':
    aas, full_names, as_counts = aa_counter(human)
    plot_aminoacid_histogram(aas)
    plt.show()
    protein_attributs_df = calc_protein_attributs(full_names, as_counts, 'human')
    longest_protein = protein_attributs_df.sort_values('length', axis=0, ascending=False).head(1)
    shortest_protein = protein_attributs_df.sort_values('length', axis=0, ascending=True).head(1)
    print('the longest protein of the organism is:\n', longest_protein, '\n',
          'the shortest protein of the oragnism is:\n', shortest_protein, )


import argparse
import plotly.io as pio
pio.renderers.default = "browser"
import csv
from collections import deque
import pandas as pd
import requests
import plotly.graph_objects as go
import re

class Protein_attributions():
    """
    """
    def __init__(self, entry_identifier):
        self.entry_identifier = entry_identifier
        self.protein_name = None
        self.seq = None
        self.hydro_dict = None

    def get_sequence(self):
        '''

        :param entry_identifier: the identity code of the protein you would like to analyse
                                (the specific position in the uniport databank)
        :param protein_name: The name of the protein you want to analyse
        :return:
        seq: The sequence of the protein

        '''
        url = 'https://www.uniprot.org/uniprot/' + self.entry_identifier + '.fasta?fil=reviewed:yes'
        r = requests.get(url)
        seq = ""
        seq_list = []
        protein_csv = '../data/' + self.entry_identifier + '.fasta'

        with open(protein_csv, 'wb') as file:
            file.write(r.content)
            file.close()
        with open(protein_csv) as file:
            file_dict = csv.reader(file, delimiter="\n")
            for line in file_dict:
                if ">" not in line[0]:
                    seq = seq + line[0]
                else:
                    protein_name = re.search('HUMAN(.+?)OS=', str(line))
                    if protein_name:
                        self.protein_name = (protein_name.group(1))
                        seq_list.append(seq)
                    seq = ""
        self.seq = seq
        return seq

    def get_hydropathy_dict(self, file):
        df = pd.read_csv(file)
        aa_names = df['1-letter code'].tolist()
        aa_pi_values = df['hydropathy index (Kyte-Doolittle method)'].to_numpy()
        hydropathy = dict(zip(aa_names, aa_pi_values))
        self.hydro_dict = hydropathy
        print(hydropathy)
        return hydropathy


    def calc_hydro_avg(self, window_length):
        '''

        :param seq: The sequence of the selected protein
        :param hyrdo_dict: dictonary with the aas single letter names and their corresponing hydropathy
        :param window_length: The length the stride should look at the as
        :return:
        avg_hydro: average hydropathy for every AS at their position
        '''
        avg_hydro = []
        window_entries = deque([], window_length)
        print(self.seq)
        for aminoacid in self.seq:
            hydropathy_value = self.hydro_dict[aminoacid]
            window_entries.append(hydropathy_value)
            avg_hydro.append(sum(window_entries) /len(window_entries))

        print(avg_hydro)
        return avg_hydro

    def plot_hydro_avg(self, hydro_avg, window_length):
        '''

        :param sequence: The sequence of the selected protein
        :param hydro_avg: average hydropathy for every AS at their position
        :param window_length: The length the stride should look at the as
        :return: the bar plot
        '''
        data = []
        data.append(
                go.Bar(
                    x = list(range(len(self.seq))),
                    y = hydro_avg
                )
            )
        fig = go.Figure(data=data)
        title = self.protein_name + ' average amino acid hydropathy over ' + str(window_length) + ' aas'
        fig.update_layout(
            title_text=title,
            title_x=0.5,
            title_font_size=42,
            font_size=18)
        fig.update_layout(
            xaxis_title="Position in the sequence",
            yaxis_title="Hydropathy",
        )
        fig.update_layout(template='plotly_dark')
        fig.show()



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--entry_identifier", help="The identity of the Protein (position in uniport ds)", type=str, default="NaN", required=True)
    parser.add_argument("--protein_name", help="Name of the protein", type=str, default="NaN", required=True)
    parser.add_argument("--window_length", help="Length of the sliding window.", type=int, required=True)
    parser.add_argument("--directory_aa_prop", help="where is the file @ bro?",type=str, required=True)
    args = parser.parse_args()

    as_prop = args.directory_aa_prop
    entry_identifier = args.entry_identifier
    window_length = args.window_length
    # protein_name = args.protein_name

    p1 = Protein_attributions("P32249")

    dict_aa_prop = p1.get_hydropathy_dict(as_prop)
    hydro_avg = p1.calc_hydro_avg(window_length)
    p1.plot_hydro_avg(hydro_avg, window_length)

    # python homework_paulkittner_02.py --entry_identifier "P32249" --protein_name "human_gpcr" --window_length 5 --directory_aa_prop "./data/amino_acid_properties.csv"
import pandas as pd
import logomaker
import json
import matplotlib.pyplot as plt
import numpy as np
import os

def create_logo(read_file, output_png):
    with open(read_file,'r') as file:
        data = json.load(file)

    logos = data["logo"]

    aa_list = set(["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"])

    aa_dict = {}

    for aa in aa_list:
        aa_dict[aa] = list()


    for num in range(len(logos) - 1):
        num += 1
        col = logos[num]

        seen = set()
        for pair in col:
            code = pair["code"]
            seen.add(code)
            bit = pair["bits"]
            aa_dict[code].append(bit)
        
        diff = aa_list - seen
        for aa_unseen in diff:
            aa_dict[aa_unseen].append(0.0)

    print(aa_dict)

    df = pd.DataFrame(aa_dict)

    logomaker.Logo(df)

    plt.savefig(output_png, format='png')


if __name__ == "__main__":
    # Set up command-line argument parsing

    for filename in os.listdir("./todo"):
        try:
            create_logo("./" + filename[0:-7]+"/params.json", "./" + filename[0:-7]+"/seq_logo.png")
        except FileNotFoundError:
            print(f"The file {filename} was not found.")



    
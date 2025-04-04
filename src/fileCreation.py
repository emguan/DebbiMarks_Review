import gzip
import sys
import argparse
import os

def decompress_gz(gz_file_path, txt_file_path_seq, txt_file_path_weight):
    try:
        print(os.getcwd())
        # Ensure the output directory exists
        os.makedirs(os.path.dirname(txt_file_path_seq), exist_ok=True)

        # Read from the .gz file and write to the .txt file
        with gzip.open(gz_file_path, 'rt') as gz_file, open(txt_file_path_seq, 'w') as txt_file, open(txt_file_path_weight, 'w') as weight_file:
            count = 0
            for line in gz_file:
                if count == 0:
                    count += 1
                    continue
                l = line.split()
                txt_file.write("> a\n" + l[0] + "\n")
                weight_file.write(l[2] + "\n")
        
    except FileNotFoundError:
        print(f"The file {gz_file_path} was not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    # Set up command-line argument parsing

    for filename in os.listdir("./todo"):
        decompress_gz("./todo/" + filename, "./" + filename[0:-7]+"/output.txt", "./" + filename[0:-7]+"/weights.txt")


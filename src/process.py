import gzip
import sys
import argparse
import os

def decompress_gz(gz_file_path, txt_file_path):
    try:
        # Ensure the output directory exists
        os.makedirs(os.path.dirname(txt_file_path), exist_ok=True)

        # Read from the .gz file and write to the .txt file
        with gzip.open(gz_file_path, 'rt') as gz_file, open(txt_file_path, 'w') as txt_file:
            for line in gz_file:
                l = line.split()
                txt_file.write("> a\n" + l[0] + "\n")
        print(f"Content successfully written to {txt_file_path}")
    except FileNotFoundError:
        print(f"The file {gz_file_path} was not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    # Set up command-line argument parsing

    gz_file = './todo copy/NGS-166.00.txt.gz'  # Replace with your .gz file path
    txt_file = "./NGS-166.00/output.txt"  # Replace with your desired .txt file path


    # Call the decompression function with the provided file paths
    decompress_gz(gz_file, txt_file)

#!/usr/bin/env python3
import argparse
import random
import math
import numpy as np
import struct

# Standard 20 amino acids (used if the file’s alphabet isn’t correct)
DEFAULT_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"

def load_plmc_model(filename):
    """
    Load the Potts model parameters from a binary model file produced by plmc
    (using OutputParametersFull). This function extracts:
      - nSites: number of sites (sequence length)
      - nCodes: alphabet size
      - alphabet: read as nCodes bytes
      - h: sitewise parameters (fields), shape (nSites, nCodes)
      - J: couplings, shape (nSites, nSites, nCodes, nCodes)
    """
    with open(filename, 'rb') as f:
        # 1. Read header values (all little-endian)
        nSites = struct.unpack('<i', f.read(4))[0]
        nCodes = struct.unpack('<i', f.read(4))[0]
        nSeqs = struct.unpack('<i', f.read(4))[0]
        nSkippedSeqs = struct.unpack('<i', f.read(4))[0]
        maxIter = struct.unpack('<i', f.read(4))[0]
        theta = struct.unpack('<f', f.read(4))[0]
        lambdaH = struct.unpack('<f', f.read(4))[0]
        lambdaE = struct.unpack('<f', f.read(4))[0]
        lambdaGroup = struct.unpack('<f', f.read(4))[0]
        nEff = struct.unpack('<f', f.read(4))[0]
        
        # 2. Read alphabet (nCodes bytes)
        alphabet_bytes = f.read(nCodes)
        try:
            alphabet = alphabet_bytes.decode('ascii')
        except UnicodeDecodeError:
            alphabet = alphabet_bytes.decode('latin1')
        
        # 3. Skip weights: (nSeqs + nSkippedSeqs) floats
        total_seq = nSeqs + nSkippedSeqs
        f.read(4 * total_seq)
        
        # 4. Skip target sequence: nSites bytes
        f.read(nSites)
        
        # 5. Skip offset map: nSites * 4 bytes
        f.read(4 * nSites)
        
        # 6. Skip sitewise marginals fi: nSites * nCodes floats
        f.read(4 * nSites * nCodes)
        
        # 7. Read sitewise parameters hi (fields): nSites * nCodes floats
        h_data = f.read(4 * nSites * nCodes)
        h = np.frombuffer(h_data, dtype='<f4').reshape((nSites, nCodes))
        
        # 8. Skip pairwise marginals fij:
        count = (nSites * (nSites - 1) // 2) * (nCodes * nCodes)
        f.read(4 * count)
        
        # 9. Read couplings eij: same count of floats
        couplings_data = f.read(4 * count)
        couplings_flat = np.frombuffer(couplings_data, dtype='<f4')
        # Reassemble into full symmetric coupling tensor of shape (nSites, nSites, nCodes, nCodes)
        J = np.zeros((nSites, nSites, nCodes, nCodes), dtype=np.float32)
        idx = 0
        for i in range(nSites - 1):
            for j in range(i + 1, nSites):
                block = couplings_flat[idx: idx + nCodes * nCodes].reshape((nCodes, nCodes))
                J[i, j, :, :] = block
                J[j, i, :, :] = block  # enforce symmetry
                idx += nCodes * nCodes
                
    return nSites, nCodes, alphabet, h, J

def energy(seq, h, J):
    """
    Compute the Potts model energy of a sequence.
    seq: list of integers (indices) of length nSites.
    h: fields array of shape (nSites, nCodes).
    J: coupling tensor of shape (nSites, nSites, nCodes, nCodes).
    """
    E = 0.0
    L = len(seq)
    # Single-site contributions
    for i in range(L):
        E += h[i, seq[i]]
    # Pairwise contributions (sum over i < j)
    for i in range(L):
        for j in range(i + 1, L):
            E += J[i, j, seq[i], seq[j]]
    return E

def metropolis_sample(h, J, alphabet, num_samples=100, steps_per_sample=1000, T=2.0):
    """
    Use Metropolis-Hastings sampling to generate synthetic sequences.
    Returns a list of sequences (as strings) using the provided alphabet.
    """
    L, q = h.shape
    idx_to_aa = {i: aa for i, aa in enumerate(alphabet)}
    
    # Initialize current sequence randomly
    current = [random.randint(0, q - 1) for _ in range(L)]
    current_E = energy(current, h, J)
    
    samples = []
    for s in range(num_samples):
        for _ in range(steps_per_sample):
            i = random.randint(0, L - 1)  # choose a random position
            new_aa = random.randint(0, q - 1)
            if new_aa == current[i]:
                continue
            proposed = current.copy()
            proposed[i] = new_aa
            E_new = energy(proposed, h, J)
            dE = E_new - current_E
            accept_prob = math.exp(-dE / T) if dE > 0 else 1.0
            if random.random() < accept_prob:
                current = proposed
                current_E = E_new
        # Map indices to letters to get the sequence string
        seq_str = ''.join(idx_to_aa[i] for i in current)
        samples.append(seq_str)
    return samples

def write_fasta(samples, filename):
    """
    Write the list of sequences to a FASTA file.
    """
    with open(filename, 'w') as f:
        for i, seq in enumerate(samples):
            f.write(f">sample_{i}\n{seq}\n")

def main():
    parser = argparse.ArgumentParser(
        description="Generate synthetic sequences from a binary plmc model file using MCMC sampling."
    )
    parser.add_argument("--model", required=True, help="Path to the binary model file (output.model produced by plmc)")
    parser.add_argument("--num_samples", type=int, default=100, help="Number of synthetic sequences to generate")
    parser.add_argument("--steps", type=int, default=1000, help="Number of MCMC steps per sample")
    parser.add_argument("--T", type=float, default=2.0, help="Temperature parameter for MCMC sampling")
    parser.add_argument("--output", required=True, help="Output FASTA file for synthetic sequences")
    args = parser.parse_args()

    print("🔍 Loading binary plmc model...")
    nSites, nCodes, alphabet, h, J = load_plmc_model(args.model)
    print(f"Model loaded: nSites = {nSites}, nCodes = {nCodes}, alphabet = {alphabet}")

    print("🧬 Sampling synthetic sequences...")
    samples = metropolis_sample(h, J, alphabet, num_samples=args.num_samples, steps_per_sample=args.steps, T=args.T)

    print(f"💾 Writing {len(samples)} sequences to {args.output}...")
    write_fasta(samples, args.output)
    print("✅ Done.")

if __name__ == "__main__":
    main()
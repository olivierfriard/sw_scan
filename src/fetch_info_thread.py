"""
Retrieve sequence info from NCBI with the Entrez utility using threads

usage:
python3 fetch_info_threaded.py FILE

FILE is a file with one Accession number by line

Must be used with python >= 3.8
"""

import sys
import time
import threading
from queue import Queue
from Bio import Entrez, GenBank

Entrez.email = "olivier.friard@unito.it"


# Funzione per elaborare una singola accessione
def fetch_info(accession, output_queue):
    try:
        handle = Entrez.efetch(
            db="nucleotide", id=accession, rettype="gb", retmode="text"
        )
        record = GenBank.read(handle)
        result = (
            f"{accession}\t{record.accession[0]}\t{record.size}\t{record.definition}"
        )
        output_queue.put(result)
    except Exception as e:
        output_queue.put(f"ERROR with {accession}: {str(e)}")


# Funzione principale
def main(input_file):
    # Coda per memorizzare i risultati
    output_queue = Queue()
    threads = []

    # Lettura delle accessioni dal file
    with open(input_file, "r") as f_in:
        accessions = [line.strip() for line in f_in if line.strip()]

    # Creazione di un numero limitato di thread
    max_threads = 3  # Limite massimo di thread simultanei
    for accession in accessions:
        while threading.active_count() > max_threads:
            time.sleep(1)  # Attendi finché non ci sono abbastanza thread disponibili
        thread = threading.Thread(target=fetch_info, args=(accession, output_queue))
        threads.append(thread)
        thread.start()

    # Attesa della terminazione di tutti i thread
    for thread in threads:
        thread.join()

    # Stampa dei risultati
    while not output_queue.empty():
        print(output_queue.get())


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 fetch_info_threaded.py FILE")
        sys.exit(1)

    main(sys.argv[1])
